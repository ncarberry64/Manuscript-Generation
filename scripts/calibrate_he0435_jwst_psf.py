"""Execute a frozen field-only calibration; never download or fit lens pixels.

This runner requires an explicit field-only execution flag and the exact frozen
protocol hash. All four native product names, this runner, its helper, and the
algorithm JSON must be bound in that protocol before any calibrated SCI array
is accessed. Failure never authorizes target-pixel measurement or error inflation.
"""
from __future__ import annotations

import argparse
from datetime import datetime, timezone
import hashlib
import json
from pathlib import Path
import sys

import numpy as np

ROOT=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT/'code'))
from strong_lensing_field_psf import (DEFAULT_ALGORITHM,field_mask_from_gwcs,
    guarded_field_arrays,select_field_stars,assign_physical_stars,
    common_pool_dither_diagnostic,build_empirical_psf,fit_heldout_field_star)
from strong_lensing_pixel_geometry import gwcs_centroid,combine_exposures

SOURCES=['scripts/calibrate_he0435_jwst_psf.py','code/strong_lensing_field_psf.py',
         'code/strong_lensing_pixel_geometry.py']


def sha(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def dump(path,value):
    Path(path).write_text(json.dumps(value,indent=2,allow_nan=False)+'\n',encoding='utf-8')


def validate_execution(protocol_path,expected_protocol_sha256,authorized):
    """Reject absent authorization before examining protocol or data paths."""
    if not authorized:
        raise PermissionError('Field-only execution must be explicitly authorized; no native data opened')
    path=Path(protocol_path)
    if sha(path)!=expected_protocol_sha256:
        raise ValueError('Frozen root protocol SHA256 mismatch')
    protocol=json.loads(path.read_text(encoding='utf-8'))
    binding=protocol['field_psf_calibration']
    algorithm_path=ROOT/binding['algorithm_file']
    if sha(algorithm_path)!=binding['algorithm_sha256']:
        raise ValueError('Frozen field-PSF algorithm SHA256 mismatch')
    cfg=json.loads(algorithm_path.read_text(encoding='utf-8'))
    if cfg!=DEFAULT_ALGORITHM:
        raise ValueError('Frozen algorithm differs from reviewed helper defaults')
    for rel in SOURCES:
        if sha(ROOT/rel)!=binding['source_sha256'][rel]:
            raise ValueError('Field calibration source changed after freeze: '+rel)
    products=binding['native_products']
    if len(products)!=4 or len(set(products))!=4 or any(Path(p).name!=p or not p.endswith('_nrcb3_cal.fits') for p in products):
        raise ValueError('Exactly four distinct native NRCB3 CAL basenames required')
    return cfg,products


def _field_stamp(frame,row,star_id):
    x0,y0,x1,y1=row['stamp_bounds_xyxy']
    if not np.all(frame['field_mask'][y0:y1,x0:x1]):
        raise ValueError('Calibration stamp intersects protected target region')
    return {'id':star_id,'data':frame['data'][y0:y1,x0:x1].copy(),
            'error':frame['error'][y0:y1,x0:x1].copy(),
            'dq':frame['dq'][y0:y1,x0:x1].copy(),
            'center':np.asarray(row['xy'])-np.array([x0,y0])}


def calibrate(protocol_path,expected_protocol_sha256,data_directory,output_directory,authorized=False):
    cfg,products=validate_execution(protocol_path,expected_protocol_sha256,authorized)
    # Imports that can touch native image datamodels occur only after the guard.
    from stdatamodels.jwst import datamodels
    from astropy.io import fits
    import photutils
    out=Path(output_directory);out.mkdir(parents=True,exist_ok=True)
    lock=out/'calibration_started.json'
    stamp={'started_utc':datetime.now(timezone.utc).isoformat(),
           'protocol_sha256':expected_protocol_sha256,'native_products':products,
           'target_pixels_authorized':False,'field_only_authorized':True}
    with lock.open('x',encoding='utf-8') as stream:
        json.dump(stamp,stream,indent=2);stream.write('\n')
    report={'status':'FIELD_CALIBRATION_STARTED_TARGETS_REMAIN_SEALED',
            'protocol_sha256':expected_protocol_sha256,'algorithm':cfg,
            'source_sha256':{rel:sha(ROOT/rel) for rel in SOURCES},
            'photutils_version':photutils.__version__,'native_products':products,
            'input_sha256':{},'selection':{},'psf_files_by_exposure':{},'file_sha256':{},
            'overall_adequacy_pass':False,'target_pixel_values_inspected':False,
            'conditional_uncertainty_note':'Known native ERR and a fixed fitted field ePSF; finite calibration-star noise, background calibration, spatial/chromatic PSF variation and absolute GWCS uncertainty are not marginalized. No residual variance rescaling, covariance floor or inflation is allowed.'}
    frames={}
    try:
        for name in products:
            path=Path(data_directory)/name
            report['input_sha256'][name]=sha(path)
            with datamodels.open(path) as model:
                wcs=model.meta.wcs
                field=field_mask_from_gwcs(model.data.shape,wcs,cfg['target_center_icrs_deg'],cfg['exclude_target_radius_arcsec'])
                data,error,dq=guarded_field_arrays(model.data,model.err,model.dq,field)
                catalogue,diagnostic=select_field_stars(data,error,dq,field,wcs,cfg)
                # Copy WCS separately; the sanitized arrays never retain target values.
                import copy
                frames[name]={'data':data,'error':error,'dq':dq,'field_mask':field,
                              'wcs':copy.deepcopy(wcs),'catalogue':catalogue,
                              'detector':str(model.meta.instrument.detector)}
                report['selection'][name]={**diagnostic,'protected_pixels':int(np.sum(~field)),
                                           'catalogue':catalogue}
        groups=assign_physical_stars({n:f['catalogue'] for n,f in frames.items()},cfg)
        report['physical_stars']=groups
        training=[g for g in groups if g['split']=='training']
        heldout=[g for g in groups if g['split']=='heldout']
        if len(training)<cfg['minimum_training_physical_stars'] or len(heldout)<cfg['minimum_heldout_physical_stars']:
            raise ValueError('Insufficient independent physical field stars for frozen training/heldout split')
        per_exposure={n:[] for n in products}
        for group in training:
            for row in group['detections']:
                per_exposure[row['exposure']].append(_field_stamp(frames[row['exposure']],row,group['physical_star_id']))
        enough=all(len(v)>=cfg['minimum_training_stamps_per_exposure'] for v in per_exposure.values())
        grids={};report['build_diagnostics']={}
        if enough:
            report['psf_sharing_mode']='per_exposure'
            for name in products:
                grid,diag=build_empirical_psf(per_exposure[name],cfg)
                grids[name]=grid;report['build_diagnostics'][name]=diag
        else:
            dither=common_pool_dither_diagnostic(training,products,cfg['maximum_common_pool_dither_native_pixel'])
            report['common_pool_dither_diagnostic']=dither
            if not cfg['allow_common_pool_for_this_visit']:
                raise ValueError('Frozen visit forbids a common pooled PSF: metadata show approximately270-pixel dithers; three training field stars per exposure are required')
            if len({f['detector'] for f in frames.values()})!=1 or not dither['small_dither_condition_pass']:
                raise ValueError('Too few per-exposure training stars and common-PSF small-dither condition fails')
            report['psf_sharing_mode']='common_same_detector_small_dither'
            grid,diag=build_empirical_psf([s for stamps in per_exposure.values() for s in stamps],cfg)
            for name in products:grids[name]=grid;report['build_diagnostics'][name]=diag
        for name,grid in grids.items():
            psf_path=out/(name.replace('_cal.fits','_field_epsf.fits'))
            hdu=fits.PrimaryHDU(grid)
            hdu.header['OVERSMP']=cfg['oversampling'];hdu.header['ORIGINX']=(grid.shape[1]-1)/2
            hdu.header['ORIGINY']=(grid.shape[0]-1)/2;hdu.header['PIXINT']=True
            hdu.header['HISTORY']='Field-only empirical pixel-integrated PSF; target circle excluded.'
            hdu.writeto(psf_path,overwrite=False)
            report['psf_files_by_exposure'][name]=psf_path.resolve().relative_to(ROOT.resolve()).as_posix()
            report['file_sha256'][report['psf_files_by_exposure'][name]]=sha(psf_path)
        report['heldout_fits']=[];sky_consistency=[]
        for group in heldout:
            positions=[];covariances=[];exposure_names=[]
            for row in group['detections']:
                name=row['exposure'];frame=frames[name]
                result=fit_heldout_field_star(frame['data'],frame['error'],frame['dq'],row['xy'],grids[name],cfg)
                position,cov,jac=gwcs_centroid(frame['wcs'],result['centroid_xy_pixel'],result['centroid_covariance_pixel2'],cfg['target_center_icrs_deg'])
                positions.append(position);covariances.append(cov);exposure_names.append(name)
                report['heldout_fits'].append({'physical_star_id':group['physical_star_id'],
                    'exposure':name,**result,'sky_tangent_arcsec':position.tolist(),
                    'sky_covariance_arcsec2':cov.tolist()})
            if len(positions)>1:
                sky_consistency.append({'physical_star_id':group['physical_star_id'],
                    'exposures':exposure_names,**combine_exposures(positions,covariances),
                    'scope':'Diagnostic absolute-GWCS consistency only; absolute pointing uncertainty is absent from ERR and is removed from the separate intrinsic lens-geometry statistic. This diagnostic does not rescale errors or tune the PSF.'})
        report['absolute_gwcs_exposure_consistency']=sky_consistency
        hold_exposures={row['exposure'] for row in report['heldout_fits']}
        report['all_exposures_have_heldout_field_star']=hold_exposures==set(products)
        report['overall_adequacy_pass']=bool(hold_exposures==set(products)
            and all(v['converged'] for v in report['build_diagnostics'].values())
            and all(row['absolute_adequacy_pass'] for row in report['heldout_fits']))
        report['status']=('FIELD_PSF_ADEQUATE_CONDITIONALLY_TARGETS_REMAIN_SEALED' if report['overall_adequacy_pass']
            else 'FIELD_PSF_INADEQUATE_TARGET_MEASUREMENT_BLOCKED')
    except Exception as error:
        report['status']='FIELD_CALIBRATION_FAILED_TARGET_MEASUREMENT_BLOCKED'
        report['failure']={'type':type(error).__name__,'message':str(error)}
    report['completed_utc']=datetime.now(timezone.utc).isoformat()
    dump(out/'field_psf_calibration.json',report)
    return report


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--protocol',type=Path,required=True)
    parser.add_argument('--expected-protocol-sha256',required=True)
    parser.add_argument('--data-directory',type=Path,default=ROOT/'data/strong_lensing_geometry_20261021/native')
    parser.add_argument('--output-directory',type=Path,default=ROOT/'artifacts/strong_lensing_geometry_20261021/field_psf')
    parser.add_argument('--execute-field-calibration',action='store_true')
    args=parser.parse_args()
    report=calibrate(args.protocol,args.expected_protocol_sha256,args.data_directory,
                     args.output_directory,args.execute_field_calibration)
    print(json.dumps({'status':report['status'],'overall_adequacy_pass':report['overall_adequacy_pass'],
                      'failure':report.get('failure')},indent=2))
    if not report['overall_adequacy_pass']:raise SystemExit(2)


if __name__=='__main__':main()
