"""Prespecified external-LMC PSF assessment on guarded Gaia FIELD-star wings."""
from __future__ import annotations

import argparse,csv,hashlib,json,ssl,subprocess,sys
from datetime import datetime,timezone
from pathlib import Path
from urllib.request import Request,urlopen
import certifi
import numpy as np

ROOT=Path(__file__).resolve().parents[1];sys.path.insert(0,str(ROOT/'code'))
from strong_lensing_standard_psf import (stellar_identity,propagate_gaia_field_position,
    associate_saturated_field_core,sample_spatial_standard_psf,registered_reference_guesses,apply_minimum_pixel_gate)
from strong_lensing_field_psf import DEFAULT_ALGORITHM,field_mask_from_gwcs,guarded_field_arrays,fit_heldout_field_star

SOURCES=['scripts/calibrate_he0435_standard_psf.py','code/strong_lensing_standard_psf.py',
         'code/strong_lensing_field_psf.py','code/strong_lensing_pixel_geometry.py']


def sha(path):return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def dump(path,value):
    def default(item):
        if isinstance(item,np.generic):return item.item()
        if isinstance(item,np.ndarray):return item.tolist()
        raise TypeError(type(item).__name__)
    Path(path).write_text(json.dumps(value,indent=2,allow_nan=False,default=default)+'\n',encoding='utf-8')


def frozen_inputs(protocol_path,expected_sha,authorized):
    if not authorized:raise PermissionError('External PSF execution is field-only and must be explicitly authorized')
    p=Path(protocol_path)
    if sha(p)!=expected_sha:raise ValueError('Frozen calibration protocol SHA mismatch')
    cfg=json.loads(p.read_text(encoding='utf-8'))
    paths={**cfg['input_sha256'],**cfg['source_sha256']}
    paths[p.resolve().relative_to(ROOT).as_posix()]=expected_sha
    for rel,expected in paths.items():
        path=ROOT/rel
        if sha(path)!=expected:raise ValueError('Frozen input hash mismatch: '+rel)
        if subprocess.check_output(['git','show','HEAD:'+rel],cwd=ROOT)!=path.read_bytes():
            raise ValueError('Calibration input is not committed: '+rel)
    for rel in SOURCES:
        if rel not in cfg['source_sha256']:raise ValueError('Missing imported source binding: '+rel)
    return cfg


def execute(protocol_path,expected_sha,authorized=False):
    cfg=frozen_inputs(protocol_path,expected_sha,authorized)
    from stdatamodels.jwst import datamodels
    from astropy.io import fits
    from astropy.time import Time
    from photutils.psf import GriddedPSFModel
    out=ROOT/cfg['output_directory'];out.mkdir(parents=True,exist_ok=True)
    lock=out/'calibration_started.json'
    with lock.open('x',encoding='utf-8') as stream:
        json.dump({'protocol_sha256':expected_sha,'started_utc':datetime.now(timezone.utc).isoformat(),
                   'field_only_authorized':True,'target_measurement_authorized':False},stream,indent=2)
    report={'status':'INDEPENDENT_STANDARD_PSF_FIELD_CALIBRATION_STARTED','overall_adequacy_pass':False,
            'protocol_sha256':expected_sha,'physical_prediction_sha256':cfg['physical_prediction_sha256'],
            'target_pixel_values_inspected':False,'PSF_trained_on_HE0435':False,'source_sha256':cfg['source_sha256'],
            'input_sha256':{},
            'psf_files_by_exposure':{},'file_sha256':{},'field_fits':[],
            'field_reference_guesses_xy_pixel_by_exposure':{},'field_registration':{},
            'conditional_uncertainty_note':cfg['uncertainty_scope']}
    try:
        rows=list(csv.DictReader((ROOT/cfg['gaia_field_catalogue']).open(encoding='utf-8')))
        selected=[r for r in rows if r['source_id']==cfg['selected_field_source_id']]
        if len(selected)!=1:raise ValueError('Exactly one prespecified independent Gaia FIELD identity required')
        row=selected[0];identity=stellar_identity(row,cfg['minimum_star_probability'],cfg['minimum_astrometric_significance'],cfg['maximum_ruwe'])
        report['selected_gaia_field_identity']={**row,'assessment':identity}
        if not identity['accepted']:raise ValueError('Prespecified catalogue source lacks independent stellar evidence')
        library=ROOT/cfg['library_local_path'];library.parent.mkdir(parents=True,exist_ok=True)
        if library.exists():raise FileExistsError('External library array must not have been opened before this new freeze')
        with urlopen(Request(cfg['library_url'],headers={'User-Agent':'Independent scientific field PSF calibration'}),
                     timeout=60,context=ssl.create_default_context(cafile=certifi.where())) as response:
            raw=response.read()
        library.write_bytes(raw)
        if len(raw)!=cfg['library_expected_bytes']:raise ValueError('Unexpected independent PSF-library byte count')
        if hashlib.sha256(raw[:2880]).hexdigest()!=cfg['library_primary_header_sha256']:
            raise ValueError('Independent library header differs from the prespecified metadata')
        report['library']={'path':cfg['library_local_path'],'url':cfg['library_url'],'bytes':len(raw),
                           'sha256':hashlib.sha256(raw).hexdigest(),'training':'LMC COM-1074/CAL-1476, independent HE0435 photons'}
        standard=GriddedPSFModel.read(library,format='stdpsf')
        root_protocol=json.loads((ROOT/cfg['reference_geometry_protocol']).read_text(encoding='utf-8'))
        anchor=np.array(root_protocol['anchor_radec_deg']);offsets=np.array(root_protocol['association_reference_relative_arcsec'])
        reference_sky=anchor+offsets/np.array([3600*np.cos(np.deg2rad(anchor[1])),3600])
        receipt=json.loads((ROOT/cfg['native_acquisition_receipt']).read_text(encoding='utf-8'))
        by_name={r['filename']:r for r in receipt['products']}
        for name in cfg['native_products']:
            path=ROOT/cfg['native_data_directory']/name
            actual_sha=sha(path)
            if actual_sha!=by_name[name]['sha256']:raise ValueError('Native acquisition hash mismatch')
            report['input_sha256'][name]=actual_sha
            with datamodels.open(path) as native:
                wcs=native.meta.wcs
                field=field_mask_from_gwcs(native.data.shape,wcs,cfg['field_center_icrs_deg'],cfg['target_guard_radius_arcsec'])
                data,error,dq=guarded_field_arrays(native.data,native.err,native.dq,field)
                epoch=Time(str(native.meta.observation.date)+'T'+str(native.meta.observation.time)).jyear
                sky=propagate_gaia_field_position(row,epoch)
                nominal=np.array(wcs.invert(*sky,with_bounding_box=False)[:2],float)
                association=associate_saturated_field_core(dq,field,nominal,cfg['association_search_radius_pixel'],cfg['minimum_saturated_component_pixels'])
                initial=np.array(association['center_xy_pixel'])
                ix,iy=np.rint(initial).astype(int);r=cfg['star_fit_stamp_radius_pixel']
                if not field[iy-r:iy+r+1,ix-r:ix+r+1].all():raise ValueError('Field-star fitting stamp intersects target guard')
                star_grid=sample_spatial_standard_psf(standard,initial,cfg['oversampling'],tuple(cfg['grid_shape']))
                settings=dict(DEFAULT_ALGORITHM);settings['fit_stamp_radius_native_pixel']=r
                settings['heldout_centroid_bound_native_pixel']=cfg['star_centroid_bound_pixel']
                settings['heldout_minimum_absolute_fit_p_value']=cfg['minimum_field_fit_p_value']
                result=fit_heldout_field_star(data,error,dq,initial,star_grid,settings)
                result=apply_minimum_pixel_gate(result,cfg['minimum_valid_native_pixels'])
                guesses,translation=registered_reference_guesses(wcs,reference_sky,nominal,result['centroid_xy_pixel'])
                report['field_fits'].append({'exposure':name,'epoch_jyear':epoch,'nominal_catalogue_xy':nominal.tolist(),
                    'saturated_component_association':association,**result,
                    'saturated_pixels_in_star_stamp':int(np.count_nonzero(dq[iy-r:iy+r+1,ix-r:ix+r+1]&2))})
                report['field_reference_guesses_xy_pixel_by_exposure'][name]=guesses.tolist()
                report['field_registration'][name]={'translation_pixel':translation.tolist(),
                    'scope':'Changes association initial guesses only; no target affine/scale/shear or physical lens parameter fitted.'}
                paths={}
                for label,position in [('FIELD',initial)]+list(zip('ABCD',guesses)):
                    grid=star_grid if label=='FIELD' else sample_spatial_standard_psf(standard,position,cfg['oversampling'],tuple(cfg['grid_shape']))
                    target=out/name.replace('_cal.fits','_'+label+'_standard_epsf.fits')
                    hdu=fits.PrimaryHDU(grid);hdu.header['OVERSMP']=cfg['oversampling'];hdu.header['ORIGINX']=(grid.shape[1]-1)/2
                    hdu.header['ORIGINY']=(grid.shape[0]-1)/2;hdu.header['PIXINT']=True;hdu.header['SOURCE']='LMC STDPSF'
                    hdu.writeto(target,overwrite=False);rel=target.resolve().relative_to(ROOT).as_posix()
                    report['file_sha256'][rel]=sha(target)
                    if label!='FIELD':paths[label]=rel
                report['psf_files_by_exposure'][name]=paths
        report['overall_adequacy_pass']=bool(len(report['field_fits'])==4 and all(r['absolute_adequacy_pass'] for r in report['field_fits']))
        report['status']='EXTERNAL_STANDARD_PSF_FIELD_ADEQUATE_CONDITIONALLY' if report['overall_adequacy_pass'] else 'EXTERNAL_STANDARD_PSF_FIELD_INADEQUATE_TARGET_MEASUREMENT_BLOCKED'
    except Exception as error:
        report['status']='EXTERNAL_STANDARD_PSF_CALIBRATION_FAILED_TARGET_MEASUREMENT_BLOCKED'
        report['failure']={'type':type(error).__name__,'message':str(error)}
    report['completed_utc']=datetime.now(timezone.utc).isoformat();dump(out/'field_psf_calibration.json',report)
    return report


def main():
    parser=argparse.ArgumentParser(description=__doc__);parser.add_argument('--protocol',type=Path,required=True)
    parser.add_argument('--expected-protocol-sha256',required=True);parser.add_argument('--execute-field-only',action='store_true')
    args=parser.parse_args();report=execute(args.protocol,args.expected_protocol_sha256,args.execute_field_only)
    print(json.dumps({'status':report['status'],'overall_adequacy_pass':report['overall_adequacy_pass'],
                     'field_p_values':[r['p_value'] for r in report['field_fits']],'failure':report.get('failure')},indent=2))
    if not report['overall_adequacy_pass']:raise SystemExit(2)


if __name__=='__main__':main()
