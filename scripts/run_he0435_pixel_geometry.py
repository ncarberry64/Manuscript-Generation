"""One-time post-freeze native JWST geometry application; never downloads targets.

Only --execute opens numerical image planes. Protocol and predictions must
already be identical to their committed bytes at --freeze-commit. Every helper
listed in the frozen protocol is hash checked before the exposure marker.
"""
from __future__ import annotations

import argparse
import hashlib
import json
import subprocess
import sys
from datetime import datetime,timezone
from pathlib import Path

import numpy as np

ROOT=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT/'code'))
from strong_lensing_pixel_geometry import (fit_moffat_point,fit_empirical_point,gwcs_centroid,
    intrinsic_geometry,intrinsic_with_covariance,combine_exposures,gaussian_score)


def sha(path):return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def dump(path,value):
    Path(path).parent.mkdir(parents=True,exist_ok=True)
    Path(path).write_text(json.dumps(value,indent=2,allow_nan=False)+'\n',encoding='utf-8')


def verify_frozen(path,commit):
    path=Path(path).resolve()
    relative=path.relative_to(ROOT).as_posix()
    committed=subprocess.check_output(['git','show',f'{commit}:{relative}'],cwd=ROOT)
    if committed!=path.read_bytes():raise RuntimeError(f'Uncommitted frozen input: {relative}')


def reference_sky(protocol):
    if 'reference_image_radec_deg' in protocol:
        r=np.asarray(protocol['reference_image_radec_deg'],float)
    else:
        anchor=np.asarray(protocol['anchor_radec_deg'],float)
        offsets=np.asarray(protocol['association_reference_relative_arcsec'],float).reshape(4,2)
        r=anchor+offsets/np.array([3600*np.cos(np.deg2rad(anchor[1])),3600])
    if r.shape!=(4,2) or not np.all(np.isfinite(r)):raise ValueError('Frozen A/B/C/D sky guesses required')
    return r


def reference_pixel_guesses(protocol,native_filename,wcs,guesses_sky):
    """Frozen field-registration association only; it never modifies geometry."""
    mapping=protocol.get('reference_guesses_xy_pixel_by_exposure')
    if mapping is not None:
        guesses=np.asarray(mapping[native_filename],float)
    else:
        guesses=np.array([wcs.invert(float(a),float(d),with_bounding_box=False)[:2]
                          for a,d in guesses_sky],float)
    if guesses.shape!=(4,2) or not np.all(np.isfinite(guesses)):
        raise ValueError('Four finite frozen initial A/B/C/D pixel associations required')
    return guesses


def score_models(observed,covariance,models,alpha,required_reduction=.25):
    scores={}
    for name,model in models.items():
        relative=np.asarray(model['relative_B_C_D_to_A_east_north_arcsec'],float).reshape(3,2)
        p=np.vstack([np.zeros((1,2)),relative])
        q=intrinsic_geometry(p)
        result=gaussian_score(np.asarray(observed)-q,covariance)
        result['predicted_intrinsic_q_arcsec']=q.tolist()
        result['intrinsic_residual_mas']=(np.asarray(result['residual'])*1000).tolist()
        result['intrinsic_coordinate_rms_mas']=float(np.sqrt(np.mean(np.asarray(result['residual'])**2))*1000)
        result['formal_geometry_compatible_at_alpha']=bool(result['p_value']>=alpha)
        result['baseline_model']=model.get('baseline_model')
        scores[name]=result
    for name,result in scores.items():
        baseline=result['baseline_model']
        result['improves_declared_baseline_raw_q_RMS']=(result['intrinsic_coordinate_rms_mas']<scores[baseline]['intrinsic_coordinate_rms_mas']
                                                     if baseline in scores else None)
        old_rms=scores[baseline]['intrinsic_coordinate_rms_mas'] if baseline in scores else None
        fraction=(1-result['intrinsic_coordinate_rms_mas']/old_rms) if old_rms is not None and old_rms>0 else None
        result['raw_q_RMS_reduction_fraction']=fraction
        result['required_raw_q_RMS_reduction_fraction']=required_reduction
        result['meets_declared_baseline_raw_q_RMS_reduction']=bool(fraction is not None and fraction>=required_reduction)
    return scores


def frozen_psf_calibration(protocol,commit):
    if protocol.get('psf_model','moffat')=='moffat':return None
    if protocol.get('psf_model')!='field_empirical':raise ValueError('Unknown frozen PSF model')
    receipt=protocol['psf_calibration'];path=ROOT/receipt['path']
    if sha(path)!=receipt['sha256']:raise RuntimeError('Field PSF calibration receipt hash changed')
    verify_frozen(path,commit)
    calibration=json.loads(path.read_text(encoding='utf-8'))
    if calibration.get('overall_adequacy_pass') is not True:
        raise RuntimeError('Field-only PSF calibration failed; quasar targets must remain sealed')
    checked=set()
    for native in calibration['psf_files_by_exposure']:
        for relative in psf_paths_by_label(calibration,native).values():
            if relative in checked:continue
            psf=ROOT/relative
            expected=calibration['file_sha256'].get(relative)
            if expected!=sha(psf):raise RuntimeError(f'Field PSF grid hash changed: {relative}')
            verify_frozen(psf,commit)
            checked.add(relative)
    return calibration


def psf_paths_by_label(calibration,native_filename):
    """Resolve flat legacy or position-dependent A/B/C/D field-frozen PSFs."""
    entry=calibration['psf_files_by_exposure'][native_filename]
    if isinstance(entry,str):return {label:entry for label in 'ABCD'}
    if not isinstance(entry,dict) or set(entry)!=set('ABCD'):
        raise ValueError('Position-dependent PSF mapping requires exactly A/B/C/D')
    if not all(isinstance(value,str) and value for value in entry.values()):
        raise ValueError('Each position-dependent PSF must name a frozen grid file')
    return dict(entry)


def verify_native_inputs(paths,protocol,calibration):
    """Check exact declared files against frozen field receipt before exposure."""
    declared=protocol['primary_observation']['products']
    basenames=[Path(path).name for path in paths]
    if len(declared)!=4 or len(set(declared))!=4 or len(basenames)!=4 or set(basenames)!=set(declared):
        raise RuntimeError('Exactly the four declared native CAL files are required')
    expected=(calibration or {}).get('input_sha256',protocol.get('native_sha256',{}))
    for path in paths:
        if expected.get(Path(path).name)!=sha(path):
            raise RuntimeError(f'Native CAL hash missing or changed before target opening: {Path(path).name}')


def archive_native_windows(model,guesses,path,out,source_sha):
    """Exact native copies, called only inside the explicitly exposed application."""
    from astropy.io import fits
    paths={}
    for label,guess in zip('ABCD',guesses):
        ix,iy=np.rint(guess).astype(int)
        x0,x1=max(0,ix-14),min(model.data.shape[1],ix+15)
        y0,y1=max(0,iy-14),min(model.data.shape[0],iy+15)
        primary=fits.PrimaryHDU()
        for key,value in [('X0',x0),('Y0',y0),('LABEL',label),('CALFILE',path.name),('CAL_SHA',source_sha)]:
            primary.header[key]=value
        primary.header['HISTORY']='Exact SCI/ERR/DQ native copies, no resampling or error modification.'
        hdus=[primary]+[fits.ImageHDU(np.array(array[y0:y1,x0:x1]),name=name)
                        for name,array in [('SCI',model.data),('ERR',model.err),('DQ',model.dq)]]
        output=out/'native_stamps'/f'{path.stem}_{label}.fits';output.parent.mkdir(parents=True,exist_ok=True)
        fits.HDUList(hdus).writeto(output,overwrite=False)
        paths[label]={'path':f'native_stamps/{output.name}','sha256':sha(output),
                      'native_bounds_xyxy':[int(v) for v in [x0,y0,x1,y1]]}
    return paths


def apply_native_files(paths,protocol,models,out,calibration=None):
    # Lazy import prevents even the metadata command from requiring JWST runtime.
    from stdatamodels.jwst import datamodels
    guesses_sky=reference_sky(protocol)
    origin=protocol.get('tangent_origin_radec_deg',guesses_sky[0].tolist())
    primary=protocol.get('primary_filter','F115W')
    expected=int(protocol.get('expected_target_exposures',4))
    planned_fits=int(protocol.get('planned_point_fits',expected*4))
    family_alpha=float(protocol.get('pixel_familywise_alpha',.016))
    pixel_threshold=family_alpha/planned_fits
    consistency_alpha=float(protocol.get('exposure_consistency_alpha',.001))
    geometry_alpha=float(protocol.get('geometry_alpha',.05))
    required_reduction=float(protocol.get('required_raw_q_rms_reduction_fraction',.25))
    radius=int(protocol.get('stamp_radius_pixels',8))
    bound=float(protocol.get('centroid_bound_pixels',5.))
    quadrature=int(protocol.get('quadrature_order',7))
    psf_model=protocol.get('psf_model','moffat')
    background_degree=int(protocol.get('background_degree',1))
    exposures=[];non_target=[]
    for path in sorted(set(Path(p).resolve() for p in paths)):
        with datamodels.open(path) as model:
            filter_name=str(model.meta.instrument.filter)
            if filter_name!=primary:
                raise ValueError(f'Unexpected filter {filter_name}; unopened reserve bands prohibited')
            wcs=model.meta.wcs
            try:
                guesses=reference_pixel_guesses(protocol,path.name,wcs,guesses_sky)
            except Exception as exc:
                non_target.append({'file':str(path),'reason':'inverse WCS failure','message':str(exc)});continue
            shape=tuple(model.data.shape)
            if len(shape)!=2:raise ValueError('Native SCI must be two-dimensional')
            margin=radius+bound
            inside=(np.all(np.isfinite(guesses)) and np.all(guesses[:,0]>=margin)
                    and np.all(guesses[:,1]>=margin) and np.all(guesses[:,0]<shape[1]-margin)
                    and np.all(guesses[:,1]<shape[0]-margin))
            if not inside:
                non_target.append({'file':str(path),'reason':'All four frozen target windows not in footprint',
                                   'shape':list(shape)});continue
            item={'file':str(path),'sha256':sha(path),'filter':filter_name,
                  'detector':str(model.meta.instrument.detector),
                  'measurement_errors':[], 'points':[]}
            item['native_stamp_archives']=archive_native_windows(model,guesses,path,out,item['sha256'])
            try:
                import asdf,copy
                wcs_path=out/'wcs'/f'{path.stem}_wcs.asdf';wcs_path.parent.mkdir(parents=True,exist_ok=True)
                asdf.AsdfFile({'wcs':copy.deepcopy(wcs),'native_CAL_SHA256':item['sha256']}).write_to(wcs_path)
                item['official_GWCS_archive']={'path':f'wcs/{wcs_path.name}','sha256':sha(wcs_path)}
            except Exception as error:
                item['official_GWCS_archive_error']={'type':type(error).__name__,'message':str(error),
                    'replay':'Recover official GWCS from public native CAL URL and retained SHA256'}
            for key in ['filename','telescope','date','observation','exposure','wcsinfo','ref_file','cal_step']:
                try:item[key]=model.meta[key].instance
                except Exception:pass
            points=[];covariance=np.zeros((8,8))
            psfs_by_label={}
            if psf_model=='field_empirical':
                from astropy.io import fits
                paths_by_label=psf_paths_by_label(calibration,path.name)
                loaded={}
                for label,relative in paths_by_label.items():
                    if relative not in loaded:
                        with fits.open(ROOT/relative) as psf:
                            grid=np.array(psf[0].data,float)
                            over=int(psf[0].header['OVERSMP'])
                            psf_origin=[float(psf[0].header['ORIGINX']),float(psf[0].header['ORIGINY'])]
                        loaded[relative]=(grid,over,psf_origin)
                    psfs_by_label[label]=loaded[relative]
                item['field_psf_paths_by_label']=paths_by_label
                item['field_psf_sha256_by_label']={label:sha(ROOT/relative) for label,relative in paths_by_label.items()}
                if len(set(paths_by_label.values()))==1:
                    relative=next(iter(paths_by_label.values()))
                    item['field_psf_path']=relative;item['field_psf_sha256']=sha(ROOT/relative)
            for i,(label,guess) in enumerate(zip('ABCD',guesses)):
                try:
                    if psf_model=='field_empirical':
                        psf_grid,oversampling,psf_origin=psfs_by_label[label]
                        fit=fit_empirical_point(model.data,model.err,model.dq,guess,psf_grid,
                                                oversampling,psf_origin,radius,bound,background_degree)
                        fit['field_psf_path']=item['field_psf_paths_by_label'][label]
                        fit['field_psf_sha256']=item['field_psf_sha256_by_label'][label]
                    else:
                        fit=fit_moffat_point(model.data,model.err,model.dq,guess,radius,bound,quadrature)
                    skyxy,skycov,jac=gwcs_centroid(wcs,fit['centroid_xy_pixel'],fit['centroid_covariance_pixel2'],origin)
                    fit['label']=label;fit['initial_guess_xy_pixel']=guess.tolist()
                    fit['tangent_east_north_arcsec']=skyxy.tolist()
                    fit['tangent_covariance_arcsec2']=skycov.tolist()
                    fit['GWCS_jacobian_arcsec_per_pixel']=jac.tolist()
                    fit['pixel_model_adequate_at_frozen_threshold']=bool(fit['valid_measurement_fit'] and fit['p_value']>=pixel_threshold)
                    item['points'].append(fit);points.append(skyxy)
                    covariance[2*i:2*i+2,2*i:2*i+2]=skycov
                except Exception as exc:
                    item['measurement_errors'].append({'label':label,'error':type(exc).__name__,'message':str(exc)})
            if len(points)==4:
                points=np.asarray(points)
                q,qcov=intrinsic_with_covariance(points,covariance)
                item['points_tangent_east_north_arcsec']=points.tolist()
                item['full_centroid_covariance_arcsec2']=covariance.tolist()
                item['intrinsic_q_arcsec']=q.tolist();item['intrinsic_covariance_arcsec2']=qcov.tolist()
                item['relative_B_C_D_to_A_east_north_arcsec']=(points[1:]-points[0]).tolist()
                item['model_scores']=score_models(q,qcov,models,geometry_alpha,required_reduction)
            dump(out/'exposures'/f'{path.stem}.json',item)
            exposures.append(item)
    complete=[x for x in exposures if 'intrinsic_q_arcsec' in x]
    count_ok=len(complete)==expected and len(exposures)==expected
    fit_ok=count_ok and all(p['pixel_model_adequate_at_frozen_threshold'] for x in exposures for p in x['points'])
    result={'measurement_family':('Field-frozen pixel-integrated empirical ePSF plus '+
                                  ('quadratic' if background_degree==2 else 'affine')+' background, fixed known ERR'
                                   if psf_model=='field_empirical' else 'Pixel-integrated elliptical Moffat plus affine background, fixed known ERR'),
            'target_filter':primary,'exposures_found':len(exposures),'complete_exposures':len(complete),
            'expected_target_exposures':expected,'exposure_count_matches_protocol':count_ok,
            'pixel_fit_familywise_alpha':family_alpha,'pixel_fit_bonferroni_threshold':pixel_threshold,
            'all_pixel_models_adequate':fit_ok,'non_target_files':non_target,
            'exposure_result_files':[f'exposures/{Path(x["file"]).stem}.json' for x in exposures],
            'conditional_covariance_scope':'Known ERR, disjoint point stamps, independent exposures; instrument distortion/reference common covariance and instrumental PSF uncertainty are not supplied by ERR. No arbitrary floor or residual scale inferred.',
            'physical_admission':False,'unchanged_physical_guard':2e-10,
            'exposure_consistency_alpha':consistency_alpha,'geometry_alpha':geometry_alpha}
    if complete:
        combined=combine_exposures([x['intrinsic_q_arcsec'] for x in complete],
                                   [x['intrinsic_covariance_arcsec2'] for x in complete])
        consistency=combined['exposure_consistency_p_value']
        consistent=bool(count_ok and consistency is not None and consistency>=consistency_alpha)
        result['combined']=combined;result['exposures_consistent']=consistent
        result['models']=score_models(combined['mean'],combined['covariance'],models,geometry_alpha,required_reduction)
        for name,s in result['models'].items():
            s['conditional_geometry_validation_pass']=bool(fit_ok and consistent and s['formal_geometry_compatible_at_alpha']
                    and s['meets_declared_baseline_raw_q_RMS_reduction'])
        result['measurement_status']='ADEQUATE_UNDER_CONDITIONAL_KNOWN_ERR_MODEL' if fit_ok and consistent else 'FAILED_MEASUREMENT_ADEQUACY_OR_EXPOSURE_CONSISTENCY'
    else:
        result['measurement_status']='NO_COMPLETE_FOUR_POINT_MEASUREMENT'
        result['exposures_consistent']=False;result['models']={}
    return result


def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--protocol',type=Path,required=True)
    p.add_argument('--predictions',type=Path,required=True)
    p.add_argument('--freeze-commit',required=True)
    p.add_argument('--cal',type=Path,nargs='*',default=[])
    p.add_argument('--output',type=Path,required=True)
    p.add_argument('--execute',action='store_true')
    args=p.parse_args()
    verify_frozen(args.protocol,args.freeze_commit);verify_frozen(args.predictions,args.freeze_commit)
    protocol=json.loads(args.protocol.read_text(encoding='utf-8'))
    prediction=json.loads(args.predictions.read_text(encoding='utf-8'))
    required=['code/strong_lensing_pixel_geometry.py','scripts/run_he0435_pixel_geometry.py',
              'tests/test_strong_lensing_pixel_geometry.py']
    recorded=protocol.get('source_sha256',protocol.get('source_files_sha256',{}))
    for relative in required:
        if recorded.get(relative)!=sha(ROOT/relative):raise RuntimeError(f'Frozen helper hash missing/mismatch: {relative}')
        verify_frozen(ROOT/relative,args.freeze_commit)
    if protocol.get('prediction_sha256')!=sha(args.predictions):raise RuntimeError('Prediction hash missing/mismatch')
    calibration=frozen_psf_calibration(protocol,args.freeze_commit)
    if not args.execute:
        print(json.dumps({'frozen_inputs_verified':True,'target_arrays_opened':False,
                          'note':'--execute required; no SCI planes opened'},indent=2));return
    if not args.cal:raise ValueError('Explicit native CAL filenames required')
    verify_native_inputs(args.cal,protocol,calibration)
    args.output.mkdir(parents=True,exist_ok=True)
    marker=args.output/'target_exposure_started.json'
    with marker.open('x',encoding='utf-8') as stream:
        json.dump({'started_utc':datetime.now(timezone.utc).isoformat(),'freeze_commit':args.freeze_commit,
                   'protocol_sha256':sha(args.protocol),'prediction_sha256':sha(args.predictions),
                   'native_filenames':[str(f) for f in args.cal],
                   'note':'Created exclusively before numerical native-array opening; targets consumed even if subsequent measurement fails'},stream,indent=2)
        stream.write('\n')
    result=apply_native_files(args.cal,protocol,prediction['models'],args.output,calibration)
    result['freeze_commit']=args.freeze_commit
    dump(args.output/'pixel_geometry_results.json',result)
    print(json.dumps({'result':str(args.output/'pixel_geometry_results.json'),
                      'measurement_status':result['measurement_status'],
                      'models':{k:{kk:vv for kk,vv in v.items() if kk in ['chi2','p_value','intrinsic_coordinate_rms_mas','conditional_geometry_validation_pass']}
                                for k,v in result['models'].items()}},indent=2))


if __name__=='__main__':main()
