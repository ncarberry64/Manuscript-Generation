"""Post-calibration geometry scores after the local family/results are fixed.

This script alone reads physical predictions. It cannot fit or select the PSF
family. Both all-image fitted response and cross-exposure fitted response are
development measurements with failed pixel response and conditional errors.
"""
from pathlib import Path
import json,sys,hashlib
import numpy as np
import asdf,gwcs
HERE=Path(__file__).resolve().parent;REPO=HERE.parent/'repo';ART=REPO/'artifacts/strong_lensing_geometry_20261021'
sys.path.insert(0,str(REPO/'code'))
from strong_lensing_pixel_geometry import gwcs_centroid,intrinsic_with_covariance,combine_exposures,gaussian_score

def main():
    path=HERE/'local_quasar_refined_results.json';raw=path.read_bytes();r=json.loads(raw)
    assert all(f['train']['success'] for f in r['folds']) and r['all_image_fit']['success']
    # Calibration is finished and immutable for this calculation before these
    # physical predictions are loaded. No numerical optimizer is imported.
    predictions=json.loads((ART/'jwst_target_diagnostic/pixel_geometry_results.json').read_text())['models']
    bymethod={'all_image_response':r['all_image_fit']['fits'],
              'held_exposure_response':[p for f in r['folds'] for p in f['test']['fits']]}
    origin=None;out={'status':'DEVELOPMENT_ONLY_FAILED_RESPONSE_CALIBRATION','physical_validation_pass':False,
       'calibration_result_sha256':hashlib.sha256(raw).hexdigest(),'scope':'No lens-fitting or mean-response tuning. Conditional centroid covariances omit learned-shape uncertainty.',
       'measurements':{}}
    for method,points in bymethod.items():
        exposures=sorted(set(p['name'].rsplit('_',1)[0] for p in points));exprows=[]
        for exp in exposures:
            pts=sorted([p for p in points if p['name'].rsplit('_',1)[0]==exp],key=lambda p:p['name'])
            row={'exposure':exp,'geometry':{}}
            with asdf.open(ART/'jwst_target_diagnostic/wcs'/f'{exp}_wcs.asdf',lazy_load=False) as af:
                w=af['wcs']
                if origin is None:origin=np.array(w(*pts[0]['centroid_raw_pixel'])[:2])
                for mode in ['raw','gauge']:
                    positions=[];cov=np.zeros((8,8))
                    for k,p in enumerate(pts):
                        xy,c,j=gwcs_centroid(w,p[f'centroid_{mode}_pixel'],p['covariance_pixel2_conditional_on_shape'],origin)
                        positions.append(xy);cov[2*k:2*k+2,2*k:2*k+2]=c
                    q,c=intrinsic_with_covariance(positions,cov)
                    row['geometry'][mode]={'q':q.tolist(),'covariance':c.tolist(),'tangent_positions':np.asarray(positions).tolist()}
            exprows.append(row)
        scores={}
        for mode in ['raw','gauge']:
            combined=combine_exposures([e['geometry'][mode]['q'] for e in exprows],[e['geometry'][mode]['covariance'] for e in exprows]);models={}
            for name,p in predictions.items():
                diff=np.array(combined['mean'])-p['predicted_intrinsic_q_arcsec'];s=gaussian_score(diff,combined['covariance'])
                s['rms_mas']=float(np.sqrt(np.mean(diff**2))*1000);s['physical_validation_pass']=False;models[name]=s
            scores[mode]={'combined':combined,'models':models}
        out['measurements'][method]={'exposures':exprows,'scores':scores}
    (HERE/'local_quasar_geometry_results.json').write_text(json.dumps(out,indent=2))
    print(json.dumps({m:{mode:{'G1_RMS_mas':s['models']['revised_G1_conventional']['rms_mas'],
       'G1_Q':s['models']['revised_G1_conventional']['chi2'],'scatter_Q':s['combined']['exposure_consistency_chi2']}
       for mode,s in row['scores'].items()} for m,row in out['measurements'].items()},indent=2))

if __name__=='__main__':main()
