"""Transfer the already-learned positive field response unchanged to targets.

All F115W target data are consumed development data. Field calibration still
fails, so every physical-validation flag remains false. The shape and all
physical predictions are fixed; no lens residual or score trains the response.
"""
from pathlib import Path
import json,sys
import numpy as np
from astropy.io import fits
import asdf,gwcs
from scipy.interpolate import RectBivariateSpline
HERE=Path(__file__).resolve().parent;REPO=HERE.parent/'repo';ART=REPO/'artifacts/strong_lensing_geometry_20261021'
sys.path.insert(0,str(HERE));from positive_field_response import get_fields,Basis,fit,serial,read_primary_fits
sys.path.insert(0,str(REPO/'code'));from strong_lensing_pixel_geometry import gwcs_centroid,intrinsic_with_covariance,combine_exposures,gaussian_score

def main():
    field=get_fields();basis=Basis(field[0]['grid']);saved=np.load(HERE/'positive_field_shape.npz');a=saved['shape']
    assert np.max(abs(saved['null']-basis.null))==0
    old=json.loads((ART/'jwst_target_diagnostic/pixel_geometry_results.json').read_text())
    origin=None;exposures=[]
    for path in sorted((ART/'jwst_target_diagnostic/exposures').glob('*.json')):
        ex=json.loads(path.read_text());name=path.stem+'.fits';records=[]
        with fits.open(REPO/'data/strong_lensing_geometry_20261021/native'/name,memmap=False) as h:
            for p in ex['points']:
                x0,y0,x1,y1=p['stamp_bounds_xyxy'];center=np.asarray(p['centroid_xy_pixel'])
                d,e,q=[np.asarray(h[k].data[y0:y1,x0:x1]) for k in ['SCI','ERR','DQ']]
                yy,xx=np.mgrid[y0:y1,x0:x1];good=np.isfinite(d)&np.isfinite(e)&(e>0)&((q.astype(np.uint64)&1)==0)
                grid,_=read_primary_fits(REPO/p['field_psf_path']);pedestal=max(0.,-float(grid.min()))*1.1+1e-8
                records.append({'name':p['label'],'x':xx[good]-center[0],'y':yy[good]-center[1],
                    'data':d[good].astype(float),'error':e[good].astype(float),'grid':grid+pedestal,
                    'signed_grid':grid,'pedestal':pedestal,'origin':center,'sat_restored':int(np.sum(good&((q&2)!=0))),
                    'old_flux':p['flux_in_SCI_pixel_units'],'background_degree':2,
                    'interpolator':RectBivariateSpline(np.arange(101)/4-12.5,np.arange(101)/4-12.5,(grid+pedestal).T)})
        result=fit(records,basis,a_fixed=a);row={'name':name,'fit':serial(result),'geometry':{}}
        with asdf.open(ART/'jwst_target_diagnostic/wcs'/name.replace('_cal.fits','_cal_wcs.asdf'),lazy_load=False) as af:
            w=af['wcs']
            if origin is None:origin=np.array(w(*ex['points'][0]['centroid_xy_pixel'])[:2])
            for mode in ['raw','gauge']:
                positions=[];cov=np.zeros((8,8))
                for k,r in enumerate(result['fits']):
                    xy,c,j=gwcs_centroid(w,r[f'centroid_{mode}_pixel'],r['covariance_pixel2_conditional_on_shape'],origin)
                    positions.append(xy);cov[2*k:2*k+2,2*k:2*k+2]=c
                q,c=intrinsic_with_covariance(positions,cov)
                row['geometry'][mode]={'intrinsic_q_arcsec':q.tolist(),'covariance_arcsec2':c.tolist(),'positions_tangent_arcsec':np.asarray(positions).tolist()}
        exposures.append(row)
        print(json.dumps({'name':name,'pixel_chi2_per_dof':[r['chi2']/r['dof_conditional_on_shape'] for r in result['fits']]}),flush=True)
    report={'status':'INCOMPLETE_CALIBRATION_DEVELOPMENT_ONLY','physical_validation_pass':False,
        'response_training':'Same positive field response fixed in positive_field_shape.npz; targets affect only point flux,centroid,quadratic background.',
        'uncertainty':'Centroid inverse information conditional on learned response; field-shape calibration uncertainty and shared systematic covariance not included. No error rescaling.',
        'exposures':exposures,'combined':{}}
    for mode in ['raw','gauge']:
        combo=combine_exposures([e['geometry'][mode]['intrinsic_q_arcsec'] for e in exposures],[e['geometry'][mode]['covariance_arcsec2'] for e in exposures]);scores={}
        for name,m in old['models'].items():
            residual=np.array(combo['mean'])-m['predicted_intrinsic_q_arcsec'];s=gaussian_score(residual,combo['covariance'])
            s['rms_mas']=float(np.sqrt(np.mean(residual**2))*1000);s['physical_validation_pass']=False;scores[name]=s
        report['combined'][mode]={'measurement':combo,'scores':scores}
    (HERE/'positive_target_transfer.json').write_text(json.dumps(report,indent=2))
    print(json.dumps({k:{'scatter':v['measurement']['exposure_consistency_chi2'],'rms_G1_mas':v['scores']['revised_G1_conventional']['rms_mas'],'Q_G1':v['scores']['revised_G1_conventional']['chi2']} for k,v in report['combined'].items()}),flush=True)

if __name__=='__main__':main()
