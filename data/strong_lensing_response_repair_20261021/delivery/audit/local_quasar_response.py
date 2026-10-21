"""Local quasar-colour PSF development using only consumed F115W stamps.

Family fixed here: positive shared log-response with RBFspacing0.5,width0.35,
same coefficient priorSD0.2, free independent centroids/fluxes/quadratic
backgrounds. Four leave-one-entire-exposure-out checks train on12 images and
evaluate4 images. Physical lens positions/scores never enter training or family
selection. All data have already been used; these are development diagnostics.
"""
from pathlib import Path
import json,sys,hashlib
import numpy as np
from astropy.io import fits
from scipy.interpolate import RectBivariateSpline
from scipy.linalg import svdvals
HERE=Path(__file__).resolve().parent;REPO=HERE.parent/'repo';ART=REPO/'artifacts/strong_lensing_geometry_20261021'
sys.path.insert(0,str(HERE))
from positive_field_response import fit,serial,read_primary_fits
from fine_field_response import FineBasis

def get_targets():
    result=[]
    for path in sorted((ART/'jwst_target_diagnostic/exposures').glob('*.json')):
        ex=json.loads(path.read_text());name=path.stem+'.fits'
        with fits.open(REPO/'data/strong_lensing_geometry_20261021/native'/name,memmap=False) as h:
            for p in ex['points']:
                x0,y0,x1,y1=p['stamp_bounds_xyxy'];center=np.asarray(p['centroid_xy_pixel'])
                d,e,q=[np.asarray(h[k].data[y0:y1,x0:x1]) for k in ['SCI','ERR','DQ']]
                yy,xx=np.mgrid[y0:y1,x0:x1];good=np.isfinite(d)&np.isfinite(e)&(e>0)&((q.astype(np.uint64)&1)==0)
                grid,_=read_primary_fits(REPO/p['field_psf_path']);pedestal=max(0.,-float(grid.min()))*1.1+1e-8
                result.append({'name':path.stem+'_'+p['label'],'exposure':path.stem,'label':p['label'],
                    'x':xx[good]-center[0],'y':yy[good]-center[1],
                    'data':d[good].astype(float),'error':e[good].astype(float),'grid':grid+pedestal,
                    'signed_grid':grid,'pedestal':pedestal,'origin':center,'sat_restored':int(np.sum(good&((q&2)!=0))),
                    'pixel_x':xx[good],'pixel_y':yy[good],'pixel_dq':q[good],
                    'old_flux':p['flux_in_SCI_pixel_units'],'background_degree':2,
                    'interpolator':RectBivariateSpline(np.arange(101)/4-12.5,np.arange(101)/4-12.5,(grid+pedestal).T)})
    return result

def profile_information(records,basis,result):
    a=result['shape'];nshape=len(a);blocks=[]
    for i,r in enumerate(records):
        p=result['parameters'][nshape+9*i:nshape+9*(i+1)]
        cx,cy,lf,bg,bx,by,bxx,bxy,byy=p;x,y=r['x']-cx,r['y']-cy
        b,dx,dy=basis.evaluate(x,y,r);corr=np.exp(b@a);base=r['interpolator'].ev(x,y)
        px=r['interpolator'].ev(x,y,dx=1);py=r['interpolator'].ev(x,y,dy=1);f=np.exp(lf);point=f*base*corr
        js=point[:,None]*b/r['error'][:,None]
        jn=np.column_stack([-f*corr*(px+base*(dx@a)),-f*corr*(py+base*(dy@a)),point,np.ones(len(x)),
            r['x'],r['y'],r['x']**2,r['x']*r['y'],r['y']**2])/r['error'][:,None]
        blocks.append(js-jn@np.linalg.lstsq(jn,js,rcond=1e-12)[0])
    singular=svdvals(np.concatenate(blocks));rank=int(np.sum(singular>singular[0]*1e-10))
    return {'native_data_rows':sum(len(r['x']) for r in records),'free_local_nuisance_parameters':9*len(records),
            'shape_coefficients':nshape,'profiled_shape_numerical_rank_relative1e-10':rank,
            'unidentified_directions_including3_gauge':nshape-rank,
            'linearized_effective_shape_df_with_prior':float(np.sum(singular*singular/(singular*singular+25))),
            'largest_smallest_singular':[float(singular[0]),float(singular[-1])],
            'scope':'Conditional point covariances exclude response-learning uncertainty; numerical rank is a local linearization, not a unique physical optical solution.'}

def residuals(records,basis,result,free_shape=True):
    a=result['shape'];offset=len(a) if free_shape else 0;rows=[]
    for i,r in enumerate(records):
        p=result['parameters'][offset+9*i:offset+9*(i+1)]
        cx,cy,lf,bg,bx,by,bxx,bxy,byy=p;x,y=r['x']-cx,r['y']-cy
        log=basis.evaluate(x,y,r)[0]@a;point=np.exp(lf)*r['interpolator'].ev(x,y)*np.exp(log)
        model=point+bg+bx*r['x']+by*r['y']+bxx*r['x']**2+bxy*r['x']*r['y']+byy*r['y']**2
        rr=(r['data']-model)/r['error'];rad=np.hypot(x,y)
        ann=[]
        for lo,hi in [(0,1),(1,2),(2,4),(4,6),(6,12)]:
            m=(rad>=lo)&(rad<hi);ann.append({'r':[lo,hi],'n':int(m.sum()),'chi2':float(np.sum(rr[m]**2))})
        rows.append({'name':r['name'],'label':r['label'],'exposure':r['exposure'],'annuli':ann,
            'pixel_x':r['pixel_x'].tolist(),'pixel_y':r['pixel_y'].tolist(),'pixel_data':r['data'].tolist(),
            'pixel_error':r['error'].tolist(),'pixel_model':model.tolist(),'standardized_residual':rr.tolist(),
            'centroid_raw_pixel':(p[:2]+r['origin']).tolist(),'flux_raw':float(np.exp(lf))})
    return rows

def main():
    records=get_targets();names=sorted(set(r['exposure'] for r in records))
    # Retain exactly the prior field trial's reference-grid gauge; this uses the
    # independent library file, not field-star SCI values or lens predictions.
    grid,_=read_primary_fits(ART/'standard_psf/jw01198007001_02101_00001_nrcb3_FIELD_standard_epsf.fits')
    basis=FineBasis(grid+max(0.,-float(grid.min()))*1.1+1e-8)
    report={'status':'CONSUMED_TARGET_DATA_LOCAL_RESPONSE_DEVELOPMENT_ONLY','native_errors_unchanged':True,
      'physical_predictions_used_for_calibration':False,'all_F115W_targets_previously_consumed':True,
      'family':{'positive_shared_log_response':True,'RBF_spacing_pixel':.5,'RBF_width_pixel':.35,'coefficient_prior_sd':.2,
          'coefficient_count':basis.null.shape[1],'local_nuisance':'independent2Dcentroid,positive flux,quadraticbackground for everypoint'},
      'uncertainty_scope':'Known native ERR conditional on learned PSF. Calibration posterior/systematic uncertainty not propagated; no covariance rescaling.',
      'validation_pass':False,'folds':[]}
    for hold,name in enumerate(names):
        trainrecs=[r for r in records if r['exposure']!=name];testrecs=[r for r in records if r['exposure']==name]
        print(json.dumps({'starting_fold':hold+1,'training_images':len(trainrecs),'held_images':len(testrecs)}),flush=True)
        train=fit(trainrecs,basis,max_nfev=90);test=fit(testrecs,basis,a_fixed=train['shape'],max_nfev=90)
        info=profile_information(trainrecs,basis,train)
        row={'held_exposure':name,'train':serial(train),'test':serial(test),'rank':info,
             'held_pixel_results':residuals(testrecs,basis,test,free_shape=False),'validation_pass':False}
        report['folds'].append(row)
        np.savez(HERE/f'local_quasar_fold{hold+1}.npz',shape=train['shape'],parameters=train['parameters'],test_parameters=test['parameters'])
        (HERE/'local_quasar_results.json').write_text(json.dumps(report,indent=2))
        print(json.dumps({'held_exposure':hold+1,'training_success':train['success'],'training_nfev':train['nfev'],
          'held_chi2_per_dof':[r['chi2']/r['dof_conditional_on_shape'] for r in test['fits']],
          'rank':info['profiled_shape_numerical_rank_relative1e-10'],'effective_df':info['linearized_effective_shape_df_with_prior']}),flush=True)
    print(json.dumps({'starting_all_image_fit':16}),flush=True)
    allfit=fit(records,basis,max_nfev=100);report['all_image_fit']=serial(allfit)
    report['all_rank']=profile_information(records,basis,allfit);report['all_pixel_results']=residuals(records,basis,allfit)
    report['held_exposure_chi2_total']=sum(r['chi2'] for f in report['folds'] for r in f['test']['fits'])
    report['held_exposure_nominal_dof_conditional_on_shape']=sum(r['dof_conditional_on_shape'] for f in report['folds'] for r in f['test']['fits'])
    np.savez(HERE/'local_quasar_shape.npz',shape=allfit['shape'],parameters=allfit['parameters'],nodes=basis.nodes,gauge=basis.gauge,width=basis.width)
    (HERE/'local_quasar_results.json').write_text(json.dumps(report,indent=2))
    print(json.dumps({'held_total_chi2':report['held_exposure_chi2_total'],'held_nominal_dof':report['held_exposure_nominal_dof_conditional_on_shape'],
         'all_training_success':allfit['success'],'all_rank':report['all_rank']}),flush=True)

if __name__=='__main__':main()
