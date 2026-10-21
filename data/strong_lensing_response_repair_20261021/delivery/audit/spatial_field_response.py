"""Bounded field-only tests of residual response variation along dither track.

No lens prediction, target score, unused filter or error rescaling is used.
Two finite alternatives: one amplitude gradient of the learned mean correction,
and one smooth low-resolution log-response gradient along detector track.
Neither extrapolates an unmeasured perpendicular spatial derivative to targets.
"""
import json,sys
from pathlib import Path
import numpy as np
from scipy.linalg import svd
HERE=Path(__file__).resolve().parent;sys.path.insert(0,str(HERE))
from positive_field_response import get_fields,Basis,fit,serial,ART

class SpatialBasis:
    def __init__(self,base,track_center,track_direction,track_scale,kind,amplitude_shape=None):
        self.base=base;self.center=track_center;self.direction=track_direction;self.scale=track_scale;self.kind=kind
        self.gx,self.gy=base.gx,base.gy;self.nbase=base.null.shape[1]
        if kind=='amplitude_gradient':
            self.amplitude_shape=amplitude_shape;self.nslope=1
        else:
            self.slope=Basis(reference_grid,spacing=3.,width=1.5);self.nslope=self.slope.null.shape[1]
        self.null=np.zeros((1,self.nbase+self.nslope))
    def t(self,rec):return float((rec['origin']-self.center)@self.direction/self.scale)
    def evaluate(self,x,y,rec=None):
        b=self.base.evaluate(x,y);t=self.t(rec)
        if self.kind=='amplitude_gradient':s=[(v@self.amplitude_shape)[:,None] for v in b]
        else:s=self.slope.evaluate(x,y)
        return tuple(np.column_stack([v,t*w]) for v,w in zip(b,s))
    def response_grid(self,rec,a):
        log=self.base.grid_basis@a[:self.nbase];t=self.t(rec)
        if self.kind=='amplitude_gradient':log=log+t*(self.base.grid_basis@self.amplitude_shape)*a[-1]
        else:log=log+t*(self.slope.grid_basis@a[self.nbase:])
        return rec['grid']*np.exp(log.reshape(101,101))

def residual_summary(recs,basis,result):
    a=result['shape'];nshape=len(a);rows=[]
    for i,r in enumerate(recs):
        params=result['parameters'][nshape+6*i:nshape+6*(i+1)]
        cx,cy,lf,bg,bx,by=params;x,y=r['x']-cx,r['y']-cy
        log=basis.evaluate(x,y,r)[0]@a;p=r['interpolator'].ev(x,y)*np.exp(log)*np.exp(lf)
        pred=p+bg+bx*r['x']+by*r['y'];res=(r['data']-pred)/r['error'];radius=np.hypot(x,y)
        sat=(r['pixel_dq']&2)!=0
        def stats(m):
            if not np.any(m):return {'n':0}
            return {'n':int(m.sum()),'chi2':float(np.sum(res[m]**2)),'rms_standardized':float(np.sqrt(np.mean(res[m]**2))),
                    'mean_standardized':float(np.mean(res[m])),'mean_SCI_residual':float(np.mean((r['data']-pred)[m]))}
        rows.append({'name':r['name'],'detector_origin_pixel':r['origin'].tolist(),
          'subpixel_phase_fitted':((r['origin']+params[:2])%1).tolist(),'fitted_flux_raw':float(np.exp(lf)),
          'background_SCI':[float(bg),float(bx),float(by)],'max_SCI':float(r['data'].max()),
          'saturated_flagged_valid_rate':stats(sat),'not_saturated_flagged':stats(~sat),
          'annuli':[{'radius':[lo,hi],**stats((radius>=lo)&(radius<hi))} for lo,hi in [(0,1),(1,2),(2,4),(4,6),(6,12)]],
          'pixel_x':r['pixel_x'].tolist(),'pixel_y':r['pixel_y'].tolist(),'standardized_residual':res.tolist(),
          'pixel_model':pred.tolist(),'pixel_error':r['error'].tolist(),'pixel_data':r['data'].tolist()})
    return rows

def track_report(recs):
    points=np.array([r['origin'] for r in recs]);center=points.mean(axis=0);u,s,vt=np.linalg.svd(points-center,full_matrices=False)
    direction=vt[0]
    if direction[0]<0:direction=-direction
    perp=np.array([-direction[1],direction[0]])
    along=(points-center)@direction;across=(points-center)@perp;scale=max(abs(along))
    design=np.column_stack([np.ones(len(points)),along/scale,across/scale]);sing=np.linalg.svd(design,compute_uv=False)
    inv=np.linalg.inv(design.T@design);targets=[]
    for path in sorted((ART/'jwst_target_diagnostic/exposures').glob('*.json')):
        ex=json.loads(path.read_text())
        for p in ex['points']:
            pos=np.array(p['centroid_xy_pixel']);ax=(pos-center)@direction;ay=(pos-center)@perp
            d=np.array([1,ax/scale,ay/scale]);targets.append({'exposure':path.name,'label':p['label'],
               'along_track_pixel':float(ax),'perpendicular_to_track_pixel':float(ay),
               'affine_calibration_prediction_leverage':float(d@inv@d)})
    rep={'calibrator_center_pixel':center.tolist(),'track_unit_vector':direction.tolist(),'perp_unit_vector':perp.tolist(),
         'coordinate_singular_values_pixel':s.tolist(),'calibrator_along_track_pixel':along.tolist(),
         'calibrator_perp_track_pixel':across.tolist(),'affine_design_singular_values':sing.tolist(),
         'affine_design_condition_number':float(sing[0]/sing[-1]),'targets':targets,
         'interpretation':'Four calibrator positions are almost collinear. A full2D affine response is numerically invertible but its transverse derivative is extremely ill-conditioned; target locations require large transverse extrapolation.'}
    return rep,center,direction,scale

def main():
    global reference_grid
    recs=get_fields();reference_grid=recs[0]['grid'];base=Basis(reference_grid)
    track,c,d,s=track_report(recs);baseline=fit(recs,base)
    report={'status':'EXPLORATORY_FIELD_ONLY_SPATIAL_RESPONSE_NO_VALIDATION','native_errors_unchanged':True,
        'lens_predictions_used':False,'target_geometry_used_to_fit':False,'spatial_track':track,
        'shared_model':serial(baseline),'shared_residuals':residual_summary(recs,base,baseline),'models':[]}
    print(json.dumps({'track_condition':track['affine_design_condition_number'],
        'field_perp':track['calibrator_perp_track_pixel'],
        'target_perp_range':[min(abs(t['perpendicular_to_track_pixel']) for t in track['targets']),max(abs(t['perpendicular_to_track_pixel']) for t in track['targets'])]}),flush=True)
    for kind in ['amplitude_gradient','smooth_gradient']:
        model={'kind':kind,'folds':[]}
        for hold in range(4):
            trainrecs=[r for i,r in enumerate(recs) if i!=hold]
            common=fit(trainrecs,base)
            basis=SpatialBasis(base,c,d,s,kind,common['shape'])
            train=fit(trainrecs,basis,prior_sd=.2,max_nfev=80)
            test=fit([recs[hold]],basis,prior_sd=.2,a_fixed=train['shape'])
            model['folds'].append({'hold':hold+1,'shape_parameters':basis.null.shape[1],
                                    'train':serial(train),'test':serial(test)})
            print(json.dumps({'kind':kind,'held':hold+1,'train_success':train['success'],
              'held_chi2':test['fits'][0]['chi2'],'held_chi2_per_dof':test['fits'][0]['chi2']/283}),flush=True)
        basis=SpatialBasis(base,c,d,s,kind,baseline['shape']);allfit=fit(recs,basis,prior_sd=.2,max_nfev=80)
        model['all_field_fit']=serial(allfit);model['all_field_residuals']=residual_summary(recs,basis,allfit)
        model['held_exposure_chi2_total']=sum(f['test']['fits'][0]['chi2'] for f in model['folds'])
        report['models'].append(model);(HERE/'spatial_field_results.json').write_text(json.dumps(report,indent=2))
    print(json.dumps({'model_results':[{k:v for k,v in r.items() if k in ['kind','held_exposure_chi2_total']} for r in report['models']]}),flush=True)

if __name__=='__main__':main()
