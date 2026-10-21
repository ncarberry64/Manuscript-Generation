"""Strict held-target solve and raw diagnostics for the converged static mean.

No target from exposure4 updates the training response or shared sky.
"""
from pathlib import Path
import hashlib,json,sys
import numpy as np
from scipy.optimize import least_squares
from astropy.io import fits

ROOT=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT/'scripts'));sys.path.insert(0,str(ROOT/'code'))
import run_he0435_joint_scene as prior
import strong_lensing_response_completion as rc
import export_he0435_joint_scene_residuals as export
js=prior.js
OUT=ROOT/'artifacts/strong_lensing_response_completion_20261021/static_diagnostics'


def maps(label,problem,p,eta,full):
    prior.set_eta(problem.records,eta)
    report=export.export(label,problem,p,eta,full,float(problem.evaluate(p)[0]@problem.evaluate(p)[0]))
    data=np.load(OUT/f'{label}.npz');stats=[]
    for i,row in enumerate(report['records']):
        mask=data['mask'][i];z=data['standardized_residual'][i];x=data['x'][i];y=data['y'][i]
        center=p[3*i:3*i+2];radius=np.hypot(x-center[0],y-center[1]);bins=[]
        for lo,hi in [(0,1.5),(1.5,3),(3,5),(5,13)]:
            keep=mask&(radius>=lo)&(radius<hi);vals=z[keep]
            bins.append({'radius_pixel':[lo,hi],'pixels':int(keep.sum()),'Q':float(vals@vals),'standardized_mean':float(vals.mean()) if len(vals) else None})
        correlations={}
        for name,first,second in [('x',z[:,:-1],z[:,1:]),('y',z[:-1],z[1:])]:
            good=np.isfinite(first)&np.isfinite(second);a=first[good];b=second[good]
            correlations[name]=float(np.corrcoef(a,b)[0,1])
        stats.append({**row,'radial_residual_bins':bins,'adjacent_standardized_residual_correlation':correlations})
    report['records']=stats
    return report


def main():
    OUT.mkdir(exist_ok=True);export.OUT=OUT
    records,fields,eta0,wcss,origin=prior.load_inputs();exposures=sorted({r['exposure'] for r in records})
    training=js.SceneProblem([r for r in records if r['exposure']!=exposures[-1]])
    saved=np.load(OUT.parent/'training/solution.npz');trainp=saved['target_parameters'];eta=saved['eta'];joint=saved['joint_covariance']
    prior.set_eta(records,eta);held=[r for r in records if r['exposure']==exposures[-1]];test=js.SceneProblem(held,scene_fixed=training.sky(trainp))
    start=np.load(OUT.parent/'training/held_target_solution.npz')['parameters']
    lo,hi=test.bounds()
    fit=least_squares(lambda p:test.evaluate(p)[0],start,jac=lambda p:test.evaluate(p)[1],bounds=(lo,hi),x_scale=test.scales,
                      ftol=1e-14,xtol=1e-14,gtol=1e-9,max_nfev=200)
    r,a,b,rows=test.evaluate(fit.x,True);g=a.T@r;pg,atlo,athi=rc.projected_gradient(fit.x,g,lo,hi,test.scales)
    scenej=np.zeros((len(r),36));offset=0
    for rec in held:
        n=len(rec['x']);scenej[offset:offset+n,9*rec['patch']:9*rec['patch']+9]=rec['scene_matrix']/rec['error'][:,None];offset+=n
    shared=np.r_[np.arange(training.size,training.size+len(eta)),np.arange(training.sky_start,training.bg_start)]
    vshared=joint[np.ix_(shared,shared)];conditional,cut,sensitivity,rank=js.propagated_cut_covariance(a,np.c_[b,scenej],vshared,scales=test.scales)
    take=np.array([3*i+d for i in range(test.n) for d in [0,1]])
    moment=np.c_[np.vstack([rec['response'].moment_derivative for rec in held]),np.zeros((2*len(held),36))]
    held_transfer=sensitivity[take]+moment;heldcov=conditional[np.ix_(take,take)]+held_transfer@vshared@held_transfer.T
    train_transform=np.zeros((2*training.n,len(joint)))
    train_take=np.array([3*i+d for i in range(training.n) for d in [0,1]])
    train_transform[np.arange(len(train_take)),train_take]=1
    train_transform[:,training.size:training.size+len(eta)]=np.vstack([rec['response'].moment_derivative for rec in training.records])
    traincov=train_transform@joint@train_transform.T;cross=train_transform@joint[:,shared]@held_transfer.T
    fullcov=np.block([[traincov,cross],[cross.T,heldcov]])
    allproblem=js.SceneProblem(records);allp=allproblem.initial();byname={rec['name']:i for i,rec in enumerate(training.records)};heldname={rec['name']:i for i,rec in enumerate(held)}
    for i,rec in enumerate(records):
        source,idx=(trainp,byname[rec['name']]) if rec['name'] in byname else (fit.x,heldname[rec['name']])
        allp[3*i:3*i+3]=source[3*idx:3*idx+3]
    allp[allproblem.sky_start:allproblem.bg_start]=training.sky(trainp)
    allp[allproblem.bg_start:]=np.r_[trainp[training.bg_start:],fit.x[test.bg_start:]]
    # Recover full native stamps, preserving invalid pixels rather than
    # replacing them with a compact mask-only vector.
    full={}
    for e in exposures:
        with fits.open(ROOT/'data/strong_lensing_geometry_20261021/native'/e,memmap=False) as h:
            for rec in records+fields:
                if rec['exposure']!=e:continue
                x0,y0,x1,y1=rec['bounds'];s,er,dq=[np.array(h[k].data[y0:y1,x0:x1]) for k in ['SCI','ERR','DQ']]
                mask,_=prior.jwst_native_rate_mask(s,er,dq,primary_header=h[0].header);full[rec['name']]=(s,er,dq,mask)
    reports={}
    reports['held_target_strict']=maps('held_target_strict',test,fit.x,eta,full)
    reports['training_target']=maps('training_target',training,trainp,eta,full)
    alld=np.load(OUT.parent/'all_data/solution.npz');reports['all_target']=maps('all_target',allproblem,alld['target_parameters'],alld['eta'],full)
    # All-data field nuisance parameters are actually retained by the new
    # optimizer, so field mean/residual maps require no additional fit.
    for label,solution,target in [('all_field',alld,allproblem),('training_field',saved,training)]:
        coupled=rc.CoupledScene(target,fields,eta0);params=solution['parameters'];et=params[coupled.eta_start:coupled.field_start]
        arrays={k:[] for k in ['name','x','y','SCI','ERR','DQ','mask','mean','raw_residual','standardized_residual']};stats=[]
        for i,rec in enumerate(fields):
            rec['response'].update(et);nu=params[coupled.field_start+6*i:coupled.field_start+6*i+6]
            rr,_,_=coupled.field_application(nu,rec);s,er,dq,mask=full[rec['name']];x0,y0,x1,y1=rec['bounds'];yy,xx=np.mgrid[y0:y1,x0:x1]
            mm=np.full(mask.shape,np.nan);raw=np.full(mask.shape,np.nan);zz=np.full(mask.shape,np.nan);zz[mask]=rr;raw[mask]=rr*rec['error'];mm[mask]=rec['data']+raw[mask]
            for key,val in [('name',rec['name']),('x',xx),('y',yy),('SCI',s),('ERR',er),('DQ',dq),('mask',mask),('mean',mm),('raw_residual',raw),('standardized_residual',zz)]:arrays[key].append(val)
            stats.append({'name':rec['name'],'Q':float(rr@rr),'pixels':len(rr),'centroid':nu[:2]})
        path=OUT/f'{label}.npz';np.savez_compressed(path,**{k:np.asarray(v) for k,v in arrays.items()});reports[label]={'file':str(path.relative_to(ROOT)).replace('\\','/'),'sha256':hashlib.sha256(path.read_bytes()).hexdigest(),'records':stats}
    prior.set_eta(records,eta)
    result={'scope':'Strict independent held-point solve after converged static training; saved target and field mean replay, raw spatial residual diagnostics. AllF115Wconsumed; no training use of heldE4target pixels.',
            'held_target':{'raw_chi2':float(r@r),'projected_KKT_scaled_gradient_max':float(abs(pg).max()),'fit_nfev':fit.nfev,'scipy_message':fit.message,'active_lower_indices':np.flatnonzero(atlo),'active_upper_indices':np.flatnonzero(athi),'records':rows},
            'covariance_scope':'Full32centroid train+held covariance: shared training eta+sky covariance and all their crossblocks propagated into held measurements; conditional nativeERR targetnoise independent. UNCONSTRAINED localcurvature, not truncatedpositiveposterior or calibratedphysicalcoverage.',
            'full32_moment_centroid_covariance_pixel2':fullcov,'minimum_eigenvalue_pixel2':float(np.linalg.eigvalsh(fullcov).min()),
            'all_exposure_geometry_training_plus_held':prior.geometry_from_centroids(allproblem,allp,fullcov,wcss,origin,'baseline_relative_moment'),
            'raw_maps':reports,'physical_validation_pass':False,'physical_admission':False,'physical_guard':2e-10,
            'source_sha256':{p:hashlib.sha256((ROOT/p).read_bytes()).hexdigest() for p in ['scripts/finalize_he0435_static_response.py','code/strong_lensing_response_completion.py']}}
    prior.dump(OUT/'result.json',result);np.savez(OUT/'strict_held_solution.npz',parameters=fit.x,eta=eta,shared_sky=training.sky(trainp),full32_moment_centroid_covariance=fullcov)
    print('held strict',float(r@r),'KKT',float(abs(pg).max()),'nfev',fit.nfev,flush=True)


if __name__=='__main__':main()
