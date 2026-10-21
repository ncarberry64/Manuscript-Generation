"""Converge the existing optical family on actual common first-two-group data."""
from pathlib import Path
import argparse,hashlib,sys,json
import numpy as np

ROOT=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT/'scripts'));sys.path.insert(0,str(ROOT/'code'))
import run_he0435_joint_scene as prior
import strong_lensing_response_completion as rc
import strong_lensing_early_rate as early
from scipy.optimize import least_squares
js=prior.js
OUT=ROOT/'artifacts/strong_lensing_response_completion_20261021/common_early'


def original_metric(problem,p):
    prior_residual=problem.evaluate(p)[0];rows=[];offset=0;total=0.
    for record in problem.records:
        n=len(record['data']);mean=record['data']+prior_residual[offset:offset+n]*record['error'];offset+=n
        r=(mean-record['original_full_native_SCI'])/record['original_full_native_ERR'];q=float(r@r);total+=q
        rows.append({'name':record['name'],'Q_original_native_ERR':q,'pixels':n})
    return {'scope':'Uncorrected low-charge optical mean compared with original full-native SCI and unchanged originalERR on matching selectedrows. This is a comparator, not the missing full nonlinear history-forward prediction.','Q':total,'records':rows}


def main():
    parser=argparse.ArgumentParser();parser.add_argument('--label',choices=['all_data','training'],default='all_data');parser.add_argument('--resume',action='store_true');args=parser.parse_args()
    dest=OUT/args.label;dest.mkdir(parents=True,exist_ok=True)
    records,fields,eta0,wcss,origin=prior.load_inputs();histories={};input_hash={}
    paths=sorted((ROOT/'data/strong_lensing_response_completion_20261021/histories').glob('*_histories.npz'))
    for path in paths:
        a=np.load(path);histories[path.name.replace('_histories.npz','_cal.fits')]=a;input_hash[str(path.relative_to(ROOT)).replace('\\','/')]=hashlib.sha256(path.read_bytes()).hexdigest()
    for r in records:
        h=histories[r['exposure']];early.adapt_record(r,h,list(h['labels']).index(r['name'][-1]))
    for r in fields:early.adapt_record(r,histories[r['exposure']],4)
    exposures=sorted(histories);train=records if args.label=='all_data' else [r for r in records if r['exposure']!=exposures[-1]]
    target=js.SceneProblem(train);coupled=rc.CoupledScene(target,fields,eta0)
    previous=OUT.parent/args.label/'solution.npz';initial=np.load(previous)['parameters']
    if args.resume and (dest/'checkpoint.npz').exists():initial=np.load(dest/'checkpoint.npz')['parameters']
    first=coupled.diagnostics(initial);prior.dump(dest/'initial.json',first)
    def callback(x,d):
        np.savez(dest/'checkpoint.npz',parameters=x);prior.dump(dest/'progress.json',d)
        print(args.label,'early iter',len(coupled.history),'objective',d['objective'],'KKT',d['scaled_projected_KKT_gradient_max'],flush=True)
    fit,final=coupled.solve(initial,callback,max_nfev=400)
    p=fit.x[:coupled.n];eta=fit.x[coupled.eta_start:coupled.field_start];r,j=coupled.evaluate(fit.x,True);cov=js.svd_fit(j,scales=coupled.scales)
    ratio=coupled.history[-1]['relative_objective_change'] if coupled.history else None
    result={'scope':'Actual matched first2GROUP temporal-projection diagnostic on consumed F115W targets and matchingfield; same318positivefamily/coherentsky/priors. Independent E4target held from training.',
      'noise_scope':'Fixed positive nativeillumination pilot gives exactκ0Poisson frame-time contraction plus independent post-correction CDS marginal hypothesis. No SIRS factor or noise multiplier added. This does not establish raw/read/reference crosscovariance, source_complete_noise=false. Objectiveconvergence is conditional mean convergence under this fixed pilot metric, not a completed parameter-dependent statistical likelihood.',
      'initial':first,'final':final,'converged':bool(final['scaled_projected_KKT_gradient_max']<=1e-3 and ratio is not None and ratio<=1e-10),
      'convergence_definition':'Projected scaledhalf-objectiveKKT<=1e-3 andlastrelativeobjective<=1e-10; sameoriginalmeanparameters/priors.','scipy_nfev':fit.nfev,'scipy_message':fit.message,
      'source_complete_noise':False,'joint_rank':cov['rank'],'joint_columns':cov['columns'],'history':coupled.history,
      'original_full_native_mean_comparator':original_metric(target,p),'early_product_receipts':[r['early_product_receipt'] for r in train+fields],
      'covariance_scope':'UNCONSTRAINED localcurvature under declaredconditionalpilotnoise, fullsharedresponse/sky/fieldnuisance correlations; no physicalcoverage or truncated-boundposterior.',
      'input_sha256':input_hash,'source_sha256':{p:hashlib.sha256((ROOT/p).read_bytes()).hexdigest() for p in ['code/strong_lensing_early_rate.py','scripts/run_he0435_common_early.py','code/strong_lensing_response_completion.py']},
      'physical_validation_pass':False,'physical_admission':False,'physical_guard':2e-10}
    prior.dump(dest/'result.json',result);np.savez(dest/'solution.npz',parameters=fit.x,target_parameters=p,eta=eta,joint_covariance=cov['covariance'],scales=coupled.scales)
    if args.label=='training':
        held=[r for r in records if r['exposure']==exposures[-1]];prior.set_eta(held,eta);test=js.SceneProblem(held,scene_fixed=target.sky(p));start=test.initial();lo,hi=test.bounds()
        heldfit=least_squares(lambda q:test.evaluate(q)[0],start,jac=lambda q:test.evaluate(q)[1],bounds=(lo,hi),x_scale=test.scales,ftol=1e-14,xtol=1e-14,gtol=1e-9,max_nfev=200)
        rr,aa,_,rows=test.evaluate(heldfit.x);pg,_,_=rc.projected_gradient(heldfit.x,aa.T@rr,lo,hi,test.scales)
        prior.dump(dest/'held_target.json',{'scope':'Alreadyconsumed E4target, trainedeta/skyfixed andindependentheldpoint/flux/plane fit. Sameconditionalpilotnoise; no unusedvalidation.','early_raw_chi2':float(rr@rr),'projected_KKT_max':float(abs(pg).max()),'records':rows,
           'original_full_native_mean_comparator':original_metric(test,heldfit.x),'physical_validation_pass':False,'physical_admission':False,'physical_guard':2e-10})
        np.savez(dest/'held_solution.npz',parameters=heldfit.x,eta=eta,shared_sky=target.sky(p))
    print(args.label,'EARLYFINAL',result['converged'],final['objective'],final['scaled_projected_KKT_gradient_max'],'originalQ',result['original_full_native_mean_comparator']['Q'],flush=True)


if __name__=='__main__':main()
