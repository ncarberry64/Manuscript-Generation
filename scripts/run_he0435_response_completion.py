"""Converge the existing joint response/scene mean on consumed F115W rates."""
from pathlib import Path
from datetime import datetime,timezone
import argparse,hashlib,json,sys,time
import numpy as np

ROOT=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT/'code'));sys.path.insert(0,str(ROOT/'scripts'))
import strong_lensing_response_completion as rc
import run_he0435_joint_scene as prior
js=prior.js
OUT=ROOT/'artifacts/strong_lensing_response_completion_20261021'


def dump(path,value):prior.dump(path,value)


def main():
    parser=argparse.ArgumentParser();parser.add_argument('--label',choices=['all_data','training'],default='all_data');parser.add_argument('--max-nfev',type=int,default=400);parser.add_argument('--resume',action='store_true');args=parser.parse_args()
    dest=OUT/args.label;dest.mkdir(parents=True,exist_ok=True)
    initial_meta={'started_utc':datetime.now(timezone.utc).isoformat(),'scope':'Consumed F115W development, exact existing mean and priors; no physical lens objective or error scaling.','physical_validation_pass':False,'physical_admission':False,'physical_guard':2e-10}
    paths=['code/strong_lensing_response_completion.py','scripts/run_he0435_response_completion.py','code/strong_lensing_joint_scene.py','scripts/run_he0435_joint_scene.py']
    initial_meta['source_sha256']={p:hashlib.sha256((ROOT/p).read_bytes()).hexdigest() for p in paths}
    if (dest/'start.json').exists() and not args.resume:raise RuntimeError('Existing application; explicit --resume required')
    if not (dest/'start.json').exists():dump(dest/'start.json',initial_meta)
    records,fields,eta0,wcss,origin=prior.load_inputs()
    exposures=sorted({r['exposure'] for r in records});selected=records if args.label=='all_data' else [r for r in records if r['exposure']!=exposures[-1]]
    target=js.SceneProblem(selected);problem=rc.CoupledScene(target,fields,eta0)
    saved=np.load(prior.OUT/('all_data_linearized_arrays.npz' if args.label=='all_data' else 'held_exposure_training_linearized_arrays.npz'))
    previous_result=None
    if args.resume and (dest/'checkpoint.npz').exists():
        initial=np.load(dest/'checkpoint.npz')['parameters']
        if (dest/'result.json').exists():
            previous_result=json.loads((dest/'result.json').read_text())
            if not (dest/'termination_before_stationarity.json').exists():dump(dest/'termination_before_stationarity.json',previous_result)
    else:initial=problem.initial(saved['parameters_updated'],saved['eta_updated'])
    first=problem.diagnostics(initial);dump(dest/'initial.json',first)
    print(args.label,'initial',first['objective'],first['scaled_projected_KKT_gradient_max'],flush=True)
    def callback(x,d):
        np.savez(dest/'checkpoint.npz',parameters=x,scales=problem.scales)
        dump(dest/'progress.json',d)
        print(args.label,'iteration',len(problem.history),'objective',d['objective'],'KKT',d['scaled_projected_KKT_gradient_max'],'seconds',round(d['seconds'],2),flush=True)
    fit,final=problem.solve(initial,callback,max_nfev=args.max_nfev)
    r,j=problem.evaluate(fit.x,True);info=js.svd_fit(j,scales=problem.scales)
    targetp=fit.x[:problem.n];eta=fit.x[problem.eta_start:problem.field_start]
    n=target.n;take=np.array([3*i+d for i in range(n) for d in [0,1]])
    transform=np.zeros((2*n,problem.size));transform[np.arange(2*n),take]=1
    transform[:,problem.eta_start:problem.field_start]=np.vstack([rec['response'].moment_derivative for rec in target.records])
    centroid_covariance=transform@info['covariance']@transform.T
    targetr,targetj,_,rows=target.evaluate(targetp)
    relative=problem.history[-1]['relative_objective_change'] if problem.history else None
    kkt_pass=final['scaled_projected_KKT_gradient_max']<=1e-3
    relative_pass=relative is not None and relative<=1e-10
    result={'scope':'Nonlinear completion of exactly the existing conditional coarse positive response/coherent-sky mean, on consumed F115W rate and field observations. No new PSF family/lens fit/error inflation.',
            'initial':first,'final':final,'scipy_success':bool(fit.success),'scipy_status':int(fit.status),'scipy_message':fit.message,'scipy_nfev':fit.nfev,'scipy_njev':fit.njev,
            'convergence_definition':'Projected KKT gradient of half-objective in declared dimensionless scales <=1e-3 AND last accepted relative full-objective change <=1e-10. Bound signs checked explicitly; native ERR unchanged.',
            'projected_KKT_converged':kkt_pass,'relative_objective_converged':relative_pass,'converged':bool(kkt_pass and relative_pass),
            'last_relative_objective_change':relative,'target_fit':js.failure_summary(targetr,targetj),'target_records':rows,
            'joint_rank':info['rank'],'joint_columns':info['columns'],'full_column_rank':info['full_column_rank'],'svd_rcond':js.SVD_RCOND,
            'covariance_scope':'UNCONSTRAINED local curvature of the specified joint target+field+declared prior likelihood. Active bounds recorded, not a truncated posterior or physical coverage claim.',
            'moment_geometry':prior.geometry_from_centroids(target,targetp,centroid_covariance,wcss,origin),
            'baseline_relative_moment_geometry':prior.geometry_from_centroids(target,targetp,centroid_covariance,wcss,origin,'baseline_relative_moment'),
            'history':([] if previous_result is None else previous_result['history'])+problem.history,
            'previous_termination_preserved':previous_result is not None,'physical_validation_pass':False,'physical_admission':False,'physical_guard':2e-10,
            'source_sha256':initial_meta['source_sha256']}
    dump(dest/'result.json',result)
    np.savez(dest/'solution.npz',parameters=fit.x,column_scales=problem.scales,joint_covariance=info['covariance'],joint_singular_values=info['singular_values'],target_moment_centroid_covariance=centroid_covariance,eta=eta,target_parameters=targetp)
    if args.label=='training':
        held=[r for r in records if r['exposure']==exposures[-1]];prior.set_eta(held,eta);test=js.SceneProblem(held,scene_fixed=target.sky(targetp));heldfit=test.fit(max_nfev=200)
        rr,aa,bb,heldrows=test.evaluate(heldfit.x,True)
        dump(dest/'held_target.json',{'scope':'Already-consumed wholeE4target exposure held from joint response/scene training. All4separatefieldstarstamps remain in training; inherited targetROI/guesses available. No unused validation.','fit_success':bool(heldfit.success),'fit_optimality':float(heldfit.optimality),'fit_nfev':heldfit.nfev,'known_ERR_fit':js.failure_summary(rr,aa),'records':heldrows,'physical_validation_pass':False,'physical_admission':False,'physical_guard':2e-10})
        np.savez(dest/'held_target_solution.npz',parameters=heldfit.x,eta=eta,shared_sky=target.sky(targetp))
    print(args.label,'FINAL',json.dumps({k:result[k] for k in ['converged','projected_KKT_converged','relative_objective_converged','scipy_nfev']}),final['objective'],final['scaled_projected_KKT_gradient_max'],flush=True)


if __name__=='__main__':main()
