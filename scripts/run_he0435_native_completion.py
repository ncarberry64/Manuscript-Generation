"""Converge consumed full-native IPC/history/optical means, then replay OLS_C."""
from pathlib import Path
import argparse,hashlib,json,sys,time
import numpy as np
from astropy.io import fits
ROOT=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT/'scripts'));sys.path.insert(0,str(ROOT/'code'))
import run_he0435_joint_scene as prior
from strong_lensing_native_completion import NativeCoupledScene
from strong_lensing_ipc_response import inverse_deconvolution_kernel,convolution_matrix
js=prior.js
BASE=ROOT/'data/strong_lensing_response_completion_20261021'
OUT=ROOT/'artifacts/strong_lensing_response_completion_20261021/native_coupled'


def build(label):
    records,fields,eta0,wcss,origin=prior.load_inputs();exposures=sorted(wcss)
    histories={}
    for path in sorted((BASE/'histories').glob('*_histories.npz')):
        with np.load(path) as f:histories[path.name.replace('_histories.npz','_cal.fits')]={k:f[k] for k in f.files}
    train=records if label=='all_data' else [r for r in records if r['exposure']!=exposures[-1]]
    target=js.SceneProblem(train);saved=np.load(OUT.parent/label/'solution.npz');initial=np.r_[saved['parameters'],0.]
    prior.set_eta(train,saved['eta'])
    for record in fields:record['response'].update(saved['eta'])
    gkernel=fits.getdata(BASE/'native_cache/jwst_nircam_ipc_0022.fits',1)
    info=inverse_deconvolution_kernel(gkernel,order=8)
    g=convolution_matrix(gkernel,(65,65));h=convolution_matrix(info['kernel'],(65,65))
    problem=NativeCoupledScene(target,fields,saved['eta'],histories,wcss,origin,g,h,initial)
    return problem,initial,records,wcss,origin,info


def replay_native(problem,x,dest):
    # Official estimator value application, distinct from its differentiable
    # fixed-branch double source-equation extension used during optimization.
    from apply_he0435_bounded_native_history import native_apply
    from stdatamodels.jwst import datamodels
    from stdatamodels.jwst.datamodels import dqflags
    from stcal.ramp_fitting import ramp_fit
    problem.evaluate(x,True);rows=[];arrays={};templates={}
    for ctx,(state,plane,weights),record in zip(problem.contexts,problem.last_states,problem.mean_records):
        rec=ctx['original'];e=rec['exposure'];stem=e.replace('_cal.fits','')
        if e not in templates:
            with datamodels.RampModel(BASE/'native_cache'/f'{stem}_uncal.fits') as source:
                templates[e]=ramp_fit.create_ramp_fit_class(source,'OLS_C',dqflags.pixel,False)
        native=native_apply(templates[e],ctx['history'],ctx['patch'],state['groups_DN'],state['zeroframe_DN']).ravel()[ctx['selected']]
        mean=native*ctx['native_SCI_factor']+plane;residual=mean-rec['data'];rr=residual/rec['error']
        rows.append({'name':rec['name'],'field':ctx['field'],'original_native_ERR_chi2':float(rr@rr),
            'raw_RMS_SCI':float(np.sqrt(np.mean(residual**2))),
            'max_abs_native_minus_continuous_mean_SCI':float(np.max(abs(mean-record['mean_SCI']))),
            'continuous_original_ERR_chi2':record['raw_chi2'],'receipt':state['diagnostics']})
        arrays[rec['name']+'_x']=rec['x'];arrays[rec['name']+'_y']=rec['y'];arrays[rec['name']+'_SCI']=rec['data'];arrays[rec['name']+'_ERR']=rec['error']
        arrays[rec['name']+'_mean']=mean;arrays[rec['name']+'_continuous_mean']=record['mean_SCI'];arrays[rec['name']+'_raw_residual']=residual
        arrays[rec['name']+'_standardized_residual']=rr;arrays[rec['name']+'_native_groups_DN']=state['groups_DN'];arrays[rec['name']+'_native_zeroframe_DN']=state['zeroframe_DN']
        arrays[rec['name']+'_effective_photon_rate_e_s']=problem.optical(x,ctx)[0].reshape(65,65)
    np.savez_compressed(dest/'native_replay_maps.npz',**arrays)
    return {'records':rows,'target_original_native_ERR_Q':sum(r['original_native_ERR_chi2'] for r in rows if not r['field']),
        'field_original_native_ERR_Q':sum(r['original_native_ERR_chi2'] for r in rows if r['field']),
        'scope':'Official OLS_C float32 value replay of saved converged continuous-source-equation mean. No independent observations or altered ERR.'}


def main():
    parser=argparse.ArgumentParser();parser.add_argument('--label',choices=['all_data','training'],default='all_data');parser.add_argument('--resume',action='store_true');parser.add_argument('--max-nfev',type=int,default=150);parser.add_argument('--benchmark',action='store_true');args=parser.parse_args()
    dest=OUT/args.label;dest.mkdir(parents=True,exist_ok=True)
    started=time.perf_counter();problem,initial,records,wcss,origin,ipc=build(args.label)
    print('BUILT seconds',round(time.perf_counter()-started,2),'columns',problem.size,flush=True)
    if args.resume:initial=np.load(dest/'checkpoint.npz')['parameters']
    first=problem.diagnostics(initial);prior.dump(dest/'initial.json',first)
    print('INITIAL',first['objective'],first['target_raw_chi2'],first['field_raw_chi2'],'KKT',first['scaled_projected_KKT_gradient_max'],'seconds',round(time.perf_counter()-started,2),flush=True)
    if args.benchmark:
        np.savez(dest/'initial_benchmark.npz',parameters=initial);return
    def callback(x,d):
        np.savez(dest/'checkpoint.npz',parameters=x);prior.dump(dest/'progress.json',d)
        print(args.label,'NATIVE_ITER',len(problem.history),'Q',d['objective'],'KKT',d['scaled_projected_KKT_gradient_max'],'kappa',x[-1],'seconds',round(d['seconds'],2),flush=True)
    fit,final=problem.solve(initial,callback,max_nfev=args.max_nfev)
    r,j=problem.evaluate(fit.x,True);info=js.svd_fit(j,scales=problem.scales)
    ratio=problem.history[-1]['relative_objective_change'] if problem.history else None
    kkt=final['scaled_projected_KKT_gradient_max']<=1e-3;relative=ratio is not None and ratio<=1e-10
    np.savez(dest/'solution.npz',parameters=fit.x,target_parameters=fit.x[:problem.n],eta=fit.x[problem.eta_start:problem.field_start],kappa_e_inverse=fit.x[-1],joint_covariance=info['covariance'],singular_values=info['singular_values'],scales=problem.scales)
    sourcepaths=['code/strong_lensing_native_completion.py','scripts/run_he0435_native_completion.py','code/strong_lensing_native_calibration.py','code/strong_lensing_native_ols.py','code/strong_lensing_calibration_order.py','code/strong_lensing_bounded_collection.py','code/strong_lensing_joint_scene.py']
    result={'scope':'Consumed F115W full-native coupled optical318/coherent36sky/point/field nuisances and positive kappa, with source-calibrated IPC/linearity ordering and actual fixed native selection. No new PSF family or lens objective.',
        'initial':first,'final':final,'kappa_e_inverse':float(fit.x[-1]),'converged':bool(kkt and relative),'projected_KKT_converged':kkt,'relative_objective_converged':relative,
        'convergence_definition':'Scaled half-objective projected KKT<=1e-3 and last accepted relative objective<=1e-10 for continuous fixed-branch native OLS source equations; original float32 native estimator replayed separately.',
        'scipy_message':fit.message,'scipy_nfev':fit.nfev,'history':problem.history,'joint_rank':info['rank'],'joint_columns':info['columns'],
        'domain_rejected_trials':problem.rejected_domain_trials,'fixed_domain_prescriptions':[{ 'name':ctx['name'],**ctx['domain_prescription']} for ctx in problem.contexts],
        'native_branch_nonunique_count':sum(row['native_branch_nonunique_count'] for row in problem.mean_records),'native_SNR_boundary_count':sum(row['native_SNR_boundary_count'] for row in problem.mean_records),
        'statistical_scope':'Original native ERR fixed conditional estimating metric. G0 electronic-intercept conditioning is nonlocal through IPC; full residual noise covariance and parameter-dependent Gaussian NLL are not established by this fit. The source-calibration electronic covariance is a separate continuation, not automatically added to native ERR.',
        'remaining_physical_relations':['Post-full-well charge/overflow law is a conditional storage-ceiling leading response, not measured closure.','Reference SAT threshold translated from initial intercept is not physical full well.','Photon point+sky versus additive electronic plane partition remains conditional.','Time-dependent optical response, source SED, and larger extended scene remain possible inadequacy sources.','Nonunique native median derivatives and adaptive branch changes remain individually reported.'],
        'covariance_scope':'Full UNCONSTRAINED local curvature of this conditional estimating problem; active bounds recorded, not truncated posterior or coverage validation.',
        'IPC_inverse_tail_l1_bound':ipc['tail_l1_bound'],'error_rescaling':False,'unused_observations_opened':False,'physical_validation_pass':False,'physical_admission':False,'physical_guard':2e-10,
        'source_sha256':{p:hashlib.sha256((ROOT/p).read_bytes()).hexdigest() for p in sourcepaths},
        'input_sha256':{str(p.relative_to(ROOT)).replace('\\','/'):hashlib.sha256(p.read_bytes()).hexdigest() for p in sorted((BASE/'histories').glob('*_histories.npz'))}}
    prior.dump(dest/'result.json',result)
    try:
        result['official_native_value_replay']=replay_native(problem,fit.x,dest)
    except ImportError as exc:result['official_native_value_replay_pending']=str(exc)
    prior.dump(dest/'result.json',result)
    print('FINAL',args.label,result['converged'],final['objective'],'KKT',final['scaled_projected_KKT_gradient_max'],'kappa',fit.x[-1],flush=True)


if __name__=='__main__':main()
