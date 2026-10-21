"""Internal E4 target check with native response/sky trained on E1--3 targets.

All four FIELD stamps remain calibration observations. Original E1 sky anchors
and target ROI/old initial guesses are inherited. No reserved data are opened.
"""
from pathlib import Path
import hashlib,sys
import numpy as np
from scipy.optimize import least_squares
ROOT=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT/'scripts'));sys.path.insert(0,str(ROOT/'code'))
import run_he0435_native_completion as run
import strong_lensing_response_completion as rc
from strong_lensing_native_completion import NativeCoupledScene
prior=run.prior;js=prior.js


def main():
    out=run.OUT/'held_target';out.mkdir(exist_ok=True)
    trained,_,records,wcss,origin,ipc=run.build('training')
    saved=np.load(run.OUT/'training/solution.npz');tp=saved['parameters'];v=saved['joint_covariance']
    fulltarget=js.SceneProblem(records);n=fulltarget.size;d=trained.d
    x=np.zeros(n+d+25);byname={r['name']:i for i,r in enumerate(trained.target.records)}
    e4=sorted(wcss)[-1];held=[r for r in records if r['exposure']==e4]
    for i,record in enumerate(records):
        x[3*i:3*i+3]=tp[3*byname[record['name']]:3*byname[record['name']]+3] if record['name'] in byname else np.r_[record['guess'],np.log(record['flux_guess'])]
    x[fulltarget.sky_start:fulltarget.bg_start]=trained.target.sky(tp[:trained.n])
    x[fulltarget.bg_start:fulltarget.bg_start+9]=tp[trained.target.bg_start:trained.n]
    x[fulltarget.bg_start+9]=np.median(np.concatenate([r['data'] for r in held]))
    x[n:n+d]=tp[trained.eta_start:trained.field_start];x[n+d:-1]=tp[trained.field_start:-1];x[-1]=tp[-1]
    prior.set_eta(records,x[n:n+d])
    for record in trained.fields:record['response'].update(x[n:n+d])
    histories={ctx['original']['exposure']:ctx['history'] for ctx in trained.contexts}
    cal=trained.contexts[0]['calibration']
    problem=NativeCoupledScene(fulltarget,trained.fields,x[n:n+d],histories,wcss,origin,cal.g,cal.h,x)
    # Preserve the original fulltarget object (and hence original E1 anchors),
    # but evaluate only E4 target rows. Frozen priors add constants only.
    problem.contexts=[ctx for ctx in problem.contexts if not ctx['field'] and ctx['original']['exposure']==e4]
    indices=np.r_[np.arange(36,48),np.arange(fulltarget.bg_start+9,fulltarget.bg_start+12)]
    lo,hi=problem.bounds();scales=problem.scales[indices]
    def expand(q):
        ans=x.copy();ans[indices]=q;return ans
    def residual(q):return problem.evaluate(expand(q))[0]
    def jacobian(q):return problem.evaluate(expand(q),True)[1][:,indices]
    history=[]
    def callback(intermediate_result):
        q=intermediate_result.x;r=residual(q);a=jacobian(q);pg,_,_=rc.projected_gradient(q,a.T@r,lo[indices],hi[indices],scales)
        history.append({'objective':float(r@r),'KKT':float(abs(pg).max())})
        print('HELD_ITER',len(history),history[-1],flush=True)
    fit=least_squares(residual,x[indices],jac=jacobian,bounds=(lo[indices],hi[indices]),x_scale=scales,
        ftol=1e-14,xtol=1e-14,gtol=1e-9,max_nfev=150,callback=callback)
    finalx=expand(fit.x);r,j=problem.evaluate(finalx,True);a=j[:,indices]
    # Map shared training variables into the full-layout derivative, retaining
    # all eta/sky/kappa crossblocks and indirect trained-point correlations.
    sharedtrain=np.r_[np.arange(trained.target.sky_start,trained.target.bg_start),np.arange(trained.eta_start,trained.field_start),trained.kappa_index]
    sharedfull=np.r_[np.arange(fulltarget.sky_start,fulltarget.bg_start),np.arange(n,n+d),problem.kappa_index]
    b=j[:,sharedfull];vshared=v[np.ix_(sharedtrain,sharedtrain)]
    cond,cut,sens,info=js.propagated_cut_covariance(a,b,vshared,scales=scales)
    train_transform=np.zeros((24,len(tp)));train_take=np.array([3*i+k for i in range(12) for k in [0,1]])
    train_transform[np.arange(24),train_take]=1
    train_transform[:,trained.eta_start:trained.field_start]=np.vstack([r['response'].moment_derivative for r in trained.target.records])
    heldtake=np.array([3*i+k for i in range(4) for k in [0,1]])
    moment=np.zeros((8,len(sharedtrain)));moment[:,36:36+d]=np.vstack([r['response'].moment_derivative for r in held])
    transfer=sens[heldtake]+moment
    traincov=train_transform@v@train_transform.T
    heldcov=cond[np.ix_(heldtake,heldtake)]+transfer@vshared@transfer.T
    cross=train_transform@v[:,sharedtrain]@transfer.T
    fullcov=np.block([[traincov,cross],[cross.T,heldcov]])
    pg,atlo,athi=rc.projected_gradient(fit.x,a.T@r,lo[indices],hi[indices],scales)
    pointrows=sum(len(ctx['original']['data']) for ctx in problem.contexts);rawq=float(r[:pointrows]@r[:pointrows])
    replay=run.replay_native(problem,finalx,out)
    np.savez(out/'solution.npz',held_parameters=fit.x,full_parameters=finalx,held_conditional_covariance=cond,held_cut_covariance=cut,
        shared_training_indices=sharedtrain,shared_transfer=sens,full32_moment_centroid_covariance_pixel2=fullcov)
    prior.dump(out/'result.json',{'scope':'Already-consumed E4 target pixels absent from nonlinear optical/sky/kappa training; all4 FIELD calibrations retained. Held positions, fluxes and exposure plane measured with frozen shared response. Original E1 anchors retained. No unused validation.',
        'raw_original_ERR_Q_continuous_mean':rawq,'projected_KKT_scaled_max':float(abs(pg).max()),'fit_nfev':fit.nfev,'scipy_message':fit.message,'history':history,
        'kappa_e_inverse':float(finalx[-1]),'records':problem.mean_records,'official_native_value_replay':replay,
        'full32_covariance_min_eigenvalue_pixel2':float(np.linalg.eigvalsh(fullcov).min()),'full32_moment_centroid_covariance_pixel2':fullcov,
        'covariance_scope':'Cut propagation of shared trained sky/eta/kappa uncertainty into held nuisances, retaining training crossblocks. Local UNCONSTRAINED curvature of original-ERR conditional estimating problem; not an exact posterior or physical precision validation.',
        'geometry_training_plus_held_baseline_relative_moment':prior.geometry_from_centroids(fulltarget,finalx[:n],fullcov,wcss,origin,'baseline_relative_moment'),
        'physical_validation_pass':False,'physical_admission':False,'physical_guard':2e-10,'source_sha256':hashlib.sha256(Path(__file__).read_bytes()).hexdigest()})
    print('HELD_FINAL',rawq,'KKT',float(abs(pg).max()),'native',replay['target_original_native_ERR_Q'],flush=True)


if __name__=='__main__':main()
