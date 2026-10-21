"""Frozen-body plane-as-carrier counterfactual; no fit, new PSF or ERR change."""
from pathlib import Path
import hashlib,sys
import numpy as np
ROOT=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT/'scripts'))
import run_he0435_native_completion as run
from diagnose_he0435_temporal_curvature import fit_curves
from apply_he0435_bounded_native_history import native_apply
from stdatamodels.jwst import datamodels
from stdatamodels.jwst.datamodels import dqflags
from stcal.ramp_fitting import ramp_fit


def full_plane(problem,x,ctx):
    rec=ctx['original'];xx=ctx['x'];yy=ctx['y']
    if ctx['field']:
        off=problem.field_start+6*ctx['field_index']+3
        return np.c_[np.ones(len(xx)),xx-rec['guess'][0],yy-rec['guess'][1]]@x[off:off+3]
    center=np.mean([r['guess'] for r in problem.target.records if r['exposure']==rec['exposure']],axis=0)
    off=problem.target.bg_start+3*problem.target.eindex[rec['exposure']]
    return np.c_[np.ones(len(xx)),(xx-center[0])/100,(yy-center[1])/100]@x[off:off+3]


def main():
    dest=run.OUT/'background_partition';dest.mkdir(exist_ok=True)
    problem,_,_,_,_,_=run.build('all_data');original_contexts=problem.contexts;rows=[];failures=[];arrays={};templates={};inputs={}
    for label in ['all_data','held_target']:
        solution=run.OUT/label/'solution.npz';inputs[str(solution.relative_to(ROOT))]=hashlib.sha256(solution.read_bytes()).hexdigest()
        x=np.load(solution)['parameters' if label=='all_data' else 'full_parameters'];e4=problem.exposure_order[-1]
        problem.contexts=original_contexts if label=='all_data' else [c for c in original_contexts if not c['field'] and c['original']['exposure']==e4]
        problem.last_x=None;problem.evaluate(x);savedmaps=np.load(run.OUT/label/'native_replay_maps.npz')
        for ctx,(baseline,offset,oldweights),record in zip(problem.contexts,problem.last_states,problem.mean_records):
            rec=ctx['original'];h=ctx['history'];p=ctx['patch'];sel=ctx['selected'];cal=ctx['calibration'];name=ctx['name'];e=rec['exposure']
            stem=e.replace('_cal.fits','');source=problem.optical(x,ctx)[0];plane=full_plane(problem,x,ctx)
            baseline_replay=float(np.max(abs(baseline['groups_DN']-savedmaps[name+'_native_groups_DN'])))
            if e not in templates:
                with datamodels.RampModel(run.BASE/'native_cache'/f'{stem}_uncal.fits') as model:
                    templates[e]=ramp_fit.create_ramp_fit_class(model,'OLS_C',dqflags.pixel,False)
            data=h['corrected_group_DN'][p].astype(float);valid=h['valid_group_DNU_SAT_only'][p].copy()
            valid[:,np.any((h['groupdq'][p]&4)!=0,axis=0)]=False
            pilot=h['native_rate_DN_s'][p].astype(float)*float(h['group_time_header_s'])/(5*float(h['frame_time_s']))
            gain=h['gain_e_per_DN'][p].astype(float);noise=h['readnoise_CDS_DN'][p].astype(float)
            cx,cy=h['prior_guess_native_xy'][p];radius=np.hypot(h['x_pixel'][p]-cx,h['y_pixel'][p]-cy)
            roi=(abs(h['x_pixel'][p]-round(cx))<=8)&(abs(h['y_pixel'][p]-round(cy))<=8)&((h['archived_cal_DQ'][p]&1)==0)
            for rho in [0.,1.]:
                try:
                    rate=source+rho*plane*ctx['conversion_e_per_s_per_SCI']
                    uncontinued=cal.gain*(cal.g@(rate/cal.gain));trusted=cal.rate_context.known
                    state=baseline if rho==0 else cal.evaluate(rate,x[-1])
                    value=native_apply(templates[e],h,p,state['groups_DN'],state['zeroframe_DN']).ravel()[sel]*ctx['native_SCI_factor']+(1-rho)*offset
                    residual=value-rec['data'];z=residual/rec['error'];_,weights,branches=problem.native_ols(ctx,state)
                    coef,_,_,linearq=fit_curves(data-state['groups_DN'],valid,pilot,gain,noise,h['native_read_times_s'])
                    use=roi&valid.all(axis=0)&(radius>=3)&(radius<6)&np.isfinite(coef[2])
                    row={'body':label,'name':name,'field':ctx['field'],'photon_plane_fraction':rho,'kappa_e_inverse':float(x[-1]),
                        'original_native_ERR_Q':float(z@z),'raw_RMS_SCI':float(np.sqrt(np.mean(residual**2))),
                        'trusted_preIPC_min_rate_e_s':float(uncontinued[trusted].min()),
                        'beta_recomputed_by_original_G0_equation':True,'max_abs_native_group_slope_Jacobian_change_per_s':float(np.max(abs(weights-oldweights))),
                        'baseline_full65_group_replay_max_abs_DN':baseline_replay,
                        'near_3_to_6_pixel_all9valid_pixels':int(use.sum()),
                        'near_curvature_raw_RMS_DN_s2':float(np.sqrt(np.mean(coef[2,use]**2))) if use.any() else None,
                        'near_temporal_shape_Q_conditional_metric':float(np.sum(linearq[use]))}
                    rows.append(row);tag=f'{label}_{name}_rho{rho:g}'
                    arrays[tag+'_native_groups_DN']=state['groups_DN'];arrays[tag+'_mean_SCI']=value;arrays[tag+'_raw_residual_SCI']=residual;arrays[tag+'_curvature_residual_DN_s2']=coef[2]
                except ValueError as exc:
                    failures.append({'body':label,'name':name,'photon_plane_fraction':rho,'error':str(exc),'no_clipping':True})
    aggregate=[]
    for body in ['all_data','held_target']:
        for rho in [0.,1.]:
            take=[r for r in rows if r['body']==body and r['photon_plane_fraction']==rho]
            aggregate.append({'body':body,'photon_plane_fraction':rho,'completed_patches':len(take),
                'target_original_ERR_Q':sum(r['original_native_ERR_Q'] for r in take if not r['field']),
                'field_original_ERR_Q':sum(r['original_native_ERR_Q'] for r in take if r['field'])})
    np.savez_compressed(dest/'counterfactual_maps.npz',**arrays)
    run.prior.dump(dest/'result.json',{'scope':'Frozen converged physical optical/sky/kappa body, alternative partition of the same fitted native plane into photon carriers versus late electronics. No new data, refit, noise inflation or new PSF family.',
        'equation':'lambda_effective(rho)=lambda_point+sky+rho*plane*physical_clock_conversion; output=native_mean(lambda_effective(rho))+(1-rho)*plane, rho=0 or1. For uniform background, L[(S+b)^2]=L[S^2]+2bL[S].',
        'boundary_scope':'Same fixed trusted preIPC domains and harmonic exterior; negative trusted carrier inputs rejected individually. Beta re-solved using original actualG0 boundary relation; exact native branch-dependent estimator applied.',
        'temporal_noise_scope':'Original conditional Poisson/CDS comparator for group-shape Q; raw curvature RMS separate. Not a completed anchored group likelihood.',
        'rows':rows,'domain_failures':failures,'aggregate':aggregate,'input_sha256':inputs,'source_sha256':hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
        'physical_admission':False,'physical_validation_pass':False,'physical_guard':2e-10})
    print(aggregate,flush=True);print('FAILURES',failures,flush=True)


if __name__=='__main__':main()
