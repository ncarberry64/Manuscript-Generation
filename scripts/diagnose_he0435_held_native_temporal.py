"""Actual E4 group-shape diagnostic after independent held-point measurement."""
from pathlib import Path
import hashlib,sys
import numpy as np
ROOT=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT/'scripts'))
import run_he0435_native_completion as run
from diagnose_he0435_temporal_curvature import fit_curves


def main():
    dest=run.OUT/'held_target';saved=dest/'solution.npz';x=np.load(saved)['full_parameters']
    problem,_,_,_,_,_=run.build('all_data');e4=sorted({ctx['original']['exposure'] for ctx in problem.contexts})[-1]
    problem.contexts=[ctx for ctx in problem.contexts if not ctx['field'] and ctx['original']['exposure']==e4]
    problem.evaluate(x);rows=[];maps={}
    for ctx,(state,_,_) in zip(problem.contexts,problem.last_states):
        h=ctx['history'];p=ctx['patch'];cal=ctx['calibration'];null=cal.evaluate(problem.optical(x,ctx)[0],0.)['groups_DN']
        data=h['corrected_group_DN'][p].astype(float);valid=h['valid_group_DNU_SAT_only'][p].copy()
        valid[:,np.any((h['groupdq'][p]&4)!=0,axis=0)]=False
        native=h['native_rate_DN_s'][p].astype(float)*float(h['group_time_header_s'])/(5*float(h['frame_time_s']))
        gain=h['gain_e_per_DN'][p].astype(float);noise=h['readnoise_CDS_DN'][p].astype(float)
        xx=h['x_pixel'][p];yy=h['y_pixel'][p];cx,cy=h['prior_guess_native_xy'][p];radius=np.hypot(xx-cx,yy-cy)
        roi=(abs(xx-round(cx))<=8)&(abs(yy-round(cy))<=8)&((h['archived_cal_DQ'][p]&1)==0)
        comparison={}
        for name,values in [('observed',data),('frozen_body_kappa0_residual',data-null),('fitted_history_residual',data-state['groups_DN'])]:
            coef,cov,q,ql=fit_curves(values,valid,native,gain,noise,h['native_read_times_s'])
            use=roi&valid.all(axis=0)&(radius>=3)&(radius<6)&np.isfinite(coef[2])
            comparison[name]={'pixels':int(use.sum()),'curvature_raw_RMS_DN_s2':float(np.sqrt(np.mean(coef[2,use]**2))),
                'temporal_shape_Q_after_fitted_static_intercept_and_rate':float(np.sum(ql[use]))}
            maps[ctx['name']+'_'+name+'_curvature_DN_s2']=coef[2]
        rows.append({'name':ctx['name'],'near_wing_3_to_6_pixel':comparison})
    run.prior.dump(dest/'temporal_residual_diagnosis.json',{'scope':'Consumed E4 groups after independently measured E4 point nuisances and frozen E1--3-trained sky/eta/kappa. Temporal curvature is an internal development diagnostic, not unused validation. Same fitted body kappa0 counterfactual included.',
        'kappa_e_inverse':float(x[-1]),'rows':rows,'noise_scope':'Conditional old Poisson/CDS-marginal comparator only; raw curvature RMS remains separate from formal Q. No error scaling, SIRS addition or claim of complete anchored group likelihood.',
        'fit_solution_sha256':hashlib.sha256(saved.read_bytes()).hexdigest(),'source_sha256':hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
        'physical_validation_pass':False,'physical_admission':False,'physical_guard':2e-10})
    np.savez_compressed(dest/'temporal_curvature_maps.npz',**maps)
    print(rows,flush=True)


if __name__=='__main__':main()
