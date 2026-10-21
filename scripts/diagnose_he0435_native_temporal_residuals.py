"""Saved fitted-history mean versus actual group curvature; no optical refit."""
from pathlib import Path
import argparse,hashlib,sys
import numpy as np
ROOT=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT/'scripts'));sys.path.insert(0,str(ROOT/'code'))
import run_he0435_native_completion as run
from diagnose_he0435_temporal_curvature import fit_curves
prior=run.prior


def main():
    parser=argparse.ArgumentParser();parser.add_argument('--label',choices=['all_data','training'],default='all_data');args=parser.parse_args()
    dest=run.OUT/args.label/'temporal_residuals';dest.mkdir(exist_ok=True)
    saved=np.load(run.OUT/args.label/'solution.npz');x=saved['parameters']
    problem,_,_,_,_,_=run.build(args.label);problem.evaluate(x)
    rows=[];arrays={}
    for ctx,(state,plane,weights) in zip(problem.contexts,problem.last_states):
        history=ctx['history'];p=ctx['patch'];cal=ctx['calibration'];rate=problem.optical(x,ctx)[0]
        null=cal.evaluate(rate,0.)['groups_DN'];model=state['groups_DN'];data=history['corrected_group_DN'][p].astype(float)
        valid=history['valid_group_DNU_SAT_only'][p].copy()
        valid[:,np.any((history['groupdq'][p]&4)!=0,axis=0)]=False
        pilot=history['native_rate_DN_s'][p].astype(float)*float(history['group_time_header_s'])/(5*float(history['frame_time_s']))
        gain=history['gain_e_per_DN'][p].astype(float);noise=history['readnoise_CDS_DN'][p].astype(float);times=history['native_read_times_s']
        xx=history['x_pixel'][p];yy=history['y_pixel'][p];center=history['prior_guess_native_xy'][p]
        radius=np.hypot(xx-center[0],yy-center[1]);roi=(abs(xx-round(center[0]))<=8)&(abs(yy-round(center[1]))<=8)
        full=valid.all(axis=0)&roi&((history['archived_cal_DQ'][p]&1)==0)
        variants={}
        for label,values in [('observed',data),('frozen_body_kappa0_residual',data-null),('fitted_history_residual',data-model)]:
            coef,cov,q,linearq=fit_curves(values,valid,pilot,gain,noise,times);bins=[]
            arrays[ctx['name']+'_'+label+'_curvature_DN_s2']=coef[2]
            for lower,upper in [(0,3),(3,6),(6,8.6)]:
                use=full&(radius>=lower)&(radius<upper)&np.isfinite(coef[2])
                bins.append({'radius_pixel':[lower,upper],'all9valid_pixels':int(use.sum()),
                    'curvature_raw_RMS_DN_s2':float(np.sqrt(np.mean(coef[2,use]**2))) if use.any() else None,
                    'temporal_shape_Q_after_fitted_intercept_and_linear_rate':float(np.nansum(linearq[use])),
                    'quadratic_shape_Q':float(np.nansum(q[use]))})
            variants[label]=bins
        rows.append({'name':ctx['name'],'field':ctx['field'],'kappa_e_inverse':float(x[-1]),'radial_comparisons':variants})
    prior.dump(dest/'result.json',{'scope':'Actual consumed9group temporal trajectories versus fitted whole-native mean and same frozen optical body with kappa=0. No refit, unused test or change of original ERR. Temporal intercept and linear rate profiled separately per pixel.',
        'mean_equations':'positive optical+coherent sky -> G/gain -> bounded leading collection -> implicit G0 intercept -> inverse calibrated L -> H -> mean of actual2read supports -> native L -> native dark subtraction. The kappa0 counterfactual retains the same fitted optical state and boundary.',
        'noise_scope':'Fixed actual-rate pilot Poisson+postcorrection CDS-marginal comparator from earlier trajectory diagnostic; no SIRS added. Conditional temporal Q is not source-complete likelihood significance.',
        'rows':rows,'physical_admission':False,'physical_validation_pass':False,'physical_guard':2e-10,
        'source_sha256':hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),'fit_solution_sha256':hashlib.sha256((run.OUT/args.label/'solution.npz').read_bytes()).hexdigest()})
    np.savez_compressed(dest/'curvature_maps.npz',**arrays)
    for row in rows:print(row['name'],[(v['radius_pixel'],v['curvature_raw_RMS_DN_s2']) for v in row['radial_comparisons']['fitted_history_residual']],flush=True)


if __name__=='__main__':main()
