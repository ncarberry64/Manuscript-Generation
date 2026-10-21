"""Actual source-calibrated IPC/linearity ordering on the consumed native ramps."""
from pathlib import Path
import argparse,hashlib,json,sys
import numpy as np
from astropy.io import fits
ROOT=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT/'code'));sys.path.insert(0,str(ROOT/'scripts'))
from strong_lensing_ipc_response import inverse_deconvolution_kernel,convolution_matrix
from strong_lensing_native_calibration import CalibrationOrderedHistory
from apply_he0435_bounded_native_history import native_apply,dump,BASE,OUT
from stdatamodels.jwst import datamodels
from stdatamodels.jwst.datamodels import dqflags
from stcal.ramp_fitting import ramp_fit


def main():
    parser=argparse.ArgumentParser();parser.add_argument('--kappa',type=float,default=0.);args=parser.parse_args()
    dest=OUT/'calibration_order'/f'k{args.kappa:g}';dest.mkdir(parents=True,exist_ok=True)
    gkernel=fits.getdata(BASE/'native_cache/jwst_nircam_ipc_0022.fits',1)
    kernel=inverse_deconvolution_kernel(gkernel,order=8)
    g=convolution_matrix(gkernel,(65,65));forward=convolution_matrix(kernel['kernel'],(65,65))
    rows=[];failures=[]
    for pilotpath in sorted((OUT/'common_early/pilot_means').glob('training_plus_held_jw*_pilot.npz')):
        stem=pilotpath.name.removeprefix('training_plus_held_').removesuffix('_pilot.npz')
        with np.load(BASE/'histories'/f'{stem}_histories.npz') as source:h={key:source[key] for key in source.files}
        with np.load(pilotpath) as source:pilot={key:source[key] for key in source.files}
        with datamodels.RampModel(BASE/'native_cache'/f'{stem}_uncal.fits') as model:
            template=ramp_fit.create_ramp_fit_class(model,'OLS_C',dqflags.pixel,False)
        for p,label in enumerate(h['labels']):
            roi=np.zeros((65,65),bool);roi[24:41,24:41]=True
            valid=roi&np.isfinite(h['archived_cal_MJy_sr'][p])&np.isfinite(h['archived_cal_ERR'][p])&(h['archived_cal_ERR'][p]>0)&((h['native_rate_DQ'][p]&1)==0)
            # Pilot amplitudes are native DN/s, whose producing estimator uses
            # rounded TGROUP. Convert to physical arrivals on unrounded times.
            clock_scale=float(h['group_time_header_s'])/(5*float(h['frame_time_s']))
            effective=pilot['point_plus_sky_e_s'][p]*clock_scale
            for cap in [1.,1.10]:
                tag=f'{stem}_{label}_cap{cap:g}'.replace('.','p')
                try:
                    history=CalibrationOrderedHistory(h,p,g,forward,effective,ceiling_multiplier=cap)
                    state=history.evaluate(effective,args.kappa)
                    groups=state['groups_DN'];zero=state['zeroframe_DN']
                    native=native_apply(template,h,p,groups,zero)
                    mean=native*float(h['photmjsr'])/h['flat'][p]+pilot['detector_plane_offset_SCI'][p]
                    residual=mean[valid]-h['archived_cal_MJy_sr'][p][valid]
                    if not np.isfinite(residual).all():raise ValueError('Nonfinite selected native mean')
                    delta=state['rate_preipc_e_s']-state['uncontinued_rate_preipc_e_s']
                    effective_change=history.gain*(forward@(delta/history.gain))
                    row={'exposure':stem,'label':str(label),'kappa_e_inverse':args.kappa,'ceiling_multiplier':cap,'pixels':int(valid.sum()),
                        'original_native_ERR_chi2':float(np.sum((residual/h['archived_cal_ERR'][p][valid])**2)),
                        'raw_RMS_SCI':float(np.sqrt(np.mean(residual**2))),
                        'max_harmonic_preIPC_context_effect_on_inner17_lowcharge_e_s':float(np.max(abs(effective_change.reshape((65,65))[roi]),initial=0.)),
                        'reported_native_to_physical_rate_clock_factor':clock_scale,
                        **state['diagnostics']}
                    np.savez_compressed(dest/f'{tag}.npz',native_groups_DN=groups,native_zeroframe_DN=zero,
                         native_model_rate_DN_s=native,mean_SCI=mean,residual_SCI=residual,valid=valid,
                         preIPC_rate_e_s=state['rate_preipc_e_s'].reshape(65,65),
                         continued_context_mask=~history.rate_context.known.reshape(65,65),
                         beta_preIPC_DN=state['beta_preipc_DN'].reshape(65,65))
                    rows.append(row);print(stem,label,cap,row['original_native_ERR_chi2'],flush=True)
                except (ValueError,RuntimeError,ArithmeticError) as exc:
                    row={'exposure':stem,'label':str(label),'ceiling_multiplier':cap,'error':str(exc)}
                    failures.append(row);print('DOMAIN_FAILURE',row,flush=True)
            dump(dest/'progress.json',{'rows':rows,'domain_failures':failures})
    aggregate=[]
    for cap in [1.,1.1]:
        selected=[r for r in rows if r['ceiling_multiplier']==cap]
        target=[r for r in selected if r['label']!='FIELD']
        aggregate.append({'ceiling_multiplier':cap,'completed_targets':len(target),'completed_fields':len(selected)-len(target),
            'target_original_native_ERR_Q':sum(r['original_native_ERR_chi2'] for r in target),
            'held_E4_original_native_ERR_Q':sum(r['original_native_ERR_chi2'] for r in target if r['exposure'].endswith('00004_nrcb3')),
            'field_original_native_ERR_Q':sum(r['original_native_ERR_chi2'] for r in selected if r['label']=='FIELD')})
    dump(dest/'result.json',{'scope':'Actual nominal source-calibrated IPC/linearity composition, converged early optical pilots and original full native ERR. Fixed observed selection, conditional caps/electronic intercept/context and photon partition; no model refit.',
        'source_sha256':{str(p.relative_to(ROOT)):hashlib.sha256(p.read_bytes()).hexdigest() for p in [Path(__file__),ROOT/'code/strong_lensing_native_calibration.py',ROOT/'code/strong_lensing_calibration_order.py',ROOT/'code/strong_lensing_ipc_response.py']},
        'inverse_kernel_tail_l1_bound':kernel['tail_l1_bound'],'reference_renormalized':False,
        'rows':rows,'domain_failures':failures,'aggregate':aggregate,'error_rescaling':False,
        'unused_observations_opened':False,'physical_admission':False,'physical_validation_pass':False,'physical_guard':2e-10})


if __name__=='__main__':main()
