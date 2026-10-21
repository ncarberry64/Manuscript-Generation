"""Carry a converged early optical mean through actual native full histories.

Conditional storage/readout ceiling trials are evaluated as failed or improved
raw means under the ORIGINAL native ERR. This is not covariance calibration,
new observations, or an optimization against physical lens positions.
Run in .venv-native-history, which contains the documented native OLS_C.
"""
from pathlib import Path
import argparse, copy, hashlib, json, sys
import numpy as np
ROOT=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT/'code'))
from strong_lensing_bounded_collection import bounded_leading_collection
from strong_lensing_stochastic_response import nearest_neighbor_edges, harmonic_rate_initializer
from strong_lensing_group_calibration import polynomial_derivatives, averaged_read_observation
from stdatamodels.jwst import datamodels
from stdatamodels.jwst.datamodels import dqflags
from stcal.ramp_fitting import ramp_fit

BASE=ROOT/'data/strong_lensing_response_completion_20261021'
OUT=ROOT/'artifacts/strong_lensing_response_completion_20261021'


def sha(path):return hashlib.sha256(path.read_bytes()).hexdigest()


def dump(path,value):
    path.write_text(json.dumps(value,indent=2,allow_nan=False)+'\n',encoding='utf-8',newline='\n')


def native_apply(template,h,p,groups,zero):
    local=copy.copy(template)
    local.set_arrays(np.array(groups[None],dtype=np.float32,copy=True),
        np.array(h['groupdq'][p,None],copy=True),np.array(h['pixeldq'][p],copy=True),
        np.array(h['average_dark_current_DN_s'][p],copy=True),
        zeroframe=np.array(zero[None],dtype=np.float32,copy=True))
    local.start_row=0;local.num_rows=groups.shape[-2]
    local.zframe_mat=None;local.zframe_locs=None;local.zframe_cnt=0
    local.one_groups_locs=None;local.one_groups_time=None
    image,_,_=ramp_fit.ramp_fit_data(local,False,
        np.array(h['readnoise_CDS_DN'][p],copy=True),np.array(h['gain_e_per_DN'][p],copy=True),
        'OLS_C','optimal','none')
    return np.array(image['slope'],copy=True)


def patch_mean(h,p,rate,kappa,prescription,ceiling_multiplier):
    gain=h['gain_e_per_DN'][p].astype(float);shape=gain.shape
    # Outside the measured optical stamp the finite response can be undefined
    # or slightly negative through spline interpolation. Keep the explicit
    # same-context harmonic continuation receipt rather than silently clip.
    rate,continuation=harmonic_rate_initializer(rate,np.isfinite(rate)&(rate>=0))
    coeff=h['linearity_coefficients'][p].copy()
    flags=h['pixeldq'][p]
    identity=((flags&1048576)!=0)|~np.isfinite(coeff).all(axis=0)
    coeff[:,identity]=0.;coeff[1,identity]=1.
    times=h['native_read_times_s'].astype(float)
    dark=h['dark_RAPID42_DN'][p][h['native_dark_RAPID_frame_indices']].astype(float)
    dark_group=h['native_subtracted_dark_group_DN'][p].astype(float)
    anchors=h['pre_linearity_group_DN'][p].astype(float)
    observed=h['corrected_group_DN'][p].astype(float)
    # An electronic intercept conditional on actual group0 is a nuisance
    # boundary datum, not an assumed physical reset charge or an independent
    # held group. It is re-applied before inverse-linearity/read averaging.
    beta=observed[0]-rate/gain*times[0].mean()
    offset=h['raw_group_DN'][p,0].astype(float)-anchors[0]
    threshold_adc=h['saturation_threshold_DN'][p].astype(float)-offset
    ceiling_linear=polynomial_derivatives(coeff,threshold_adc)[0]
    caps=gain*(ceiling_linear-beta-dark[0].mean(axis=0))*ceiling_multiplier
    caps,cap_continuation=harmonic_rate_initializer(caps,np.isfinite(caps)&(caps>0))
    cap_continuation['meaning']='Provisional harmonic positive-ceiling continuation where the reference/translation does not define a positive value; not measured full well.'
    # ZEROFRAME precedes group0's second read. Preserve its actual native
    # branch/absence flags. It is not an extra independent observation.
    alltimes=np.r_[float(h['frame_time_s']),times.ravel()]
    order=np.argsort(alltimes,kind='stable');inverse=np.argsort(order)
    physical=bounded_leading_collection(rate.ravel(),alltimes[order],kappa,caps.ravel(),
         nearest_neighbor_edges(shape),prescription=prescription,
         ceiling_source=f'{ceiling_multiplier:g} times the initial-intercept group0-translated SAT calibration cutoff, frozen before intercept refinement; conditional trial, not measured full well')
    means=physical['observed_read_mean_e'][inverse].reshape((-1,*shape))
    read=means[1:].reshape((9,2,*shape))/gain
    used=((h['groupdq'][p]&3)==0)&((flags[None]&1)==0)
    # Only groups used by the native estimator need an inverse. Retain native
    # values on excluded groups; their means never contribute to the score.
    anchor_q=polynomial_derivatives(coeff,anchors)[0]
    anchor_reads=np.repeat(anchors[:,None],2,axis=1)
    for _ in range(5):
        q=read+dark+beta[None,None]
        q=np.where(used[:,None],q,anchor_q[:,None])
        mapped=averaged_read_observation(q,coeff,anchor_reads,dark_group_dn=dark_group)
        correction=np.where(used[0],observed[0]-mapped['mean_DN'][0],0.)
        beta+=correction/np.sum(mapped['read_to_group_jacobian'][0],axis=0)
        if np.max(abs(correction[used[0]]),initial=0.)<1e-7:break
    q=read+dark+beta[None,None]
    q=np.where(used[:,None],q,anchor_q[:,None])
    mapped=averaged_read_observation(q,coeff,anchor_reads,dark_group_dn=dark_group)
    groups=np.where(used,mapped['mean_DN'],observed)
    # Native DarkCurrentStep subtracts group darks, but not ZEROFRAME dark.
    zero=beta+means[0]/gain+dark[0,0]
    zero=np.where(h['native_corrected_zeroframe_DN'][p]!=0,zero,0.)
    return groups,zero,{'context_continuation':continuation,
       'ceiling_context_continuation':cap_continuation,
       'identity_linearity_pixels':int(identity.sum()),
       'intercept_boundary':'Actual group0 conditional per-pixel electronic intercept; no independent prediction of group0.',
       'group0_intercept_match_max_abs_DN':float(np.max(abs(groups[0][used[0]]-observed[0][used[0]]),initial=0.)),
       'ceiling_multiplier':ceiling_multiplier,'ceiling_not_physical_full_well':True,
       'ceiling_intercept_convention':'Frozen from initial group0 intercept, before nonlinear intercept refinement.',
       'first_order_collection_drift_sum_max_abs_e2':float(np.max(abs(np.array(physical['diagnostics']['sum_internal_first_order_drift_e2'])))),
       'conditional_lost_charge_last_read_e':float(physical['lost_mean_e'][-1]),
       'observed_group_selection_held_fixed':True}


def main():
    parser=argparse.ArgumentParser();parser.add_argument('--pilot-label',default='training_plus_held',choices=['all_data','training_plus_held']);args=parser.parse_args()
    dest=OUT/'bounded_native'/args.pilot_label;dest.mkdir(parents=True,exist_ok=True)
    pilots=sorted((OUT/'common_early/pilot_means').glob(args.pilot_label+'_jw*_pilot.npz'))
    if len(pilots)!=4:raise RuntimeError('Four matching converged early pilot files required')
    rows=[];failures=[];contexts={}
    for pilotpath in pilots:
        stem=pilotpath.name.removeprefix(args.pilot_label+'_').removesuffix('_pilot.npz')
        hpath=BASE/'histories'/f'{stem}_histories.npz'
        with np.load(hpath) as source:h={key:source[key] for key in source.files}
        with np.load(pilotpath) as source:pilot={key:source[key] for key in source.files}
        raw=BASE/'native_cache'/f'{stem}_uncal.fits'
        with datamodels.RampModel(raw) as model:
            template=ramp_fit.create_ramp_fit_class(model,'OLS_C',dqflags.pixel,False)
        for p,label in enumerate(h['labels']):
            roi=np.zeros(h['gain_e_per_DN'][p].shape,bool);roi[24:41,24:41]=True
            valid=roi&np.isfinite(h['archived_cal_MJy_sr'][p])&np.isfinite(h['archived_cal_ERR'][p])&(h['archived_cal_ERR'][p]>0)&((h['native_rate_DQ'][p]&1)==0)
            observed=h['archived_cal_MJy_sr'][p][valid];error=h['archived_cal_ERR'][p][valid]
            rate=pilot['point_plus_sky_e_s'][p]
            offset=pilot['detector_plane_offset_SCI'][p]
            canonical={}
            for cap in [1.,1.10]:
                for kind in ['storage_ceiling','readout_ceiling']:
                    for kappa in [0.,1e-7,2.3665e-7]:
                        tag=f'{stem}_{label}_{kind}_cap{cap:g}_k{kappa:g}'.replace('.','p')
                        try:
                            groups,zero,receipt=patch_mean(h,p,rate,kappa,kind,cap)
                            native=native_apply(template,h,p,groups,zero)
                            mean=native*float(h['photmjsr'])/h['flat'][p]+offset
                            residual=mean[valid]-observed
                            context_key=f'{stem}_{label}_cap{cap:g}'
                            contexts[context_key]={key:receipt.pop(key) for key in ['context_continuation','ceiling_context_continuation']}
                            row={'exposure':stem,'label':str(label),'prescription':kind,'ceiling_multiplier':cap,'kappa_e_inverse':kappa,
                                'pixels':int(valid.sum()),'raw_RMS_SCI':float(np.sqrt(np.mean(residual**2))),
                                'original_native_ERR_chi2':float(np.sum((residual/error)**2)),
                                'context_receipt_key':context_key,**receipt}
                            if not np.isfinite(residual).all():raise ValueError('Nonfinite native selected mean')
                            if cap==1. and kind=='storage_ceiling':
                                np.savez_compressed(dest/f'{tag}.npz',native_groups_DN=groups,native_zeroframe_DN=zero,
                                    native_model_rate_DN_s=native,mean_SCI=mean,valid=valid,residual_SCI=residual,
                                    original_native_ERR_SCI=error)
                                canonical[kappa]=(tag,groups,zero,native)
                            reference=canonical[kappa]
                            row['canonical_array_file']=reference[0]+'.npz'
                            row['max_abs_group_mean_difference_from_canonical_DN']=float(np.max(abs(groups-reference[1])))
                            row['float32_native_input_equal_canonical']=bool(np.array_equal(groups.astype(np.float32),reference[1].astype(np.float32)) and np.array_equal(zero.astype(np.float32),reference[2].astype(np.float32)))
                            row['native_rate_equal_canonical']=bool(np.array_equal(native,reference[3],equal_nan=True))
                            rows.append(row)
                            print(stem,label,kind,cap,kappa,row['original_native_ERR_chi2'],flush=True)
                        except (ValueError,RuntimeError,ArithmeticError) as exc:
                            failure={'exposure':stem,'label':str(label),'prescription':kind,'ceiling_multiplier':cap,'kappa_e_inverse':kappa,'error':str(exc)}
                            failures.append(failure);print('DOMAIN_FAILURE',failure,flush=True)
            dump(dest/'progress.json',{'rows':rows,'domain_failures':failures})
            dump(dest/'context_receipts.json',contexts)
    aggregate=[]
    for cap in [1.,1.10]:
        for kind in ['storage_ceiling','readout_ceiling']:
            for k in [0.,1e-7,2.3665e-7]:
                selected=[r for r in rows if r['ceiling_multiplier']==cap and r['prescription']==kind and r['kappa_e_inverse']==k]
                target=[r for r in selected if r['label']!='FIELD'];held=[r for r in target if r['exposure'].endswith('00004_nrcb3')]
                aggregate.append({'prescription':kind,'ceiling_multiplier':cap,'kappa_e_inverse':k,'completed_target_patches':len(target),
                    'target_original_ERR_Q':sum(r['original_native_ERR_chi2'] for r in target),
                    'held_E4_original_ERR_Q':sum(r['original_native_ERR_chi2'] for r in held),
                    'field_original_ERR_Q':sum(r['original_native_ERR_chi2'] for r in selected if r['label']=='FIELD')})
    dump(dest/'result.json',{'scope':'Conditional first-order native-history forward applications of converged common-early optical means, original full native ERR unchanged. No refit, physical validation, or independent held-group prediction.',
        'photon_partition_assumption':'Positive point+coherent sky assigned to photon arrival rate; detector plane remains an additive electronic rate offset. This partition is not identified by broad-band rates.',
        'source_sha256':{str(p.relative_to(ROOT)):sha(p) for p in [Path(__file__),ROOT/'code/strong_lensing_bounded_collection.py',ROOT/'code/strong_lensing_group_calibration.py']},
        'pilots':{str(p.relative_to(ROOT)):sha(p) for p in pilots},'rows':rows,'domain_failures':failures,'aggregate':aggregate,
        'error_rescaling':False,'unused_observations_opened':False,'physical_validation_pass':False,'physical_admission':False,'physical_guard':2e-10})


if __name__=='__main__':main()
