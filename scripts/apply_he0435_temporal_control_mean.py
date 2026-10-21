"""One frozen outside-ROI covariance-conditioned E4A noise-mean experiment."""
from pathlib import Path
import hashlib,json,sys,time
import numpy as np
from scipy import sparse,linalg
from astropy.io import fits
ROOT=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT/'code'))
from strong_lensing_dark_covariance import BinnedSpectralModel,native_to_clock_coordinates
from strong_lensing_temporal_controls import affine_free_temporal_contrasts,condition_noise_mean
from strong_lensing_ipc_response import inverse_deconvolution_kernel,convolution_matrix
from strong_lensing_native_calibration import CalibrationOrderedHistory

BASE=ROOT/'data/strong_lensing_response_completion_20261021'
ART=ROOT/'artifacts/strong_lensing_response_completion_20261021'
STEM='jw01198007001_02101_00004_nrcb3'


def sha(p):return hashlib.sha256(Path(p).read_bytes()).hexdigest()


def dump(path,value):
    def default(x):
        if isinstance(x,np.ndarray):return x.tolist()
        if isinstance(x,np.generic):return x.item()
        raise TypeError(type(x).__name__)
    path.write_text(json.dumps(value,indent=2,allow_nan=False,default=default)+'\n',encoding='utf-8',newline='\n')


def nominal_control_photon_covariance(cal,state,control,contrast):
    """Fixed-beta independent Poisson source, actual inverseL/H, no ADC scaling.

    This is a separate source-count contamination term in the control channel,
    not another electronic variance or a complete target counting likelihood.
    It is admitted only on unsaturated H-support in this nominal state.
    """
    hs=cal.h[control];support=np.unique(hs.indices)
    active=state['physical']['clipped_read_mask'][:,support]
    if active.any():return None,{'used':False,'reason':'Declared cap active on control photon-transfer support','clipped_read_support_pairs':int(active.sum())}
    rate=state['rate_preipc_e_s']
    ng,nc=len(contrast[0]),len(control);nl=len(contrast)
    maps=[];times=[]
    for j in range(ng):
        expand=sparse.coo_matrix((contrast[:,j][None,:].repeat(nc,axis=0).ravel(),
            (np.arange(nc*nl),np.repeat(np.arange(nc),nl))),shape=(nc*nl,nc)).tocsr()
        for r in range(2):
            inv=state['observation']['inverse_linearity_derivative'][j,r]
            maps.append(expand@hs@sparse.diags(inv/(2*cal.gain)))
            times.append(cal.times[j,r])
    times=np.asarray(times);unique=np.unique(times)
    answer=np.zeros((nc*nl,)*2);previous=0.
    for t in unique:
        remaining=sum((m for m,s in zip(maps,times) if s>=t),sparse.csr_matrix((nc*nl,cal.n)))
        weighted=remaining@sparse.diags(rate*(t-previous))
        answer+=(weighted@remaining.T).toarray();previous=t
    return (answer+answer.T)/2,{'used':True,'law':'Independent incident Poisson electrons at fixed positive preIPC pilot rates; CovN(t,s)=diag(lambda min(t,s)); actual mean inverseL and H before two-read group average; fixed beta and SIRS mean/operator; finite context.',
        'source_rate_e_s_range_on_transfer_support':[float(rate[support].min()),float(rate[support].max())],
        'context_pixels_in_transfer_support':len(support),'additional_limits':'Does not define unknown state-dependent collection/full-well noise or photon-driven changes of SIRS control coefficients. Supplies control-channel shot contamination only, not a completed target covariance.'}


def main():
    started=time.perf_counter();dest=ART/'dark_calibration/E4A_temporal_control_mean';dest.mkdir(parents=True,exist_ok=True)
    protocolpath=ROOT/'preregistration/he0435_E4A_temporal_control_20261021.json'
    protocol=json.loads(protocolpath.read_text())
    inputs=[protocolpath,BASE/'histories'/f'{STEM}_histories.npz',
        ART/'dark_calibration/nonlocal_G0_E4A/nonlocal_READ_DARK_component.npz',
        BASE/'dark_calibration/earliest_epoch_binned_cross_spectrum_512.npz',
        ART/'common_early/pilot_means'/f'training_plus_held_{STEM}_pilot.npz',
        BASE/'native_cache/jwst_nircam_ipc_0022.fits']
    with np.load(inputs[1]) as f:h={k:f[k] for k in f.files}
    with np.load(inputs[2]) as f:anchor={k:f[k] for k in f.files}
    offsets=np.asarray(protocol['offsets_xy'],int);points=offsets+np.asarray(protocol['center_context_xy'],int)
    if np.any(np.all(abs(offsets)<=8,axis=1)):raise ValueError('Frozen control is inside target ROI')
    original=points[:,1]*65+points[:,0]
    flags=h['groupdq'][0].reshape(9,-1)[:,original]
    values=h['pre_linearity_group_DN'][0].reshape(9,-1)[:,original]
    accepted=np.all((flags&7)==0,axis=0)&np.all(np.isfinite(values),axis=0)
    control=original[accepted];values=values[:,accepted]
    if len(control)<1:raise ValueError('No eligible frozen controls; no replacement is authorized')
    contrast=affine_free_temporal_contrasts(np.asarray(h['native_read_times_s']).mean(axis=1))
    z=(contrast@values).T.ravel()
    native=np.c_[h['x_pixel'][0].ravel(),h['y_pixel'][0].ravel()]
    selected=anchor['selected_context_index'];context=anchor['anchor_context_index'];nc=len(control);nl=len(contrast);m=len(selected)
    pointsjoint=np.vstack((native[selected],native[context],np.repeat(native[control],nl,axis=0)))
    amp,clock=native_to_clock_coordinates(pointsjoint)
    weights=np.r_[anchor['observed_native_ADC_to_SCI_weights'],np.c_[np.ones(len(context)),np.zeros((len(context),8))],np.tile(contrast,(nc,1))]
    with np.load(inputs[3]) as f:model=BinnedSpectralModel(f['spectral_matrices'],f['y_frequency_bin'],f['x_frequency_bin'],f['spatial_shape'],int(f['namp']),int(f['ngroup']))
    print('FIXED_CONTROL_DOMAIN',nc,'contrasts',len(z),'target',m,'joint',len(amp),flush=True)
    joint=model.project([{'amp':amp,'clock_xy':clock,'group_weights':weights}])[0]
    count=m+len(context);targettransform=np.c_[np.eye(m),-anchor['nonlocal_prediction_G0_ADC_to_SCI']]
    ctt=targettransform@joint[:count,:count]@targettransform.T;ctt=(ctt+ctt.T)/2
    replay=float(np.max(abs(ctt-anchor['residual_READ_DARK_covariance_SCI2'])))
    if replay>1e-12:raise ValueError('Target electronic contraction no longer matches preserved source component')
    ctc=targettransform@joint[:count,count:];ccc=joint[count:,count:]
    mean,conditional,transfer=condition_noise_mean(ctt,ctc,ccc,z)
    before=anchor['unchanged_raw_residual_SCI'];error=anchor['unchanged_original_native_ERR'];after=before-mean
    with np.load(inputs[4]) as f:pilot={k:f[k] for k in f.files}
    gkernel=fits.getdata(inputs[5],1);inverse=inverse_deconvolution_kernel(gkernel,order=8)
    g=convolution_matrix(gkernel,(65,65));forward=convolution_matrix(inverse['kernel'],(65,65))
    effective=pilot['point_plus_sky_e_s'][0]*(float(h['group_time_header_s'])/(5*float(h['frame_time_s'])))
    cal=CalibrationOrderedHistory(h,0,g,forward,effective,ceiling_multiplier=1.)
    state=cal.evaluate(effective,0.,derivatives=True)
    curvature=(contrast@state['observation']['coupled_ADC_group_DN'][:,control]).T.ravel()
    curvature_mean=transfer@curvature # Diagnostic only; never subtracted/tuned.
    photon,photonmeta=nominal_control_photon_covariance(cal,state,control,contrast)
    variants={}
    arrays={'frozen_control_offsets_xy':offsets,'fixed_controls_accepted':accepted,'accepted_control_context_index':control,
        'accepted_control_native_xy':native[control],'temporal_contrast_matrix':contrast,
        'observed_control_contrasts_DN':z,'nominal_constant_light_control_curvature_DN':curvature,
        'target_READ_DARK_covariance_SCI2':ctt,'target_control_READ_DARK_cross_covariance':ctc,
        'control_READ_DARK_covariance_DN2':ccc,'conditional_noise_mean_SCI':mean,
        'conditional_READ_DARK_covariance_SCI2':conditional,'unchanged_original_native_ERR':error,
        'original_raw_residual_SCI':before,'noise_mean_corrected_raw_residual_SCI':after,
        'predicted_control_curvature_contamination_in_noise_mean_SCI':curvature_mean}
    if photon is not None:
        pm,pc,pt=condition_noise_mean(ctt,ctc,ccc+photon,z)
        arrays.update(control_nominal_Poisson_covariance_DN2=photon,shot_aware_electronic_noise_mean_SCI=pm,
            shot_aware_linear_prediction_error_READ_DARK_covariance_SCI2=pc,shot_aware_raw_residual_SCI=before-pm)
        variants['nominal_control_Poisson_contamination']={
            'interpretation':'Best linear prediction of the electronic target residual and its error covariance from electronics plus independent Poisson control moments; not generally the exact conditional expectation after adding non-Gaussian Poisson noise.',
            'raw_RMS_SCI':float(np.sqrt(np.mean((before-pm)**2))),
            'original_ERR_Q':float(np.sum(((before-pm)/error)**2)),
            'conditional_mean_RMS_SCI':float(np.sqrt(np.mean(pm**2))),
            'control_Poisson_variance_DN2_range':[float(np.diag(photon).min()),float(np.diag(photon).max())],
            'full_target_counting_covariance_supplied':False}
    output=dest/'actual_conditional_noise_mean.npz';np.savez_compressed(output,**arrays)
    report={
        'scope':'One frozen outside-ROI temporal-control source-calibrated electronic noise mean application. No target ROI data enters z_control; existing saved pilot contributes covariance/state and the separately reported curvature diagnostic only. No fit or PSD retuning.',
        'exposure':STEM,'patch':'A','frozen_controls':len(original),'accepted_controls':nc,'contrast_dimension':len(z),
        'accepted_offsets_xy':offsets[accepted],'rejected_offsets_xy':offsets[~accepted],
        'protocol_sha256':sha(protocolpath),'equations':{
            'control':'z_c = P [ADC_control_group0,...,ADC_control_group8]; P*1=P*t_mean=0, P P^T=I7',
            'noise_mean':'m_noise=C_target,control C_control,control^-1 z_control',
            'conditional_component':'C_target|control=C_target-C_target,control C_control^-1 C_control,target',
            'unchanged_error_comparison':'r_before=data-mean; r_after=r_before-m_noise; original ERR is not changed'},
        'target_covariance_replay_max_abs_SCI2':replay,
        'before':{'original_ERR_Q':float(np.sum((before/error)**2)),'raw_RMS_SCI':float(np.sqrt(np.mean(before**2)))},
        'after_READ_DARK_only_control_law':{'original_ERR_Q':float(np.sum((after/error)**2)),
            'raw_RMS_SCI':float(np.sqrt(np.mean(after**2))),
            'conditional_noise_mean_RMS_SCI':float(np.sqrt(np.mean(mean**2))),
            'conditional_noise_mean_max_abs_SCI':float(np.max(abs(mean))),
            'raw_residual_dot_noise_mean_SCI2':float(before@mean)},
        'source_READ_DARK_variance_before_SCI2_range':[float(np.diag(ctt).min()),float(np.diag(ctt).max())],
        'source_READ_DARK_variance_conditioned_SCI2_range':[float(np.diag(conditional).min()),float(np.diag(conditional).max())],
        'minimum_conditional_covariance_eigenvalue_SCI2':float(np.linalg.eigvalsh(conditional).min()),
        'conditional_variance_reduction_trace_fraction':float(np.trace(ctt-conditional)/np.trace(ctt)),
        'observed_control_contrast_RMS_DN':float(np.sqrt(np.mean(z*z))),
        'nominal_constant_light_control_curvature_RMS_DN':float(np.sqrt(np.mean(curvature**2))),
        'curvature_only_prediction_RMS_SCI':float(np.sqrt(np.mean(curvature_mean**2))),
        'optional_control_photon_projection':photonmeta,'separate_declared_variants':variants,
        'limits':[
            'READ+DARK-only conditioning assumes all retained affine-free ADC control fluctuations are electronics/dark. Actual control source photon noise and nonlinear/variable light may contaminate this assumption; nominal source-shot variant is separate.',
            'Affine temporal contrasts remove any static light only under an affine ADC response. The saved inverse-linearity/cap optical pilot predicts the reported curvature; it is not fitted to control agreement and not subtracted from observed controls.',
            'The direct beta target-noise map holds physical caps/charge history fixed. Full translated cap/raw/post-SIRS/reference stochastic paths and target photon counting covariance remain unresolved in this component.',
            'Stationary masked periodic spectral transfer and one independent calibration difference realization per epoch remain conditional. No target-error rescaling or spectral refit occurs.',
            'This is previously consumed E4A, not unused validation, and no BHSM state kernel, Wick factorization, physical preparation or admission is inferred.'],
        'input_sha256':{str(p.relative_to(ROOT)).replace('\\','/'):sha(p) for p in inputs},
        'source_sha256':{str(p.relative_to(ROOT)).replace('\\','/'):sha(p) for p in [Path(__file__),ROOT/'code/strong_lensing_temporal_controls.py',ROOT/'code/strong_lensing_dark_covariance.py',ROOT/'code/strong_lensing_native_calibration.py']},
        'result_file':str(output.relative_to(ROOT)).replace('\\','/'),'result_sha256':sha(output),
        'seconds':time.perf_counter()-started,'noise_rescaling':False,'physical_validation_pass':False,
        'physical_admission':False,'source_complete_covariance':False,'reserved_observations_opened':False,'physical_guard':2e-10}
    dump(dest/'result.json',report)
    print('CONTROL_MEAN_RESULT',report['before'],report['after_READ_DARK_only_control_law'],flush=True)


if __name__=='__main__':main()
