"""One actual nonlocal G0 electronic-intercept READ+DARK contraction.

Uses already consumed E4A, saved kappa=0/cap1 optical pilot, and the unchanged
independent earliest-dark spectrum. No fit, reserved observation, variance
floor, photon covariance substitution, or physical validation occurs here.
"""
from pathlib import Path
import hashlib,json,sys,time
import numpy as np
from scipy import sparse,linalg
from astropy.io import fits
ROOT=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT/'code'))
from strong_lensing_ipc_response import inverse_deconvolution_kernel,convolution_matrix
from strong_lensing_native_calibration import CalibrationOrderedHistory
from strong_lensing_native_ols import piecewise_ols_arrays
from strong_lensing_anchor_covariance import electronic_anchor_rate_map,anchored_noise_contraction
from strong_lensing_dark_covariance import BinnedSpectralModel,native_to_clock_coordinates

BASE=ROOT/'data/strong_lensing_response_completion_20261021'
ART=ROOT/'artifacts/strong_lensing_response_completion_20261021'
STEM='jw01198007001_02101_00004_nrcb3'


def sha(path):return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def dump(path,value):
    def default(x):
        if isinstance(x,np.ndarray):return x.tolist()
        if isinstance(x,np.generic):return x.item()
        raise TypeError(type(x).__name__)
    path.write_text(json.dumps(value,indent=2,allow_nan=False,default=default)+'\n',encoding='utf-8',newline='\n')


def ols(h,groups,zero,selected,cast):
    p=0;s=selected
    return piecewise_ols_arrays(np.asarray(groups).reshape(9,-1)[:,s],h['groupdq'][p].reshape(9,-1)[:,s],
        gain=h['gain_e_per_DN'][p].ravel()[s],readnoise_cds_dn=h['readnoise_CDS_DN'][p].ravel()[s],
        group_time=float(h['group_time_header_s']),frame_time=float(h['frame_time_s']),nframes=2,
        dark_dn_s=h['average_dark_current_DN_s'][p].ravel()[s],zeroframe=np.asarray(zero).ravel()[s],
        pixeldq=h['pixeldq'][p].ravel()[s],cast_input=cast)


def main():
    started=time.perf_counter()
    dest=ART/'dark_calibration/nonlocal_G0_E4A';dest.mkdir(parents=True,exist_ok=True)
    inputs={
        'history':BASE/'histories'/f'{STEM}_histories.npz',
        'pilot':ART/'common_early/pilot_means'/f'training_plus_held_{STEM}_pilot.npz',
        'saved_application':ART/'calibration_order/k0'/f'{STEM}_A_cap1.npz',
        'spectrum':BASE/'dark_calibration/earliest_epoch_binned_cross_spectrum_512.npz',
        'IPC':BASE/'native_cache/jwst_nircam_ipc_0022.fits'}
    with np.load(inputs['history']) as f:h={k:f[k] for k in f.files}
    with np.load(inputs['pilot']) as f:pilot={k:f[k] for k in f.files}
    with np.load(inputs['saved_application']) as f:saved={k:f[k] for k in f.files}
    gkernel=fits.getdata(inputs['IPC'],1)
    inverse=inverse_deconvolution_kernel(gkernel,order=8)
    g=convolution_matrix(gkernel,(65,65));forward=convolution_matrix(inverse['kernel'],(65,65))
    clock=float(h['group_time_header_s'])/(5*float(h['frame_time_s']))
    effective=pilot['point_plus_sky_e_s'][0]*clock
    cal=CalibrationOrderedHistory(h,0,g,forward,effective,ceiling_multiplier=1.)
    state=cal.evaluate(effective,0.,derivatives=True)
    if not np.array_equal(state['groups_DN'],saved['native_groups_DN']):
        raise ValueError('Reevaluated conditional state differs from preserved original application')
    selected=np.flatnonzero(saved['valid'].ravel())
    observed,wo,bo=ols(h,h['corrected_group_DN'][0],h['native_corrected_zeroframe_DN'][0],selected,True)
    predicted,wp,bp=ols(h,state['groups_DN'],state['zeroframe_DN'],selected,False)
    actual=h['native_rate_pre_gain_scale_DN_s'][0].ravel()[selected]
    if not np.array_equal(observed.astype(np.float32),actual.astype(np.float32)):
        raise ValueError('Observed continuous native extension does not reproduce original C values')
    if np.any(wo[9]!=0) or np.any(wp[9]!=0):
        raise ValueError('Actual source consumes ZEROFRAME; dark9group-only projection is insufficient')
    b,k=electronic_anchor_rate_map(forward,g,selected,wp[:9],
        state['observation']['native_linearity_derivative'],state['observation']['inverse_linearity_derivative'],
        native_zero_weights=wp[9])
    nativefactor=float(h['photmjsr'])/h['flat'][0].ravel()[selected]
    gain_scale=h['native_rate_DN_s'][0].ravel()[selected]/actual
    if not np.all(gain_scale==1):raise ValueError('Unaccounted native gain_scale')
    k=sparse.diags(nativefactor)@k
    context=np.unique(k.indices)
    # Only the direct intercept contraction is considered: conditional charge
    # history/caps, trusted continuation mask, SIRS response and references fixed.
    kk=k[:,context].toarray()
    observed_L=np.asarray(h['linearity_polynomial_first_derivative'][0],float).reshape(9,-1).copy()
    observed_L[:,cal.identity]=1.
    weights=wo[:9].T*observed_L[:,selected].T*nativefactor[:,None]
    native=np.c_[h['x_pixel'][0].ravel(),h['y_pixel'][0].ravel()]
    amp,xy=native_to_clock_coordinates(np.vstack((native[selected],native[context])))
    request={'amp':amp,'clock_xy':xy,'group_weights':np.r_[weights,np.c_[np.ones(len(context)),np.zeros((len(context),8))]]}
    with np.load(inputs['spectrum']) as f:
        model=BinnedSpectralModel(f['spectral_matrices'],f['y_frequency_bin'],f['x_frequency_bin'],f['spatial_shape'],int(f['namp']),int(f['ngroup']))
    print('NONLOCAL_SOURCE_DOMAIN',len(selected),len(context),'joint',len(amp),flush=True)
    joint=model.project([request])[0]
    covariance,transform=anchored_noise_contraction(joint,kk,len(selected))
    unanchored=joint[:len(selected),:len(selected)]
    cross=-joint[:len(selected),len(selected):]@kk.T
    added=kk@joint[len(selected):,len(selected):]@kk.T
    residual=h['archived_cal_MJy_sr'][0].ravel()[selected]-saved['mean_SCI'].ravel()[selected]
    error=h['archived_cal_ERR'][0].ravel()[selected]
    chol=linalg.cholesky(covariance,lower=True)
    white=linalg.solve_triangular(chol,residual,lower=True)
    q=float(white@white);logdet=float(2*np.log(np.diag(chol)).sum())
    diagonal=np.zeros_like(kk)
    for i,p in enumerate(selected):
        j=np.searchsorted(context,p)
        if j<len(context) and context[j]==p:diagonal[i,j]=kk[i,j]
    diagonal_c,_=anchored_noise_contraction(joint,diagonal,len(selected))
    output=dest/'nonlocal_READ_DARK_component.npz'
    np.savez_compressed(output,selected_context_index=selected,anchor_context_index=context,
        selected_native_xy=native[selected],anchor_native_xy=native[context],
        source_joint_READ_DARK_covariance=joint,nonlocal_prediction_G0_ADC_to_SCI=kk,
        observed_native_ADC_to_SCI_weights=weights,
        residual_READ_DARK_covariance_SCI2=covariance,
        unanchored_READ_DARK_covariance_SCI2=unanchored,
        anchor_cross_term_SCI2=cross+cross.T,anchor_self_term_SCI2=added,
        diagonal_anchor_counterfactual_covariance_SCI2=diagonal_c,
        unchanged_raw_residual_SCI=residual,unchanged_original_native_ERR=error,
        conditional_READ_DARK_only_whitened_residual=white)
    readcap=state['physical']['observed_cap_jacobian']
    allused=np.unique(forward[selected].indices)
    capactive=[int(np.count_nonzero(readcap[i][allused].data)) for i in cal.inverse]
    receipt={
        'scope':'Actual E4A direct electronic-intercept noisy-G0 READ+DARK contraction at preserved kappa=0/cap1 mean. Conditional charge history/caps, references, SIRS mapping and masks are held fixed. No mean or source calibration fit.',
        'exposure':STEM,'patch':'A','kappa_e_inverse':0.,'ceiling_multiplier':1.,
        'equations':{
            'B':'sum_j diag(Wpred_j Lnativeprime_j) H diag(mean_read invLprime_j), plus actual ZEROFRAME term when consumed',
            'anchor':'d beta/d ADCobserved0 = diag(1/mean_read invLprime_0) G, at fixed incident charge history and caps',
            'K':'PHOTMJSR/flat times B diag(1/mean_read invLprime_0) G',
            'residual':'R=R_observed-R_predicted(theta,ADCobserved0); J0=diag(Wobs0 Lobsprime0 PHOTMJSR/flat)-K; Jj=diag(Wobsj Lobsprimej PHOTMJSR/flat), j>0',
            'covariance':'C_R = [I,-K] C_same_realization(observed_native_rate,ADC0_context) [I,-K]^T; shared cross terms retained'},
        'selected_outputs':len(selected),'G0_context_pixels':len(context),'combined_source_channels':len(amp),
        'G0_native_bounds_xyxy':[int(native[context,0].min()),int(native[context,1].min()),int(native[context,0].max()+1),int(native[context,1].max()+1)],
        'observed_native_values_bitwise_float32':True,
        'observed_nonunique_branch_count':sum(not v.get('locally_unique_derivative',False) for v in bo),
        'predicted_nonunique_branch_count':sum(not v.get('locally_unique_derivative',False) for v in bp),
        'zero_weights_proven_zero':True,'replayed_state_groups_bitwise':True,
        'original_native_ERR_Q':float(np.sum((residual/error)**2)),
        'raw_RMS_SCI':float(np.sqrt(np.mean(residual**2))),
        'conditional_READ_DARK_only_Q':q,'conditional_READ_DARK_only_Q_per_pixel':q/len(selected),
        'conditional_READ_DARK_only_logdet':logdet,
        'conditional_READ_DARK_only_Gaussian_NLL':.5*(q+logdet+len(selected)*np.log(2*np.pi)),
        'unanchored_variance_SCI2_range':[float(np.diag(unanchored).min()),float(np.diag(unanchored).max())],
        'anchored_variance_SCI2_range':[float(np.diag(covariance).min()),float(np.diag(covariance).max())],
        'minimum_covariance_eigenvalue_SCI2':float(np.linalg.eigvalsh(covariance).min()),
        'anchor_K_nonlocal_max_abs_SCI_per_DN':float(np.max(abs(kk-diagonal))),
        'full_vs_diagonal_anchor_relative_covariance_Frobenius':float(np.linalg.norm(covariance-diagonal_c)/np.linalg.norm(covariance)),
        'anchor_self_term_max_abs_SCI2':float(np.max(abs(added))),
        'anchor_same_realization_cross_term_max_abs_SCI2':float(np.max(abs(cross+cross.T))),
        'component_identity_max_abs_SCI2':float(np.max(abs(covariance-(unanchored+cross+cross.T+added)))),
        'cap_derivative_active_nonzeros_on_H_output_support_by_read':capactive,
        'limits':[
            'This closes the direct beta-intercept path at a FIXED conditional cap/incident charge state, not the unconditional cap/threshold noise path. Translated caps depend on group0 and reference-control correction; their joint noise and physical full-well stochastic law remain separate.',
            'Measured source is raw READ+DARK after same-current-channel SIRS, including spatial/time/reference cross terms under the frozen periodic stationary spectral model. Finite epoch, masked spectral-window and calibration uncertainty remain conditional.',
            'No source photon/IPC/saturation counting covariance, flat/linearity shared calibration covariance, or optical dressed-state statistical kernel is supplied by this component.',
            'No native ERR or native VAR_RNOISE is added to this component; original ERR scores and raw residuals remain separate. A READ+DARK-only Gaussian score is a component diagnostic, not a full statistical likelihood or validation.',
            'The optical pilot and E4 native data were already consumed; this is an internal held-channel contraction, not a historically unused science test.'],
        'input_sha256':{str(p.relative_to(ROOT)).replace('\\','/'):sha(p) for p in inputs.values()},
        'source_sha256':{str(p.relative_to(ROOT)).replace('\\','/'):sha(p) for p in [Path(__file__),ROOT/'code/strong_lensing_anchor_covariance.py',ROOT/'code/strong_lensing_dark_covariance.py',ROOT/'code/strong_lensing_native_calibration.py',ROOT/'code/strong_lensing_native_ols.py']},
        'result_file':str(output.relative_to(ROOT)).replace('\\','/'),'result_sha256':sha(output),
        'seconds':time.perf_counter()-started,'noise_rescaling':False,'physical_validation':False,
        'source_complete_covariance':False,'reserved_observations_opened':False,'physical_guard':2e-10}
    dump(dest/'result.json',receipt)
    print('NONLOCAL_COMPONENT_RESULT',q,'original',receipt['original_native_ERR_Q'],'RMS',receipt['raw_RMS_SCI'],flush=True)


if __name__=='__main__':main()
