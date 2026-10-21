"""Calibration-ordered IPC/linearity read map on an explicit padded context.

Linearity0054 was fitted after IPC correction, whereas the consumed native
science channel skips IPC. G is the exact recorded deconvolution and H=G^-1.
Thus a latent preIPC linearized mean is inverse-linearized, coupled by H,
averaged, then passed through the ACTUAL native output L. This is a nominal
signal-independent calibration composition, not an identified nonlinear IPC
or post-full-well law. An effective optical PSF must NOT be convolved twice.
"""
from __future__ import annotations
import numpy as np
from scipy import sparse
from strong_lensing_group_calibration import inverse_linearity, polynomial_derivatives


def ipc_averaged_read_observation(linear_read_dn, coefficients,
                                 initial_deconvolved_adc, forward_ipc, *,
                                 dark_group_dn=0., adc_lower=-np.inf,
                                 adc_upper=np.inf, derivatives=False):
    """L_native[ H mean_r L_cal^-1(q_jr) ] - D_native,j.

    Arrays have group,read,pixel shape; coefficients have power,pixel shape.
    H acts on raw DN BEFORE native output linearity, not on electron counts
    or already measured dark/electronic covariance. The external context and
    inverse-kernel tail are defined by the caller. The returned sparse read
    Jacobians retain interpixel coupling and belong on the LEFT of counting
    covariance maps. NO_LIN_CORR context treatment must be supplied explicitly.
    """
    q=np.asarray(linear_read_dn,float);h=sparse.csr_matrix(forward_ipc)
    if (q.ndim!=3 or q.shape[1]<1 or h.shape!=(q.shape[2],q.shape[2])
            or not np.isfinite(h.data).all()):
        raise ValueError('Explicit group,read,pixel support and finite matching IPC map required')
    adc=inverse_linearity(q,coefficients,initial_deconvolved_adc,
                          adc_lower=adc_lower,adc_upper=adc_upper)
    _,pre_d,_=polynomial_derivatives(coefficients,adc)
    mean_deconv=adc.mean(axis=1)
    coupled=np.asarray((h@mean_deconv.T).T)
    mean,native_d,native_dd=polynomial_derivatives(coefficients,coupled)
    result={'mean_DN':mean-np.asarray(dark_group_dn),
            'inverse_deconvolved_ADC_reads_DN':adc,
            'coupled_ADC_group_DN':coupled,
            'inverse_linearity_derivative':1/pre_d,
            'native_linearity_derivative':native_d,
            'native_linearity_second_derivative':native_dd,
            'IPC_after_inverse_linearity':True,'extra_PSF_convolution':False}
    if derivatives:
        result['read_to_native_group_jacobian']=tuple(tuple(
            sparse.diags(native_d[j])@h@sparse.diags(1/pre_d[j,r]/q.shape[1])
            for r in range(q.shape[1])) for j in range(q.shape[0]))
    return result


def first_group_preipc_intercept(read_signal_dn, coefficients,
                                 target_deconvolved_adc, initial_beta, *,
                                 initial_deconvolved_adc=None,
                                 atol_dn=1e-8,maxiter=12):
    """Condition an electronic beta on the ACTUAL first averaged ADC group.

    Solve mean_r L^-1(beta+signal_r)=G*ADC_group0. This fixes an electronic
    boundary nuisance, not physical reset charge or a predicted independent
    group0. The derivative d beta/d signal_r is returned for the same implicit
    relation, so a coupled mean Jacobian does not freeze this fitted nuisance.
    Readwise dark belongs in signal under its separately declared convention.
    """
    signal=np.asarray(read_signal_dn,float)
    target=np.asarray(target_deconvolved_adc,float)
    beta=np.asarray(initial_beta,float).copy()
    if (signal.ndim!=2 or signal.shape[1:]!=target.shape or beta.shape!=target.shape
            or not np.isfinite(signal).all() or not np.isfinite(target).all()):
        raise ValueError('Read,pixel signal and matching finite first-group target required')
    anchor=np.broadcast_to(target if initial_deconvolved_adc is None else initial_deconvolved_adc,signal.shape)
    for _ in range(maxiter):
        adc=inverse_linearity(beta[None]+signal,coefficients,anchor)
        _,d,_=polynomial_derivatives(coefficients,adc)
        inverse_d=1/d
        error=adc.mean(axis=0)-target
        if np.max(abs(error),initial=0.)<=atol_dn:
            weight=inverse_d/inverse_d.sum(axis=0)
            return {'beta_preipc_DN':beta,'inverse_ADC_reads_DN':adc,
                    'beta_signal_derivative':-weight,
                    'first_group_inverse_ADC_match_max_abs_DN':float(np.max(abs(error),initial=0.))}
        beta-=error/inverse_d.mean(axis=0)
    raise RuntimeError('First-group preIPC electronic intercept did not converge')
