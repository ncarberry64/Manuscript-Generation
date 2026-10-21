"""Fixed outside-ROI temporal controls for a conditional detector-noise mean.

This is a same-realization conditional READ+DARK experiment. It changes the
conditional mean by an externally calibrated cross covariance, not by widening
ERR. A static affine-ADC light ramp is removed; nonlinear signal contamination
and uncalibrated collection/reference uncertainty remain separate.
"""
import numpy as np
from scipy import linalg


def affine_free_temporal_contrasts(group_mean_times):
    times=np.asarray(group_mean_times,float)
    if times.ndim!=1 or len(times)<3 or not np.isfinite(times).all() or np.ptp(times)==0:
        raise ValueError('At least three distinct finite actual group mean times required')
    design=np.c_[np.ones(len(times)),(times-times.mean())/np.ptp(times)]
    q,_=linalg.qr(design,mode='full')
    contrast=q[:,2:].T
    for row in contrast:
        if row[np.argmax(abs(row))]<0:row*=-1
    return contrast


def condition_noise_mean(target_covariance,cross_covariance,control_covariance,control_values):
    """Return E[epsilon_target|z_control] and Cov, with no variance inflation.

    Cross and marginal blocks must refer to the SAME source noise realization.
    For a declared Gaussian electronics law these are conditional mean/Cov.
    For general specified moments they are the best LINEAR predictor and its
    error covariance, not generally the exact conditional expectation/Cov.
    No factorization/closure of the BHSM optical state or PSF refit occurs.
    """
    t,x,c,z=(np.asarray(v,float) for v in (target_covariance,cross_covariance,control_covariance,control_values))
    if (t.ndim!=2 or t.shape[0]!=t.shape[1] or c.shape!=(len(z),len(z))
            or x.shape!=(len(t),len(z)) or z.ndim!=1
            or any(not np.isfinite(v).all() for v in (t,x,c,z))):
        raise ValueError('Same-channel target/control covariance blocks and finite control vector required')
    factor=linalg.cho_factor(c,lower=True) # No jitter, floor, or scale.
    transfer=linalg.cho_solve(factor,x.T).T
    mean=transfer@z
    conditional=t-transfer@x.T
    conditional=(conditional+conditional.T)/2
    return mean,conditional,transfer
