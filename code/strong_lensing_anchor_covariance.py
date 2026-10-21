"""Nonlocal electronic-G0 anchor contraction on a fixed conditional state.

This differentiates the fitted electronic intercept at FIXED incident charge
history and cap state. Dependence of a translated physical cap on observed G0
is a separate path; it is not silently included by this direct contraction.
"""
from __future__ import annotations
import numpy as np
from scipy import sparse


def electronic_anchor_rate_map(forward_ipc, deconvolution, selected,
                               native_group_weights, native_group_Lprime,
                               inverse_read_Lprime, *, native_zero_weights=None,
                               native_zero_Lprime=None, inverse_zero_Lprime=None):
    """Return B=dR/dbeta and K=dR/dADC_observed,0, retaining H diag(r) G.

    W is group,selected; native L' is group,context; inverse L' is
    group,read,context. Slopes are caller units (normally native DN/s).
    ZEROFRAME enters B only if its actual estimator weight is nonzero.
    """
    h, g = sparse.csr_matrix(forward_ipc), sparse.csr_matrix(deconvolution)
    sel = np.asarray(selected, int)
    w = np.asarray(native_group_weights, float)
    ld = np.asarray(native_group_Lprime, float)
    inv = np.asarray(inverse_read_Lprime, float)
    if (h.shape != g.shape or h.shape[0] != h.shape[1]
            or ld.ndim != 2 or inv.ndim != 3 or inv.shape[0] != ld.shape[0]
            or inv.shape[2] != h.shape[0] or ld.shape[1] != h.shape[0]
            or w.shape != (ld.shape[0], len(sel))
            or np.any(sel < 0) or np.any(sel >= h.shape[0])
            or any(not np.isfinite(v).all() for v in (w, ld, inv))
            or np.any(inv <= 0)):
        raise ValueError("Actual matching monotone read/group maps and selected native weights required")
    hs = h[sel]
    b = sparse.csr_matrix((len(sel), h.shape[0]))
    for j in range(len(w)):
        b += sparse.diags(w[j]*ld[j, sel])@hs@sparse.diags(inv[j].mean(axis=0))
    if native_zero_weights is not None:
        wz = np.asarray(native_zero_weights, float)
        if wz.shape != (len(sel),) or not np.isfinite(wz).all():
            raise ValueError("Matching actual ZEROFRAME weights required")
        if np.any(wz != 0):
            zl, zi = np.asarray(native_zero_Lprime, float), np.asarray(inverse_zero_Lprime, float)
            if (zl.shape != (h.shape[0],) or zi.shape != zl.shape
                    or not np.isfinite(zl).all() or not np.isfinite(zi).all() or np.any(zi <= 0)):
                raise ValueError("Actual ZEROFRAME inverse/native linearity derivatives required")
            b += sparse.diags(wz*zl[sel])@hs@sparse.diags(zi)
    k = b@sparse.diags(1/inv[0].mean(axis=0))@g
    b.eliminate_zeros(); k.eliminate_zeros()
    return b, k


def anchored_noise_contraction(joint_covariance, prediction_map, output_count):
    """Cov([native observed rate, ADC0context]) -> Cov(rate-K ADC0).

    The source block MUST retain observed/context cross terms from the SAME
    calibration realization. No independent anchor-noise covariance is added.
    """
    c = np.asarray(joint_covariance, float)
    k = np.asarray(prediction_map, float)
    m = int(output_count)
    if (m < 1 or k.ndim != 2 or k.shape[0] != m
            or c.shape != (m+k.shape[1],)*2 or not np.isfinite(c).all()
            or not np.isfinite(k).all()):
        raise ValueError("Same-realization joint covariance and nonlocal prediction map required")
    transform = np.c_[np.eye(m), -k]
    result = transform@c@transform.T
    return (result+result.T)/2, transform
