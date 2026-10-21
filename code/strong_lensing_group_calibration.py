"""Native polynomial/averaged-read observation map, distinct from detector physics.

The pipeline applies its linearity polynomial to an ADC group average. A latent
linearized read mean must therefore be inverse-mapped BEFORE read averaging.
Units here are DN at the recorded pre-flat calibration stage. No fitted noise
scale or physically selected charge-transfer law is introduced by this module.
"""
from __future__ import annotations

import numpy as np


def polynomial_derivatives(coefficients, adc):
    """Evaluate L, L', L''; coefficient axis zero, ascending powers."""
    c = np.asarray(coefficients, dtype=float)
    x = np.asarray(adc, dtype=float)
    if c.ndim < 1 or c.shape[0] < 2 or not np.isfinite(c).all():
        raise ValueError("Finite native coefficients with degree at least one required")
    value = np.zeros_like(x) + c[-1]
    first = np.zeros_like(value)
    second = np.zeros_like(value)
    for coefficient in c[-2::-1]:
        second = second*x + 2*first
        first = first*x + value
        value = value*x + coefficient
    return value, first, second


def inverse_linearity(linear_dn, coefficients, initial_adc, *,
                      adc_lower=-np.inf, adc_upper=np.inf,
                      atol=2e-9, rtol=2e-13, maxiter=40):
    """Invert the monotone branch anchored by actual pre-linearity native data.

    An anchor is required rather than silently selecting a polynomial branch.
    Bounds must describe the valid calibration interval when supplied. A failed
    monotone inverse is a domain error, not permission to extrapolate a PSF.
    """
    target = np.asarray(linear_dn, dtype=float)
    x = np.broadcast_to(np.asarray(initial_adc, float), target.shape).copy()
    lower, upper = np.broadcast_arrays(adc_lower, adc_upper)
    if (not np.isfinite(target).all() or not np.isfinite(x).all()
            or np.any(lower >= upper)):
        raise ValueError("Finite read means/anchors and ordered ADC bounds required")
    x = np.clip(x, lower, upper)
    tolerance = atol + rtol*np.abs(target)
    for _ in range(maxiter):
        value, first, _ = polynomial_derivatives(coefficients, x)
        if np.any(first <= 0) or not np.isfinite(first).all():
            raise ValueError("Native linearity branch is not monotone at the anchor")
        residual = value-target
        if np.all(np.abs(residual) <= tolerance):
            return x
        step = residual/first
        factor = 1.
        for _ in range(30):
            trial = np.clip(x-factor*step, lower, upper)
            vt, dt, _ = polynomial_derivatives(coefficients, trial)
            if (np.all(dt > 0) and np.isfinite(vt).all()
                    and np.all(np.abs(vt-target) <= np.abs(residual)+tolerance)):
                break
            factor *= .5
        else:
            raise ValueError("No valid monotone inverse-linearity step")
        if np.array_equal(trial, x):
            raise ValueError("Latent read mean lies outside the native calibration interval")
        x = trial
    raise RuntimeError("Native inverse-linearity map did not converge")


def averaged_read_observation(linear_read_dn, coefficients, initial_adc, *,
                              dark_group_dn=0., adc_lower=-np.inf,
                              adc_upper=np.inf, read_axis=1):
    """Return L(mean_r L^-1(q_r))-dark and exact local read derivatives.

    Input layout is group, read, ...pixel by default. The latent q_r must include
    the actual readwise dark contribution before this map; the pipeline's group
    dark is subtracted only afterwards. Segment intercepts are part of q, not an
    extra ADC term added after the nonlinear map. Outputs include the Hessian
    needed for a declared second-order mean correction, not an exact probability
    distribution for transformed counts.
    """
    q = np.asarray(linear_read_dn, float)
    axis = read_axis % q.ndim
    if q.ndim < 2 or q.shape[axis] < 1:
        raise ValueError("Explicit group/read support required")
    adc = inverse_linearity(q, coefficients, initial_adc,
                           adc_lower=adc_lower, adc_upper=adc_upper)
    _, d, dd = polynomial_derivatives(coefficients, adc)
    average = np.mean(adc, axis=axis)
    mean, group_d, group_dd = polynomial_derivatives(coefficients, average)
    n = q.shape[axis]
    inverse_d = 1/d
    inverse_dd = -dd/d**3
    gd = np.expand_dims(group_d, axis)
    jacobian = gd*inverse_d/n
    # Place the two read indices immediately after the group index.
    inv = np.moveaxis(inverse_d, axis, 1)
    inv2 = np.moveaxis(inverse_dd, axis, 1)
    hessian = (group_dd[:, None, None]*inv[:, :, None]*inv[:, None, :]/n**2)
    for r in range(n):
        hessian[:, r, r] += group_d*inv2[:, r]/n
    return {"mean_DN": mean-np.asarray(dark_group_dn),
            "inverse_ADC_reads_DN": adc,
            "read_to_group_jacobian": jacobian,
            "read_to_group_hessian": hessian,
            "native_group_linearity_derivative": group_d,
            "native_group_linearity_second_derivative": group_dd,
            "unit_jacobian_approximation": False}


def second_order_mean_correction(hessian, within_group_read_covariance):
    """1/2 Hessian:C for declared latent-DN read covariance (no ERR scaling)."""
    h = np.asarray(hessian, float)
    c = np.asarray(within_group_read_covariance, float)
    if h.shape != c.shape or h.ndim < 3:
        raise ValueError("Matching group,read,read,...pixel tensors required")
    return .5*np.sum(h*c, axis=(1, 2))


def shared_reference_factor(reference_derivative, reference_error):
    """J diag(sigma) for a supplied shared reference, retaining signed effects.

    Rows sharing a physical reference pixel must share its factor column across
    groups/exposures. A reference ERR map supplies marginal errors, not missing
    spatial/temporal reference correlations; callers must not invent independent
    group errors from repeated uses of the same reference pixel.
    """
    j = np.asarray(reference_derivative, float)
    sigma = np.asarray(reference_error, float)
    if j.ndim != 2 or sigma.shape != (j.shape[1],):
        raise ValueError("Observation by shared-reference Jacobian and ERR required")
    if not np.isfinite(j).all() or not np.isfinite(sigma).all() or np.any(sigma < 0):
        raise ValueError("Finite reference derivatives and nonnegative ERR required")
    return j*sigma
