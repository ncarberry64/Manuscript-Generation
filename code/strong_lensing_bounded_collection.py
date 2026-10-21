"""Explicit first-order bounded collection continuations for consumed ramps.

The leading unsaturated source is lambda*t + kappa*t**2/2*L(lambda**2).
It does not determine post-saturation charge storage. Two alternatives are
therefore evaluated separately: a finite storage ceiling with lost overflow,
and an electronic readout ceiling with continued stored charge. Neither is an
identified NIRCam full-well law. No threshold reference is relabeled full well.
"""
from __future__ import annotations

import numpy as np
from scipy import sparse
from strong_lensing_stochastic_response import nearest_neighbor_edges


def _inputs(rate_e_s, times_s, kappa_e_inverse, charge_ceiling_e, edges):
    rate = np.asarray(rate_e_s, float)
    times = np.asarray(times_s, float)
    caps = np.broadcast_to(np.asarray(charge_ceiling_e, float), rate.shape)
    kappa = float(kappa_e_inverse)
    edge = np.asarray(edges, int).reshape(-1, 2)
    if (rate.ndim != 1 or not len(rate) or times.ndim != 1 or not len(times)
            or not np.isfinite(rate).all() or np.any(rate < 0)
            or not np.isfinite(times).all() or np.any(times < 0)
            or np.any(np.diff(times) < 0) or not np.isfinite(caps).all()
            or np.any(caps <= 0) or not np.isfinite(kappa) or kappa < 0):
        raise ValueError("Finite nonnegative rates/times/kappa and positive finite ceilings required")
    if (edge.size and (edge.min() < 0 or edge.max() >= len(rate)
                      or np.any(edge[:, 0] == edge[:, 1]))):
        raise ValueError("Different in-domain neighbor endpoints required")
    if len(set(map(tuple, np.sort(edge, axis=1)))) != len(edge):
        raise ValueError("Duplicate physical collection edge")
    return rate, times, kappa, caps, edge


def bounded_leading_collection(rate_e_s, times_s, kappa_e_inverse,
                               charge_ceiling_e, edges, *, prescription,
                               ceiling_source, derivatives=False):
    """Evaluate a declared conditional first-order same-domain continuation.

    Let H_p(t)=integral_0^t min(lambda_p*s,Q_p) ds. For storage_ceiling,
      D_p(t)=sum_(p,q) (lambda_p+lambda_q)*(H_q-H_p),
      m_p(t)=min(lambda_p*t+kappa*D_p(t), Q_p).
    Excess carriers are LOST in this conditional completion, not assigned to
    an unmodeled neighbor. In readout_ceiling, H=lambda*t**2/2 and the stored
    mean remains unbounded; only the reported per-read mean is clipped.
    Both reproduce the original leading mean exactly before any ceiling.
    Collection depends on zeroth-order stored means: this is first order in
    kappa, not a resummed nonlinear saturation or stochastic overflow solve.

    Charge ceilings are supplied physical/electronic model parameters. Their
    source and meaning are mandatory. Group SAT flags alone do not fix them.
    Reset/electronic intercepts and the native inverse-linearity/group map are
    separate caller-supplied applications. Piecewise analytic Jacobians hold
    on the selected active branch; exact transition kinks are reported.
    """
    if prescription not in ("storage_ceiling", "readout_ceiling"):
        raise ValueError("Choose storage_ceiling or readout_ceiling explicitly")
    if not str(ceiling_source).strip():
        raise ValueError("Ceiling defining relation/source is required")
    rate, times, kappa, caps, edge = _inputs(
        rate_e_s, times_s, kappa_e_inverse, charge_ceiling_e, edges)
    t = times[:, None]
    linear = t*rate[None, :]
    early = linear < caps
    if prescription == "storage_ceiling":
        invrate = np.divide(1., rate, out=np.zeros_like(rate), where=rate > 0)
        late_h = caps[None, :]*t-.5*(caps*caps*invrate)[None, :]
        h = np.where(early, .5*t*t*rate[None, :], late_h)
        dh_rate = np.where(early, .5*t*t, .5*(caps*caps*invrate**2)[None, :])
        dh_cap = np.where(early, 0., t-(caps*invrate)[None, :])
    else:
        h = .5*t*t*rate[None, :]
        dh_rate = np.broadcast_to(.5*t*t, linear.shape)
        dh_cap = np.zeros_like(linear)
    p, q = edge.T
    summed_rate = rate[p]+rate[q]
    difference = h[:, q]-h[:, p]
    drift = np.zeros_like(linear)
    for j in range(len(times)):
        np.add.at(drift[j], p, summed_rate*difference[j])
        np.add.at(drift[j], q, -summed_rate*difference[j])
    before = linear+kappa*drift
    if not np.isfinite(before).all() or np.any(before < 0):
        raise ValueError("First-order continuation produced negative/invalid stored charge")
    active = before < caps
    observed = np.minimum(before, caps)
    stored = observed if prescription == "storage_ceiling" else before
    lost = np.sum(linear-stored, axis=1) if prescription == "storage_ceiling" else np.zeros(len(times))
    if np.any(lost < -1e-9*np.maximum(1., linear.sum(axis=1))):
        raise ArithmeticError("Conditional overflow bookkeeping creates charge")
    eps = 32*np.finfo(float).eps
    kink = (np.abs(before-caps) <= eps*np.maximum(1., caps))
    if prescription == "storage_ceiling":
        kink |= np.abs(linear-caps) <= eps*np.maximum(1., caps)
    answer = {
        "stored_mean_e": stored, "observed_read_mean_e": observed,
        "pre_ceiling_mean_e": before, "first_order_drift_e2": drift,
        "lost_mean_e": lost, "clipped_read_mask": ~active,
        "nondifferentiable_transition_mask": kink,
        "diagnostics": {
            "prescription": prescription, "ceiling_source": str(ceiling_source),
            "order": "first order in kappa; zeroth-order stored-charge boundary response",
            "overflow": "lost carriers" if prescription == "storage_ceiling" else "none; readout clipping only",
            "no_exterior_edges": True, "finite_context_pixels": len(rate),
            "sum_internal_first_order_drift_e2": drift.sum(axis=1).tolist(),
            "full_well_identified": False, "stochastic_overflow_law": "not defined by this mean continuation",
            "native_group_map_included": False, "physical_validation": False,
            "error_rescaling": False}}
    if derivatives:
        rate_jac = []; cap_jac = []; stored_rate_jac = []
        for j, tj in enumerate(times):
            dp = difference[j]-summed_rate*dh_rate[j, p]
            dq = difference[j]+summed_rate*dh_rate[j, q]
            rr = np.r_[p, p, q, q]
            cc = np.r_[p, q, p, q]
            d = sparse.coo_matrix((kappa*np.r_[dp, dq, -dp, -dq], (rr, cc)),
                                  shape=(len(rate), len(rate))).tocsr()+sparse.eye(len(rate), format="csr")*tj
            dc = sparse.coo_matrix((kappa*np.r_[-summed_rate*dh_cap[j, p], summed_rate*dh_cap[j, q],
                                               summed_rate*dh_cap[j, p], -summed_rate*dh_cap[j, q]], (rr, cc)),
                                   shape=(len(rate), len(rate))).tocsr()
            out_d = d.multiply(active[j, :, None]).tocsr()
            out_dc = dc.multiply(active[j, :, None]).tocsr()+sparse.diags((~active[j]).astype(float), format="csr")
            rate_jac.append(out_d); cap_jac.append(out_dc)
            stored_rate_jac.append(out_d if prescription == "storage_ceiling" else d)
        answer.update({"observed_rate_jacobian_s": tuple(rate_jac),
                       "observed_cap_jacobian": tuple(cap_jac),
                       "stored_rate_jacobian_s": tuple(stored_rate_jac),
                       "observed_kappa_jacobian_e2": active*drift,
                       "stored_kappa_jacobian_e2": active*drift if prescription == "storage_ceiling" else drift})
    return answer
