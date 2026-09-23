"""Prospective two-gate statistics; the CLI produces SYNTHETIC validation only.

Gate A assumes a fixed, full-rank N x 2 design and known Gaussian covariance.
Gate B supplies an exact Monte Carlo test only for an independently fixed
Gaussian null covariance and iid sampling. Survey applications must freeze
their selection-aware null simulator before inspecting targets. No clipping.
"""
from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path

import numpy as np
from scipy.linalg import solve_triangular
from scipy.stats import chi2

ROOT = Path(__file__).resolve().parents[1]
PROTOCOL = ROOT / "preregistration/bhsm_geometry_first_make_or_break_v1.json"
ALPHA = 0.01


def _real(value, name):
    if np.iscomplexobj(value):
        raise ValueError(f"{name} must be real for this sampling model")
    value = np.asarray(value, dtype=float)
    if not np.all(np.isfinite(value)):
        raise ValueError(f"{name} must be finite")
    return value


def _symmetric(value, name):
    value = _real(value, name)
    if value.ndim != 2 or value.shape[0] != value.shape[1] or not value.size:
        raise ValueError(f"{name} must be a nonempty square matrix")
    if not np.allclose(value, value.T, rtol=1e-12, atol=0):
        raise ValueError(f"{name} must be symmetric")
    return value


def shared_state_test(data, kernels, covariance):
    """GLS projection, with df=N-2. No fitted nuisance parameters in this API."""
    d = _real(data, "data")
    F = _real(kernels, "kernels")
    S = _symmetric(covariance, "covariance")
    if d.ndim != 1 or d.size <= 2 or F.shape != (d.size, 2):
        raise ValueError("require N>2 measurements and exactly two state columns")
    if S.shape != (d.size, d.size):
        raise ValueError("covariance shape does not match data")
    try:
        L = np.linalg.cholesky(S)
    except np.linalg.LinAlgError as exc:
        raise ValueError("covariance must be positive definite") from exc
    dw = solve_triangular(L, d, lower=True)
    Fw = solve_triangular(L, F, lower=True)
    state, _, rank, singular = np.linalg.lstsq(Fw, dw, rcond=None)
    if rank != 2:
        raise ValueError("both state directions must be identifiable")
    residual = dw - Fw @ state
    statistic = float(residual @ residual)
    p = float(chi2.sf(statistic, d.size - 2))
    return {"state": state.tolist(), "chi2_perp": statistic,
            "dof": int(d.size - 2), "p_value": p, "reject": p < ALPHA,
            "whitened_design_condition": float(singular[0] / singular[-1]),
            "scope": "specific two-state deterministic realization"}


def covariance_diagnostics(residual, scales):
    """Signed eigenvalues of a fixed-unit-scaled residual; no PSD projection.

    Scales must be positive and fixed independently of target data. Ratios and
    minors are diagnostics, not additional uncorrected significance tests.
    """
    C = _symmetric(residual, "residual")
    scale = _real(scales, "scales")
    if scale.shape != (C.shape[0],) or np.any(scale <= 0) or C.shape[0] < 2:
        raise ValueError("require >=2 observables and fixed positive scales")
    C = C / np.outer(scale, scale)
    eig = np.linalg.eigvalsh(C)[::-1]
    minors = [float(C[i, i] * C[j, j] - C[i, j] ** 2)
              for i in range(len(C)) for j in range(i + 1, len(C))]
    return {"eigenvalues_descending": eig.tolist(),
            "lambda2": float(eig[1]),
            "rho21": float(eig[1] / eig[0]) if eig[0] > 0 else None,
            "principal_minors": minors,
            "negative_eigenvalues_retained": True}


def monte_carlo_pvalue(observed, null_statistics):
    """Upper-tail Monte Carlo p=(1+#(T_sim>=T_obs))/(B+1)."""
    null = _real(null_statistics, "null statistics")
    if null.ndim != 1 or null.size < 999 or not np.isfinite(observed):
        raise ValueError("require a finite statistic and at least 999 null draws")
    return float((1 + np.count_nonzero(null >= observed)) / (len(null) + 1))


def fixed_gaussian_covariance_test(samples, standard_noise, loading, scales,
                                   *, seed=20260919, draws=9999):
    """Gate B benchmark for an independently fixed rank <=1 Gaussian null.

    loading includes sqrt(channel variance); it is NEVER fitted here. The
    total null is standard_noise + loading loading^T. This tests that fixed
    null, not the composite union of all possible rank-one covariances.
    A survey with estimated nuisances/loading needs its own calibrated test.
    """
    y = _real(samples, "samples")
    S = _symmetric(standard_noise, "standard plus noise covariance")
    v = _real(loading, "loading")
    if y.ndim != 2 or y.shape[0] < 3 or S.shape != (y.shape[1], y.shape[1]):
        raise ValueError("samples must have >=3 rows and compatible covariance")
    if v.shape != (y.shape[1],):
        raise ValueError("loading must contain one fixed channel")
    if not isinstance(draws, int) or draws < 999:
        raise ValueError("draws must be an integer >=999")
    try:
        np.linalg.cholesky(S)
        L = np.linalg.cholesky(S + np.outer(v, v))
    except np.linalg.LinAlgError as exc:
        raise ValueError("benchmark standard/noise covariance must be positive definite") from exc
    diagnostic = covariance_diagnostics(np.cov(y, rowvar=False, ddof=1) - S, scales)
    rng = np.random.default_rng(seed)
    null = np.empty(draws)
    for i in range(draws):
        simulated = rng.standard_normal(y.shape) @ L.T
        null[i] = covariance_diagnostics(
            np.cov(simulated, rowvar=False, ddof=1) - S, scales)["lambda2"]
    p = monte_carlo_pvalue(diagnostic["lambda2"], null)
    return {**diagnostic, "p_value": p, "reject": p < ALPHA,
            "draws": draws, "seed": seed,
            "primary_statistic": "unprojected scaled lambda2 (upper tail)",
            "scope": "independently fixed one-channel iid Gaussian null",
            "composite_rank_null_calibrated": False}


def synthetic_report():
    F = np.array([[1., 0.], [0., 1.], [1., 1.], [1., -1.]])
    state = np.array([0.2, -0.1])
    v = np.array([1., 0.5, -0.3])
    S = np.eye(3) * 0.1
    rng = np.random.default_rng(731)
    samples = rng.standard_normal((300, 3)) @ np.linalg.cholesky(S + np.outer(v, v)).T
    return {"schema": "bhsm-make-or-break-validation-v1",
            "data_kind": "SYNTHETIC_ONLY_NO_TARGET_DATA",
            "protocol_sha256": hashlib.sha256(PROTOCOL.read_bytes()).hexdigest(),
            "gate_a_consistent": shared_state_test(F @ state, F, np.eye(4)),
            "gate_a_inconsistent": shared_state_test(F @ state + np.array([10., 10., -10., 0.]), F, np.eye(4)),
            "gate_b_fixed_null_example": fixed_gaussian_covariance_test(samples, S, v, np.ones(3)),
            "target_readiness": "BLOCKED_PENDING_TARGET_DESIGN_AND_INDEPENDENT_NUISANCE_FREEZE"}


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(synthetic_report(), indent=2, allow_nan=False) + "\n", encoding="utf-8")
