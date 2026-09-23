from __future__ import annotations

import numpy as np


def static_seam_export(
    r,
    r_star: float,
    q_sigma: float,
    q_coupling: float = 1.0,
    orientation: float = 1.0,
):
    """Static spherical seam export.

    E0 = orientation * (q_sigma/q_coupling) * Theta(r-r_star)/r^2.
    """
    if r_star <= 0.0:
        raise ValueError("r_star must be positive")
    if q_coupling == 0.0:
        raise ValueError("q_coupling must be nonzero")

    rr = np.asarray(r, dtype=float)
    if np.any(rr <= 0.0):
        raise ValueError("r must be positive")

    return (
        float(orientation)
        * float(q_sigma)
        / float(q_coupling)
        * (rr >= r_star).astype(float)
        / (rr * rr)
    )


def dynamic_seam_export_components(
    r_star: float,
    q_sigma: float,
    q_sigma_dot: float,
    r_star_dot: float,
    q_coupling: float = 1.0,
    orientation: float = 1.0,
):
    """Distributional coefficients of d/dt of the static seam export.

    Returns:
      bulk_prefactor multiplying Theta(r-r_star)/r^2
      seam_delta_coefficient multiplying delta(r-r_star)

    The latter already includes the 1/r_star^2 factor.
    """
    if r_star <= 0.0:
        raise ValueError("r_star must be positive")
    if q_coupling == 0.0:
        raise ValueError("q_coupling must be nonzero")

    s = float(orientation) / float(q_coupling)
    bulk = s * float(q_sigma_dot)
    seam_delta = (
        -s
        * float(q_sigma)
        * float(r_star_dot)
        / (float(r_star) ** 2)
    )
    return bulk, seam_delta


def local_metric_response(
    a: float,
    b: float,
    c: float,
    u: float,
    v: float,
):
    """Exact -A_g^{-1} B_gpsi for a two-coordinate physical metric basis."""
    det = float(a) * float(c) - float(b) ** 2
    if abs(det) < 1.0e-15:
        raise ValueError("local physical metric Hessian is singular")

    r_phi = (float(b) * float(v) - float(c) * float(u)) / det
    r_psi = (float(b) * float(u) - float(a) * float(v)) / det
    return np.array([r_phi, r_psi], dtype=float)


def weyl_response_coefficient(
    a: float,
    b: float,
    c: float,
    u: float,
    v: float,
) -> float:
    """Coefficient of (Phi+Psi)/2 per unit local topographic field."""
    r = local_metric_response(a, b, c, u, v)
    return float(0.5 * (r[0] + r[1]))


def slip_response_coefficient(
    a: float,
    b: float,
    c: float,
    u: float,
    v: float,
) -> float:
    """Coefficient of (Phi-Psi)/2 per unit local topographic field."""
    r = local_metric_response(a, b, c, u, v)
    return float(0.5 * (r[0] - r[1]))


def pushforward_covariance(
    optical_map: np.ndarray,
    metric_response: np.ndarray,
    green: np.ndarray,
    source_covariance: np.ndarray,
) -> np.ndarray:
    """Q_opt = O R G C_J G^T R^T O^T for finite discretizations."""
    O = np.asarray(optical_map, dtype=float)
    R = np.asarray(metric_response, dtype=float)
    G = np.asarray(green, dtype=float)
    C = np.asarray(source_covariance, dtype=float)

    M = O @ R @ G
    out = M @ C @ M.T
    return 0.5 * (out + out.T)


def compound_poisson_source_covariance(
    rate: float,
    amplitude_second_moment: np.ndarray,
) -> np.ndarray:
    """White shot-noise covariance lambda * E[q q^T]."""
    if rate < 0.0:
        raise ValueError("rate must be nonnegative")
    M2 = np.asarray(amplitude_second_moment, dtype=float)
    if M2.ndim != 2 or M2.shape[0] != M2.shape[1]:
        raise ValueError("amplitude_second_moment must be square")
    return float(rate) * M2
