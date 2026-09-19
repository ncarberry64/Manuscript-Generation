from __future__ import annotations

import numpy as np


def pushforward_channel_covariance(
    response: np.ndarray,
    channel_covariance: np.ndarray,
) -> np.ndarray:
    """C_O = M C_q M^T."""
    M = np.asarray(response, dtype=float)
    Cq = np.asarray(channel_covariance, dtype=float)
    if M.ndim != 2:
        raise ValueError("response must be a matrix")
    if Cq.shape != (M.shape[1], M.shape[1]):
        raise ValueError("channel covariance has incompatible shape")
    C = M @ Cq @ M.T
    return 0.5 * (C + C.T)


def covariance_rank(
    covariance: np.ndarray,
    tol: float = 1e-10,
) -> int:
    C = np.asarray(covariance, dtype=float)
    s = np.linalg.svd(C, compute_uv=False)
    if s.size == 0:
        return 0
    threshold = float(tol) * max(float(s[0]), 1.0)
    return int(np.sum(s > threshold))


def rank_diagnostics(
    covariance: np.ndarray,
) -> dict:
    C = np.asarray(covariance, dtype=float)
    if C.ndim != 2 or C.shape[0] != C.shape[1]:
        raise ValueError("covariance must be square")
    C = 0.5 * (C + C.T)
    eig = np.linalg.eigvalsh(C)[::-1]
    pos = np.maximum(eig, 0.0)
    total = float(np.sum(pos))
    lam1 = float(pos[0]) if len(pos) else 0.0
    lam2 = float(pos[1]) if len(pos) > 1 else 0.0
    return {
        "eigenvalues_descending": eig.tolist(),
        "positive_trace": total,
        "R1": lam1 / total if total > 0 else 0.0,
        "rho21": lam2 / lam1 if lam1 > 0 else 0.0,
    }


def principal_minor_residuals(
    covariance: np.ndarray,
) -> np.ndarray:
    """All Cii*Cjj-Cij^2 for i<j."""
    C = np.asarray(covariance, dtype=float)
    if C.ndim != 2 or C.shape[0] != C.shape[1]:
        raise ValueError("covariance must be square")
    out = []
    for i in range(C.shape[0]):
        for j in range(i + 1, C.shape[0]):
            out.append(C[i,i] * C[j,j] - C[i,j] ** 2)
    return np.asarray(out, dtype=float)


def nearest_psd_eigenclip(covariance: np.ndarray) -> np.ndarray:
    """Predeclared simple PSD projection: clip negative eigenvalues to zero."""
    C = np.asarray(covariance, dtype=float)
    C = 0.5 * (C + C.T)
    vals, vecs = np.linalg.eigh(C)
    vals = np.maximum(vals, 0.0)
    return (vecs * vals[None,:]) @ vecs.T


def residual_covariance(
    observed: np.ndarray,
    standard: np.ndarray,
    noise: np.ndarray,
) -> np.ndarray:
    C = (
        np.asarray(observed, dtype=float)
        - np.asarray(standard, dtype=float)
        - np.asarray(noise, dtype=float)
    )
    return 0.5 * (C + C.T)
