from __future__ import annotations

import numpy as np


def spatial_mixed_source_tensor(
    fp: np.ndarray,
    fpp: np.ndarray,
    j: np.ndarray,
    grad_phi: np.ndarray,
    A: np.ndarray,
    g: np.ndarray,
) -> np.ndarray:
    """Pointwise static spatial source tensor S_ij.

    Arrays:
      fp, fpp: (...,)
      j, grad_phi: (..., d)
      A, g: (..., d, d)
    """
    fp = np.asarray(fp, dtype=float)
    fpp = np.asarray(fpp, dtype=float)
    j = np.asarray(j, dtype=float)
    grad_phi = np.asarray(grad_phi, dtype=float)
    A = np.asarray(A, dtype=float)
    g = np.asarray(g, dtype=float)

    if j.shape != grad_phi.shape:
        raise ValueError("j and grad_phi must have same shape")
    if A.shape != g.shape:
        raise ValueError("A and g must have same shape")
    if A.shape[:-2] != j.shape[:-1]:
        raise ValueError("tensor/vector batch shapes are incompatible")
    if A.shape[-1] != A.shape[-2] or A.shape[-1] != j.shape[-1]:
        raise ValueError("spatial dimension mismatch")
    if fp.shape != j.shape[:-1] or fpp.shape != j.shape[:-1]:
        raise ValueError("fp/fpp batch shape mismatch")

    s = np.sum(j * grad_phi, axis=-1)
    sym = 0.5 * (
        j[..., :, None] * grad_phi[..., None, :]
        + grad_phi[..., :, None] * j[..., None, :]
    )

    return (
        2.0 * fpp[..., None, None] * s[..., None, None] * A
        + 2.0 * fp[..., None, None] * sym
        - fp[..., None, None] * s[..., None, None] * g
    )


def project_mixed_coefficients(
    weights: np.ndarray,
    metric_basis_A: np.ndarray,
    metric_basis_psi: np.ndarray,
    source_tensor: np.ndarray,
) -> np.ndarray:
    """Discrete/quadrature projection giving (u,v).

    weights has shape (N,), basis/source tensors have shape (N,d,d).
    """
    w = np.asarray(weights, dtype=float)
    eA = np.asarray(metric_basis_A, dtype=float)
    eP = np.asarray(metric_basis_psi, dtype=float)
    S = np.asarray(source_tensor, dtype=float)

    if w.ndim != 1:
        raise ValueError("weights must be one-dimensional")
    if eA.shape != eP.shape or eA.shape != S.shape:
        raise ValueError("metric basis tensors and source tensor must match")
    if eA.shape[0] != w.shape[0]:
        raise ValueError("quadrature length mismatch")

    u = np.sum(w * np.einsum("nij,nij->n", eA, S))
    v = np.sum(w * np.einsum("nij,nij->n", eP, S))
    return np.array([u, v], dtype=float)


def support_diagnostics(
    weights: np.ndarray,
    j: np.ndarray,
    grad_phi: np.ndarray,
) -> dict:
    w = np.asarray(weights, dtype=float)
    j = np.asarray(j, dtype=float)
    gp = np.asarray(grad_phi, dtype=float)
    if j.shape != gp.shape or j.shape[0] != w.shape[0]:
        raise ValueError("shape mismatch")

    j2 = np.sum(j*j, axis=-1)
    g2 = np.sum(gp*gp, axis=-1)
    overlap = np.sum(j*gp, axis=-1)

    return {
        "j_norm": float(np.sqrt(np.sum(w * j2))),
        "grad_phi_norm": float(np.sqrt(np.sum(w * g2))),
        "j_dot_grad_phi_integral": float(np.sum(w * overlap)),
        "support_possible": bool(
            np.any(np.linalg.norm(j, axis=-1) > 0.0)
            and np.any(np.linalg.norm(gp, axis=-1) > 0.0)
        ),
    }


def rank_one_source_covariance(
    mixed_column: np.ndarray,
    amplitude_variance: float,
) -> np.ndarray:
    b = np.asarray(mixed_column, dtype=float)
    if b.ndim != 1:
        raise ValueError("mixed_column must be a vector")
    if amplitude_variance < 0.0:
        raise ValueError("amplitude_variance must be nonnegative")
    return float(amplitude_variance) * np.outer(b, b)


def multichannel_source_covariance(
    mixed_columns: np.ndarray,
    amplitude_covariance: np.ndarray,
) -> np.ndarray:
    """B C_q B^T, with columns of B corresponding to texture channels."""
    B = np.asarray(mixed_columns, dtype=float)
    C = np.asarray(amplitude_covariance, dtype=float)
    if B.ndim != 2:
        raise ValueError("mixed_columns must be a matrix")
    if C.shape != (B.shape[1], B.shape[1]):
        raise ValueError("amplitude covariance has incompatible shape")
    out = B @ C @ B.T
    return 0.5 * (out + out.T)
