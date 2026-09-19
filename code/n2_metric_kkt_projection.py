from __future__ import annotations

import json
import numpy as np


def metric_matrix(
    a_aa: float,
    a_ap: float,
    a_pp: float,
) -> np.ndarray:
    """Real symmetric physical scalar metric block in the (A,psi) basis."""
    return np.array(
        [
            [a_aa, a_ap],
            [a_ap, a_pp],
        ],
        dtype=float,
    )


def mixed_matrix(
    b_aq: float,
    b_api: float,
    b_pq: float,
    b_ppi: float,
) -> np.ndarray:
    """Map X2=(q2,Pi2) into metric sources (A,psi)."""
    return np.array(
        [
            [b_aq, b_api],
            [b_pq, b_ppi],
        ],
        dtype=float,
    )


def metric_response(A: np.ndarray, B: np.ndarray) -> np.ndarray:
    det = float(np.linalg.det(A))
    if abs(det) < 1e-14:
        raise ValueError("metric block is singular")
    return -np.linalg.solve(A, B)


def schur_scalar(C: np.ndarray, A: np.ndarray, B: np.ndarray) -> np.ndarray:
    """Second-order Schur complement C-B^T A^-1 B."""
    return C - B.T @ np.linalg.solve(A, B)


def diagnostics() -> dict:
    # Pure algebra sanity witness. These are not physical coefficients.
    A = metric_matrix(2.0, 0.25, 3.0)
    B = mixed_matrix(0.1, 0.02, -0.05, 0.03)
    C = np.array([[5.0, -1.0], [-1.0, 13.0]])

    R = metric_response(A, B)
    K = schur_scalar(C, A, B)

    return {
        "status": "N2_METRIC_KKT_PROJECTION_ALGEBRA_VALID",
        "basis_metric": ["A", "psi"],
        "basis_mode": ["q2", "Pi2"],
        "metric_det": float(np.linalg.det(A)),
        "response": R.tolist(),
        "schur": K.tolist(),
        "note": (
            "numbers are algebra witnesses only; physical matrix elements "
            "must come from the retained action on the closed-FLRW n=2 domain"
        ),
    }


if __name__ == "__main__":
    print(json.dumps(diagnostics(), indent=2, sort_keys=True))
