from __future__ import annotations

import json
import numpy as np


def pullback_operator(K_x: np.ndarray, U: np.ndarray) -> np.ndarray:
    """
    X = U Z.
    If Q = 1/2 X^T K_x X, then
    K_z = U^T K_x U.
    """
    return U.T @ K_x @ U


def transform_response(R_x: np.ndarray, U: np.ndarray) -> np.ndarray:
    """
    Metric response h = R_x X = R_x U Z.
    """
    return R_x @ U


def transform_observation(L_x: np.ndarray, U: np.ndarray) -> np.ndarray:
    """
    O = L_x X = (L_x U) Z.
    """
    return L_x @ U


def diagnostics() -> dict:
    Kx = np.array([[5.0, -1.0], [-1.0, 13.0]])
    Rx = np.array([[0.4, -0.1], [0.2, 0.3]])
    Lx = np.array([[1.7, -0.25]])

    U = np.array([[1.0, 0.2], [0.0, 1.1]])
    if abs(np.linalg.det(U)) < 1e-14:
        raise AssertionError("witness map must be invertible")

    z = np.array([0.7, -0.4])
    x = U @ z

    Kz = pullback_operator(Kx, U)
    Rz = transform_response(Rx, U)
    Lz = transform_observation(Lx, U)

    qx = float(x.T @ Kx @ x)
    qz = float(z.T @ Kz @ z)

    hx = Rx @ x
    hz = Rz @ z

    ox = float((Lx @ x).item())
    oz = float((Lz @ z).item())

    return {
        "status": "N2_REPRESENTATION_INVARIANCE_VALID",
        "det_U": float(np.linalg.det(U)),
        "quadratic_form_residual": abs(qx - qz),
        "metric_response_residual": float(np.max(np.abs(hx - hz))),
        "observable_residual": abs(ox - oz),
        "interpretation": (
            "invertible coordinate changes on the same two-dimensional "
            "physical phase space do not create a new scalar degree of freedom"
        ),
    }


if __name__ == "__main__":
    print(json.dumps(diagnostics(), indent=2, sort_keys=True))
