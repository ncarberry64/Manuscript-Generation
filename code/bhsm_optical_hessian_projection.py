from __future__ import annotations

from dataclasses import dataclass
import json
from pathlib import Path

import numpy as np


@dataclass(frozen=True)
class LocalOpticalCoefficients:
    a: float
    b: float
    c: float
    u: float
    v: float
    c_phi: float | None = None

    @property
    def metric_matrix(self) -> np.ndarray:
        return np.array(
            [[self.a, self.b], [self.b, self.c]],
            dtype=float,
        )

    @property
    def mixed_column(self) -> np.ndarray:
        return np.array([self.u, self.v], dtype=float)

    @property
    def determinant(self) -> float:
        return float(self.a * self.c - self.b * self.b)

    def response(self) -> np.ndarray:
        A = self.metric_matrix
        B = self.mixed_column
        if abs(self.determinant) < 1.0e-14:
            raise ValueError("singular physical metric block")
        return -np.linalg.solve(A, B)

    def response_residual(self) -> float:
        R = self.response()
        return float(
            np.linalg.norm(
                self.metric_matrix @ R + self.mixed_column
            )
        )

    def weyl_response(self) -> float:
        R = self.response()
        return float(0.5 * (R[0] + R[1]))

    def slip_response(self) -> float:
        R = self.response()
        return float(0.5 * (R[0] - R[1]))

    def schur_curvature(self) -> float | None:
        if self.c_phi is None:
            return None
        B = self.mixed_column
        return float(
            self.c_phi
            - B.T @ np.linalg.solve(self.metric_matrix, B)
        )


def project_symmetric_block(
    hessian: np.ndarray,
    basis: np.ndarray,
) -> np.ndarray:
    """Project H to columns of basis: basis.T H basis."""
    H = np.asarray(hessian, dtype=float)
    P = np.asarray(basis, dtype=float)
    if H.ndim != 2 or H.shape[0] != H.shape[1]:
        raise ValueError("hessian must be square")
    if P.ndim != 2 or P.shape[0] != H.shape[0]:
        raise ValueError("basis has incompatible shape")
    out = P.T @ H @ P
    return 0.5 * (out + out.T)


def project_mixed_column(
    mixed_operator: np.ndarray,
    metric_basis: np.ndarray,
    topographic_profile: np.ndarray,
) -> np.ndarray:
    """P_g.T H_gphi phi for finite discretizations."""
    H = np.asarray(mixed_operator, dtype=float)
    P = np.asarray(metric_basis, dtype=float)
    phi = np.asarray(topographic_profile, dtype=float)
    if H.ndim != 2:
        raise ValueError("mixed_operator must be a matrix")
    if P.ndim != 2 or P.shape[0] != H.shape[0]:
        raise ValueError("metric basis has incompatible shape")
    if phi.ndim != 1 or phi.shape[0] != H.shape[1]:
        raise ValueError("topographic profile has incompatible shape")
    return P.T @ H @ phi


def coefficients_from_projected(
    A: np.ndarray,
    B: np.ndarray,
    c_phi: float | None = None,
) -> LocalOpticalCoefficients:
    A = np.asarray(A, dtype=float)
    B = np.asarray(B, dtype=float)
    if A.shape != (2, 2):
        raise ValueError("A must be 2x2")
    if B.shape != (2,):
        raise ValueError("B must be length 2")
    if not np.allclose(A, A.T, rtol=0.0, atol=1e-12):
        raise ValueError("A must be symmetric")
    return LocalOpticalCoefficients(
        a=float(A[0,0]),
        b=float(A[0,1]),
        c=float(A[1,1]),
        u=float(B[0]),
        v=float(B[1]),
        c_phi=None if c_phi is None else float(c_phi),
    )


def basis_transform(
    coeff: LocalOpticalCoefficients,
    S: np.ndarray,
) -> tuple[np.ndarray, np.ndarray]:
    """Coordinate change P -> P S."""
    S = np.asarray(S, dtype=float)
    if S.shape != (2,2):
        raise ValueError("S must be 2x2")
    if abs(np.linalg.det(S)) < 1e-14:
        raise ValueError("S must be invertible")
    A2 = S.T @ coeff.metric_matrix @ S
    B2 = S.T @ coeff.mixed_column
    return A2, B2


def validate_basis_covariance(
    coeff: LocalOpticalCoefficients,
    S: np.ndarray,
) -> float:
    """Return norm of R_new - S^-1 R_old."""
    A2, B2 = basis_transform(coeff, S)
    c2 = coefficients_from_projected(A2, B2, coeff.c_phi)
    lhs = c2.response()
    rhs = np.linalg.solve(np.asarray(S, dtype=float), coeff.response())
    return float(np.linalg.norm(lhs - rhs))


def write_contract_template(path: str | Path) -> None:
    template = {
        "schema": "bhsm-local-optical-hessian-coefficients-v1",
        "status": "AWAITING_ACTION_EVALUATION",
        "physical_domain": {
            "background": "closed-FLRW R1 physical background",
            "metric_basis": ["A", "psi"],
            "constraint_reduction": "momentum constraint solved",
            "scalar_longitudinal_gauge": "E=0",
            "endpoint_domain": "action-derived physical/KKT domain",
        },
        "normalizations": {
            "metric_basis_A": None,
            "metric_basis_psi": None,
            "topographic_profile": None,
        },
        "coefficients": {
            "a_AA": None,
            "b_Apsi": None,
            "c_psipsi": None,
            "u_Aphi": None,
            "v_psiphi": None,
            "c_phi": None,
        },
        "required_diagnostics": {
            "Delta_g": None,
            "metric_response_residual": None,
            "schur_curvature": None,
            "legendre_gate": None,
            "common_domain": None,
        },
        "forbidden": [
            "fit any coefficient to SN residuals",
            "substitute algebra-witness coefficients",
            "mix coefficients from different backgrounds or domains",
        ],
    }
    Path(path).write_text(
        json.dumps(template, indent=2, sort_keys=True) + "\n",
        encoding="utf-8",
    )
