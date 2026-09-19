from __future__ import annotations

import json
import math
import numpy as np


def candidate_matrix(
    lambda_rho: float = 0.0,
    I2: complex = 1.0,
    L: float = 1.0,
) -> np.ndarray:
    """Canonical harmonic-preserving n=2 candidate.

    Assumptions:
    - lowest Robin zero-flux collar radial mode -> lambda_rho = 0;
    - normalized harmonic-preserving trace/extension candidate -> I2 = 1;
    - physical n=2 S^3 eigenvalue -> -Delta_B = 8.
    """
    if L <= 0:
        raise ValueError("L must be positive")
    return np.array(
        [
            [lambda_rho + 5.0, -I2],
            [-np.conjugate(I2), 13.0],
        ],
        dtype=complex,
    ) / (L * L)


def positivity_margin(
    lambda_rho: float = 0.0,
    I2: complex = 1.0,
) -> float:
    return 13.0 * (lambda_rho + 5.0) - abs(I2) ** 2


def diagnostics() -> dict:
    M = candidate_matrix()
    eig = np.linalg.eigvalsh(M).real

    expected = np.array(
        [9.0 - math.sqrt(17.0), 9.0 + math.sqrt(17.0)]
    )

    return {
        "status": "CONDITIONAL_HARMONIC_PRESERVING_N2_ST_CANDIDATE",
        "lambda_rho": 0.0,
        "I2": 1.0,
        "matrix_L1": M.real.tolist(),
        "positivity_margin": positivity_margin(),
        "eigenvalues_L1": eig.tolist(),
        "expected_eigenvalues_L1": expected.tolist(),
        "positive": bool(np.min(eig) > 0.0),
        "claim_boundary": (
            "lambda_rho=0 and I2=1 are the canonical lowest-Robin/"
            "harmonic-preserving extension candidate, not yet a retained-"
            "action proof of the full nonhomogeneous trace/extension domain"
        ),
    }


def main() -> None:
    print(json.dumps(diagnostics(), indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
