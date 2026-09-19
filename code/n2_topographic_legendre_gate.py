from __future__ import annotations

import json


def transverse_eigenvalue(
    kappa1: float,
    X: float,
    w: float = 1.0,
) -> float:
    return w * (kappa1 + X**3)


def parallel_eigenvalue(
    kappa1: float,
    X: float,
    D0eta2: float,
    w: float = 1.0,
) -> float:
    return w * (
        kappa1
        + X**3
        - 6.0 * X**2 * D0eta2
    )


def diagnostics() -> dict:
    # Representative algebra witnesses only.
    cases = {
        "stationary": dict(kappa1=1.0, X=1.0, D0eta2=0.0),
        "positive_dynamic": dict(
            kappa1=1.0,
            X=1.0,
            D0eta2=0.05,
        ),
        "negative_dynamic": dict(
            kappa1=1.0,
            X=1.0,
            D0eta2=0.40,
        ),
    }

    out = {}
    for name, p in cases.items():
        l_perp = transverse_eigenvalue(p["kappa1"], p["X"])
        l_par = parallel_eigenvalue(
            p["kappa1"],
            p["X"],
            p["D0eta2"],
        )
        out[name] = {
            **p,
            "lambda_perp": l_perp,
            "lambda_parallel": l_par,
            "legendre_positive": (
                l_perp > 0.0 and l_par > 0.0
            ),
        }

    return {
        "status": "BHSM_TOPOGRAPHIC_LEGENDRE_GATE_RECORDED",
        "positive_condition": (
            "kappa1+X^3-6 X^2 |D0 eta|^2 > 0"
        ),
        "stationary_zero_momentum_consequence": (
            "D0 eta=0 and J_i^eta=0 on retained "
            "X>=0,kappa1>0,w>0 branch"
        ),
        "cases": out,
    }


if __name__ == "__main__":
    print(json.dumps(diagnostics(), indent=2, sort_keys=True))
