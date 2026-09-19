from __future__ import annotations

import math
import sys
from pathlib import Path

import numpy as np
import pytest

ROOT = Path(__file__).resolve().parents[1]
CODE = ROOT / "code"
if str(CODE) not in sys.path:
    sys.path.insert(0, str(CODE))

import background
import bianchi_dictionary_n2 as bd


def _background():
    v0 = background.shoot_v0()
    return background.integrate(v0)


def test_adm_newtonian_dictionary_round_trip():
    values = dict(
        delta_n=0.013,
        chi=-0.027,
        chi_dot=0.008,
        zeta=0.021,
        H=1.17,
    )
    Phi, Psi, pi = bd.adm_to_newtonian(**values)
    recovered = bd.newtonian_to_adm(
        Phi=Phi,
        Psi=Psi,
        pi=pi,
        pi_dot=-values["chi_dot"],
        H=values["H"],
    )
    target = (
        values["delta_n"],
        values["chi"],
        values["chi_dot"],
        values["zeta"],
    )
    np.testing.assert_allclose(
        recovered, target, rtol=0.0, atol=2e-15
    )


def test_acceleration_matrix_determinant_is_2H2Dkin():
    sol = _background()

    for z in [2.1, 1.5, 1.0, 0.8, 0.5, 0.3, 0.0]:
        N = math.log(1.0 / (1.0 + z))
        direct, analytic = bd.acceleration_determinant_identity(
            sol, N
        )
        assert math.isclose(
            direct, analytic, rel_tol=3e-11, abs_tol=3e-12
        )
        assert analytic > 0.0


def test_initial_state_is_constraint_and_tangency_consistent():
    sol = _background()
    N0, y0 = bd.constraint_consistent_initial_state(sol)

    rn = bd.normalized_constraint_residuals(sol, N0, y0)
    D = bd.constraint_tangency_operator(sol, N0)
    tangent = D @ y0

    assert rn[0] < 2e-10
    assert rn[1] < 2e-10
    assert np.linalg.norm(tangent) / np.linalg.norm(y0) < 2e-7


@pytest.mark.xfail(strict=True, reason="Known curvature-dependent EFT compatibility defect; BHSM defining transport uses the action-reduced physical n=2 operator.")
def test_dae_preserves_00_0i_constraints_from_index_consistent_seed():
    sol = _background()
    pert = bd.integrate_dae(sol)
    report = bd.constraint_preservation_report(sol, pert)

    assert report["max_constraint_00_normalized"] < 2e-5
    assert report["max_constraint_0i_normalized"] < 2e-5


def test_topographic_response_remains_off():
    text = Path(bd.__file__).read_text(encoding="utf-8-sig")
    assert "topographic_response" in text
    assert '"off"' in text


def test_gp_eq12_alphaB_derivative_source_term():
    """Regression for GP Eq. (12), arXiv:2412.01781v1.

    The published equation contains
        H * alpha_B_dot * (3 H^2 + Hdot)
        + H * d/dt(H * alpha_B_dot).

    This test guards the literal source transcription.
    """
    sol = _background()
    N = -0.3
    bg = bd.background_derivatives(sol, N)

    H = bg["H"]
    Hdot = bg["Hdot"]

    expected_middle = H * bg["aBdot"] * (3.0 * H * H + Hdot)

    h = 2.0e-5

    def h_abdot(n):
        return (
            bd._H(sol, n)
            * (
                bd._H(sol, n)
                * (
                    bd._aB(sol, n + h)
                    - bd._aB(sol, n - h)
                )
                / (2.0 * h)
            )
        )

    direct = H * (h_abdot(N + h) - h_abdot(N - h)) / (2.0 * h)

    assert np.isclose(
        bg["H_aBdot_dot"],
        direct,
        rtol=2e-6,
        atol=2e-9,
    )

    assert np.isfinite(expected_middle)
