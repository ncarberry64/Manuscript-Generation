from __future__ import annotations

import math
import sys
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
CODE = ROOT / "code"
if str(CODE) not in sys.path:
    sys.path.insert(0, str(CODE))

import background
import closed_horndeski_n2 as ch
import matter_sourced_n2 as ms


def _background():
    v0 = background.shoot_v0()
    return background.integrate(v0)


def test_constraint_determinant_identity_and_invertibility():
    sol = _background()

    for z in [10.0, 3.0, 2.1, 1.0, 0.5, 0.0]:
        N = math.log(1.0 / (1.0 + z))
        y = np.asarray(sol.sol(N), dtype=float)
        _, v, H = [float(x) for x in y]
        direct, analytic = ms.determinant_identity(v, H)
        q = ch.alpha_quantities(v, H)

        assert math.isclose(
            direct, analytic, rel_tol=2e-12, abs_tol=2e-12
        )
        assert q["D_kin"] > 0.0
        assert abs(direct) > 0.0


def test_sourced_constraint_inverse_closes_exactly():
    sol = _background()
    probe = ms.SourceState(
        Phi=-1.7e-5,
        pi=2.1e-4,
        delta_m=3.4e-5,
        v_m=-8.0e-5,
        delta_r=4.6e-5,
        v_r=3.0e-5,
    )

    for z in [3.0, 2.1, 1.0, 0.3, 0.0]:
        N = math.log(1.0 / (1.0 + z))
        y = np.asarray(sol.sol(N), dtype=float)
        pdot, pidot = ms.solve_metric_scalar_rates(N, y, probe)
        residual = ms.constraint_residuals(
            N, y, probe, pdot, pidot
        )
        np.testing.assert_allclose(
            residual, np.zeros(2), rtol=0.0, atol=2e-11
        )


def test_no_slip_for_ideal_sources():
    assert ms.no_slip_residual(0.123, 0.123, 0.0) == 0.0
    assert math.isclose(
        ms.no_slip_residual(0.12, 0.10, -0.02),
        0.0,
        abs_tol=1e-15,
    )


def test_local_transfer_matrix_is_linear():
    sol = _background()
    N = math.log(1.0 / 2.0)
    y = np.asarray(sol.sol(N), dtype=float)
    T = ms.local_transfer_matrix(N, y)

    x1 = np.array([1e-5, -2e-4, 3e-5, 4e-5, -5e-5, 6e-5])
    x2 = np.array([-7e-6, 1e-4, 2e-5, -3e-5, 4e-5, 2e-5])

    def rates(x):
        s = ms.SourceState(*x)
        Phi_dot, pi_dot = ms.solve_metric_scalar_rates(N, y, s)
        f = ms.ideal_fluid_rates(N, y, s, Phi_dot)
        return np.array([
            Phi_dot,
            pi_dot,
            f["delta_m_dot"],
            f["v_m_dot"],
            f["delta_r_dot"],
            f["v_r_dot"],
        ])

    np.testing.assert_allclose(
        T @ x1, rates(x1), rtol=1e-12, atol=1e-13
    )
    np.testing.assert_allclose(
        T @ (x1 + x2),
        rates(x1) + rates(x2),
        rtol=1e-12,
        atol=1e-13,
    )


def test_topographic_operator_is_not_silently_present():
    text = Path(ms.__file__).read_text(encoding="utf-8-sig")
    assert "topographic reduced-response operator is OFF" in text
