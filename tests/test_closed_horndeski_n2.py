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


def test_theta_alpha_identity():
    v0 = background.shoot_v0()
    sol = background.integrate(v0)

    for z in [10.0, 3.0, 2.1, 1.0, 0.5, 0.0]:
        N = math.log(1.0 / (1.0 + z))
        _, v, E = [float(x) for x in sol.sol(N)]
        theta, _ = ch.theta_sigma(v, E)
        alpha_b = ch.alpha_quantities(v, E)["alpha_B"]
        expected = E * (1.0 - 0.5 * alpha_b)
        assert math.isclose(theta, expected, rel_tol=0.0, abs_tol=2e-13)


def test_closed_n2_kinetic_identity_and_positivity():
    v0 = background.shoot_v0()
    sol = background.integrate(v0)

    for z in [10.0, 3.0, 2.1, 1.5, 1.0, 0.8, 0.5, 0.3, 0.0]:
        N = math.log(1.0 / (1.0 + z))
        _, v, E = [float(x) for x in sol.sol(N)]
        g1 = ch.gs_n2(v, E)
        g2 = ch.gs_n2_from_dkin(v, E)
        assert math.isclose(g1, g2, rel_tol=2e-12, abs_tol=2e-13)
        assert g1 > 0.0


def test_n2_constraint_solution_closes_original_equations():
    v0 = background.shoot_v0()
    sol = background.integrate(v0)

    probes = [
        (2.1, 0.17, -0.031),
        (1.0, -0.08, 0.047),
        (0.3, 0.11, 0.019),
        (0.0, -0.03, -0.021),
    ]

    for z, zeta, zeta_dot in probes:
        N = math.log(1.0 / (1.0 + z))
        _, v, E = [float(x) for x in sol.sol(N)]
        a = math.exp(N)
        r1, r2 = ch.constraint_residuals_n2(
            a, v, E, zeta, zeta_dot
        )
        scale = max(1.0, abs(zeta), abs(zeta_dot))
        assert abs(r1) < 2e-11 * scale
        assert abs(r2) < 2e-11 * scale


def test_formal_gravity_sector_fs_is_finite_on_reference_samples():
    v0 = background.shoot_v0()
    sol = background.integrate(v0)

    for z in [3.0, 2.1, 1.5, 1.0, 0.8, 0.5, 0.3, 0.0]:
        N = math.log(1.0 / (1.0 + z))
        y = np.asarray(sol.sol(N), dtype=float)
        fs = ch.formal_fs_n2_gravity_sector(N, y)
        assert math.isfinite(fs)
