from pathlib import Path
import sys

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
CODE = ROOT / "code"
if str(CODE) not in sys.path:
    sys.path.insert(0, str(CODE))

import r1_n2_reduced_propagator as rp


def test_anchor_is_identity():
    sol = rp.integrate_propagator()
    U = rp.transfer_matrix(sol, rp.Z_START)
    assert np.allclose(U, np.eye(2), atol=1e-12)


def test_reduced_coefficients_positive_on_prediction_interval():
    scan = rp.stability_scan(120)
    assert scan["min_G_S_n2"] > 0.0
    assert scan["min_F_S_n2_gravity_scalar"] > 0.0
    assert scan["min_omega2"] > 0.0


def test_closed_n2_constraints_remain_exact():
    sol = rp.integrate_propagator()
    assert rp.constraint_residual_max(sol, 30) < 1e-10


def test_liouville_identity():
    sol = rp.integrate_propagator()
    assert rp.liouville_residual(sol) < 1e-8
