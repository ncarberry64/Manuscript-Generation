"""Physical static action, exact vacuum constraints, and finite Dirichlet domain."""
from pathlib import Path
import sys

import numpy as np
import pytest
from scipy.integrate import quad

sys.path.insert(0, str(Path(__file__).resolve().parents[1]/"code"))
from strong_lensing_local_hessian import (dirichlet_profile, evaluate_local_hessian,
                                         exterior_reference)


def test_flat_limit_has_derived_Einstein_block_and_positive_phase_curvature():
    row = evaluate_local_hessian(0.0)
    expected = np.pi**2
    np.testing.assert_allclose(row["metric_hessian_over_KR"],
                               [[0.0, 2*expected], [2*expected, -2*expected]], atol=1e-12)
    assert row["c_phi_over_k_eta_R"] == pytest.approx(expected, rel=1e-13)
    assert row["profile_norm"] == pytest.approx(1.0, abs=5e-14)


@pytest.mark.parametrize("compactness", [1e-6, 0.02, 0.15])
def test_curved_original_euler_constraints_null_response_and_positive_schur(compactness):
    row = evaluate_local_hessian(compactness)
    assert row["original_metric_Euler_residual_max_in_R"] < 1e-12
    assert abs(row["a_over_KR"]) < 1e-12
    assert row["determinant_over_KR_squared"] < 0.0
    assert row["schur_over_k_eta_R"] > 0.0
    assert row["metric_response_residual"] == 0.0
    assert row["vacuum_simplified_block_residual_over_KR"] < 2e-12
    assert row["metric_response"] == [0.0, 0.0]
    assert row["eta_boundary_conormal_residual"] == 0.0


def test_metric_block_is_second_variation_of_original_exact_static_energy():
    mu = 0.15  # appreciable curvature makes this more than a flat-limit check
    row = evaluate_local_hessian(mu)
    inner = row["domain_in_R"]["isotropic_inner"]
    outer = row["domain_in_R"]["isotropic_outer"]
    direction = np.asarray([0.31, -0.27])

    def energy(t):
        def density(r):
            A, P, Ap, Pp, _, _ = exterior_reference(r, mu)
            f, fp = dirichlet_profile(r, inner, outer)
            A += t*direction[0]*f
            P += t*direction[1]*f
            Ap += t*direction[0]*fp
            Pp += t*direction[1]*fp
            return 4*np.pi*r*r*np.exp(A-P)*(2*Ap*Pp-Pp*Pp)
        return quad(density, inner, outer, epsabs=1e-13, epsrel=1e-13)[0]

    epsilon = 3e-4
    finite = (energy(epsilon)-2*energy(0)+energy(-epsilon))/epsilon**2
    exact = direction@np.asarray(row["metric_hessian_over_KR"])@direction
    assert finite == pytest.approx(exact, rel=2e-8, abs=2e-8)


def test_no_horizon_or_empty_projection_domain_is_accepted():
    with pytest.raises(ValueError):
        evaluate_local_hessian(0.5)
    with pytest.raises(ValueError):
        evaluate_local_hessian(0.1, 0.9)
