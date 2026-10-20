"""Checks against the original action equations, not optical witness fits."""
from pathlib import Path
import sys

import numpy as np
import pytest

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT/"code"))

import action_native_matter_n2 as an
import strong_lensing_envelope as env


@pytest.fixture(scope="module")
def background():
    return an.frozen_background()


@pytest.mark.parametrize("z", [0., .3, 1., 2.1])
def test_sequential_fluid_elimination_matches_original_constraints(background, z):
    row = env.adm_projection_snapshot(z, background)
    assert row["original_constraints_max"] < 2e-12
    assert row["simultaneous_aux_reconstruction_max_norm"] < 2e-12
    assert row["analytic_vs_projected_max"] < 2e-12
    assert row["response_residual_norm"] < 2e-12
    assert row["explicit_inverse_vs_solve_norm"] < 2e-12
    assert row["schur_vs_native_reduction_norm"] < 2e-12
    assert row["legendre_gate"]


def test_projected_hessian_is_same_polynomial_action(background):
    N = -np.log1p(.3)
    y = background.sol(N)
    blocks = env.adm_projected_blocks(N, y)
    # Nonzero, mixed-sign phase-space probe; fluid auxiliaries are solved,
    # while metric auxiliaries are left arbitrary to check the whole Hessian.
    q = np.array([.13, -.7, .4, .2, -.31, .19])
    metric = np.array([.11, -.23])
    _, _, k, _, _ = an.scales(N)
    u = np.r_[metric, metric[1]+q[1:3]/k]
    direct = an.unreduced_lagrangian(N, y, u, q[3:], q[:3])
    projected = (.5*metric@blocks["metric"]@metric
                 +metric@blocks["mixed"]@q+.5*q@blocks["state"]@q)
    np.testing.assert_allclose(direct, projected, rtol=1e-13, atol=1e-13)


def test_static_spatial_phase_null_does_not_imply_full_lensing_null():
    row = env.static_spatial_phase_null()
    for application in row["applications"].values():
        assert application["u"] == application["v"] == application["source_norm"] == 0.
    assert row["schur_curvature"].startswith("c_phi")
    assert "alone do not imply" in row["not_inferred"]


def test_snapshot_keeps_native_and_phase_domains_distinct(background):
    row = env.checkpoint_snapshot(background=background)
    assert row["adm_native_projection"]["metric_coordinates"] == ["alpha_ADM_lapse", "chi_ADM_shift"]
    assert not row["muon_input_consumed"]
    assert not row["realized_physical_envelope_amplitude_claimed"]
    with pytest.raises(ValueError, match="audited R1 domain"):
        env.adm_projection_snapshot(2.2, background)
