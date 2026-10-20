"""Focused physical-interface checks for the axial first comparison."""
import math

import numpy as np
import pytest
from scipy.integrate import quad

import action_native_matter_n2 as an
import closed_horndeski_n2 as owner
import strong_lensing_transport as slt


@pytest.fixture(scope="module")
def transport():
    return slt.ActionNativeAxialTransport()


@pytest.fixture(scope="module")
def unit_response(transport):
    return transport.response()


def test_dimensional_s3_profile_normalization_and_eigenvalue(transport):
    K, norm = transport.K, transport.normalization
    root = math.sqrt(K)
    # Angular integral cos(theta)^2 dOmega = 4pi/3. This integral uses
    # dimensional dV=S_K(chi)^2 dchi dOmega, not the unit-radius volume.
    integral = (4*math.pi/3)*quad(
        lambda c: (norm*math.sin(2*root*c))**2*transport.S(c)**2,
        0., math.pi/root, epsabs=1e-12)[0]
    assert abs(integral-1) < 2e-13
    for c in (.2, 1., 4.):
        q = transport.q(c)
        qp = transport.q_prime(c)
        # Radial scalar Laplacian for angular ell=1.
        lap = -4*K*q+2*math.cos(root*c)/transport.S(c)*qp-2*q/transport.S(c)**2
        assert abs(lap+8*K*q) < 1e-14


def test_all_six_physical_anchor_rows_and_original_constraints(transport):
    N = an.N_START
    y = transport.background.sol(N)
    a = math.exp(N)
    vel, aux = an.canonical_maps(N, y)
    for i in range(6):
        state = transport.anchor_to_canonical[:, i]
        xd, u = vel@state, aux@state
        actual = [state[0], 2*a**3*owner.gs_n2(y[1], y[2])*xd[0],
                  state[1]-3*state[0], (4/3)*(state[2]-3*state[0]), u[2], u[3]]
        np.testing.assert_allclose(actual, np.eye(6)[i], atol=2e-12, rtol=0.)
        assert max(abs(an.constraint_residuals(N, y, state[:3], xd, u))) < 3e-11
    # Historical two columns are recovered, rather than confused with the
    # full six-dimensional matter/radiation state.
    np.testing.assert_allclose(transport.anchor_to_canonical[:, :2],
                               an.topographic_seed(transport.background), atol=1e-12, rtol=0.)
    # Although the coefficient of Pi2 happens to equal one on this slice,
    # p_zeta also depends on the other five physical coordinates.
    assert np.linalg.norm(transport.anchor_to_canonical[3]-np.eye(6)[1]) > .1


def test_background_three_maps_match_closed_flrw_and_not_subtraction(transport):
    maps = transport.maps()
    cl, cs = maps.lens["chi"], maps.source["chi"]
    expected = np.array([transport.S(cl)/1.3, transport.S(cs)/2,
                         transport.S(cs-cl)/2])
    np.testing.assert_allclose(maps.scalar_maps, expected, rtol=2e-10, atol=2e-12)
    assert abs(maps.D_ls[0, 0]-(maps.D_os[0, 0]-maps.D_ol[0, 0])) > .05
    assert abs(maps.lens["observed_z"]-.3) < 2e-12
    assert abs(maps.source["observed_z"]-1.) < 2e-12


def test_reverse_reciprocity_all_segments_and_t_symmetry(transport):
    # This includes nonzero observer velocity and moving measured-z events.
    maps = transport.maps(np.array([1e-3, 3e-5, -2e-4, 1e-3, 4e-5, -1e-3]))
    for key in ("relative_residual", "relative_residual_ol", "relative_residual_os"):
        assert abs(maps.reciprocity[key]) < 3e-10
    T = 1.3*maps.D_ol.T@np.linalg.solve(maps.D_ls, maps.D_os)
    np.testing.assert_array_equal(T, T.T)
    assert abs(maps.observer["omega"]-1) < 2e-14
    assert abs(maps.lens["omega"]-1.3) < 2e-12
    assert abs(maps.source["omega"]-2) < 2e-12
    assert abs(maps.observer["peculiar_velocity"]) > 1e-8
    # Holding measured z does not hold either the background source N or
    # the lens's affine coordinate fixed.
    assert abs(maps.lens["N"]+math.log(1.3)) > 1e-9
    assert abs(maps.source["N"]+math.log(2.)) > 1e-9


def test_six_unit_responses_stable_under_independent_half_step(transport, unit_response):
    full, half = unit_response, transport.response(step=1e-3)
    first, second = np.asarray(full["d_maps_d_anchor"]), np.asarray(half["d_maps_d_anchor"])
    assert first.shape == (3, 6)
    assert np.all(np.linalg.norm(first, axis=0) > 1e-9)
    np.testing.assert_allclose(first, second, rtol=3e-5, atol=1e-10)
    events = np.asarray(full["d_events_d_anchor"])
    assert events.shape == (4, 6)
    assert np.linalg.norm(events[0]) > .01
    for trial in full["finite_trials"]:
        for sign in ("plus", "minus"):
            assert abs(trial[sign+"_endpoints"]["lens"]["observed_z"]-.3) < 2e-12
            assert abs(trial[sign+"_endpoints"]["source"]["observed_z"]-1) < 2e-12
            assert abs(trial[sign+"_reciprocity"]["relative_residual"]) < 1e-10
            assert trial[sign+"_perturbation_diagnostics"]["max_relative_metric_perturbation"] < 1e-3


def test_time_interpolation_grid_refinement_for_mixed_six_coordinate_trial(transport):
    refined = slt.ActionNativeAxialTransport(transport.background, transport.propagator,
                                             grid_size=1601)
    anchor = np.array([.001, 2e-5, .0004, .001, -2e-5, -.001])
    coarse, fine = transport.maps(anchor), refined.maps(anchor)
    np.testing.assert_allclose(coarse.scalar_maps, fine.scalar_maps, rtol=0., atol=2e-11)
    assert coarse.perturbation_diagnostics["max_relative_metric_perturbation"] > 1e-5
    assert coarse.perturbation_diagnostics["max_relative_metric_perturbation"] < 1e-3


def test_inherited_seven_state_jacobi_control_matches_full_metric_maps(transport, unit_response):
    control = transport.seven_state_control(unit_response)
    assert control["max_absolute_difference"] < 3e-10
    np.testing.assert_allclose(control["d_maps_d_anchor"],
                               control["full_metric_d_maps_d_anchor"], rtol=2e-5, atol=2e-10)
    for name in ("lens_seven_state_unit_columns", "source_seven_state_unit_columns"):
        states = np.array(control[name])
        assert states.shape == (7, 6)
        # These slots are unused frequency and axial shear controls, not a
        # statement that the physical frequency/redshift response vanishes.
        np.testing.assert_array_equal(states[[0, 3, 4, 5, 6]], 0.)
    assert max(abs(np.array(control["initial_kappa"]))) > .01


def test_supported_domain_is_explicit(transport):
    with pytest.raises(ValueError, match="six"):
        transport.maps(np.ones(2))
    with pytest.raises(ValueError, match="domain"):
        transport.maps(z_l=.3, z_s=2.1)
