from pathlib import Path
import sys
import hashlib
import numpy as np
import pytest
from scipy.stats import chi2

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "code"))
import bhsm_make_or_break_prediction as mb


def test_shared_state_known_projection_and_basis_invariance():
    F = np.array([[1., 0.], [0., 1.], [0., 0.], [0., 0.]])
    S = np.diag([2., 3., 4., 5.])
    d = np.array([2., -3., 4., 5.])
    result = mb.shared_state_test(d, F, S)
    assert result["chi2_perp"] == pytest.approx(9.)
    assert result["p_value"] == pytest.approx(chi2.sf(9., 2))
    transform = np.array([[2., 1.], [1., 3.]])
    rotated = mb.shared_state_test(d, F @ transform, S)
    assert rotated["chi2_perp"] == pytest.approx(result["chi2_perp"])
    np.testing.assert_allclose(F @ result["state"], F @ transform @ rotated["state"])


def test_correlated_covariance_and_incompatible_extra_observable():
    F = np.array([[1., 0.], [0., 1.], [1., 1.]])
    S = np.array([[2., .2, .3], [.2, 1., .1], [.3, .1, 3.]])
    assert mb.shared_state_test(F @ [2., 3.], F, S)["chi2_perp"] < 1e-25
    assert mb.shared_state_test([2., 3., 100.], F, S)["reject"]


@pytest.mark.parametrize("F,S", [
    (np.ones((3, 2)), np.eye(3)),
    (np.ones((3, 3)), np.eye(3)),
    (np.array([[1., 0.], [0., 1.], [1., 1.]]), -np.eye(3)),
    (np.array([[1., 0.], [0., 1.], [1., 1.]]), np.ones((3, 3))),
])
def test_invalid_design_or_covariance_fails_closed(F, S):
    with pytest.raises(ValueError):
        mb.shared_state_test(np.ones(3), F, S)


def test_unprojected_negative_modes_and_zero_signal():
    d = mb.covariance_diagnostics(np.diag([1., -.2, -.3]), np.ones(3))
    assert d["rho21"] == -.2
    assert d["principal_minors"][0] == -.2
    assert mb.covariance_diagnostics(np.zeros((2, 2)), np.ones(2))["rho21"] is None
    v = np.array([1., 2., -3.])
    assert max(abs(x) for x in mb.covariance_diagnostics(np.outer(v, v), np.ones(3))["principal_minors"]) == 0


def test_unit_scaling_invariance():
    C = np.array([[2., .3], [.3, -.1]])
    units = np.array([100., .01])
    first = mb.covariance_diagnostics(C, np.ones(2))
    second = mb.covariance_diagnostics(C * np.outer(units, units), units)
    np.testing.assert_allclose(first["eigenvalues_descending"], second["eigenvalues_descending"])


def test_mc_tail_never_zero_and_strict_threshold():
    assert mb.monte_carlo_pvalue(2., np.zeros(999)) == .001
    assert mb.monte_carlo_pvalue(0., np.zeros(999)) == 1.
    with pytest.raises(ValueError):
        mb.monte_carlo_pvalue(1., np.zeros(10))


def test_fixed_null_simulator_detects_strong_second_mode_reproducibly():
    rng = np.random.default_rng(78)
    samples = rng.standard_normal((100, 3)) * [1., 5., 4.]
    result = mb.fixed_gaussian_covariance_test(samples, np.eye(3)*.01, [1., 0., 0.], np.ones(3), draws=999)
    assert result["reject"]
    assert result["p_value"] == .001
    assert not result["composite_rank_null_calibrated"]


def test_one_local_noise_type_can_generate_rank_two_after_path_integration():
    # Independent increments at two distances, same noise type, different response.
    response1, response2 = np.array([1., 0.]), np.array([0., 1.])
    integrated = np.outer(response1, response1) + np.outer(response2, response2)
    assert np.linalg.matrix_rank(integrated) == 2


def test_nonfinite_and_asymmetric_inputs_rejected():
    with pytest.raises(ValueError):
        mb.covariance_diagnostics([[1., .2], [.3, 1.]], [1., 1.])
    with pytest.raises(ValueError):
        mb.covariance_diagnostics([[1., 0.], [0., np.nan]], [1., 1.])
    with pytest.raises(ValueError):
        mb.covariance_diagnostics(np.eye(2), [1., 0.])


def test_structural_protocol_hash():
    expected = (ROOT / "preregistration/bhsm_geometry_first_make_or_break_v1.sha256").read_text().split()[0]
    assert hashlib.sha256(mb.PROTOCOL.read_bytes()).hexdigest() == expected
