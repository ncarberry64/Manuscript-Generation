import sys
from pathlib import Path

import numpy as np
import pytest

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT/"code"))
from strong_lensing_degeneracy import (
    factorization_check, saved_array_receipt, scalar_map_coefficients,
)


@pytest.fixture(scope="module")
def saved():
    directory = ROOT/"artifacts/strong_lensing_first_comparison_20261020"
    with np.load(directory/"response_arrays.npz") as arrays:
        values = {key: arrays[key].copy() for key in arrays.files}
    return saved_array_receipt(directory), values


@pytest.mark.parametrize("family,expected", [("SIS", 2.452944603865752e-12),
                                             ("smooth_core", 6.377724303772925e-12)])
def test_saved_observable_responses_factor_with_minimal_physical_nuisances(saved, family, expected):
    receipt, _ = saved
    check = receipt["families"][family]
    assert check["relative_whitened_factorization_discrepancy"] == pytest.approx(expected, rel=1e-4)
    assert check["relative_whitened_projected_response"] < 1e-10
    assert check["minimal_nuisance_rank"] == (4 if family == "smooth_core" else 3)
    assert "log_H0_fixed_physical_source_and_galaxy" not in check["nuisance_names"]
    assert "independent_external_convergence_nuisance_only" not in check["nuisance_names"]


@pytest.mark.parametrize("family,prefix", [("SIS", ""), ("smooth_core", "smooth_")])
def test_null_survives_correlated_observation_coordinates_but_needs_delay_scale(saved, family, prefix):
    receipt, arrays = saved
    check = receipt["families"][family]
    B = arrays[prefix+"B"]
    N = arrays[prefix+"N"][:, check["nuisance_indices"]]
    R = np.array(check["coefficient_matrix"])
    size = len(B)
    # Off-diagonal SPD covariance verifies that the identity is independent
    # of the illustrative diagonal errors in the original checkpoint.
    mixing = np.eye(size)+.07*np.tril(np.ones((size, size)), -1)
    covariance = mixing@mixing.T
    result = factorization_check(B, N, R, covariance)
    assert result["relative_whitened_projected_response"] < 1e-10
    removed = factorization_check(B, N[:, :-1], R[:-1], covariance)
    assert removed["relative_whitened_projected_response"] > 1e-4


def test_endpoint_receipt_retains_actual_radial_and_affine_displacements(saved):
    receipt, _ = saved
    endpoint = receipt["source_endpoint_unit_zeta"]
    assert endpoint["delta_chi_event"] == pytest.approx(-.005611965495588844)
    assert endpoint["delta_lambda_event"] == pytest.approx(.00002815169979075716)
    assert endpoint["converted_affine_displacement"] == pytest.approx(.00011260679916302863)
    assert abs(endpoint["delta_chi_event"]-endpoint["converted_affine_displacement"]) > .005
    assert endpoint["inferred_delta_chi_ray_fixed_lambda"] == pytest.approx(-.005724572294751873)


def test_invalid_coefficient_domain_and_nonpositive_covariance_are_rejected():
    with pytest.raises(ValueError, match="positive finite"):
        scalar_map_coefficients([1, 0, 2], np.zeros((3, 6)), .2, core=False)
    with pytest.raises(ValueError, match="three scalar distances"):
        scalar_map_coefficients(np.eye(3), np.zeros((3, 6)), .2, core=False)
    with pytest.raises(np.linalg.LinAlgError):
        factorization_check(np.ones((2, 1)), np.ones((2, 1)), np.ones((1, 1)), np.diag([1, -1]))
