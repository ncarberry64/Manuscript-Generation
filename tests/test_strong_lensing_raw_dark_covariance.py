"""Source-stage and science/reference cross-term regressions on controlled fields."""
import importlib.util
from pathlib import Path
import numpy as np
import pytest

spec = importlib.util.spec_from_file_location("dark_calibration",
    Path(__file__).resolve().parents[1]/"scripts/measure_he0435_raw_dark_covariance.py")
module = importlib.util.module_from_spec(spec)
spec.loader.exec_module(module)


def test_pair_difference_removes_fixed_mean_and_has_correct_read_variance():
    rng = np.random.default_rng(63)
    mean = np.linspace(100., 500., 100000)
    first = mean+rng.normal(0, 7., (3, 100000))
    second = mean+rng.normal(0, 7., (3, 100000))
    delta = module.paired_integration_difference(first, second)
    np.testing.assert_allclose(np.var(delta, axis=1), 49., rtol=.013)
    assert abs(delta.mean()) < .03


def test_shared_reference_cross_terms_are_not_independent_variance_additions():
    rng = np.random.default_rng(4)
    reference = rng.normal(size=(3, 20, 30))
    science = .7*reference+rng.normal(size=(3, 20, 30))
    correction = .5*reference
    post = science-correction
    good = np.ones((20, 30), bool)
    s = module.pooled_time_covariance(science, good)["zero_mean_second_moment"]
    c = module.pooled_time_covariance(correction, good)["zero_mean_second_moment"]
    p = module.pooled_time_covariance(post, good)["zero_mean_second_moment"]
    cross = science[:, good]@correction[:, good].T/good.sum()
    np.testing.assert_allclose(p, s+c-cross-cross.T, atol=2e-15)
    assert np.max(abs(p-(s+c))) > .5
    assert np.linalg.eigvalsh(p).min() > 0


def test_lag_is_nonwrapped_and_preserves_physical_noise_filter_definition():
    a = np.arange(30., dtype=float).reshape(5, 6)
    good = np.ones_like(a, bool)
    good[0, 0] = False
    result = module.spatial_lag_second_moment(a, good, 0, 1)
    assert result["pairs"] == 24
    assert result["adjacent_difference_filtered_variance_DN2"] == .5
    with pytest.raises(ValueError):
        module.spatial_lag_second_moment(a, good, 0, 0)


def test_pooled_spatial_mean_is_kept_separately_from_local_covariance():
    a = np.array([np.full((4, 5), 2.), np.full((4, 5), -3.)])
    result = module.pooled_time_covariance(a, np.ones((4, 5), bool))
    np.testing.assert_allclose(result["zero_mean_second_moment"], [[4., -6.], [-6., 9.]])
    np.testing.assert_array_equal(result["spatially_centered_covariance"], 0.)
    np.testing.assert_array_equal(result["spatial_mean"], [2., -3.])
