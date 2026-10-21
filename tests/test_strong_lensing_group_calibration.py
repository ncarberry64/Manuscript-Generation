"""Independent derivative, native-average and shared-reference regressions."""
import sys
from pathlib import Path
import numpy as np
import pytest

sys.path.insert(0, str(Path(__file__).resolve().parents[1]/"code"))
from strong_lensing_group_calibration import (
    polynomial_derivatives, inverse_linearity, averaged_read_observation,
    second_order_mean_correction, shared_reference_factor,
)


def test_native_adc_average_is_not_average_of_corrected_reads():
    c = np.array([0., 1., 1e-4])[:, None]
    adc = np.array([[[10.], [30.]], [[50.], [90.]]])
    q = polynomial_derivatives(c, adc)[0]
    answer = averaged_read_observation(q, c, adc)
    expected = polynomial_derivatives(c, adc.mean(axis=1))[0]
    np.testing.assert_allclose(answer["mean_DN"], expected, rtol=1e-13)
    assert not np.allclose(expected, q.mean(axis=1))
    np.testing.assert_allclose(inverse_linearity(q, c, adc*.97), adc, atol=2e-9)


def test_read_jacobian_and_hessian_match_numerical_derivatives():
    c = np.array([.2, 1.01, 3e-5, 2e-9])[:, None]
    adc = np.array([[[100.], [145.]]])
    q = polynomial_derivatives(c, adc)[0]
    result = averaged_read_observation(q, c, adc)
    h = .01
    for r in range(2):
        delta = np.zeros_like(q); delta[:, r] = h
        plus = averaged_read_observation(q+delta, c, adc)
        minus = averaged_read_observation(q-delta, c, adc)
        np.testing.assert_allclose((plus["mean_DN"]-minus["mean_DN"])/(2*h),
                                   result["read_to_group_jacobian"][:, r], rtol=1e-8)
        np.testing.assert_allclose((plus["read_to_group_jacobian"]-minus["read_to_group_jacobian"])/(2*h),
                                   result["read_to_group_hessian"][:, :, r], rtol=2e-6, atol=2e-10)


def test_identity_linearity_gives_exact_average_and_zero_mean_correction():
    c = np.array([0., 1.])[:, None, None]
    q = np.arange(24., dtype=float).reshape(3, 2, 2, 2)+10
    result = averaged_read_observation(q, c, q, dark_group_dn=2.)
    np.testing.assert_allclose(result["mean_DN"], q.mean(axis=1)-2.)
    np.testing.assert_allclose(result["read_to_group_jacobian"], .5)
    np.testing.assert_array_equal(second_order_mean_correction(
        result["read_to_group_hessian"], np.ones((3, 2, 2, 2, 2))), 0.)


def test_inverse_rejects_invalid_branch_and_out_of_domain():
    with pytest.raises(ValueError, match="not monotone"):
        inverse_linearity(np.array([2.]), np.array([0., -1.]), np.array([1.]))
    with pytest.raises(ValueError, match="outside"):
        inverse_linearity(np.array([200.]), np.array([0., 1.]), np.array([20.]), adc_lower=0., adc_upper=100.)


def test_shared_reference_factor_preserves_signed_cross_group_covariance():
    derivative = np.array([[1., 0.], [-1., 0.], [2., 1.]])
    factor = shared_reference_factor(derivative, np.array([.5, 2.]))
    covariance = factor@factor.T
    assert covariance[0, 1] == -.25
    assert covariance[0, 2] == .5
    np.testing.assert_allclose(np.linalg.eigvalsh(covariance).min(), 0., atol=1e-15)
