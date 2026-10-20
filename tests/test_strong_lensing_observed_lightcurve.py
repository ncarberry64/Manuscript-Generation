"""Meaningful exact-GP, holdout leakage, driver conditioning and delay checks."""
from dataclasses import replace
from pathlib import Path
import sys

import numpy as np
from scipy.linalg import cho_factor, cho_solve
from scipy.optimize import minimize_scalar

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT/'code'))
from strong_lensing_observed_lightcurve import (
    LightCurves, OUParameters, concatenate, dense_covariance, log_likelihood,
    predict_heldout, score_predictions, training_delay_curvature,
)


def fixture():
    parameters = OUParameters((0., 2.7, -1.3), source_sigma_mag=.18,
                              source_tau_days=17., extrinsic_sigma_mag=.04,
                              extrinsic_tau_days=190., jitter_mag=.002,
                              offset_sigma_mag=.8, magnitude_zero=-10.)
    rng = np.random.default_rng(2081)
    times = np.array([0., 2., 5., 8., 11., 15., 18.])
    time = np.repeat(times, 3)
    image = np.tile(np.arange(3), len(times))
    placeholder = LightCurves(time, image, np.full(len(time), -10.),
                              rng.uniform(.005, .02, len(time)), 3)
    y = -10.+np.linalg.cholesky(dense_covariance(placeholder, parameters))@rng.normal(size=len(time))
    data = LightCurves(time, image, y, placeholder.error_mag, 3)
    # Mixed chronological points deliberately interleave after time shifting.
    training_mask = time < 11.
    return data.subset(training_mask), data.subset(~training_mask), parameters


def dense_conditional(train, test, parameters):
    data = concatenate(train, test)
    covariance = dense_covariance(data, parameters)
    size = train.time_days.size
    factor = cho_factor(covariance[:size, :size], lower=True)
    cross = covariance[size:, :size]
    mean = parameters.magnitude_zero+cross@cho_solve(
        factor, train.magnitude-parameters.magnitude_zero)
    conditional = covariance[size:, size:]-cross@cho_solve(factor, cross.T)
    return mean, conditional


def test_kalman_normalized_likelihood_matches_independent_dense_gp():
    train, test, parameters = fixture()
    data = concatenate(train, test)
    factor = cho_factor(dense_covariance(data, parameters), lower=True)
    centered = data.magnitude-parameters.magnitude_zero
    reference = -.5*(len(centered)*np.log(2*np.pi)
                      +2*np.log(np.diag(factor[0])).sum()
                      +centered@cho_solve(factor, centered))
    np.testing.assert_allclose(log_likelihood(data, parameters), reference, rtol=0, atol=1e-9)


def test_train_only_rts_predictions_match_dense_conditioning_at_interleaved_times():
    train, test, parameters = fixture()
    # Reference image observations are allowed driving data; other images stay hidden.
    conditioning = concatenate(train, test.subset(test.image == 1))
    targets = test.subset(test.image != 1)
    reference_mean, reference_covariance = dense_conditional(conditioning, targets, parameters)
    prediction = predict_heldout(conditioning, targets, parameters)
    np.testing.assert_allclose(prediction['mean_magnitude'], reference_mean, rtol=0, atol=1e-9)
    np.testing.assert_allclose(prediction['variance_mag2'], np.diag(reference_covariance), rtol=0, atol=1e-10)


def test_withheld_magnitudes_cannot_change_prediction_but_change_scoring():
    train, test, parameters = fixture()
    altered = LightCurves(test.time_days, test.image, test.magnitude+2., test.error_mag, test.nimage)
    first = predict_heldout(train, test, parameters)
    second = predict_heldout(train, altered, parameters)
    np.testing.assert_array_equal(first['mean_magnitude'], second['mean_magnitude'])
    np.testing.assert_array_equal(first['variance_mag2'], second['variance_mag2'])
    assert (score_predictions(train, altered, parameters)['joint_conditional_test_log_density']
            < score_predictions(train, test, parameters)['joint_conditional_test_log_density'])


def test_exact_correlated_conditional_test_density_matches_dense_formula():
    train, test, parameters = fixture()
    mean, covariance = dense_conditional(train, test, parameters)
    factor = cho_factor(covariance, lower=True)
    residual = test.magnitude-mean
    reference = -.5*(len(residual)*np.log(2*np.pi)
                     +2*np.log(np.diag(factor[0])).sum()
                     +residual@cho_solve(factor, residual))
    result = score_predictions(train, test, parameters)
    np.testing.assert_allclose(result['joint_conditional_test_log_density'], reference, rtol=0, atol=1e-9)
    # The product of marginal densities is deliberately not used as joint density.
    marginal_sum = -result['overall']['mean_marginal_negative_log_density']*len(residual)
    assert abs(marginal_sum-reference) > .1


def test_synthetic_time_delay_recovers_sign_and_value_without_any_test_values():
    rng = np.random.default_rng(2083)
    time = np.sort(rng.uniform(0., 600., 190))
    parameters = OUParameters((0., 6.5), source_sigma_mag=.2, source_tau_days=40.,
                              extrinsic_sigma_mag=.025, extrinsic_tau_days=1000.,
                              jitter_mag=.001, offset_sigma_mag=.8, magnitude_zero=18.)
    times = np.repeat(time, 2)
    images = np.tile(np.arange(2), len(time))
    placeholder = LightCurves(times, images, np.full(len(times), 18.),
                              np.full(len(times), .006), 2)
    sample = 18.+np.linalg.cholesky(dense_covariance(placeholder, parameters))@rng.normal(size=len(times))
    data = LightCurves(times, images, sample, placeholder.error_mag, 2)
    fitted = minimize_scalar(lambda delay: -log_likelihood(
        data, replace(parameters, delays_days=(0., float(delay)))),
        bounds=(2., 11.), method='bounded', options={'xatol': 1e-7})
    assert fitted.success
    assert abs(fitted.x-6.5) < .3
    assert log_likelihood(data, replace(parameters, delays_days=(0., fitted.x))) > log_likelihood(
        data, replace(parameters, delays_days=(0., -6.5)))+50.
    local = training_delay_curvature(data, replace(parameters, delays_days=(0., fitted.x)))
    assert local['eigenvalues_per_day2'][0] > 0


def test_common_time_translation_and_delay_zero_point_do_not_change_likelihood():
    train, _, parameters = fixture()
    moved = LightCurves(train.time_days+2.4e6, train.image, train.magnitude,
                        train.error_mag, train.nimage)
    np.testing.assert_allclose(log_likelihood(train, parameters), log_likelihood(moved, parameters),
                               rtol=0, atol=2e-8)


def test_observation_inputs_are_copied_and_immutable():
    values = np.array([18., 18.1])
    data = LightCurves(np.array([0., 1.]), np.array([0, 0]), values, np.array([.01, .01]), 1)
    values[:] = 100.
    np.testing.assert_array_equal(data.magnitude, [18., 18.1])
    assert not data.magnitude.flags.writeable


def test_json_delay_list_cannot_mutate_a_frozen_model():
    delays = [0., 2.]
    parameters = OUParameters(delays)
    delays[1] = 100.
    assert parameters.delays_days == (0., 2.)
