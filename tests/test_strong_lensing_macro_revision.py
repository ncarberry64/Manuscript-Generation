"""Physical-source consistency and central-branch discovery of the revision."""
import sys
from pathlib import Path

import numpy as np
from scipy.integrate import quad

sys.path.insert(0, str(Path(__file__).resolve().parents[1]/'code'))

import strong_lensing_macro_revision as revised
import strong_lensing_observed_macro as old
from strong_lensing_robustness import ProfileR1Lens


P = np.array([11.72, .74, -.01, .052, .09, .173, .091, 0., 0.])


def test_two_plummer_density_integrates_to_total_enclosed_mass():
    body = revised.body_from_parameters(P, 1e9)
    assert body.density_Msun_per_Mpc3(0.) > 0
    assert body.mass_and_slope(0.)[1] == 3.
    for r in (3e-5, .001, .006):
        numerical = quad(lambda x: 4*np.pi*x*x*body.density_Msun_per_Mpc3(x),
                         0., r, points=[3e-5] if r > 3e-5 else None,
                         epsabs=1e-4, epsrel=2e-11)[0]
        expected = body.mass_Msun*body.mass_and_slope(r)[0]
        np.testing.assert_allclose(numerical, expected, rtol=3e-11)


def test_zero_compact_mass_recovers_original_body_and_scalar():
    old_body = old.body_from_parameters(P)
    new_body = revised.body_from_parameters(P, 0.)
    coefficients = old.frozen_coefficients()
    old_scalar = ProfileR1Lens(old_body, coefficients)
    new_scalar = ProfileR1Lens(new_body, coefficients)
    for r in (0., 3e-5, .001, .02):
        old_app = old_body.application([r, .3*r])
        new_app = new_body.application([r, .3*r])
        for actual, expected in zip(new_app, old_app):
            np.testing.assert_allclose(actual, expected, rtol=2e-15, atol=1e-23)
        np.testing.assert_allclose(new_scalar.scalar_gradient_natural(r),
                                   old_scalar.scalar_gradient_natural(r), rtol=2e-15)


def test_total_source_obeys_original_flux_and_gauss_not_sum_of_scalar_fields():
    body = revised.body_from_parameters(P, 1e9)
    coefficients = old.frozen_coefficients()
    screened = ProfileR1Lens(body, coefficients)
    sections = revised.source_section_checks(screened, [3e-5, .001, .006])
    assert max(row['scalar_relative_residual'] for row in sections) < 2e-14
    assert max(row['metric_relative_residual'] for row in sections) < 2e-14
    assert all(row['consumed_reduced_source_J_eff'] < 0 for row in sections)
    separate = sum(ProfileR1Lens(revised.TwoPlummerLens(m, a), coefficients).scalar_gradient_natural(.006)
                   for m, a in body.components if m > 0)
    actual = screened.scalar_gradient_natural(.006)
    assert abs(separate-actual)/abs(actual) > .01


def test_two_plummer_potential_gradient_and_hessian_are_consistent():
    body = revised.body_from_parameters(P, 1e9)
    x = np.array([4e-5, -2e-5])
    _, alpha, hessian = body.application(x)
    step = 1e-8
    gradient = []
    derivative = []
    for direction in np.eye(2):
        plus = body.application(x+step*direction)
        minus = body.application(x-step*direction)
        plus2 = body.application(x+2*step*direction)
        minus2 = body.application(x-2*step*direction)
        gradient.append((-plus2[0]+8*plus[0]-8*minus[0]+minus2[0])/(12*step))
        derivative.append((-plus2[1]+8*plus[1]-8*minus[1]+minus2[1])/(12*step))
    np.testing.assert_allclose(alpha, gradient, rtol=2e-8, atol=1e-17)
    np.testing.assert_allclose(hessian, np.asarray(derivative).T, rtol=2e-8, atol=1e-15)


def test_global_polynomial_discovers_three_inner_images_when_present():
    original = revised.conventional_census(P, compact_mass_Msun=0.)
    compact = revised.conventional_census(P, compact_mass_Msun=1e8)
    assert original['image_count'] == 5
    assert compact['polynomial_degree'] == 9
    assert compact['image_count'] == 7
    assert compact['polynomial_positive_candidates'] == 7
    radii = np.linalg.norm(compact['images_light_centered_arcsec'], axis=1)
    assert sum(radii < .05) == 3
    assert compact['max_lens_residual_rad'] < 2e-19
    assert compact['min_jacobian_singular_value'] > .01


def test_original_r1_projection_verifies_continued_central_and_bright_roots():
    census = revised.image_census(P, compact_mass_Msun=1e8, model='r1')
    assert census['image_count'] == census['conventional_image_count'] == 7
    assert census['max_lens_residual_rad'] < 2e-19
    assert census['max_arrival_stationarity_Mpc'] < 1e-15


def test_injected_revised_macro_and_training_callback_recovered():
    actual = revised.macro_application(P, compact_mass_Msun=1e9)
    observed = dict(observed_relative_arcsec=actual['relative_astrometry_arcsec'],
                    sigma_arcsec=np.full(8, .002),
                    observed_images_light_centered_arcsec=actual['images_light_centered_arcsec'],
                    assumptions='synthetic injected source, independent observational application')
    calls = []
    def likelihood(delay):
        calls.append(np.array(delay))
        return 37.-.5*np.sum(((delay-actual['delays_days'])/.1)**2)
    fit = revised.fit_macro(actual['delays_days'], np.eye(3)*.1**2,
                            compact_mass_Msun=1e9, initial=P, astrometry=observed,
                            delay_log_likelihood=likelihood, exact_maxiter=2,
                            include_census=False)
    assert len(calls) > 18
    assert fit['chi2_astrometry'] < 1e-10
    np.testing.assert_allclose(fit['predicted_delays_days'], actual['delays_days'], atol=1e-6, rtol=0.)
    assert fit['max_lens_residual_rad'] < 2e-19
    assert max(abs(np.array(fit['two_plane_minus_effective_delay_days']))) < 1e-10
    assert fit['compact_mass_Msun'] == 1e9
    assert fit['physical_admission'] is False


def test_exact_only_warm_start_keeps_best_original_application_endpoint():
    actual = revised.macro_application(P, compact_mass_Msun=1e9)
    observed = dict(observed_relative_arcsec=actual['relative_astrometry_arcsec'],
                    sigma_arcsec=np.full(8, .002),
                    observed_images_light_centered_arcsec=actual['images_light_centered_arcsec'],
                    assumptions='synthetic exact-only warm-start application')
    likelihood = lambda delay: 37.-.5*np.sum(((delay-actual['delays_days'])/.1)**2)
    fit = revised.fit_macro(actual['delays_days'], np.eye(3)*.1**2,
                            compact_mass_Msun=1e9, initial=P, astrometry=observed,
                            delay_log_likelihood=likelihood, exact_maxiter=2,
                            include_census=False, skip_quadratic_initializer=True,
                            retain_best_original_endpoint=True)
    assert fit['optimizer_trace'] == []
    endpoints = fit['original_objective_endpoint_selection']
    assert endpoints['original_projection_objective'] <= endpoints['checked_endpoints'][0]['objective']
    assert fit['exact_training_objective'] == endpoints['original_projection_objective']
    assert fit['skip_quadratic_initializer'] is True
