from dataclasses import replace
import json
from pathlib import Path
import sys

import numpy as np
from scipy.integrate import quad

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT/'code'))
from strong_lensing_r1_screened import R1Coefficients
from strong_lensing_robustness import BASE_HALF_MASS_MPC, ProfileLens, ProfileR1Lens, projected_coefficient_sensitivities


def coefficients():
    saved = json.loads((ROOT/'artifacts/strong_lensing_r1_screened_20261020/results.json').read_text())
    return R1Coefficients(**saved['coefficients'])


def test_total_mass_half_mass_and_regular_center_are_common_physical_inputs():
    for name in ('tolman_vii', 'plummer'):
        base = ProfileLens(name)
        a = base.radius_Mpc
        upper = 1. if name == 'tolman_vii' else np.inf
        normalized = quad(lambda t: 4*np.pi*(a*t)**2*base.density_Msun_per_Mpc3(a*t)*a/base.mass_Msun,
                          0., upper, epsabs=1e-12, epsrel=1e-12)[0]
        np.testing.assert_allclose(normalized, 1., rtol=2e-12)
        np.testing.assert_allclose(base.mass_and_slope(BASE_HALF_MASS_MPC)[0], .5, rtol=2e-14)
        assert np.isfinite(base.density_Msun_per_Mpc3(0.))
        assert np.isfinite(base.application([0., 0.])[2]).all()


def test_plummer_projected_density_mass_and_deflection_are_the_same_profile():
    base = ProfileLens('plummer')
    a, M = base.radius_Mpc, base.mass_Msun
    for b in (.2*a, a, 3*a):
        projected_density = 2*quad(lambda t: base.density_Msun_per_Mpc3(np.sqrt(b*b+(a*t)**2))*a,
                                  0., np.inf, epsabs=1e-2, epsrel=1e-12)[0]
        expected_density = M*a*a/(np.pi*(b*b+a*a)**2)
        np.testing.assert_allclose(projected_density, expected_density, rtol=2e-12)
        projected_mass = quad(lambda s: 2*np.pi*s*M*a*a/(np.pi*(s*s+a*a)**2), 0., b,
                              epsabs=1e-3, epsrel=1e-12)[0]
        alpha = base.application([b, 0.])[1][0]
        np.testing.assert_allclose(alpha, 4*base.mass_length_Mpc/b*projected_mass/M, rtol=2e-12)


def test_both_profile_gradients_hessians_and_size_responses_are_consistent():
    step = 1e-5
    for name in ('tolman_vii', 'plummer'):
        base = ProfileLens(name)
        x = np.array([.7*base.radius_Mpc, .2*base.radius_Mpc])
        p, alpha, hessian = base.application(x)
        dx = base.radius_Mpc*1e-5
        for i in range(2):
            direction = np.eye(2)[i]*dx
            plus = base.application(x+direction)
            minus = base.application(x-direction)
            np.testing.assert_allclose((plus[0]-minus[0])/(2*dx), alpha[i], rtol=2e-9)
            np.testing.assert_allclose((plus[1]-minus[1])/(2*dx), hessian[:, i], rtol=2e-9)
        for name, parameter in (('log_mass', 'mass_Msun'), ('log_radius', 'half_mass_Mpc')):
            value = getattr(base, parameter)
            plus = replace(base, **{parameter:value*np.exp(step)}).application(x)
            minus = replace(base, **{parameter:value*np.exp(-step)}).application(x)
            dp, da = base.intrinsic(name)(x)
            np.testing.assert_allclose((plus[0]-minus[0])/(2*step), dp, rtol=2e-9)
            np.testing.assert_allclose((plus[1]-minus[1])/(2*step), da, rtol=2e-9)


def test_plummer_scalar_flux_uses_the_integrated_matter_source_and_positive_halo():
    base = ProfileLens('plummer')
    lens = ProfileR1Lens(base, coefficients())
    c = lens.coefficients
    for r in (.1*base.radius_Mpc, base.half_mass_Mpc, 10*base.radius_Mpc):
        enclosed = quad(lambda s: 4*np.pi*s*s*base.density_Msun_per_Mpc3(s), 0., r,
                        epsabs=1e-3, epsrel=2e-12)[0]
        x = r/c.hubble_length_Mpc
        p = lens.scalar_gradient_natural(r)
        flux = c.C*x*x*p-4/c.lambda_over_MplH02*x*p*p
        source = -c.b*base.mass_length_Mpc/c.hubble_length_Mpc*enclosed/base.mass_Msun
        np.testing.assert_allclose(flux, source, rtol=2e-12)
        assert lens.effective_density_Msun_per_Mpc3(r) > 0


def test_eight_recomputed_cases_preserve_state_and_have_valid_original_residuals():
    path = ROOT/'artifacts/strong_lensing_acceptance_20261020/robustness.json'
    result = json.loads(path.read_text())
    assert len(result['cases']) == 8
    assert not result['action_state']['changed'] and not result['action_state']['new_scalar_amplitude']
    assert not result['physical_admission']
    assert {x['profile'] for x in result['cases']} == {'tolman_vii', 'plummer'}
    for case in result['cases']:
        assert len(case['conventional']['images_rad']) == len(case['screened']['images_rad']) == 3
        assert case['original_scalar_flux_relative_residual'] < 1e-12
        assert case['maximum_image_residual_rad'] < 1e-18
        assert len(case['reference_relative_offsets']['delays_seconds']) == 2
        assert min(case['reference_relative_offsets']['delays_seconds']) > 0


def test_effective_mass_equivalence_and_nuisance_scope_survive_profile_change():
    result = json.loads((ROOT/'artifacts/strong_lensing_acceptance_20261020/robustness.json').read_text())
    cases = [x for x in result['cases'] if x['case'] == 'reference']
    for case in cases:
        equivalence = case['effective_mass_equivalence']
        assert equivalence['nonnegative_density']
        assert max(x['relative_residual'] for x in equivalence['enclosed_mass_checks']) < 1e-9
        assert equivalence['independent_direct_z_vs_R1_alpha_relative_error'] < 1e-9
        assert equivalence['maximum_image_difference_arcsec'] < 1e-12
        assert equivalence['maximum_delay_difference_seconds'] < 1e-6
        assert case['linearized_nuisance_check']['fixed_center']['signal_rank'] == 1
        assert case['linearized_nuisance_check']['free_center']['signal_rank'] == 0
        assert case['leading_matter_support']['central_pressure_fractional_adjustment'] > 0
    # The profile change is a physical input control, not a new fitted amplitude.
    assert result['cross_profile_baseline']['ordinary_profile_change_to_max_R1_offset_ratio'] > 1000
    np.testing.assert_allclose(cases[0]['reference_relative_offsets']['delays_seconds'],
                               [25.959258149566722, 8.649019375383205], atol=2e-8, rtol=0)


def test_projected_state_coefficient_partials_are_weighted_and_match_finite_variations():
    c = coefficients()
    for profile in ('tolman_vii', 'plummer'):
        base = ProfileLens(profile)
        lens = ProfileR1Lens(base, c)
        r = base.half_mass_Mpc
        values = projected_coefficient_sensitivities(lens, r)
        assert abs(values['identity_residual_2_db_plus_dC_minus_3']) < 2e-12
        assert 1.5 < values['d_log_extra_alpha']['ln_b_T'] < 2.
        assert -1. < values['d_log_extra_alpha']['ln_C'] < 0.
        h = 1e-4
        for key, attribute in (('ln_b_T', 'b'), ('ln_C', 'C')):
            original = getattr(c, attribute)
            plus = ProfileR1Lens(base, replace(c, **{attribute:original*np.exp(h)})).extra_alpha(r)
            minus = ProfileR1Lens(base, replace(c, **{attribute:original*np.exp(-h)})).extra_alpha(r)
            finite = (np.log(plus)-np.log(minus))/(2*h)
            np.testing.assert_allclose(finite, values['d_log_extra_alpha'][key], rtol=2e-7, atol=2e-10)
