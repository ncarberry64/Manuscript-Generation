"""Original-equation checks for the finite inverse application."""
from pathlib import Path
import json
import sys

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT/'code'))
from strong_lensing_comparison import ARCSEC
from strong_lensing_local_family import saved_geometry
from strong_lensing_refit import fit_conventional, independent_forward_solution, original_rows


def inputs():
    data = json.loads((ROOT/'artifacts/strong_lensing_r1_screened_20261020/results.json').read_text())
    _, geometry, source = saved_geometry(ROOT/'artifacts/strong_lensing_first_comparison_20261020/results.json')
    return data, geometry, source


def test_exact_mass_sheet_symmetry_preserves_images_and_scales_relative_times():
    data, geometry, source = inputs()
    body = data['body']
    seeds = np.array(data['conventional']['images_rad'])[:, 0]
    base = independent_forward_solution(np.zeros(5), geometry, body['mass_Msun'], body['radius_Mpc'], source, seeds)
    lam = .97
    beta_arcsec = source[0]/geometry.D_os[0, 0]/ARCSEC
    transformed = np.array([np.log(lam), 0., (lam-1)*beta_arcsec, 0., 1-lam])
    other = independent_forward_solution(transformed, geometry, body['mass_Msun'], body['radius_Mpc'], source, seeds)
    np.testing.assert_allclose(other['images_rad'], base['images_rad'], rtol=0, atol=2e-19)
    np.testing.assert_allclose(other['relative_delays_days'], lam*base['relative_delays_days'], rtol=2e-13)


def test_recovery_of_a_finite_conventional_injection_from_five_original_rows():
    data, geometry, source = inputs()
    mass, radius = data['body']['mass_Msun'], data['body']['radius_Mpc']
    seeds = np.array(data['conventional']['images_rad'])[:, 0]
    truth = np.array([.001, -.002, .00005, .00004, .002])
    generated = independent_forward_solution(truth, geometry, mass, radius, source, seeds)
    fit = fit_conventional(geometry, mass, radius, source, generated['images_rad'], generated['relative_delays_days'])
    np.testing.assert_allclose(fit['coordinates'], truth, rtol=0, atol=1e-10)
    assert np.max(abs(fit['image_errors_microarcsec'])) < 1e-5
    assert np.max(abs(fit['delay_errors_seconds'])) < 1e-5


def test_R1_injection_has_a_finite_free_center_conventional_fit():
    data, geometry, source = inputs()
    fit = fit_conventional(geometry, data['body']['mass_Msun'], data['body']['radius_Mpc'], source,
                          np.array(data['screened']['images_rad'])[:, 0],
                          np.array(data['screened']['relative_delays_days']))
    assert fit['optimizer_success']
    assert fit['scaled_residual_norm'] < 1e-9
    assert fit['numerical_Jacobian_singular_values'][-1] > 1.
    assert fit['forward']['max_original_lens_residual_rad'] < 1e-17
    assert np.max(abs(fit['delay_errors_seconds'])) < 1e-5
    assert np.max(abs(fit['coordinates'])) < 1e-4


def test_fixed_center_fit_retains_a_finite_residual_in_the_declared_profile_family():
    data, geometry, source = inputs()
    mass, radius = data['body']['mass_Msun'], data['body']['radius_Mpc']
    images = np.array(data['screened']['images_rad'])[:, 0]
    delays = np.array(data['screened']['relative_delays_days'])
    fit = fit_conventional(geometry, mass, radius, source, images, delays, free_center=False)
    assert fit['coordinates'][3] == 0.
    assert fit['scaled_residual_norm'] > 1e-8
    assert .001 < np.max(abs(fit['delay_errors_seconds'])) < .01
    np.testing.assert_allclose(original_rows(fit['coordinates'], geometry, mass, radius, source, images, delays),
                               fit['original_rows'], atol=0, rtol=0)
