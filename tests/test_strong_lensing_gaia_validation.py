"""Independent coordinate/covariance/density checks for the sealed Gaia test."""
import sys
from pathlib import Path
import math

import numpy as np
import pytest

sys.path.insert(0, str(Path(__file__).resolve().parents[1]/'code'))
import strong_lensing_gaia_validation as gaia


def test_relative_icrs_east_north_and_ra_wrap():
    ra = np.array([359.9999, .0001, 359.9998])
    dec = np.array([60., 60.+1/3600, 60.-2/3600])
    offsets = gaia.relative_icrs_offsets_arcsec(ra, dec)
    np.testing.assert_allclose(offsets, [[0, 0], [.36, 1.], [-.18, -2.]], atol=3e-10, rtol=0.)


def test_common_reference_errors_appear_in_all_off_diagonal_blocks():
    source = np.array([gaia.tangent_covariance(i, i+1, .2) for i in (1., 2., 3., 4.)])
    relative = gaia.relative_astrometry_covariance_mas2(source)
    np.testing.assert_allclose(relative[:2, :2], source[0]+source[1])
    np.testing.assert_allclose(relative[:2, 2:4], source[0])
    np.testing.assert_allclose(relative[2:4, 4:6], source[0])
    np.testing.assert_allclose(gaia.relative_astrometry_covariance_mas2(source, coordinate_floor_mas=2.)-relative, 4*np.eye(6))


def test_gaia_ra_error_already_is_alpha_star_no_second_cosine():
    covariance = gaia.tangent_covariance(2., 3., -.5)
    np.testing.assert_allclose(covariance, [[4., -3.], [-3., 9.]])


def test_all24_assignments_recover_unique_fixed_frame_labels():
    relative = np.array([[-1.5, .55], [-2.4, -.6], [-.9, -1.6]])
    ra = 70.+np.r_[0., relative[:, 0]]/(3600*np.cos(np.deg2rad(-12.)))
    dec = -12.+np.r_[0., relative[:, 1]]/3600
    shuffle = [2, 0, 3, 1]
    result = gaia.match_four_sources(ra[shuffle], dec[shuffle], relative)
    assert result['permutation_abcd'] == [1, 3, 0, 2]
    assert result['permutations_evaluated'] == 24
    assert result['rss_arcsec2'] < 1e-17
    assert result['runner_up_gap_arcsec2'] > 1.
    assert result['continuous_frame_parameters_fitted'] == 0


def test_normalized_joint_gaussian_matches_known_six_coordinate_density():
    observed = np.arange(1., 7.)/1000
    result = gaia.gaussian_astrometry_stats(observed, np.zeros(6), np.eye(6)*4)
    expected_chi2 = float(np.sum(np.arange(1., 7.)**2)/4)
    expected_logpdf = -.5*(6*math.log(2*math.pi)+6*math.log(4.)+expected_chi2)
    np.testing.assert_allclose(result['formal_chi2'], expected_chi2, rtol=0., atol=1e-14)
    np.testing.assert_allclose(result['joint_log_density_per_mas6'], expected_logpdf, rtol=0., atol=1e-14)
    controls = gaia.gaussian_comparison_by_floor(observed, np.zeros(6), np.eye(6)*4)
    assert controls['5.0']['formal_chi2'] < controls['2.0']['formal_chi2'] < controls['0.0']['formal_chi2']


def test_missing_correlation_or_invalid_uncertainty_cannot_become_zero():
    for args in ((1., 2., None), (0., 2., .1), (1., 2., np.nan), (1., 2., 1.)):
        with pytest.raises(ValueError):
            gaia.tangent_covariance(*args)
    with pytest.raises(np.linalg.LinAlgError):
        gaia.gaussian_astrometry_stats(np.zeros(6), np.zeros(6), np.zeros((6, 6)))
