"""Pure fixed-frame Gaia astrometry checks; no catalogue query or file access.

ICRS offsets are east/north relative to A. Gaia ``ra_error`` is already the
alpha-star (RA*cos(delta)) uncertainty and must not receive a second cosine.
Relative covariances use one common local tangent frame (the small-field
approximation) and retain the common reference-A errors in every off-diagonal
block. Optional floors are independent errors per *relative coordinate*.
"""
from __future__ import annotations

from itertools import permutations
import math

import numpy as np


def _finite_array(value, shape, name):
    if value is None:
        raise ValueError(f'{name} is required')
    array = np.asarray(value, dtype=float)
    if array.shape != shape or not np.isfinite(array).all():
        raise ValueError(f'{name} must have finite shape {shape}')
    return array


def tangent_covariance(ra_error_mas, dec_error_mas, ra_dec_corr):
    """Return the Gaia alpha-star/declination covariance in square mas."""
    if any(value is None for value in (ra_error_mas, dec_error_mas, ra_dec_corr)):
        raise ValueError('Gaia ra_error, dec_error and ra_dec_corr are required')
    errors = np.asarray([ra_error_mas, dec_error_mas], dtype=float)
    correlation = float(ra_dec_corr)
    if not np.isfinite(errors).all() or np.min(errors) <= 0:
        raise ValueError('Gaia coordinate errors must be finite and positive')
    if not np.isfinite(correlation) or abs(correlation) >= 1:
        raise ValueError('Gaia ra_dec_corr must be finite and strictly between -1 and 1')
    return np.array([[errors[0]**2, correlation*errors[0]*errors[1]],
                     [correlation*errors[0]*errors[1], errors[1]**2]])


def relative_icrs_offsets_arcsec(ra_deg, dec_deg, reference_index=0):
    """Return [delta RA*cos(delta_A), delta Dec], arcsec, in the fixed ICRS frame.

    This is the requested small-field coordinate map, with a wrapped RA
    difference. No rotation, reflection, plate scale or continuous fit is used.
    """
    ra = np.asarray(ra_deg, dtype=float)
    dec = np.asarray(dec_deg, dtype=float)
    if ra.ndim != 1 or ra.shape != dec.shape or not len(ra) or not np.isfinite(ra).all() or not np.isfinite(dec).all():
        raise ValueError('Finite equal-length RA and declination arrays required')
    if np.max(abs(dec)) >= 90 or not 0 <= reference_index < len(ra):
        raise ValueError('Nonpolar declinations and valid reference index required')
    dra = (ra-ra[reference_index]+180.) % 360.-180.
    east = dra*np.cos(np.deg2rad(dec[reference_index]))*3600.
    north = (dec-dec[reference_index])*3600.
    return np.c_[east, north]


def relative_astrometry_covariance_mas2(source_covariances, reference_index=0,
                                       other_indices=None, coordinate_floor_mas=0.):
    """Construct relative-coordinate covariance, including common-A blocks.

    For independent source positions, Cov(i-A,j-A)=delta_ij*C_i+C_A.
    The optional floor adds floor^2*I after differencing; it is explicitly a
    per-relative-coordinate floor, rather than independent noise on each star.
    Individual Gaia source cross-covariances are not supplied by this helper.
    """
    covariances = np.asarray(source_covariances, dtype=float)
    if covariances.ndim != 3 or covariances.shape[1:] != (2, 2) or not np.isfinite(covariances).all():
        raise ValueError('Finite per-source 2x2 tangent covariances required')
    if not 0 <= reference_index < len(covariances):
        raise ValueError('Invalid reference index')
    for covariance in covariances:
        if not np.allclose(covariance, covariance.T, rtol=0., atol=1e-12):
            raise ValueError('Source covariance must be symmetric')
        np.linalg.cholesky(covariance)
    indices = ([i for i in range(len(covariances)) if i != reference_index]
               if other_indices is None else list(other_indices))
    if not indices or len(set(indices)) != len(indices) or reference_index in indices or min(indices) < 0 or max(indices) >= len(covariances):
        raise ValueError('Distinct non-reference source indices required')
    floor = float(coordinate_floor_mas)
    if not np.isfinite(floor) or floor < 0:
        raise ValueError('Nonnegative finite relative-coordinate floor required')
    covariance = np.tile(covariances[reference_index], (len(indices), len(indices)))
    for block, index in enumerate(indices):
        covariance[2*block:2*block+2, 2*block:2*block+2] += covariances[index]
    covariance += floor**2*np.eye(2*len(indices))
    return covariance


def match_four_sources(ra_deg, dec_deg, conventional_relative_arcsec):
    """Assign A/B/C/D by all24 fixed-frame permutations; no continuous fit.

    Use this one conventional-reference assignment unchanged for both models.
    Its RSS and runner-up gap describe catalogue identification, not a new
    independent astrometric goodness-of-fit statistic after matching.
    """
    ra = _finite_array(ra_deg, (4,), 'four source RAs')
    dec = _finite_array(dec_deg, (4,), 'four source declinations')
    target = np.asarray(conventional_relative_arcsec, dtype=float)
    if target.shape == (3, 2):
        target = target.ravel()
    target = _finite_array(target, (6,), 'conventional B/C/D relative coordinates')
    trials = []
    for assignment in permutations(range(4)):
        positions = relative_icrs_offsets_arcsec(ra[list(assignment)], dec[list(assignment)])
        relative = positions[1:].ravel()
        rss = float(np.sum((relative-target)**2))
        trials.append({'permutation_abcd': list(assignment),
                       'relative_arcsec': relative.tolist(), 'rss_arcsec2': rss})
    trials.sort(key=lambda trial: trial['rss_arcsec2'])
    return {**trials[0], 'runner_up_permutation_abcd': trials[1]['permutation_abcd'],
            'runner_up_rss_arcsec2': trials[1]['rss_arcsec2'],
            'runner_up_gap_arcsec2': trials[1]['rss_arcsec2']-trials[0]['rss_arcsec2'],
            'permutations_evaluated': 24, 'continuous_frame_parameters_fitted': 0,
            'all_permutations': trials}


def gaussian_astrometry_stats(observed_relative_arcsec, predicted_relative_arcsec,
                              covariance_mas2, coordinate_floor_mas=0.):
    """Normalized six-coordinate Gaussian density in mas^-6 and formal chi2.

    Formal Gaia errors plus a declared sensitivity floor define this density;
    this is not a calibrated physical likelihood for all catalogue systematics.
    """
    observed = np.asarray(observed_relative_arcsec, dtype=float).ravel()
    predicted = np.asarray(predicted_relative_arcsec, dtype=float).ravel()
    observed = _finite_array(observed, (6,), 'observed relative coordinates')
    predicted = _finite_array(predicted, (6,), 'predicted relative coordinates')
    covariance = _finite_array(covariance_mas2, (6, 6), 'relative covariance')
    if not np.allclose(covariance, covariance.T, rtol=0., atol=1e-12):
        raise ValueError('Relative covariance must be symmetric')
    floor = float(coordinate_floor_mas)
    if not np.isfinite(floor) or floor < 0:
        raise ValueError('Nonnegative finite relative-coordinate floor required')
    covariance = covariance+floor**2*np.eye(6)
    cholesky = np.linalg.cholesky(covariance)
    residual = (observed-predicted)*1000.
    whitened = np.linalg.solve(cholesky, residual)
    chi2 = float(whitened@whitened)
    logdet = float(2*np.log(np.diag(cholesky)).sum())
    return {'residual_mas': residual.tolist(), 'rms_coordinate_residual_mas': float(np.sqrt(np.mean(residual**2))),
            'maximum_absolute_coordinate_residual_mas': float(np.max(abs(residual))),
            'formal_chi2': chi2, 'coordinate_count': 6,
            'joint_log_density_per_mas6': -.5*(6*math.log(2*math.pi)+logdet+chi2),
            'coordinate_floor_mas': floor, 'covariance_mas2': covariance.tolist(),
            'floor_convention': 'independent added variance per relative B/C/D coordinate, after retaining shared-A formal covariance',
            'continuous_parameters_fitted': 0}


def gaussian_comparison_by_floor(observed_relative_arcsec, predicted_relative_arcsec,
                                 covariance_mas2, floors_mas=(0., 2., 5.)):
    return {str(float(floor)): gaussian_astrometry_stats(observed_relative_arcsec,
              predicted_relative_arcsec, covariance_mas2, floor) for floor in floors_mas}
