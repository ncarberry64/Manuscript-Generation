"""Independent Gram/coordinate checks for the source-calibrated spectral model."""
import sys
from pathlib import Path
import numpy as np
import pytest
sys.path.insert(0, str(Path(__file__).resolve().parents[1]/"code"))
from strong_lensing_dark_covariance import (
    template_spectra, SpectralTemplateModel, native_to_clock_coordinates,
    pool_spectral_matrices, BinnedSpectralModel, anchored_group0_residual_jacobian,
    full_anchored_residual_maps,
)


def test_fft_covariance_exactly_matches_shared_innovation_dense_gram():
    rng = np.random.default_rng(7)
    z = rng.normal(size=(2, 3, 8, 6))
    spectra, _ = template_spectra(z, np.ones((2, 8, 6), bool))
    amp = np.array([0, 0, 1, 1])
    xy = np.array([[0, 0], [2, 3], [1, 4], [4, 7]])
    weights = rng.normal(size=(4, 3))
    actual = SpectralTemplateModel(spectra, (8, 6)).project([
        {"amp": amp, "clock_xy": xy, "group_weights": weights}])[0]
    yy, xx = np.indices((8, 6))
    rows = []
    for a, (x, y), w in zip(amp, xy, weights):
        feature = sum(w[g]*z[a, g, (y-yy) % 8, (x-xx) % 6] for g in range(3))/np.sqrt(48)
        rows.append(feature.ravel())
    rows = np.asarray(rows)
    np.testing.assert_allclose(actual, rows@rows.T, rtol=3e-15, atol=3e-15)
    assert actual[0, 2] != 0  # Cross-amplifier covariance retained.
    assert np.linalg.eigvalsh(actual).min() > 0


def test_window_normalization_matches_actual_valid_pixel_variance():
    z = np.arange(32., dtype=float).reshape(1, 1, 8, 4)
    good = np.ones((1, 8, 4), bool); good[0, 0] = False
    z[0, 0, 0] = np.nan
    f, coverage = template_spectra(z, good)
    model = SpectralTemplateModel(f, (8, 4))
    c = model.project([{"amp": np.array([0]), "clock_xy": np.array([[0, 0]]),
                         "group_weights": np.ones((1, 1))}])[0]
    assert coverage[0] == .875
    assert c[0, 0] == pytest.approx(np.mean(z[0, 0][good[0]]**2), rel=2e-15)


def test_native_clock_parity_and_row_direction_are_explicit():
    amp, xy = native_to_clock_coordinates(np.array([[4, 2043], [512+4, 2043],
        [1024+507, 4], [1536+507, 4]]))
    np.testing.assert_array_equal(amp, [0, 1, 2, 3])
    np.testing.assert_array_equal(xy, [[4, 0], [507, 0], [507, 2039], [4, 2039]])
    # Internal amplifier edges are science pixels, even when a native17stamp
    # crosses the boundary. They must not be dropped to fit a calibration cell.
    amp2, clock2 = native_to_clock_coordinates(np.array([[1023, 30], [1024, 30]]))
    np.testing.assert_array_equal(amp2, [1, 2])
    np.testing.assert_array_equal(clock2[:, 0], [0, 0])
    with pytest.raises(ValueError, match="outside"):
        native_to_clock_coordinates(np.array([[0, 0]]))


def test_frequency_pool_preserves_variance_psd_and_real_cross_symmetry():
    rng = np.random.default_rng(19)
    z = rng.normal(size=(2, 2, 16, 12))
    f, _ = template_spectra(z, np.ones((2, 16, 12), bool))
    p = pool_spectral_matrices(f)
    model = BinnedSpectralModel(p["spectral_matrices"], p["y_frequency_bin"],
        p["x_frequency_bin"], p["spatial_shape"], 2, 2)
    for a in range(2):
        for g in range(2):
            assert model.cross_kernel(a, g, a, g)[0, 0] == pytest.approx(np.mean(z[a, g]**2), rel=2e-15)
    for matrix in p["spectral_matrices"].reshape(-1, 4, 4):
        assert np.linalg.eigvalsh(matrix).min() > -1e-14
    req = {"amp": np.array([0, 0, 1]), "clock_xy": np.array([[0, 0], [2, 3], [1, 4]]),
           "group_weights": np.array([[1., -1.], [.2, .4], [-.3, 2.]])}
    c = model.project([req])[0]
    assert np.linalg.eigvalsh(c).min() > 0


def test_noisy_group0_anchor_is_an_estimating_equation_derivative():
    observed = np.array([[-1., 1.]])
    predicted = observed.copy()
    beta = np.array([[1., 1.2]])
    answer = anchored_group0_residual_jacobian(observed, predicted, beta)
    np.testing.assert_allclose(answer, [[-1.2, 1.]])
    # Identity group map makes the ordinary slope independent of a constant
    # intercept; its noisyG0 anchoring adds no spurious extra variance.
    np.testing.assert_allclose(anchored_group0_residual_jacobian(observed, predicted,
        np.ones((1, 2))), observed)


def test_full_boundary_map_retains_noncommuting_spatial_response():
    h = np.array([[.9, .1], [.2, .8]])
    g = np.linalg.inv(h)
    ratio = np.diag([1., 1.3])
    prediction = h@ratio
    observed = [np.eye(2)*-.4, np.eye(2)*.6]
    boundary = [g, np.zeros((2, 2))]
    result = full_anchored_residual_maps(observed, prediction, boundary)
    np.testing.assert_allclose(result[0], observed[0]-h@ratio@g)
    assert result[0][0, 1] != 0
    np.testing.assert_array_equal(result[1], observed[1])
