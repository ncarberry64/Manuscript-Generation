"""Independent limits, conserved moments, native stages and likelihood tests."""
from pathlib import Path
import sys
import numpy as np
import pytest
from scipy import linalg, sparse
from scipy.integrate import solve_ivp

sys.path.insert(0, str(Path(__file__).resolve().parents[1]/"code"))
from strong_lensing_stochastic_response import (
    BoundaryMomentModel, collection_boundary_generator, nearest_neighbor_edges,
    positive_collection_rates, pixel_read_projection, project_collection_covariance,
    independent_adc_read_covariance, shared_calibration_covariance,
    fixed_residual_covariance_application, gaussian_nll_and_score,
    harmonic_rate_initializer,
    within_pixel_read_covariance_first_order,
)


def project(model, times, pixels, weights, **kwargs):
    return project_collection_covariance(model, times,
        pixel_read_projection(pixels, weights, len(model.rate_e_s)),
        provenance="synthetic exact read map; output DN", **kwargs)


def test_original_t2_mean_nonuniform_gain_and_total_charge_conservation():
    gain = np.array([1.5, 2., 3.1, 1.9])
    dn_rate = np.array([5., 1., 2., 4.])
    rate = gain*dn_rate
    edges = nearest_neighbor_edges((2, 2))
    model = BoundaryMomentModel(rate, 2e-3, edges)
    lap = collection_boundary_generator(np.ones(4), .5, edges).toarray()
    np.testing.assert_allclose(model.generator@rate, 2e-3*lap@rate**2, atol=1e-16)
    lead_dn = model.original_leading_mean(.2)/gain
    np.testing.assert_allclose(lead_dn,
        dn_rate*.2+.5*.2**2*2e-3*(lap@(gain*dn_rate)**2)/gain)
    for t in (0., .01, .2, 20.):
        assert model.mean(t).sum() == pytest.approx(rate.sum()*t, rel=5e-14, abs=1e-14)
        assert model.birth_rate(t).sum() == pytest.approx(rate.sum(), rel=5e-14)
        assert model.birth_rate(t).min() >= 0
    assert not np.isclose(np.sum(lap@rate**2/gain), 0.)


def test_positive_completion_matches_local_jacobian_but_never_negative():
    rate = np.array([2., 4., 8.])
    model = BoundaryMomentModel(rate, .01, nearest_neighbor_edges((1, 3)))
    np.testing.assert_allclose(positive_collection_rates(np.zeros(3), rate, model.generator), rate)
    # Forward difference is valid at the nonnegative count boundary.
    numerical = np.column_stack([(positive_collection_rates(np.eye(3)[i]*1e-5,
        rate, model.generator)-rate)/1e-5 for i in range(3)])
    np.testing.assert_allclose(numerical, model.generator.toarray(), atol=3e-8)
    for counts in (np.array([0., 1e5, 0.]), np.array([1e6, 0., 1e6])):
        answer = positive_collection_rates(counts, rate, model.generator)
        assert answer.min() >= 0.
        assert answer.sum() == pytest.approx(rate.sum(), rel=1e-11)
    with pytest.raises(ValueError, match="Positive incident"):
        positive_collection_rates(np.zeros(3), np.array([0., 1., 2.]), model.generator)


def test_zero_response_exact_poisson_min_time_and_group_read_factors():
    rate = np.array([6., 10.])
    model = BoundaryMomentModel(rate, 0., nearest_neighbor_edges((1, 2)))
    times = np.array([1., 2., 5., 6.])
    # Two group averages of the same pixel, electron -> DN at gain two.
    weights = np.array([[.25, .25, 0., 0.], [0., 0., .25, .25]])
    answer = project(model, times, np.array([0, 0]), weights)
    kernel = np.minimum(times[:, None], times[None, :])
    np.testing.assert_allclose(answer["covariance"], 6*weights@kernel@weights.T,
                               rtol=3e-15, atol=3e-15)
    np.testing.assert_allclose(answer["mean"], weights@(6*times))
    assert answer["diagnostics"]["exponential_degree"] == 0
    # A CDS map has sigma_single^2=CDS^2/2: averaged two reads gives CDS^2/4.
    maps = pixel_read_projection([0, 0], np.array([[.5, .5, 0., 0.],
        [0., 0., .5, .5]]), 2)
    read = independent_adc_read_covariance(maps, np.array([4., 9.]),
                                           stage_assumption="synthetic raw ADC independent reads")
    np.testing.assert_allclose(read, np.eye(2)*4.)


def test_projected_covariance_matches_independent_moment_ode_and_cross_time():
    rate = np.array([2., 5., 3.])
    model = BoundaryMomentModel(rate, .003, nearest_neighbor_edges((1, 3)))
    a = model.generator.toarray()
    times = np.array([.3, 1.2, 2.])
    # All three pixels at three times, preserving temporal off-diagonal blocks.
    pixels = np.tile(np.arange(3), 3)
    weights = np.repeat(np.eye(3), 3, axis=0)
    answer = project(model, times, pixels, weights)
    def derivative(t, y):
        mu, c = y[:3], y[3:].reshape(3, 3)
        birth = rate+a@mu
        return np.r_[birth, (a@c+c@a.T+np.diag(birth)).ravel()]
    solved = solve_ivp(derivative, (0., 2.), np.zeros(12),
                       t_eval=times, rtol=2e-12, atol=2e-13)
    expected = np.zeros((9, 9))
    for j, tj in enumerate(times):
        for k, tk in enumerate(times):
            earlier = min(j, k)
            cs = solved.y[3:, earlier].reshape(3, 3)
            block = linalg.expm(a*abs(tj-tk))@cs
            expected[j*3:(j+1)*3, k*3:(k+1)*3] = block if j >= k else block.T
    np.testing.assert_allclose(answer["covariance"], expected, rtol=3e-11, atol=3e-12)
    assert np.linalg.eigvalsh(answer["covariance"]).min() > 0
    assert answer["diagnostics"]["quadrature_check_frobenius_difference"] < 1e-12
    total_weights = np.ones((1, 3))
    total = project_collection_covariance(model, [2.], [sparse.csr_matrix(total_weights)],
        provenance="all physical electrons at t=2")
    assert total["covariance"][0, 0] == pytest.approx(rate.sum()*2, rel=1e-12)


def test_flat_field_neighbor_covariance_and_mean_do_not_identify_noise_law():
    lam, kappa, t, s = 5., 1e-7, 2., 1.
    model = BoundaryMomentModel(np.full(2, lam), kappa, nearest_neighbor_edges((1, 2)))
    answer = project(model, [s, t], [0, 1], [[1., 0.], [0., 1.]])
    assert answer["covariance"][0, 1] == pytest.approx(2*kappa*lam**2*t*s, rel=5e-6)
    assert answer["covariance"][0, 0] < lam*s
    # An equally mean-compatible donor-routing response has half this flat A.
    # This explicitly demonstrates non-identification; it is not used instead.
    other = BoundaryMomentModel(model.rate_e_s, kappa/2, model.edges)
    np.testing.assert_allclose(model.original_leading_mean(t), other.original_leading_mean(t))
    alt = project(other, [s, t], [0, 1], [[1., 0.], [0., 1.]])
    assert answer["covariance"][0, 1]/alt["covariance"][0, 1] == pytest.approx(2., rel=3e-6)
    # Nonuniform fields also have distinct, conserved mean-compatible drifts.
    rate = np.array([2., 4., 7.])
    edges = nearest_neighbor_edges((1, 3))
    boundary = collection_boundary_generator(rate, .001, edges).toarray()
    lap = collection_boundary_generator(np.ones(3), .5, edges).toarray()
    donor = .001*lap@np.diag(rate)
    np.testing.assert_allclose(boundary@rate, donor@rate)
    np.testing.assert_allclose(donor.sum(axis=0), 0., atol=1e-18)
    assert not np.allclose(boundary@np.diag(rate)+np.diag(rate)@boundary.T,
                            donor@np.diag(rate)+np.diag(rate)@donor.T)


def test_shared_reference_full_cross_exposure_terms_and_no_floor():
    j = np.array([[1., .5], [-2., 1.], [1., -.5]])
    c = np.array([[4., .3], [.3, 1.]])
    actual = shared_calibration_covariance(j, c, source="synthetic full shared bias covariance")
    np.testing.assert_allclose(actual, j@c@j.T)
    assert actual[0, 1] != 0
    residual = np.array([1., -2., 3.])
    result = fixed_residual_covariance_application(residual, actual+np.eye(3),
        np.ones(3), covariance_scope="conditional synthetic counting and shared calibration")
    assert result["native_ERR_chi2"] == 14.
    assert result["raw_residual_rms"] == pytest.approx(np.sqrt(14/3))
    assert result["physical_validation"] is False
    with pytest.raises(linalg.LinAlgError):
        fixed_residual_covariance_application(residual, np.zeros((3, 3)),
            np.ones(3), covariance_scope="singular; no arbitrary covariance jitter")


def test_native_photon_chain_and_adc_read_chain_are_distinct():
    from strong_lensing_group_calibration import averaged_read_observation
    adc = np.array([[[100.], [200.]]])
    coefficients = np.array([0., 1., 1e-3])[:, None]
    q = adc+1e-3*adc**2
    native = averaged_read_observation(q, coefficients, adc)
    gain = 2.
    photon_weights = native["read_to_group_jacobian"].reshape(1, 2)/gain
    adc_weights = np.full((1, 2), native["native_group_linearity_derivative"].item()/2)
    assert not np.allclose(photon_weights*gain, adc_weights)
    photon = project(BoundaryMomentModel(np.array([3.]), 0., np.empty((0, 2), int)),
                     [1., 2.], [0], photon_weights)
    read = independent_adc_read_covariance(pixel_read_projection([0], adc_weights, 1), [4.],
                                           stage_assumption="synthetic raw ADC independent reads")
    expected = 3*photon_weights@np.array([[1., 1.], [1., 2.]])@photon_weights.T
    np.testing.assert_allclose(photon["covariance"], expected)
    assert read[0, 0] == pytest.approx(4**2/2*np.sum(adc_weights**2))


def test_parameter_dependent_gaussian_score_includes_covariance_logdet():
    data = np.array([2., -1.])
    def state(theta):
        mean = np.array([theta, theta**2])
        cov = np.array([[2.+theta, .2*theta], [.2*theta, 1.+theta**2]])
        return mean, cov
    theta = .7
    mean, covariance = state(theta)
    jacobian = np.array([[1.], [2*theta]])
    dc = np.array([[[1., .2], [.2, 2*theta]]])
    value, score = gaussian_nll_and_score(data-mean, covariance, jacobian, dc)
    step = 1e-5
    def nll(x):
        m, c = state(x)
        return gaussian_nll_and_score(data-m, c, np.zeros((2, 0)), np.zeros((0, 2, 2)))[0]
    assert np.isfinite(value)
    assert score[0] == pytest.approx((nll(theta+step)-nll(theta-step))/(2*step), rel=2e-9)


def test_missing_provenance_or_excessive_response_fails_explicitly():
    model = BoundaryMomentModel(np.array([2., 4.]), .1, nearest_neighbor_edges((1, 2)))
    with pytest.raises(ValueError, match="provenance"):
        project_collection_covariance(model, [1.], [sparse.eye(2)], provenance="")
    with pytest.raises(ValueError, match="accuracy/domain"):
        project(model, [100.], [0], [[1.]], maximum_degree=2)


def test_full_context_saturation_domain_blocks_unseen_core_extrapolation():
    rate = np.array([100., 1.])
    edges = nearest_neighbor_edges((1, 2))
    time_limited = BoundaryMomentModel(rate, 1e-6, edges, max_time_s=2.,
        domain_source="synthetic last unflagged individual core read")
    with pytest.raises(ValueError, match="time domain"):
        project(time_limited, [3.], [1], [[1.]])
    charge_limited = BoundaryMomentModel(rate, 1e-6, edges, charge_limit_e=np.array([200., 50.]),
        domain_source="synthetic full-context measured charge ceiling")
    # Only the faint wing is observed, but its charge graph includes the core.
    with pytest.raises(ValueError, match="full-well"):
        project(charge_limited, [3.], [1], [[1.]])
    valid = project(charge_limited, [1.], [1], [[1.]])
    assert valid["diagnostics"]["charge_ceiling_supplied"] is True
    assert valid["diagnostics"]["post_saturation_collection_or_overflow_law"] == "not supplied"


def test_explicit_harmonic_initializer_preserves_boundary_and_maximum_principle():
    rates = np.array([[1., 2., 4.], [2., np.nan, 8.], [4., -1., 16.]])
    good = np.isfinite(rates) & (rates >= 0)
    continued, receipt = harmonic_rate_initializer(rates, good)
    np.testing.assert_array_equal(continued[good], rates[good])
    assert continued.min() >= rates[good].min()
    assert continued.max() <= rates[good].max()
    lap = collection_boundary_generator(np.ones(9), .5, nearest_neighbor_edges((3, 3)))
    np.testing.assert_allclose((lap@continued.ravel())[~good.ravel()], 0., atol=2e-15)
    assert receipt["count"] == 2
    assert receipt["physical_validation"] is False


def test_first_order_native_hessian_input_agrees_with_resummed_small_response():
    rate = np.array([5., 9.])
    times = np.array([1., 2.])
    model = BoundaryMomentModel(rate, 1e-7, nearest_neighbor_edges((1, 2)))
    moments = within_pixel_read_covariance_first_order(model, times, [0])[:, :, 0]
    resummed = project(model, times, [0, 0], np.eye(2))["covariance"]
    np.testing.assert_allclose(moments, resummed, rtol=3e-11, atol=3e-11)
    # The first-order tensor's domain is checked instead of curing negative
    # variance with an arbitrary error floor.
    extreme = BoundaryMomentModel(rate, .2, model.edges)
    with pytest.raises(ValueError, match="PSD small-signal"):
        within_pixel_read_covariance_first_order(extreme, times, [0])
