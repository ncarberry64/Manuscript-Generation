"""Conditional collection-boundary counting covariance, not an optical state kernel.

The original electron mean is lambda*t + kappa*L(lambda**2)*t**2/2.
It does not determine a unique stochastic law.  This module makes one additional
law explicit: incident collection areas respond to neighboring stored counts.
The affine local drift A has A[p,q]=kappa*(lambda[p]+lambda[q]).  A*lambda
therefore has exactly the original leading mean.  A globally positive softmax
completion is supplied, but the implemented covariance is the local linear-noise
moment system, NOT an exact solution of that nonlinear counting process.

No observed ERR scaling, eigenvalue floor, arbitrary correlation kernel, Wick
closure, microscopic BHSM-state assumption, or instrument identification occurs.
All time and spatial projection maps must be supplied by their actual producer.
"""
from __future__ import annotations

from dataclasses import dataclass
import math
import numpy as np
from numpy.polynomial.legendre import leggauss
from scipy import linalg, sparse
from scipy.sparse.linalg import expm_multiply
from scipy.sparse.linalg import spsolve
from scipy.special import logsumexp


SCOPE = "conditional detector counting/readout law; not BHSM dressed-state statistics"


def nearest_neighbor_edges(shape):
    """Undirected edges on the declared finite no-exterior-edge context."""
    if len(shape) != 2 or min(shape) < 1 or any(int(v) != v for v in shape):
        raise ValueError("A positive integer two-dimensional context is required")
    ny, nx = map(int, shape)
    ids = np.arange(ny*nx).reshape(ny, nx)
    return np.concatenate((np.c_[ids[:, :-1].ravel(), ids[:, 1:].ravel()],
                           np.c_[ids[:-1, :].ravel(), ids[1:, :].ravel()]), axis=0)


def harmonic_rate_initializer(rate_e_s, valid):
    """Explicit conditional continuation for undefined context illumination.

    Solve L_UU*lambda_U=-L_UK*lambda_K on the same finite charge graph.
    This is a trial Dirichlet harmonic initializer, NOT a measured rate or a
    calibration law.  A discrete maximum principle preserves positive valid
    boundary data.  No unavailable slope is silently converted to zero.
    """
    rate = np.asarray(rate_e_s, float)
    known = np.asarray(valid, bool)
    if (rate.ndim != 2 or known.shape != rate.shape or not known.any()
            or not np.isfinite(rate[known]).all() or np.any(rate[known] < 0)):
        raise ValueError("At least one finite nonnegative context boundary rate required")
    answer = rate.ravel().copy()
    missing = np.flatnonzero(~known.ravel())
    if len(missing):
        fixed = np.flatnonzero(known.ravel())
        lap = collection_boundary_generator(np.ones(rate.size), .5,
                                            nearest_neighbor_edges(rate.shape))
        answer[missing] = spsolve(lap[missing][:, missing],
                                 -lap[missing][:, fixed]@answer[fixed])
        if not np.isfinite(answer).all() or np.any(answer < 0):
            raise ArithmeticError("Harmonic incident-rate continuation left its domain")
    return answer.reshape(rate.shape), {
        "meaning": "provisional harmonic context initializer, not measured flux",
        "continued_yx": np.argwhere(~known).tolist(),
        "continued_values_e_s": answer[missing].tolist(),
        "count": len(missing), "physical_validation": False}


def collection_boundary_generator(rate_e_s, kappa_e_inverse, edges):
    """Symmetric conserved affine boundary response, with units s^-1.

    These edge coefficients are an explicit conditional boundary prescription.
    They cannot be inferred uniquely from the original mean alone.
    """
    rate = np.asarray(rate_e_s, float)
    edge = np.asarray(edges, int).reshape(-1, 2)
    kappa = float(kappa_e_inverse)
    if (rate.ndim != 1 or not len(rate) or not np.isfinite(rate).all() or np.any(rate < 0)
            or not np.isfinite(kappa) or kappa < 0):
        raise ValueError("Finite nonnegative electron rates and kappa required")
    if (edge.size and (edge.min() < 0 or edge.max() >= len(rate)
                      or np.any(edge[:, 0] == edge[:, 1]))):
        raise ValueError("Edges must connect different in-domain pixels")
    if len(set(map(tuple, np.sort(edge, axis=1)))) != len(edge):
        raise ValueError("Repeated physical boundary")
    p, q = edge.T
    weight = kappa*(rate[p]+rate[q])
    rows = np.r_[p, q, p, q]
    cols = np.r_[q, p, p, q]
    return sparse.coo_matrix((np.r_[weight, weight, -weight, -weight],
                              (rows, cols)), shape=(len(rate), len(rate))).tocsr()


def positive_collection_rates(counts_e, incident_rate_e_s, generator):
    """One globally positive, total-rate conserving completion of the local law.

    a_p(N)=sum(lambda)*softmax(log(lambda)+(A*N)/lambda)_p.
    The value and first count derivative at zero are lambda and A.  Strictly
    positive incident rates are a domain requirement, never a fitted floor.
    Higher-order physics is NOT fixed by the leading-mean equation.
    """
    rate = np.asarray(incident_rate_e_s, float)
    counts = np.asarray(counts_e, float)
    a = sparse.csr_matrix(generator)
    if (counts.shape != rate.shape or rate.ndim != 1 or np.any(rate <= 0)
            or not np.isfinite(rate).all() or not np.isfinite(counts).all()
            or np.any(counts < 0) or a.shape != (len(rate), len(rate))):
        raise ValueError("Positive incident rates and nonnegative stored counts required")
    logits = np.log(rate)+np.asarray(a@counts)/rate
    return np.exp(np.log(rate.sum())+logits-logsumexp(logits))


@dataclass
class BoundaryMomentModel:
    """Local affine/LNA moments from zero reset, with conserved Poisson births.

    mu'=lambda+A*mu; C'=A*C+C*A.T+diag(lambda+A*mu).
    The deterministic birth rate exp(A*t)*lambda is nonnegative.  The resulting
    covariance is positive semidefinite, but the affine a(N)=lambda+A*N is not a
    globally valid birth propensity.  It is a local moment closure to the
    explicitly supplied positive completion.  Initial/reset/reference charge
    uncertainty must be supplied separately through its actual Jacobian.
    """
    rate_e_s: np.ndarray
    kappa_e_inverse: float
    edges: np.ndarray
    max_time_s: float | None = None
    charge_limit_e: np.ndarray | None = None
    domain_source: str = "unspecified small-signal domain; not a native full-well calibration"

    def __post_init__(self):
        self.rate_e_s = np.asarray(self.rate_e_s, float).copy()
        self.generator = collection_boundary_generator(
            self.rate_e_s, self.kappa_e_inverse, self.edges)
        if self.max_time_s is not None and (not np.isfinite(self.max_time_s) or self.max_time_s < 0):
            raise ValueError("Finite nonnegative small-signal maximum time required")
        if self.charge_limit_e is not None:
            self.charge_limit_e = np.broadcast_to(np.asarray(self.charge_limit_e, float),
                                                self.rate_e_s.shape).copy()
            if not np.isfinite(self.charge_limit_e).all() or np.any(self.charge_limit_e <= 0):
                raise ValueError("Finite positive declared charge ceilings required")
        if not str(self.domain_source).strip():
            raise ValueError("Small-signal domain/calibration source must be explicit")

    def _check_time(self, time_s):
        t = float(time_s)
        if not np.isfinite(t) or t < 0:
            raise ValueError("Nonnegative finite time required")
        if self.max_time_s is not None and t > self.max_time_s:
            raise ValueError("Read exceeds the declared unsaturated/small-signal time domain")
        return t

    def birth_rate(self, time_s):
        t = self._check_time(time_s)
        return np.asarray(expm_multiply(self.generator*t, self.rate_e_s))

    def mean(self, time_s):
        """Resummed local mean; no inverse of the singular conserved A."""
        t = self._check_time(time_s)
        n = len(self.rate_e_s)
        augmented = sparse.bmat([[self.generator, self.rate_e_s[:, None]],
                                 [None, sparse.csr_matrix((1, 1))]], format="csr")
        initial = np.r_[np.zeros(n), 1.]
        mean = np.asarray(expm_multiply(augmented*t, initial))[:n]
        if self.charge_limit_e is not None and np.any(mean >= self.charge_limit_e):
            raise ValueError("Mean stored charge reaches the declared full-well/domain ceiling")
        return mean

    def original_leading_mean(self, time_s):
        t = self._check_time(time_s)
        return self.rate_e_s*t+.5*t*t*(self.generator@self.rate_e_s)


def pixel_read_projection(selected_pixels, count_to_output_by_read, n_context):
    """Sparse per-read maps H_t, preserving caller-supplied physical units.

    count_to_output_by_read has observation-by-individual-read shape.  It may
    include native group Jacobians, electronic gain, flat/PHOTMJSR conversion,
    and the ACTUAL producing rate-estimator derivative.  Their provenance is
    required by project_collection_covariance.  Repeated selected pixels permit
    a full group/read observation vector without discarding cross-time terms.
    """
    pixels = np.asarray(selected_pixels, int)
    weights = np.asarray(count_to_output_by_read, float)
    if (pixels.ndim != 1 or weights.ndim != 2 or len(pixels) != weights.shape[0]
            or not np.isfinite(weights).all() or np.any(pixels < 0)
            or np.any(pixels >= n_context)):
        raise ValueError("Observation pixels and finite read projection required")
    rows = np.arange(len(pixels))
    return tuple(sparse.csr_matrix((weights[:, j], (rows, pixels)),
                                   shape=(len(pixels), n_context))
                 for j in range(weights.shape[1]))


def _projection_inputs(model, read_times_s, read_maps, provenance):
    times = np.asarray(read_times_s, float)
    maps = tuple(sparse.csr_matrix(h, dtype=float) for h in read_maps)
    if (times.ndim != 1 or not len(times) or np.any(times < 0)
            or not np.isfinite(times).all() or len(maps) != len(times)
            or np.any(np.diff(times) < 0)):
        raise ValueError("Sorted finite individual read times and one map each required")
    if not provenance or not str(provenance).strip():
        raise ValueError("Actual projection/units provenance is required")
    shape = maps[0].shape
    if (shape[1] != len(model.rate_e_s) or any(h.shape != shape for h in maps)
            or any(not np.isfinite(h.data).all() for h in maps)):
        raise ValueError("All read maps must share observation/context dimensions")
    return times, maps


def _operator_norm_bound(matrix):
    if matrix.nnz == 0:
        return 0.
    return math.sqrt(float(abs(matrix).sum(axis=0).max())*
                     float(abs(matrix).sum(axis=1).max()))


def _taylor_degree(a_norm_t, tolerance, maximum):
    """Conservative exponential-series remainder, in induced 2-norm."""
    if a_norm_t == 0:
        return 0, 0.
    for degree in range(maximum+1):
        log_bound = a_norm_t+(degree+1)*math.log(a_norm_t)-math.lgamma(degree+2)
        bound = math.exp(log_bound) if log_bound < 700 else math.inf
        if bound <= tolerance:
            return degree, bound
    raise ValueError("Response exceeds the declared sparse exponential accuracy/domain")


def project_collection_covariance(model, read_times_s, read_maps, *, provenance,
                                  quadrature_order=6, exponential_tolerance=1e-11,
                                  maximum_degree=30, check_quadrature=True):
    """Project interpixel AND cross-time covariance without a full read tensor.

    Z=sum_t H_t*N(t), B(u)=sum_{t>=u}H_t*exp(A*(t-u));
    Cov(Z)=integral B(u)*diag(exp(A*u)*lambda)*B(u).T du.
    This is a PSD integral of the declared local moment/LNA system.  It is not
    a diagonal-ERR replacement.  It includes every source of Poisson births on
    the declared full finite context, including those outside selected pixels.
    Quadrature is split at each actual read time.  Sparse Taylor error has a
    conservative absolute bound; a second quadrature order is an independent
    convergence diagnostic, not a fitted noise change.
    """
    times, maps = _projection_inputs(model, read_times_s, read_maps, provenance)
    # Check the full charge context, not merely observed unsaturated wing pixels.
    # A saturated unobserved core can still change its neighbors' collection.
    for t in np.unique(times):
        model.mean(t)
    if (int(quadrature_order) != quadrature_order or quadrature_order < 2
            or not 0 < exponential_tolerance < 1 or maximum_degree < 0):
        raise ValueError("Declared positive integration tolerances required")
    a = model.generator
    last = float(times[-1])
    norm_t = float(abs(a).sum(axis=1).max())*last
    degree, tail = _taylor_degree(norm_t, exponential_tolerance, maximum_degree)
    # Native pixelwise calibration has one physical electron-count pixel per
    # output row.  Factor the time weights before taking graph powers: keeping
    # 18 identical spatial power stacks would waste substantial context memory.
    single_pixel_rows = all(np.all(np.diff(h.indptr) == 1) for h in maps)
    if single_pixel_rows:
        columns = maps[0].indices
        single_pixel_rows = all(np.array_equal(h.indices, columns) for h in maps)
    if single_pixel_rows:
        read_weights = np.stack([h.data for h in maps])
        base = sparse.csr_matrix((np.ones(len(columns)),
                                 (np.arange(len(columns)), columns)), shape=maps[0].shape)
        base_powers = [base]
        for k in range(1, degree+1):
            base_powers.append((base_powers[-1]@a/k).tocsr())
        powers = None
    else:
        powers = []
        for h in maps:
            h_powers = [h]
            for k in range(1, degree+1):
                h_powers.append((h_powers[-1]@a/k).tocsr())
            powers.append(h_powers)
    output_size = maps[0].shape[0]

    def integrate(order):
        nodes, weights = leggauss(order)
        covariance = np.zeros((output_size, output_size))
        endpoints = np.unique(np.r_[0., times])
        min_birth = math.inf
        for lower, upper in zip(endpoints[:-1], endpoints[1:]):
            if upper == lower:
                continue
            for node, weight in zip(nodes, weights):
                u = (lower+upper)/2+(upper-lower)*node/2
                birth = model.birth_rate(u)
                min_birth = min(min_birth, float(birth.min(initial=math.inf)))
                if np.any(birth < 0) or not np.isfinite(birth).all():
                    raise ArithmeticError("Nonpositive/invalid conditional diffusion; no floor applied")
                b = sparse.csr_matrix((output_size, len(model.rate_e_s)))
                if single_pixel_rows:
                    active = times >= u
                    delta = times[active]-u
                    for k, power in enumerate(base_powers):
                        coefficient = np.sum(read_weights[active]*delta[:, None]**k, axis=0)
                        b = b+power.multiply(coefficient[:, None])
                else:
                    for t, h_powers in zip(times, powers):
                        if t < u:
                            continue
                        delta = t-u
                        for k, power in enumerate(h_powers):
                            if power.nnz:
                                b = b+power*delta**k
                weighted = b.multiply(np.sqrt(birth)[None, :])
                covariance += (upper-lower)*weight/2*(weighted@weighted.T).toarray()
        return (covariance+covariance.T)/2, min_birth

    covariance, min_birth = integrate(int(quadrature_order))
    quadrature_difference = None
    if check_quadrature:
        checked, _ = integrate(int(quadrature_order)+2)
        quadrature_difference = float(np.linalg.norm(checked-covariance, ord="fro"))
    h_norm = sum(_operator_norm_bound(h) for h in maps)
    dbound = h_norm*tail
    error_bound = last*float(model.rate_e_s.max(initial=0.))*(
        2*h_norm*dbound+dbound**2)
    mean = sum((h@model.mean(t) for t, h in zip(times, maps)),
               np.zeros(output_size))
    original_mean = sum((h@model.original_leading_mean(t) for t, h in zip(times, maps)),
                        np.zeros(output_size))
    return {"covariance": covariance, "mean": mean,
            "original_leading_mean": original_mean,
            "resummed_minus_original_mean": mean-original_mean,
            "diagnostics": {"scope": SCOPE, "projection_provenance": str(provenance),
                "closure": "local affine drift / linear-noise; not exact nonlinear softmax law",
                "counts_unit": "electron", "time_unit": "s", "kappa_unit": "electron^-1",
                "finite_context_pixels": len(model.rate_e_s),
                "no_exterior_edges": True, "initial_random_charge": "not included",
                "small_signal_domain_source": model.domain_source,
                "declared_maximum_time_s": model.max_time_s,
                "charge_ceiling_check": "necessary mean-domain check only; no random overflow law",
                "charge_ceiling_supplied": model.charge_limit_e is not None,
                "post_saturation_collection_or_overflow_law": "not supplied",
                "exponential_degree": degree, "generator_norm_times_max_time": norm_t,
                "shared_spatial_powers_across_reads": single_pixel_rows,
                "exponential_remainder_norm_bound": tail,
                "covariance_exponential_absolute_2norm_error_bound": error_bound,
                "quadrature_order": int(quadrature_order),
                "quadrature_check_frobenius_difference": quadrature_difference,
                "minimum_deterministic_birth_rate_e_s": min_birth if last else 0.,
                "noise_rescaling": False, "physical_validation": False}}


def independent_adc_read_covariance(read_projection, readnoise_cds_DN, *, stage_assumption):
    """Project an EXPLICIT conditional independent-single-read ADC noise law.

    read_projection[t] maps raw ADC-DN noise at individual read t to the output.
    Therefore native L'(mean ADC)/N, estimator weights, flat and units belong in
    this map; the photon-count inverse-linearity Jacobian must NOT be reused.
    sigma_single^2=CDS^2/2 requires the declared stage identification. The actual
    NRCB3 reference 0282 was made after reference AND IPC correction. It does
    NOT establish raw ADC independence or the variance before SIRS. Its marginal
    map cannot by itself determine raw covariance or justify adding a second
    copy of reference noise. stage_assumption must identify this approximation
    (or cite an independently supplied raw-ADC calibration). Shared/reference
    components are separate and must be checked against this stage convention.
    """
    projections = tuple(sparse.csr_matrix(h, dtype=float) for h in read_projection)
    sigma = np.asarray(readnoise_cds_DN, float)
    if (not str(stage_assumption).strip() or not projections or sigma.ndim != 1 or np.any(sigma < 0)
            or not np.isfinite(sigma).all()
            or any(h.shape != projections[0].shape or h.shape[1] != len(sigma)
                   or not np.isfinite(h.data).all() for h in projections)):
        raise ValueError("Actual ADC-noise projections and nonnegative CDS-DN map required")
    answer = np.zeros((projections[0].shape[0],)*2)
    for h in projections:
        factor = h.multiply((sigma/math.sqrt(2))[None, :])
        answer += (factor@factor.T).toarray()
    return (answer+answer.T)/2


def within_pixel_read_covariance_first_order(model, read_times_s, selected_pixels):
    """Explicit O(kappa) physical-count covariance at one pixel, unequal reads.

    For t>=s, C_pp(t,s)=lambda_p*s+A_pp*lambda_p*t*s
    +.5*(A*lambda)_p*s^2+O(kappa^2).  This supplies a stated leading-order
    Hessian contraction for the nonlinear native group mean. It is distinct
    from the resummed PSD covariance projector; no negative moment is floored.
    """
    times = np.asarray(read_times_s, float)
    pixels = np.asarray(selected_pixels, int)
    if (times.ndim != 1 or pixels.ndim != 1 or np.any(pixels < 0)
            or np.any(pixels >= len(model.rate_e_s))):
        raise ValueError("Explicit read times and in-context selected pixels required")
    for t in times:
        model.mean(t)
    rate = model.rate_e_s[pixels]
    derivative = model.generator.diagonal()[pixels]*rate
    drift = np.asarray(model.generator@model.rate_e_s)[pixels]
    earlier = np.minimum(times[:, None], times[None, :])
    product = times[:, None]*times[None, :]
    covariance = earlier[..., None]*rate+product[..., None]*derivative+.5*earlier[..., None]**2*drift
    for p in range(len(pixels)):
        eigen = np.linalg.eigvalsh(covariance[:, :, p])
        if eigen.min(initial=0.) < -1e-12*max(1., eigen.max(initial=0.)):
            raise ValueError("First-order read moment leaves the PSD small-signal domain")
    return covariance


def shared_calibration_covariance(jacobian, parameter_covariance, *, source):
    """Full shared J*C*J.T; repeated uses retain the SAME latent columns.

    A reference ERR map alone supplies marginal variances, not its missing
    off-diagonal covariance.  A caller choosing diagonal C must identify that
    additional approximation in source.  Never replace it by new per-use noise.
    """
    j = np.asarray(jacobian, float)
    c = np.asarray(parameter_covariance, float)
    if (j.ndim != 2 or c.shape != (j.shape[1],)*2 or not str(source).strip()
            or not np.isfinite(j).all() or not np.isfinite(c).all()
            or not np.allclose(c, c.T, rtol=1e-12, atol=0)):
        raise ValueError("Named shared-reference Jacobian and symmetric covariance required")
    if len(c) and linalg.eigvalsh(c, subset_by_index=[0, 0])[0] < 0:
        raise ValueError("Shared calibration covariance must be positive semidefinite")
    return j@c@j.T


def fixed_residual_covariance_application(residual, covariance, native_error, *,
                                         covariance_scope):
    """Preserve original residual evidence beside the conditional noise score.

    No p-value/admission flag is manufactured from a covariance change.  The
    Gaussian score is a declared approximation to counts/readout, unrelated to
    factorization of the BHSM dressed optical state.
    """
    r = np.asarray(residual, float)
    c = np.asarray(covariance, float)
    error = np.asarray(native_error, float)
    if (r.ndim != 1 or c.shape != (len(r),)*2 or error.shape != r.shape
            or np.any(error <= 0) or not np.isfinite(r).all()
            or not np.isfinite(c).all() or not np.isfinite(error).all()
            or not np.allclose(c, c.T, rtol=1e-12, atol=0)
            or not str(covariance_scope).strip()):
        raise ValueError("Finite residual/original errors and declared symmetric covariance required")
    chol = linalg.cholesky(c, lower=True)  # Deliberately no jitter/eigenvalue floor.
    whitened = linalg.solve_triangular(chol, r, lower=True)
    return {"number_of_observations": len(r),
            "raw_residual_rms": float(np.sqrt(np.mean(r*r))),
            "native_ERR_chi2": float(np.sum((r/error)**2)),
            "conditional_stochastic_chi2": float(whitened@whitened),
            "conditional_gaussian_logdet": float(2*np.log(np.diag(chol)).sum()),
            "covariance_scope": str(covariance_scope),
            "mean_was_changed": False, "error_rescaling": False,
            "physical_validation": False}


def gaussian_nll_and_score(residual, covariance, mean_jacobian, covariance_jacobian):
    """Full parameter-dependent local Gaussian likelihood score, including logdet.

    residual=data-mean; dC/dtheta must include rates, boundary response and every
    propagated calibration map that depends on theta.  Frozen-C least squares
    is not this likelihood.  Derivative tensors are explicit to avoid a silent
    frozen-noise optimization when a new physical mean changes shot noise.
    """
    r = np.asarray(residual, float)
    c = np.asarray(covariance, float)
    j = np.asarray(mean_jacobian, float)
    dc = np.asarray(covariance_jacobian, float)
    if (r.ndim != 1 or c.shape != (len(r),)*2 or j.ndim != 2
            or j.shape[0] != len(r) or dc.shape != (j.shape[1], len(r), len(r))
            or any(not np.isfinite(v).all() for v in (r, c, j, dc))
            or not np.allclose(c, c.T, rtol=1e-12, atol=0)
            or not np.allclose(dc, dc.transpose(0, 2, 1), rtol=1e-12, atol=0)):
        raise ValueError("Finite mean AND covariance parameter derivatives required")
    chol = linalg.cho_factor(c, lower=True)
    z = linalg.cho_solve(chol, r)
    precision = linalg.cho_solve(chol, np.eye(len(r)))
    score = -j.T@z+.5*np.einsum("ij,kji->k", precision, dc)-.5*np.einsum("i,kij,j->k", z, dc, z)
    nll = .5*(r@z+2*np.log(np.diag(chol[0])).sum()+len(r)*math.log(2*math.pi))
    return float(nll), score
