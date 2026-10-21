"""Training-only classical light-curve likelihood and genuinely held-out forecasts.

For image i, observed magnitude at observer time t is

    m_i(t) = m0 + b_i + q(t - d_i) + e_i(t) + noise_i(t).

q is one stationary OU/DRW process with covariance sigma_q**2 exp(-|dt|/tau_q).
The independent slow e_i have covariance sigma_e**2 exp(-|dt|/tau_e), and
b_i ~ N(0, sigma_b**2).  Quoted independent photometric error and a common
training-calibrated jitter add in quadrature.  This is a classical empirical
variability likelihood, not a BHSM statistical quantum-state contraction.

At shifted times u=t-d_i, stationary e_i(u+d_i) has the same within-image
covariance.  A (1+2*nimage)-dimensional diagonal-transition Kalman filter
therefore evaluates the exact normalized GP likelihood without forming the
nobservation-squared covariance.  A train-only RTS smoother gives marginal
predictions at the known held-out timestamps, including any points that
interleave in source time.  The full conditional test density is exactly
log p(train,test|frozen parameters)-log p(train|frozen parameters).

Delay sign convention: positive d_i means image i is observed later than
reference image 0 for the same intrinsic source feature.  All times are days.
Every fit, initializer, bound and mean level consumes training values only.
Test values enter only score_predictions; no fitting occurs there.
"""
from __future__ import annotations

from dataclasses import asdict, dataclass, replace

import numpy as np
from scipy.linalg import solve
from scipy.optimize import minimize


LOG2PI = float(np.log(2*np.pi))


@dataclass(frozen=True)
class LightCurves:
    """Flat observation arrays; image identifiers are consecutive 0,...,nimage-1."""

    time_days: np.ndarray
    image: np.ndarray
    magnitude: np.ndarray
    error_mag: np.ndarray
    nimage: int

    def __post_init__(self):
        arrays = (np.asarray(self.time_days, dtype=float),
                  np.asarray(self.image, dtype=int),
                  np.asarray(self.magnitude, dtype=float),
                  np.asarray(self.error_mag, dtype=float))
        if any(a.ndim != 1 for a in arrays) or len({a.size for a in arrays}) != 1:
            raise ValueError('Observation arrays must be equally sized vectors')
        if self.nimage < 1 or arrays[0].size == 0:
            raise ValueError('At least one image and observation are required')
        if any(not np.all(np.isfinite(a)) for a in arrays):
            raise ValueError('Every observation and uncertainty must be finite')
        if np.any(arrays[3] <= 0) or np.any(arrays[1] < 0) or np.any(arrays[1] >= self.nimage):
            raise ValueError('Invalid image identifier or nonpositive uncertainty')
        for name, value in zip(('time_days', 'image', 'magnitude', 'error_mag'), arrays):
            # Copies prevent changes to source arrays from changing a frozen fit.
            copied = value.copy()
            copied.setflags(write=False)
            object.__setattr__(self, name, copied)

    def subset(self, mask):
        return LightCurves(self.time_days[mask], self.image[mask],
                           self.magnitude[mask], self.error_mag[mask], self.nimage)


@dataclass(frozen=True)
class OUParameters:
    delays_days: tuple[float, ...]
    source_sigma_mag: float = .15
    source_tau_days: float = 200.
    extrinsic_sigma_mag: float = .05
    extrinsic_tau_days: float = 1500.
    jitter_mag: float = .005
    offset_sigma_mag: float = 3.
    magnitude_zero: float = 20.

    def __post_init__(self):
        # JSON callers provide lists; copy to a tuple so a frozen prediction
        # cannot be changed through an external mutable delay list.
        object.__setattr__(self, 'delays_days', tuple(float(x) for x in self.delays_days))
        if not self.delays_days or abs(self.delays_days[0]) > 1e-12:
            raise ValueError('Delay of reference image 0 must be zero')
        positive = (self.source_sigma_mag, self.source_tau_days,
                    self.extrinsic_sigma_mag, self.extrinsic_tau_days,
                    self.offset_sigma_mag)
        if any(not np.isfinite(x) or x <= 0 for x in positive):
            raise ValueError('OU scales and offset prior must be positive and finite')
        if self.jitter_mag < 0 or not np.isfinite(self.jitter_mag):
            raise ValueError('Jitter must be nonnegative and finite')
        if not np.isfinite(self.magnitude_zero) or not np.all(np.isfinite(self.delays_days)):
            raise ValueError('Mean level and delays must be finite')


def concatenate(first: LightCurves, second: LightCurves) -> LightCurves:
    if first.nimage != second.nimage:
        raise ValueError('Image domains must agree')
    return LightCurves(*(np.r_[getattr(first, key), getattr(second, key)]
                         for key in ('time_days', 'image', 'magnitude', 'error_mag')),
                       first.nimage)


def _state_setup(nimage, parameters):
    if len(parameters.delays_days) != nimage:
        raise ValueError('One delay per image is required')
    variance = np.r_[parameters.source_sigma_mag**2,
                     np.full(nimage, parameters.extrinsic_sigma_mag**2),
                     np.full(nimage, parameters.offset_sigma_mag**2)]
    rate = np.r_[1/parameters.source_tau_days,
                 np.full(nimage, 1/parameters.extrinsic_tau_days), np.zeros(nimage)]
    return variance, rate


def _filter(data, parameters, assimilate=None, save_states=False):
    """Exact filter; unassimilated values are never read, even in the smoother."""
    variance, rate = _state_setup(data.nimage, parameters)
    if assimilate is None:
        assimilate = np.ones(data.time_days.size, dtype=bool)
    assimilate = np.asarray(assimilate, dtype=bool)
    if assimilate.shape != data.time_days.shape:
        raise ValueError('Assimilation mask must match observation domain')
    source_time = data.time_days-np.asarray(parameters.delays_days)[data.image]
    order = np.argsort(source_time, kind='stable')
    time = source_time[order]
    state_size, count = variance.size, order.size
    mean = np.zeros(state_size)
    covariance = np.diag(variance)
    log_likelihood = 0.
    if save_states:
        predicted_mean = np.empty((count, state_size))
        filtered_mean = np.empty_like(predicted_mean)
        predicted_covariance = np.empty((count, state_size, state_size))
        filtered_covariance = np.empty_like(predicted_covariance)
        transitions = np.empty((count, state_size))
    for position, original in enumerate(order):
        dt = 0. if position == 0 else time[position]-time[position-1]
        transition = np.exp(-rate*dt)
        mean = transition*mean
        covariance = transition[:, None]*covariance*transition[None, :]
        covariance.flat[::state_size+1] += variance*(-np.expm1(-2*rate*dt))
        if save_states:
            predicted_mean[position] = mean
            predicted_covariance[position] = covariance
            transitions[position] = transition
        if assimilate[original]:
            image = data.image[original]
            selected = [0, 1+image, 1+data.nimage+image]
            projected = covariance[:, selected].sum(axis=1)
            noise = data.error_mag[original]**2+parameters.jitter_mag**2
            innovation_variance = projected[selected].sum()+noise
            residual = (data.magnitude[original]-parameters.magnitude_zero
                        -mean[selected].sum())
            log_likelihood -= .5*(LOG2PI+np.log(innovation_variance)
                                   +residual**2/innovation_variance)
            mean = mean+projected*(residual/innovation_variance)
            covariance = covariance-np.outer(projected, projected)/innovation_variance
            covariance = .5*(covariance+covariance.T)
        if save_states:
            filtered_mean[position] = mean
            filtered_covariance[position] = covariance
    result = {'log_likelihood': float(log_likelihood), 'order': order}
    if save_states:
        result.update(predicted_mean=predicted_mean, filtered_mean=filtered_mean,
                      predicted_covariance=predicted_covariance,
                      filtered_covariance=filtered_covariance,
                      transitions=transitions)
    return result


def log_likelihood(data: LightCurves, parameters: OUParameters) -> float:
    """Normalized, analytically marginalized training GP log likelihood."""
    return _filter(data, parameters)['log_likelihood']


def dense_covariance(data: LightCurves, parameters: OUParameters):
    """Independent O(N**2) reference, for small scientific validation only."""
    _state_setup(data.nimage, parameters)
    source_time = data.time_days-np.asarray(parameters.delays_days)[data.image]
    covariance = parameters.source_sigma_mag**2*np.exp(
        -np.abs(source_time[:, None]-source_time[None, :])/parameters.source_tau_days)
    same_image = data.image[:, None] == data.image[None, :]
    covariance += same_image*(parameters.offset_sigma_mag**2
        +parameters.extrinsic_sigma_mag**2*np.exp(
            -np.abs(data.time_days[:, None]-data.time_days[None, :])/parameters.extrinsic_tau_days))
    covariance.flat[::len(covariance)+1] += data.error_mag**2+parameters.jitter_mag**2
    return covariance


def predict_heldout(train: LightCurves, test: LightCurves, parameters: OUParameters):
    """Train-conditioned marginal prediction; test.magnitude is not consumed.

    Test timestamps, image labels and uncertainties are observed design data.
    Parameters must already be frozen.  Marginal variances include quoted
    error, jitter, source/extrinsic uncertainty and offset uncertainty.
    """
    all_data = concatenate(train, test)
    train_mask = np.arange(all_data.time_days.size) < train.time_days.size
    trace = _filter(all_data, parameters, assimilate=train_mask, save_states=True)
    mean = trace['filtered_mean'].copy()
    covariance = trace['filtered_covariance'].copy()
    for position in range(len(trace['order'])-2, -1, -1):
        # C=P_filtered A^T (P_predicted_next)^-1. Static offsets and stationary
        # OU states make P_predicted positive definite despite zero time steps.
        numerator = trace['filtered_covariance'][position]*trace['transitions'][position+1][None, :]
        gain = solve(trace['predicted_covariance'][position+1], numerator.T,
                     assume_a='pos', check_finite=False).T
        mean[position] += gain@(mean[position+1]-trace['predicted_mean'][position+1])
        covariance[position] += gain@(covariance[position+1]
            -trace['predicted_covariance'][position+1])@gain.T
        covariance[position] = .5*(covariance[position]+covariance[position].T)
    inverse_order = np.argsort(trace['order'])
    original_positions = inverse_order[train.time_days.size:]
    prediction_mean, prediction_variance = np.empty(test.time_days.size), np.empty(test.time_days.size)
    for j, position in enumerate(original_positions):
        image = test.image[j]
        selected = [0, 1+image, 1+test.nimage+image]
        prediction_mean[j] = parameters.magnitude_zero+mean[position, selected].sum()
        prediction_variance[j] = (covariance[position][np.ix_(selected, selected)].sum()
                                 +test.error_mag[j]**2+parameters.jitter_mag**2)
    if np.any(prediction_variance <= 0):
        raise ArithmeticError('Smoother returned a nonpositive predictive variance')
    return {'mean_magnitude': prediction_mean, 'variance_mag2': prediction_variance,
            'training_log_likelihood': trace['log_likelihood']}


def score_predictions(train: LightCurves, test: LightCurves, parameters: OUParameters):
    """One frozen application, including exact correlated conditional test density."""
    prediction = predict_heldout(train, test, parameters)
    residual = test.magnitude-prediction['mean_magnitude']
    variance = prediction['variance_mag2']
    standardized = residual/np.sqrt(variance)
    conditional_log_density = (log_likelihood(concatenate(train, test), parameters)
                               -prediction['training_log_likelihood'])
    marginal_log_density = -.5*(LOG2PI+np.log(variance)+standardized**2)

    def summary(mask):
        return {'count': int(np.count_nonzero(mask)),
                'rmse_mag': float(np.sqrt(np.mean(residual[mask]**2))),
                'mean_residual_mag': float(np.mean(residual[mask])),
                'mean_marginal_negative_log_density': float(-np.mean(marginal_log_density[mask])),
                'coverage_68_27_percent': float(np.mean(np.abs(standardized[mask]) <= 1)),
                'coverage_95_percent': float(np.mean(np.abs(standardized[mask]) <= 1.959963984540054)),
                'mean_predictive_sigma_mag': float(np.mean(np.sqrt(variance[mask])))}

    result = {'parameters': asdict(parameters),
              'training_log_likelihood': prediction['training_log_likelihood'],
              'joint_conditional_test_log_density': float(conditional_log_density),
              'joint_conditional_negative_log_density_per_observation': float(-conditional_log_density/len(residual)),
              'overall': summary(np.ones(len(residual), dtype=bool)),
              'by_image': {str(i): summary(test.image == i) for i in range(test.nimage)
                           if np.any(test.image == i)},
              'predicted_magnitude': prediction['mean_magnitude'],
              'predictive_sigma_mag': np.sqrt(variance),
              'residual_mag': residual}
    return result


def _decode(coordinates, template):
    nimage = len(template.delays_days)
    source_sigma, source_tau, extrinsic_sigma, extrinsic_tau, jitter = np.exp(coordinates[nimage-1:])
    return replace(template, delays_days=tuple(np.r_[0., coordinates[:nimage-1]]),
                   source_sigma_mag=float(source_sigma), source_tau_days=float(source_tau),
                   extrinsic_sigma_mag=float(extrinsic_sigma),
                   extrinsic_tau_days=float(extrinsic_tau), jitter_mag=float(jitter))


def _encode(parameters):
    return np.r_[parameters.delays_days[1:], np.log([parameters.source_sigma_mag,
        parameters.source_tau_days, parameters.extrinsic_sigma_mag,
        parameters.extrinsic_tau_days, parameters.jitter_mag])]


def fit_training(train: LightCurves, delay_starts, *, delay_bounds=(-30., 30.),
                 maxiter=200, magnitude_zero=20., offset_sigma_mag=3.,
                 hyperparameter_bounds=((.02, .5), (10., 1500.), (.005, .5),
                                        (200., 10000.), (.0001, .03)),
                 initial_hyperparameters=(.15, 200., .05, 1500., .005)):
    """Bounded training-only empirical Bayes fit with declared deterministic starts.

    delay_starts is an explicit list of full vectors including reference zero.
    No published full-data delay or held-out value is used by this function.
    The return contains every optimizer start/result, chosen training maximum,
    finite-difference gradient diagnostic, and active-bound flags. Point-fit
    predictive intervals condition on these fitted parameters; they do not
    integrate hyperparameter or macro-model uncertainty.
    """
    starts = [tuple(float(x) for x in row) for row in delay_starts]
    if not starts or any(len(row) != train.nimage or abs(row[0]) > 1e-12 for row in starts):
        raise ValueError('Provide at least one explicit full delay initializer with reference zero')
    if len(hyperparameter_bounds) != 5 or len(initial_hyperparameters) != 5:
        raise ValueError('Five ordered OU/jitter bounds and initial values are required')
    if initial_hyperparameters[-1] <= 0:
        raise ValueError('Fit logarithmically parameterizes a positive jitter')
    bounds = [tuple(delay_bounds)]*(train.nimage-1)+[
        tuple(np.log(row)) for row in hyperparameter_bounds]
    template = OUParameters(starts[0], *initial_hyperparameters,
                            offset_sigma_mag=offset_sigma_mag, magnitude_zero=magnitude_zero)

    def objective(coordinates):
        return -log_likelihood(train, _decode(coordinates, template))

    fits = []
    for delays in starts:
        initial = _encode(replace(template, delays_days=delays))
        optimum = minimize(objective, initial, method='L-BFGS-B', bounds=bounds,
                           options={'maxiter': int(maxiter), 'ftol': 1e-10,
                                    'gtol': 1e-5, 'eps': 1e-5,
                                    'maxls': 30})
        fits.append({'initial_delays_days': delays, 'coordinates': optimum.x.tolist(),
                     'training_log_likelihood': float(-optimum.fun),
                     'optimizer_success': bool(optimum.success),
                     'optimizer_message': str(optimum.message),
                     'iterations': int(optimum.nit), 'function_evaluations': int(optimum.nfev),
                     'maximum_absolute_optimizer_gradient': float(np.max(np.abs(optimum.jac)))})
    best = max(fits, key=lambda fit: fit['training_log_likelihood'])
    coordinates = np.asarray(best['coordinates'])
    parameters = _decode(coordinates, template)
    active = [bool(min(abs(x-lo), abs(x-hi)) < 1e-5*(hi-lo))
              for x, (lo, hi) in zip(coordinates, bounds)]
    return parameters, {'parameters': asdict(parameters), 'starts': fits,
                        'selected_start_index': fits.index(best),
                        'training_observation_count': int(train.time_days.size),
                        'coordinate_order': [f'delay_image_{i}_days' for i in range(1, train.nimage)]
                            +['log_source_sigma', 'log_source_tau', 'log_extrinsic_sigma',
                              'log_extrinsic_tau', 'log_jitter'],
                        'coordinate_bounds': bounds, 'active_bound_flags': active,
                        'training_log_likelihood': best['training_log_likelihood'],
                        'predictive_uncertainty': 'Conditional on fitted delays/hyperparameters; no parameter-posterior integration'}


def training_delay_curvature(train: LightCurves, parameters: OUParameters, step_days=.02):
    """Local fixed-hyperparameter delay curvature, not calibrated confidence limits.

    OU kernels have derivative cusps when shifted samples coincide, so the
    chosen step and eigenvalues are returned and must be inspected. Profiles
    or posterior marginalization are needed for definitive delay uncertainty.
    """
    delays = np.asarray(parameters.delays_days)
    dimensions = len(delays)-1
    matrix = np.empty((dimensions, dimensions))
    base = -log_likelihood(train, parameters)

    def nll(offset):
        return -log_likelihood(train, replace(parameters, delays_days=tuple(delays+offset)))

    for i in range(dimensions):
        ei = np.zeros_like(delays)
        ei[i+1] = step_days
        matrix[i, i] = (nll(ei)-2*base+nll(-ei))/step_days**2
        for j in range(i):
            ej = np.zeros_like(delays)
            ej[j+1] = step_days
            matrix[i, j] = matrix[j, i] = (nll(ei+ej)-nll(ei-ej)-nll(-ei+ej)+nll(-ei-ej))/(4*step_days**2)
    eigenvalues = np.linalg.eigvalsh(matrix)
    result = {'step_days': float(step_days), 'fixed_hyperparameter_hessian_per_day2': matrix.tolist(),
              'eigenvalues_per_day2': eigenvalues.tolist(),
              'interpretation': 'Local fixed-hyperparameter diagnostic; not calibrated delay confidence limits'}
    if np.all(eigenvalues > 0):
        covariance = np.linalg.inv(matrix)
        result['local_inverse_curvature_covariance_days2'] = covariance.tolist()
        result['local_inverse_curvature_sigma_days'] = np.sqrt(np.diag(covariance)).tolist()
    return result
