"""Regular compact-matter revisions of the conditional HE 0435 macro model.

The local matter is a positive sum of two finite-mass Plummer components.
The original nonlinear R1 scalar equation consumes their TOTAL enclosed mass.
The nuclear component is a prescribed regular matter surrogate, not a black
hole inference or a completed equation-of-state solution. No unused light-curve
target is read by this module. The old macro implementation remains unchanged.
"""
from __future__ import annotations

from dataclasses import dataclass, asdict
import math

import numpy as np
from numpy.polynomial import Polynomial
from scipy.optimize import least_squares, minimize, root

import strong_lensing_observed_macro as old
from strong_lensing_comparison import ARCSEC, C_KM_S, MPC_OVER_C_DAYS, lens_application, source_angle
from strong_lensing_local_family import G_MPC_KMS2_MSUN
from strong_lensing_robustness import ProfileR1Lens, DirectEffectiveHaloLens

PARAMETER_NAMES = old.PARAMETER_NAMES
LOWER = old.LOWER.copy()
UPPER = old.UPPER.copy()
COMPACT_SCALE_KPC = .03
COMPACT_MASS_CANDIDATES_MSUN = (0., 1e8, 10**8.5, 1e9, 10**9.5)


@dataclass(frozen=True)
class TwoPlummerLens:
    outer_mass_Msun: float
    outer_scale_Mpc: float
    compact_mass_Msun: float = 0.
    compact_scale_Mpc: float = COMPACT_SCALE_KPC/1000
    centre_Mpc: tuple = (0., 0.)

    def __post_init__(self):
        if min(self.outer_mass_Msun, self.outer_scale_Mpc, self.compact_scale_Mpc) <= 0 or self.compact_mass_Msun < 0:
            raise ValueError('positive outer mass/scales and nonnegative regular compact mass required')
        if not np.isfinite([self.outer_mass_Msun, self.outer_scale_Mpc, self.compact_mass_Msun,
                            self.compact_scale_Mpc, *self.centre_Mpc]).all():
            raise ValueError('finite prescribed matter parameters required')

    @property
    def mass_Msun(self):
        return self.outer_mass_Msun+self.compact_mass_Msun

    @property
    def radius_Mpc(self):
        # A projection reference scale, not a material edge.
        return self.outer_scale_Mpc

    @property
    def mass_length_Mpc(self):
        return G_MPC_KMS2_MSUN*self.mass_Msun/C_KM_S**2

    @property
    def b(self):
        return 4*self.mass_length_Mpc/self.outer_scale_Mpc

    @property
    def components(self):
        return ((self.outer_mass_Msun, self.outer_scale_Mpc),
                (self.compact_mass_Msun, self.compact_scale_Mpc))

    def mass_and_slope(self, r):
        r = float(r)
        if r < 0:
            raise ValueError('nonnegative physical radius required')
        mr3 = sum(m/self.mass_Msun/(r*r+a*a)**1.5 for m, a in self.components)
        mu = r**3*mr3
        if r == 0:
            slope = 3.
        else:
            weighted = sum(m/self.mass_Msun/(r*r+a*a)**1.5*3*a*a/(r*r+a*a)
                           for m, a in self.components)
            slope = weighted/mr3
        return mu, slope, mr3

    def density_Msun_per_Mpc3(self, r):
        return sum(3*m*a*a/(4*np.pi*(r*r+a*a)**2.5) for m, a in self.components)

    def application(self, x):
        y = np.asarray(x, float)-np.asarray(self.centre_Mpc)
        r2 = float(y@y)
        potential = 0.
        alpha = np.zeros(2)
        hessian = np.zeros((2, 2))
        for mass, scale in self.components:
            m = G_MPC_KMS2_MSUN*mass/C_KM_S**2
            denominator = r2+scale*scale
            potential += 2*m*math.log(denominator/(2*scale*scale))
            alpha += 4*m*y/denominator
            hessian += 4*m*(np.eye(2)/denominator-2*np.outer(y, y)/denominator**2)
        return potential, alpha, hessian


def body_from_parameters(parameters, compact_mass_Msun=0., *, compact_scale_kpc=COMPACT_SCALE_KPC):
    p = np.asarray(parameters, float)
    d = old.frozen_distances()
    return TwoPlummerLens(10**p[0], 10**p[1]/1000, float(compact_mass_Msun),
                          compact_scale_kpc/1000, tuple(d['ol']*ARCSEC*p[7:9]))


def macro_application(parameters, *, compact_mass_Msun=0., compact_scale_kpc=COMPACT_SCALE_KPC,
                      model='conventional', frozen_extra=None, astrometry=None,
                      require_success=True):
    p = np.asarray(parameters, float)
    d = old.frozen_distances()
    observed = old.load_astrometry() if astrometry is None else astrometry
    body = body_from_parameters(p, compact_mass_Msun, compact_scale_kpc=compact_scale_kpc)
    if model == 'conventional':
        lens = body
    elif model == 'r1':
        lens = (ProfileR1Lens(body, old.frozen_coefficients()) if frozen_extra is None
                else old.CompositeContinuationLens(body, frozen_extra))
    else:
        raise ValueError('model must be conventional or r1')
    geometry, U = old.screen_geometry(*p[4:7], distances=d)
    source = d['os']*ARCSEC*p[2:4]
    images = []
    applications = []
    failures = []
    for seed in observed['observed_images_light_centered_arcsec']:
        def fun(q):
            return lens_application(q*ARCSEC, geometry, lens, source)['residual']/ARCSEC
        def jac(q):
            return lens_application(q*ARCSEC, geometry, lens, source)['jacobian']
        solved = root(fun, seed, jac=jac, tol=2e-10)
        defect = float(np.max(abs(fun(solved.x))))
        if defect > 2e-7 or not np.isfinite(solved.x).all():
            failures.append(defect)
        images.append(solved.x)
        applications.append(lens_application(solved.x*ARCSEC, geometry, lens, source))
    images = np.asarray(images)
    separation = min(np.linalg.norm(images[i]-images[j]) for i in range(4) for j in range(i))
    if separation < .1:
        failures.append(float(.1-separation))
    if failures and require_success:
        raise ArithmeticError('requested bright image branches failed or merged')
    arrival = np.array([app['arrival_days'] for app in applications])
    relative = np.r_[(images[1:]-images[0]).ravel(), -images[0]]
    return dict(images_light_centered_arcsec=images, relative_astrometry_arcsec=relative,
                delays_days=arrival[1:]-arrival[0], arrival_days=arrival,
                max_lens_residual_rad=float(max(np.max(abs(app['residual'])) for app in applications)),
                max_arrival_stationarity_Mpc=float(max(np.max(abs(app['stationarity'])) for app in applications)),
                minimum_image_separation_arcsec=float(separation), branch_failures=failures,
                geometry=geometry, source=source, lens=lens, U=U)


def conventional_census(parameters, *, compact_mass_Msun=0., compact_scale_kpc=COMPACT_SCALE_KPC):
    """Enumerate the degree-nine radial polynomial; verify original 2D roots.

    Angular coordinates are arcseconds, z=|theta-center|^2. In the eigenframe
    of G, k(z)=A1/(z+s1^2)+A2/(z+s2^2), Q=(z+s1^2)(z+s2^2),
    P=A1(z+s2^2)+A2(z+s1^2), Zj=Q-gj P. The original equations imply
    z Z1^2 Z2^2-Q^2(beta1^2 Z2^2+beta2^2 Z1^2)=0.
    Multiplication introduces poles; original-equation checks remove them.
    """
    p = np.asarray(parameters, float)
    d = old.frozen_distances()
    body = body_from_parameters(p, compact_mass_Msun, compact_scale_kpc=compact_scale_kpc)
    geometry, _ = old.screen_geometry(*p[4:7], distances=d)
    source = d['os']*ARCSEC*p[2:4]
    eigenvalues, eigenvectors = np.linalg.eigh(geometry.G)
    center = p[7:9]
    beta = eigenvectors.T@(source_angle(geometry, source)/ARCSEC-center)
    amplitudes = np.array([4*G_MPC_KMS2_MSUN*m/C_KM_S**2/(d['ol']*ARCSEC**2)
                           for m, a in body.components])
    scales2 = np.array([(a/(d['ol']*ARCSEC))**2 for m, a in body.components])
    x = Polynomial([0., 1.])
    q = (x+scales2[0])*(x+scales2[1])
    numerator = amplitudes[0]*(x+scales2[1])+amplitudes[1]*(x+scales2[0])
    z1, z2 = [q-g*numerator for g in eigenvalues]
    polynomial = x*z1*z1*z2*z2-q*q*(beta[0]**2*z2*z2+beta[1]**2*z1*z1)
    candidate_roots = polynomial.roots()
    seeds = []
    for value in candidate_roots:
        if abs(value.imag) > 2e-6*max(1., abs(value.real)) or value.real < 0:
            continue
        z = float(value.real)
        k = float(np.sum(amplitudes/(z+scales2)))
        denominators = 1-eigenvalues*k
        if np.min(abs(denominators)) > 1e-10:
            seeds.append(center+eigenvectors@(beta/denominators))
        else:
            # A source lying on an eigenaxis may have valid roots exactly at
            # a canceled pole; reconstruct the missing component explicitly.
            for axis in range(2):
                other = 1-axis
                if abs(beta[axis]) < 1e-10 and abs(denominators[other]) > 1e-10:
                    component = beta[other]/denominators[other]
                    remaining = z-component*component
                    if remaining >= 0:
                        for sign in (-1., 1.):
                            u = np.zeros(2)
                            u[other] = component
                            u[axis] = sign*math.sqrt(remaining)
                            seeds.append(center+eigenvectors@u)
    # Add central/eigenaxis logarithmic seeds as a numerical audit of the
    # polynomial reconstruction, not as a replacement for the enumeration.
    angles = np.linspace(0., 2*np.pi, 16, endpoint=False)
    for radius in np.geomspace(max(math.sqrt(min(scales2))*1e-4, 1e-7), 6., 18):
        seeds.extend(center+radius*np.c_[np.cos(angles), np.sin(angles)])
    seeds.append(center)
    images, applications = _verified_roots(seeds, geometry, body, source)
    return _census_result(images, applications, body, polynomial_degree=len(polynomial.coef)-1,
                          polynomial_positive_candidates=int(sum(abs(z.imag) < 2e-6*max(1., abs(z.real)) and z.real >= 0 for z in candidate_roots)),
                          scope='degree-nine conventional radial polynomial plus logarithmic-seed audit; all roots checked in original 2D equations')


def _verified_roots(seeds, geometry, lens, source):
    images = []
    for seed in seeds:
        def fun(q):
            return lens_application(q*ARCSEC, geometry, lens, source)['residual']/ARCSEC
        def jac(q):
            return lens_application(q*ARCSEC, geometry, lens, source)['jacobian']
        try:
            solved = root(fun, seed, jac=jac, tol=2e-11)
            if np.isfinite(solved.x).all() and max(abs(fun(solved.x))) < 2e-8 and not any(np.linalg.norm(solved.x-other) < 1e-7 for other in images):
                images.append(solved.x)
        except (ArithmeticError, ValueError, FloatingPointError):
            continue
    images = np.asarray(sorted(images, key=lambda image: float(np.linalg.norm(image-np.asarray(lens.centre_Mpc)/geometry.D_ol[0, 0]/ARCSEC))))
    applications = [lens_application(image*ARCSEC, geometry, lens, source) for image in images]
    return images, applications


def _census_result(images, applications, lens, **extra):
    if not len(images):
        raise ArithmeticError('no verified roots in declared global image census')
    determinants = np.array([np.linalg.det(app['jacobian']) for app in applications])
    return dict(images_light_centered_arcsec=images.tolist(), image_count=len(images),
                absolute_magnifications=(1/abs(determinants)).tolist(),
                signed_magnifications=(1/determinants).tolist(),
                arrival_days=[app['arrival_days'] for app in applications],
                max_lens_residual_rad=float(max(np.max(abs(app['residual'])) for app in applications)),
                max_arrival_stationarity_Mpc=float(max(np.max(abs(app['stationarity'])) for app in applications)),
                min_jacobian_singular_value=float(min(np.linalg.svd(app['jacobian'], compute_uv=False)[-1] for app in applications)),
                **extra)


def image_census(parameters, *, compact_mass_Msun=0., compact_scale_kpc=COMPACT_SCALE_KPC,
                 model='conventional', projection_count=260):
    conventional = conventional_census(parameters, compact_mass_Msun=compact_mass_Msun,
                                        compact_scale_kpc=compact_scale_kpc)
    if model == 'conventional':
        return conventional
    if model != 'r1':
        raise ValueError('model must be conventional or r1')
    p = np.asarray(parameters, float)
    d = old.frozen_distances()
    body = body_from_parameters(p, compact_mass_Msun, compact_scale_kpc=compact_scale_kpc)
    screened = ProfileR1Lens(body, old.frozen_coefficients())
    projection = old.FrozenExtraProjection(screened, count=projection_count)
    fast = old.CompositeContinuationLens(body, projection)
    geometry, _ = old.screen_geometry(*p[4:7], distances=d)
    source = d['os']*ARCSEC*p[2:4]
    seeds = list(np.asarray(conventional['images_light_centered_arcsec']))
    angles = np.linspace(0., 2*np.pi, 12, endpoint=False)
    for radius in np.geomspace(max(compact_scale_kpc/1000/(d['ol']*ARCSEC)*1e-4, 1e-7), 6., 15):
        seeds.extend(p[7:9]+radius*np.c_[np.cos(angles), np.sin(angles)])
    images, _ = _verified_roots(seeds, geometry, fast, source)
    # Only final roots use the expensive original total-source R1 projection.
    images, applications = _verified_roots(images, geometry, screened, source)
    return _census_result(images, applications, screened,
                         conventional_image_count=conventional['image_count'],
                         scope='continuation of all conventional roots plus logarithmic-seed audit, verified against original total-source R1 projection; no theorem excluding remote caustic root creation')


def source_section_checks(screened, radii_Mpc):
    """Evaluate original integrated scalar flux and metric Gauss applications."""
    c = screened.coefficients
    length = c.hubble_length_Mpc
    sections = []
    for r in np.asarray(radii_Mpc, float):
        mu, _, _ = screened.mass_and_slope(float(r))
        x = r/length
        gradient = screened.scalar_gradient_natural(float(r))
        mass = screened.mass_length_Mpc/length*mu
        density = screened.base.density_Msun_per_Mpc3(float(r))
        density_natural = 8*np.pi*G_MPC_KMS2_MSUN/C_KM_S**2*density*length**2
        scalar_lhs = c.C*x*x*gradient-4/c.lambda_over_MplH02*x*gradient**2
        phi_prime = length*screened.mass_length_Mpc*mu/r**2*(1+screened.extra_force(float(r)))
        metric_lhs = 2*x*x*phi_prime+c.b*x*x*gradient
        sections.append(dict(radius_kpc=float(r*1000), density_over_Mpl2H02=float(density_natural),
                             consumed_reduced_source_J_eff=float(-.5*c.b*density_natural),
                             scalar_gradient_over_MplH0=float(gradient), scalar_integrated_lhs=float(scalar_lhs),
                             scalar_integrated_source=float(-c.b*mass), metric_Gauss_lhs=float(metric_lhs),
                             metric_Gauss_source=float(2*mass),
                             scalar_relative_residual=float(abs(scalar_lhs+c.b*mass)/(c.b*mass)),
                             metric_relative_residual=float(abs(metric_lhs-2*mass)/(2*mass))))
    return sections


def fit_macro(train_delays_days, train_delay_covariance_days2, *, compact_mass_Msun=0.,
              compact_scale_kpc=COMPACT_SCALE_KPC, model='conventional', initial=None,
              astrometry=None, max_nfev=400, refresh_steps=2, delay_log_likelihood=None,
              exact_maxiter=60, include_census=True, skip_quadratic_initializer=False,
              retain_best_original_endpoint=False):
    """Fit only old development observations and the prescribed compact mass.

    The optional likelihood callback takes the three B/C/D-minus-A delays.
    It must close over development data only. Nuclear mass/scale are fixed in
    each candidate; source, screen and center bounds are the original nine.
    """
    delay = np.asarray(train_delays_days, float)
    covariance = np.asarray(train_delay_covariance_days2, float)
    if delay.shape != (3,) or covariance.shape != (3, 3):
        raise ValueError('three historical delays and 3x3 covariance required')
    if not np.allclose(covariance, covariance.T, rtol=1e-10, atol=1e-12):
        raise ValueError('symmetric development delay covariance required')
    whiten = np.linalg.solve(np.linalg.cholesky(covariance), np.eye(3))
    observed = old.load_astrometry() if astrometry is None else astrometry
    p = (np.array([11.72, .74, -.01, .052, .09, .173, .091, 0., 0.])
         if initial is None else np.asarray(initial, float).copy())
    p = np.clip(p, LOWER+1e-9, UPPER-1e-9)
    initial_parameters = p.copy()
    extra = None
    trace = []
    exact_trace = []
    if skip_quadratic_initializer and delay_log_likelihood is None:
        raise ValueError('exact-only warm start requires a development likelihood callback')
    for iteration in range(0 if skip_quadratic_initializer else (1 if model == 'conventional' else refresh_steps+1)):
        if model == 'r1':
            extra = old.FrozenExtraProjection(ProfileR1Lens(body_from_parameters(p, compact_mass_Msun, compact_scale_kpc=compact_scale_kpc), old.frozen_coefficients()), count=260)
        def apply(parameters):
            return macro_application(parameters, compact_mass_Msun=compact_mass_Msun,
                                     compact_scale_kpc=compact_scale_kpc, model=model,
                                     frozen_extra=extra, astrometry=observed)
        def residual(parameters):
            try:
                app = apply(parameters)
                return np.r_[(app['relative_astrometry_arcsec']-observed['observed_relative_arcsec'])/observed['sigma_arcsec'],
                             whiten@(app['delays_days']-delay), parameters[7:9]/.003]
            except (ArithmeticError, ValueError):
                return np.r_[np.full(11, 1e6), parameters[7:9]/.003]
        fitted = least_squares(residual, p, bounds=(LOWER, UPPER), max_nfev=max_nfev,
                               ftol=2e-11, xtol=2e-11, gtol=2e-7, diff_step=2e-5, x_scale='jac')
        p = fitted.x
        trace.append(dict(iteration=iteration, chi2=float(fitted.fun@fitted.fun),
                          success=bool(fitted.success), message=str(fitted.message),
                          evaluations=int(fitted.nfev), parameters=p.tolist()))
    quadratic_parameters = p.copy()
    endpoints = [p.copy()]
    if delay_log_likelihood is not None:
        for iteration in range(1 if model == 'conventional' else 2):
            if model == 'r1':
                extra = old.FrozenExtraProjection(ProfileR1Lens(body_from_parameters(p, compact_mass_Msun, compact_scale_kpc=compact_scale_kpc), old.frozen_coefficients()), count=260)
            def objective(parameters):
                try:
                    app = macro_application(parameters, compact_mass_Msun=compact_mass_Msun,
                                            compact_scale_kpc=compact_scale_kpc, model=model,
                                            frozen_extra=extra, astrometry=observed)
                    astro = (app['relative_astrometry_arcsec']-observed['observed_relative_arcsec'])/observed['sigma_arcsec']
                    return float(astro@astro+np.sum((parameters[7:9]/.003)**2)-2*delay_log_likelihood(app['delays_days']))
                except (ArithmeticError, ValueError):
                    return 1e15
            # These are numerical coordinates, not additional priors/bounds.
            # A raw arcsecond centroid and log-mass have very different
            # curvature; scaling avoids a premature line-search plateau.
            origin = p.copy()
            scales = np.array([.03, .03, .01, .01, .1, .05, .05, .003, .003])
            qlower, qupper = (LOWER-origin)/scales, (UPPER-origin)/scales
            def scaled_objective(q):
                # A harmless common constant scale also keeps the first
                # quasi-Newton trial on the existing quartet, instead of
                # jumping into the rejected-branch penalty plateau.
                return objective(origin+scales*q)/1000.
            def scaled_gradient(q):
                gradient = []
                for axis in range(9):
                    plus = q.copy(); minus = q.copy()
                    plus[axis] = min(qupper[axis], q[axis]+2e-5)
                    minus[axis] = max(qlower[axis], q[axis]-2e-5)
                    gradient.append((scaled_objective(plus)-scaled_objective(minus))/(plus[axis]-minus[axis]))
                return np.array(gradient)
            optimized = minimize(scaled_objective, np.zeros(9), method='L-BFGS-B', jac=scaled_gradient,
                                 bounds=list(zip(qlower, qupper)),
                                 options=dict(maxiter=exact_maxiter, ftol=2e-12, gtol=1e-7, maxls=30))
            p = origin+scales*optimized.x
            endpoints.append(p.copy())
            exact_trace.append(dict(iteration=iteration, objective=float(optimized.fun*1000),
                                    success=bool(optimized.success), message=str(optimized.message),
                                    evaluations=int(optimized.nfev), parameters=p.tolist()))
    endpoint_checks = []
    if retain_best_original_endpoint:
        if delay_log_likelihood is None:
            raise ValueError('endpoint comparison requires the exact development likelihood')
        for index, endpoint in enumerate(endpoints):
            try:
                app = macro_application(endpoint, compact_mass_Msun=compact_mass_Msun,
                                        compact_scale_kpc=compact_scale_kpc, model=model,
                                        astrometry=observed)
                astro_endpoint = (app['relative_astrometry_arcsec']-observed['observed_relative_arcsec'])/observed['sigma_arcsec']
                value = float(astro_endpoint@astro_endpoint+np.sum((endpoint[7:9]/.003)**2)-2*delay_log_likelihood(app['delays_days']))
                endpoint_checks.append(dict(index=index, objective=value, parameters=endpoint.tolist(),
                                            max_lens_residual_rad=app['max_lens_residual_rad']))
            except (ArithmeticError, ValueError):
                continue
        if not endpoint_checks:
            raise ArithmeticError('no exact original-projection endpoint admitted')
        best = min(endpoint_checks, key=lambda item: item['objective'])
        p = np.array(best['parameters'])
    final = macro_application(p, compact_mass_Msun=compact_mass_Msun, compact_scale_kpc=compact_scale_kpc,
                              model=model, astrometry=observed)
    astro = (final['relative_astrometry_arcsec']-observed['observed_relative_arcsec'])/observed['sigma_arcsec']
    delta = whiten@(final['delays_days']-delay)
    result = dict(model=model, status='CONDITIONAL_TWO_PLUMMER_DEVELOPMENT_FIT',
                  compact_mass_Msun=float(compact_mass_Msun), compact_scale_kpc=float(compact_scale_kpc),
                  compact_matter_scope='positive finite regular matter surrogate; not an inferred black hole or solved material equation of state',
                  candidate_compact_masses_Msun=list(COMPACT_MASS_CANDIDATES_MSUN),
                  compact_mass_prior='fixed discrete development candidate, no continuous fitted prior; select using HST and old development data only',
                  parameters=dict(zip(PARAMETER_NAMES, p.tolist())), parameter_vector=p.tolist(),
                  parameter_lower_bounds=LOWER.tolist(), parameter_upper_bounds=UPPER.tolist(),
                  predicted_delays_days=final['delays_days'].tolist(),
                  predicted_images_light_centered_arcsec=final['images_light_centered_arcsec'].tolist(),
                  predicted_relative_astrometry_arcsec=final['relative_astrometry_arcsec'].tolist(),
                  observed_relative_astrometry_arcsec=observed['observed_relative_arcsec'].tolist(),
                  chi2_astrometry=float(astro@astro), chi2_delay_initializer=float(delta@delta),
                  chi2_center_prior=float(np.sum((p[7:9]/.003)**2)),
                  standardized_astrometry_residuals=astro.tolist(),
                  optimizer_trace=trace, quadratic_initializer_parameters=quadratic_parameters.tolist(),
                  initial_parameters=initial_parameters.tolist(),
                  skip_quadratic_initializer=bool(skip_quadratic_initializer),
                  exact_training_refinement=exact_trace,
                  optimizer_success=bool(exact_trace[-1]['success'] if exact_trace else trace[-1]['success']),
                  max_lens_residual_rad=final['max_lens_residual_rad'],
                  max_arrival_stationarity_Mpc=final['max_arrival_stationarity_Mpc'],
                  astrometry_covariance=observed['assumptions'],
                  total_local_mass_Msun=final['lens'].base.mass_Msun if model == 'r1' else final['lens'].mass_Msun,
                  physical_admission=False)
    if endpoint_checks:
        result['original_objective_endpoint_selection'] = dict(
            selected_index=best['index'], original_projection_objective=best['objective'],
            checked_endpoints=endpoint_checks,
            interpretation='best start/optimizer endpoint selected by original exact total-source projection and original development likelihood; temporary spline objective never determines endpoint acceptance')
    if delay_log_likelihood is not None:
        result['exact_training_log_likelihood'] = float(delay_log_likelihood(final['delays_days']))
        result['exact_training_objective'] = result['chi2_astrometry']+result['chi2_center_prior']-2*result['exact_training_log_likelihood']
    result['active_parameter_bounds'] = {name: ('lower' if abs(value-lo) < 1e-5 else 'upper')
        for name, value, lo, hi in zip(PARAMETER_NAMES, p, LOWER, UPPER)
        if min(abs(value-lo), abs(value-hi)) < 1e-5}
    full = np.array([old.two_plane_arrival_length(theta*ARCSEC, final['source'], final['lens'], final['U'])[0]
                     for theta in final['images_light_centered_arcsec']])
    result['two_plane_minus_effective_delay_days'] = ((full[1:]-full[0])*MPC_OVER_C_DAYS-final['delays_days']).tolist()
    if model == 'r1':
        screened = final['lens']
        impacts = np.linalg.norm(final['images_light_centered_arcsec']*ARCSEC*old.frozen_distances()['ol']-np.asarray(screened.centre_Mpc), axis=1)
        radii = np.unique(np.r_[impacts, compact_scale_kpc/1000, screened.radius_Mpc])
        result['coefficients'] = asdict(old.frozen_coefficients())
        result['rV_kpc'] = screened.rV_Mpc*1000
        result['original_source_checks'] = screened.checks(impacts)
        result['original_section_source_applications'] = source_section_checks(screened, radii)
        approximate = macro_application(p, compact_mass_Msun=compact_mass_Msun, compact_scale_kpc=compact_scale_kpc,
                                        model=model, frozen_extra=extra, astrometry=observed)
        result['continuation_projection_verification'] = dict(
            max_image_difference_arcsec=float(np.max(abs(approximate['images_light_centered_arcsec']-final['images_light_centered_arcsec']))),
            max_delay_difference_seconds=float(np.max(abs(approximate['delays_days']-final['delays_days']))*86400))
        halo = DirectEffectiveHaloLens(screened)
        differences = []
        for theta in final['images_light_centered_arcsec']:
            actual = lens_application(theta*ARCSEC, final['geometry'], screened, final['source'])
            matched = lens_application(theta*ARCSEC, final['geometry'], halo, final['source'])
            differences.append(dict(potential_Mpc=float(matched['potential']-actual['potential']),
                                    lens_equation_rad=float(np.max(abs(matched['residual']-actual['residual']))),
                                    arrival_days=float(matched['arrival_days']-actual['arrival_days']),
                                    scalar_extra_alpha_rad=float(screened.extra_alpha(float(np.linalg.norm(theta*ARCSEC*old.frozen_distances()['ol']-np.asarray(screened.centre_Mpc))))),
                                    direct_halo_extra_alpha_rad=float(halo.extra_alpha(float(np.linalg.norm(theta*ARCSEC*old.frozen_distances()['ol']-np.asarray(screened.centre_Mpc)))))))
        result['matched_effective_halo_control'] = dict(pointwise_differences=differences,
            interpretation='independent direct positive-effective-mass projection of the SAME nonlinear total-source scalar field; no halo added on top of the scalar')
    if include_census:
        census = image_census(p, compact_mass_Msun=compact_mass_Msun, compact_scale_kpc=compact_scale_kpc, model=model)
        quartet = np.asarray(final['images_light_centered_arcsec'])
        central = []
        for image, magnification in zip(census['images_light_centered_arcsec'], census['absolute_magnifications']):
            if np.min(np.linalg.norm(quartet-np.asarray(image), axis=1)) > 1e-5:
                central.append(dict(position_arcsec=image, absolute_magnification=magnification))
        magnification_a = 1/abs(np.linalg.det(lens_application(quartet[0]*ARCSEC, final['geometry'], final['lens'], final['source'])['jacobian']))
        for image in central:
            image['flux_ratio_to_image_A_achromatic_no_differential_extinction_or_microlensing'] = image['absolute_magnification']/magnification_a
        result['global_image_census'] = census
        result['additional_images'] = central
    return result
