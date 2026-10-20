"""Bounded lens-profile controls for the same state-conditioned light transfer.

The alternative is the standard Plummer sphere, with finite total mass and a
regular center. Both profiles are matched by total mass and 3D half-mass
radius, not by a fitted deflection. The R1 action/state remains unchanged.
"""
from __future__ import annotations

from dataclasses import dataclass
from functools import lru_cache
import math

import numpy as np
from scipy.integrate import quad
from scipy.optimize import brentq

from strong_lensing_comparison import (
    ARCSEC, C_KM_S, complete_map_response, nuisance_projected_information,
    response_vector, solve_images,
)
from strong_lensing_local_family import (
    G_MPC_KMS2_MSUN, TolmanVIILens, local_nuisance_columns,
)
from strong_lensing_r1_screened import R1Coefficients, R1ScreenedLens, solve_axis_images

TOLMAN_HALF_X = brentq(lambda x: .5*(5*x**3-3*x**5)-.5, .5, .8, xtol=1e-15)
BASE_HALF_MASS_MPC = .003*TOLMAN_HALF_X


@dataclass(frozen=True)
class ProfileLens:
    profile: str
    mass_Msun: float = 1e11
    half_mass_Mpc: float = BASE_HALF_MASS_MPC
    centre_Mpc: tuple = (0., 0.)

    def __post_init__(self):
        if self.profile not in ('tolman_vii', 'plummer'):
            raise ValueError('unsupported declared profile')
        if min(self.mass_Msun, self.half_mass_Mpc) <= 0:
            raise ValueError('positive total mass and half-mass radius required')

    @property
    def radius_Mpc(self):
        # Material support for Tolman VII; Plummer scale a has no hard edge.
        if self.profile == 'tolman_vii':
            return self.half_mass_Mpc/TOLMAN_HALF_X
        return self.half_mass_Mpc*math.sqrt(2**(2/3)-1)

    @property
    def mass_length_Mpc(self):
        return G_MPC_KMS2_MSUN*self.mass_Msun/C_KM_S**2

    @property
    def b(self):
        return 4*self.mass_length_Mpc/self.radius_Mpc

    def mass_and_slope(self, r):
        a = self.radius_Mpc
        x = r/a
        if self.profile == 'plummer':
            return x**3/(1+x*x)**1.5, 3/(1+x*x), 1/(r*r+a*a)**1.5
        if x >= 1:
            return 1., 0., 1/r**3
        return x**3*(2.5-1.5*x*x), (7.5-7.5*x*x)/(2.5-1.5*x*x), (2.5-1.5*x*x)/a**3

    def density_Msun_per_Mpc3(self, r):
        a = self.radius_Mpc
        if self.profile == 'plummer':
            return 3*self.mass_Msun/(4*np.pi*a**3)*(1+(r/a)**2)**-2.5
        return 15*self.mass_Msun/(8*np.pi*a**3)*max(1-(r/a)**2, 0.)

    def application(self, x):
        if self.profile == 'tolman_vii':
            return TolmanVIILens(self.mass_Msun, self.radius_Mpc, self.centre_Mpc).application(x)
        displacement = np.asarray(x, float)-np.asarray(self.centre_Mpc)
        r2 = float(displacement@displacement)
        a2 = self.radius_Mpc**2
        m = self.mass_length_Mpc
        potential = 2*m*math.log((r2+a2)/(2*a2))
        alpha = 4*m*displacement/(r2+a2)
        hessian = 4*m*(np.eye(2)/(r2+a2)-2*np.outer(displacement, displacement)/(r2+a2)**2)
        return potential, alpha, hessian

    def intrinsic(self, coordinate):
        if coordinate not in ('log_mass', 'log_radius'):
            raise ValueError('only declared mass and homologous size variations are supported')
        def partial(x):
            p, alpha, hessian = self.application(x)
            if coordinate == 'log_mass':
                return p, alpha
            v = np.asarray(x)-np.asarray(self.centre_Mpc)
            return -float(v@alpha), -alpha-hessian@v
        return partial


class ProfileR1Lens(R1ScreenedLens):
    """Same original scalar branch with the declared profile's enclosed mass."""

    def mass_and_slope(self, r):
        return self.base.mass_and_slope(r)

    def effective_density_Msun_per_Mpc3(self, r):
        if r == 0:
            return self.base.density_Msun_per_Mpc3(0.)*self.extra_force(0.)
        mu, slope, mr3 = self.mass_and_slope(r)
        S = math.sqrt(1+self.rV3_Mpc3*mr3)
        log_slope = .5*((1+1/S)*slope+3*(1-1/S))
        return self.mass_Msun*mu*self.extra_force(r)*log_slope/(4*np.pi*r**3)


def projected_coefficient_sensitivities(lens, impact):
    """Fixed-M,lambda partials of the extra deflection, not selected states.

    For S=sqrt(1+rV³ mu/r³), dln(delta Phi')/dln b_T=3/2+1/(2S)
    and dln(delta Phi')/dln C=-1/S. Projection takes positive weighted
    averages. An actual envelope-state derivative still needs its b_T,C
    Jacobian; varying these independently here does not select a new state.
    """
    denominator = lens._integral(impact)
    points = [math.acos(impact/r) for r in (lens.radius_Mpc, lens.rV_Mpc) if r > impact]
    def weighted(theta, coordinate):
        r = impact/math.cos(theta)
        S = math.sqrt(1+lens.rV3_Mpc3*lens.mass_and_slope(r)[2])
        factor = 1.5+.5/S if coordinate == 'ln_b_T' else -1/S
        return lens._theta_integrand(theta, impact, False)*factor
    values = {coordinate:quad(lambda t:weighted(t, coordinate), 0., math.pi/2,
                              points=points, epsabs=1e-18, epsrel=2e-11, limit=180)[0]/denominator
              for coordinate in ('ln_b_T', 'ln_C')}
    base_alpha = float(lens.base.application([impact, 0.])[1][0])
    extra_alpha = lens.extra_alpha(float(impact))
    return dict(impact_kpc=impact*1000, extra_alpha_rad=extra_alpha,
                d_log_extra_alpha=values,
                d_total_alpha_divided_by_total_alpha={key:value*extra_alpha/(base_alpha+extra_alpha) for key,value in values.items()},
                identity_residual_2_db_plus_dC_minus_3=2*values['ln_b_T']+values['ln_C']-3,
                scope='local coefficient partials at fixed M,lambda,profile and physical impact; actual state dependence requires the state-to-coefficient Jacobian')


class DirectEffectiveHaloLens:
    """Independent z/b projection of the positive effective mass representation.

    It does not reuse the R1 angular projection or its projected potential.
    The source mass mu*F is the analytic equivalence map being tested.
    """
    def __init__(self, screened):
        self.screened = screened
        self.base = screened.base
        self.radius_Mpc = screened.radius_Mpc
        self.mass_length_Mpc = screened.mass_length_Mpc
        self.b = screened.b
        self.centre_Mpc = screened.centre_Mpc

    def _integral(self, impact, derivative=False):
        scales = (self.radius_Mpc, self.screened.rV_Mpc)
        cuts = sorted(set([0., 1.] + [math.sqrt((r/impact)**2-1) for r in scales if r > impact]))
        cuts.append(np.inf)
        def integrand(t):
            r = impact*math.sqrt(1+t*t)
            mu, slope, mr3 = self.base.mass_and_slope(r)
            S = math.sqrt(1+self.screened.rV3_Mpc3*mr3)
            F = 2*self.screened.coefficients.unscreened_extra/(1+S)
            value = mu*F/(1+t*t)**1.5
            if derivative:
                value *= .5*((1+1/S)*slope+3*(1-1/S))-1
            return value
        return sum(quad(integrand, lo, hi, epsabs=1e-18, epsrel=2e-11, limit=180)[0]
                   for lo, hi in zip(cuts[:-1], cuts[1:]))

    @lru_cache(maxsize=4096)
    def extra_alpha(self, impact):
        if impact == 0:
            return 0.
        return 4*self.mass_length_Mpc/impact*self._integral(impact)

    def extra_alpha_prime(self, impact):
        return 4*self.mass_length_Mpc/impact**2*self._integral(impact, True)

    @lru_cache(maxsize=4096)
    def extra_potential(self, impact):
        return quad(self.extra_alpha, self.radius_Mpc, float(impact),
                    epsabs=1e-24, epsrel=2e-11, limit=180)[0]

    def application(self, x):
        p, alpha, hessian = self.base.application(x)
        v = np.asarray(x)-np.asarray(self.centre_Mpc)
        r = float(np.linalg.norm(v))
        if r == 0:
            raise ValueError('this comparison evaluates noncentral image impacts')
        unit = v/r
        extra = self.extra_alpha(r)
        extra_hessian = extra/r*np.eye(2)+(self.extra_alpha_prime(r)-extra/r)*np.outer(unit, unit)
        return p+self.extra_potential(r), alpha+extra*unit, hessian+extra_hessian


def profile_nuisance_check(geometry, base, source, conventional, screened):
    zero = np.zeros((3, 2, 2))
    response = complete_map_response(conventional['images_rad'], geometry, base, source, zero,
                                     intrinsic_encounter=screened.intrinsic_extra)
    B = response_vector(response)[:, None]
    N, names = local_nuisance_columns(geometry, base, source, conventional)
    center_columns = []
    for direction in np.eye(2):
        shift = geometry.D_ol@(direction*ARCSEC)
        def center_partial(x, shift=shift):
            _, a, h = base.application(x)
            return -float(a@shift), -h@shift
        center_columns.append(response_vector(complete_map_response(
            conventional['images_rad'], geometry, base, source, zero, intrinsic_encounter=center_partial)))
    scales = np.r_[np.full(conventional['images_rad'].size, .001),
                   np.full(len(conventional['relative_delays_days']), .1)]
    covariance = np.diag(scales**2)
    return dict(nuisance_names=names+['center_x_per_arcsec', 'center_y_per_arcsec'],
                fixed_center=nuisance_projected_information(B, N, covariance),
                free_center=nuisance_projected_information(B, np.column_stack([N, *center_columns]), covariance),
                scope='one axial source, this profile family, linearized response; illustrative .001 arcsec/.1 day weighting, not measured precision or detection significance')


def evaluate_case(geometry, coefficients, source, profile, mass_factor=1., size_factor=1., source_factor=1., *, detailed=False):
    base = ProfileLens(profile, 1e11*mass_factor, BASE_HALF_MASS_MPC*size_factor)
    screened = ProfileR1Lens(base, coefficients)
    physical_source = np.asarray(source)*source_factor
    conventional = solve_images(geometry, base, physical_source)
    response = solve_axis_images(geometry, screened, physical_source, conventional['images_rad'])
    impacts = np.linalg.norm(conventional['images_rad']@geometry.D_ol.T, axis=1)
    x = impacts/coefficients.hubble_length_Mpc
    p = np.array([screened.scalar_gradient_natural(float(r)) for r in impacts])
    mn = base.mass_length_Mpc/coefficients.hubble_length_Mpc
    mu = np.array([base.mass_and_slope(float(r))[0] for r in impacts])
    flux = coefficients.C*x*x*p-4/coefficients.lambda_over_MplH02*x*p*p
    flux_relative = np.max(np.abs(flux+coefficients.b*mn*mu)/(coefficients.b*mn*mu))
    case = dict(profile=profile, mass_factor=mass_factor, half_mass_radius_factor=size_factor,
                physical_source_factor=source_factor, mass_Msun=base.mass_Msun,
                half_mass_radius_kpc=1000*base.half_mass_Mpc, profile_scale_kpc=1000*base.radius_Mpc,
                hard_material_edge=profile == 'tolman_vii',
                physical_source_Mpc=physical_source, screening_radius_kpc=1000*screened.rV_Mpc,
                conventional=conventional, screened=response,
                reference_relative_offsets=dict(images_microarcsec=(response['images_arcsec']-conventional['images_arcsec'])*1e6,
                                                delays_seconds=(response['relative_delays_days']-conventional['relative_delays_days'])*86400),
                original_scalar_flux_relative_residual=float(flux_relative),
                maximum_image_residual_rad=max(conventional['max_lens_residual_rad'], response['max_lens_residual_rad']))
    if detailed:
        case['projected_coefficient_sensitivities'] = [projected_coefficient_sensitivities(screened, float(r)) for r in impacts]
        halo = DirectEffectiveHaloLens(screened)
        matched = solve_axis_images(geometry, halo, physical_source, response['images_rad'])
        radial = np.geomspace(base.radius_Mpc*1e-5, 20*screened.rV_Mpc, 121)
        density = np.array([screened.effective_density_Msun_per_Mpc3(float(r)) for r in radial])
        mass_checks = []
        for radius in (base.half_mass_Mpc, base.radius_Mpc*3, screened.rV_Mpc):
            enclosed = quad(lambda t: 4*np.pi*t*t*screened.effective_density_Msun_per_Mpc3(t),
                            0., radius, points=[base.radius_Mpc] if radius > base.radius_Mpc else None,
                            epsabs=1e-6, epsrel=2e-10, limit=200)[0]
            target = base.mass_Msun*base.mass_and_slope(radius)[0]*screened.extra_force(radius)
            mass_checks.append(dict(radius_kpc=radius*1000, integrated_effective_mass_Msun=enclosed,
                                    required_effective_mass_Msun=target, relative_residual=abs(enclosed-target)/target))
        case['effective_mass_equivalence'] = dict(
            same_profile_and_state=True, nonnegative_density=bool(np.all(density >= 0)),
            minimum_sampled_density_Msun_per_Mpc3=float(density.min()), enclosed_mass_checks=mass_checks,
            independent_direct_z_vs_R1_alpha_relative_error=float(max(abs(halo.extra_alpha(float(r))-screened.extra_alpha(float(r)))/screened.extra_alpha(float(r)) for r in impacts)),
            recomputed_effective_mass_images_arcsec=matched['images_arcsec'],
            recomputed_effective_mass_relative_delays_days=matched['relative_delays_days'],
            maximum_image_difference_arcsec=float(np.max(abs(matched['images_arcsec']-response['images_arcsec']))),
            maximum_delay_difference_seconds=float(np.max(abs(matched['relative_delays_days']-response['relative_delays_days']))*86400))
        case['linearized_nuisance_check'] = profile_nuisance_check(geometry, base, physical_source, conventional, screened)
        upper = 1. if profile == 'tolman_vii' else np.inf
        def pressure_integrand(t, extra=False):
            if t == 0:
                return 0.
            shape = 15/(8*np.pi)*max(1-t*t, 0.) if profile == 'tolman_vii' else 3/(4*np.pi)*(1+t*t)**-2.5
            r = base.radius_Mpc*t
            value = shape*base.mass_and_slope(r)[0]/(t*t)
            return value*(screened.extra_force(r) if extra else 1.)
        pressure_N = quad(pressure_integrand, 0., upper, epsabs=1e-13, epsrel=2e-11)[0]
        pressure_extra = quad(lambda t: pressure_integrand(t, True), 0., upper, epsabs=1e-18, epsrel=2e-11)[0]
        case['leading_matter_support'] = dict(
            equation='p(r)=integral_r^outer rho(s) Phi_prime(s) ds, outer=R for Tolman VII and infinity for Plummer',
            central_newtonian_pbar=pressure_N, central_extra_pbar=pressure_extra,
            central_pressure_fractional_adjustment=pressure_extra/pressure_N,
            pbar_units='p/c²=(M m_g/a^4) pbar, a=declared profile scale',
            scope='positive leading isotropic support for a prescribed density; not a common fixed material EOS, exact modified TOV, or distribution-function stability proof')
    return case
