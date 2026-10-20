"""Conditional HE0435 macro forecasts from the existing spherical R1 source.

An explicit remote quadratic optical plane supplies asymmetric Jacobi maps.
It is not an elliptical local matter source or a solved scalar environment.
The historical delay summary initializes the fit. An optional exact historical
light likelihood then supplies the objective, with its hyperparameters frozen.
"""
from __future__ import annotations

from dataclasses import asdict
from functools import lru_cache
import json
from pathlib import Path

import numpy as np
from scipy.integrate import quad
from scipy.interpolate import CubicSpline
from scipy.optimize import least_squares, root, minimize

import action_native_matter_n2 as an
import background as bg
from strong_lensing_comparison import ARCSEC, C_KM_S, LensGeometry, lens_application
from strong_lensing_r1_screened import R1Coefficients
from strong_lensing_robustness import ProfileLens, ProfileR1Lens, DirectEffectiveHaloLens
from strong_lensing_local_family import G_MPC_KMS2_MSUN

ROOT = Path(__file__).resolve().parents[1]
ASTROMETRY = ROOT/'data/strong_lensing_observed_20261020/candidates/HE0435-1223/astrometry_redshifts.json'
PARAMETER_NAMES = ('log10_mass_Msun','log10_plummer_a_kpc','source_x_arcsec',
                   'source_y_arcsec','screen_kappa','screen_gamma1','screen_gamma2',
                   'mass_center_x_arcsec','mass_center_y_arcsec')
LOWER = np.array([10., -.7, -.8, -.8, -.4, -.4, -.4, -.03, -.03])
UPPER = np.array([13., 1.3, .8, .8, .4, .4, .4, .03, .03])


@lru_cache(maxsize=1)
def _frozen_history():
    return an.frozen_background()


@lru_cache(maxsize=1)
def frozen_distances():
    """Evaluate, without re-shooting, the adopted closed-R1 FLRW reference."""
    history = _frozen_history()
    length = C_KM_S/67.4
    curvature = bg.KCURV/length**2
    redshifts = {'l':.4546, 'p':.8, 's':1.693}
    chi = {name:length*quad(lambda z:1/float(history.sol(-np.log1p(z))[2]),
                            0., z, epsabs=2e-12, epsrel=2e-12)[0]
           for name,z in redshifts.items()}
    def distance(a,b):
        separation=chi[b]-(0. if a=='o' else chi[a])
        return np.sin(np.sqrt(curvature)*separation)/np.sqrt(curvature)/(1+redshifts[b])
    result={name:distance(name[0],name[1]) for name in ('ol','op','os','lp','ls','ps')}
    result.update(z_l=redshifts['l'],z_p=redshifts['p'],z_s=redshifts['s'],
                  H0_km_s_Mpc=67.4,curvature_per_Mpc2=curvature)
    identity=((1+redshifts['l'])*result['ol']*result['ps']
              +(1+redshifts['p'])*result['os']*result['lp']
              -(1+redshifts['p'])*result['op']*result['ls'])
    result['distance_sum_relative_residual']=identity/((1+redshifts['p'])*result['op']*result['ls'])
    return result


@lru_cache(maxsize=1)
def frozen_coefficients():
    N=-np.log1p(.4546); y=_frozen_history().sol(N)
    v,H=float(y[1]),float(y[2]); d=H*float(bg.rhs(N,y)[1])
    b=2*v*v/bg.LAMBDA
    C=1+4*(d+2*H*v)/bg.LAMBDA-2*v**4/bg.LAMBDA**2
    if C<=0: raise ValueError('fixed observed-lens reference is not on the regular branch')
    return R1Coefficients(.4546,67.4,H,v,d,bg.LAMBDA,b,C)


def screen_geometry(kappa, gamma1, gamma2, distances=None):
    """Exact Jacobi maps after eliminating a single rear quadratic screen."""
    d=frozen_distances() if distances is None else distances
    gamma=np.array([[kappa+gamma1,gamma2],[gamma2,kappa-gamma1]],float)
    beta_lp=d['os']*d['lp']/(d['op']*d['ls'])
    if min(np.linalg.eigvalsh(np.eye(2)-gamma))<=0 or min(np.linalg.eigvalsh(np.eye(2)-beta_lp*gamma))<=0:
        raise ValueError('remote screen leaves the prescribed pre-conjugate branch')
    U=d['os']/(d['ps']*d['op'])*gamma
    geometry=LensGeometry(d['ol']*np.eye(2),d['os']*np.eye(2)-d['ps']*d['op']*U,
                          d['ls']*np.eye(2)-d['ps']*d['lp']*U,d['z_l'])
    return geometry,U


def two_plane_arrival_length(theta, source_Mpc, lens, U, distances=None):
    """Independent stationary rear-plane application of the full arrival action."""
    d=frozen_distances() if distances is None else distances
    theta=np.asarray(theta,float); beta=np.asarray(source_Mpc,float)/d['os']
    kl=(1+d['z_l'])*d['ol']*d['op']/d['lp']
    kp=(1+d['z_p'])*d['op']*d['os']/d['ps']
    W=(1+d['z_p'])*d['op']**2*np.asarray(U)
    theta_p=np.linalg.solve((kl+kp)*np.eye(2)-W,kl*theta+kp*beta)
    potential=lens.application(d['ol']*theta)[0]
    value=(.5*kl*np.sum((theta-theta_p)**2)+.5*kp*np.sum((theta_p-beta)**2)
           -(1+d['z_l'])*potential-.5*theta_p@W@theta_p)
    gradient=(-kl*(theta-theta_p)+kp*(theta_p-beta)-W@theta_p)
    return float(value),theta_p,gradient


def load_astrometry(path=ASTROMETRY):
    data=json.loads(Path(path).read_text(encoding='utf8'))['astrometry']
    rows=data['data']; positions=np.array([[row['delta_alpha'],row['delta_delta']] for row in rows])
    internal=np.array([[row['internal_sigma_alpha'],row['internal_sigma_delta']] for row in rows[1:]],float)
    sigma=np.sqrt(internal**2+float(data['additional_systematic_sigma_arcsec'])**2).ravel()
    return dict(observed_relative_arcsec=positions[1:].ravel(),sigma_arcsec=sigma,
                observed_images_light_centered_arcsec=positions[:4]-positions[4],
                assumptions='Published internal errors plus required2mas systematic in quadrature; diagonal relative-position covariance approximation; imageA fixes origin.')


def body_from_parameters(parameters):
    p=np.asarray(parameters,float); d=frozen_distances()
    mass=10**p[0]; scale=10**p[1]/1000
    # ProfileLens uses the3D half-mass radius, not the Plummer a parameter.
    half_mass=scale/np.sqrt(2**(2/3)-1)
    return ProfileLens('plummer',mass,half_mass,tuple(d['ol']*ARCSEC*p[7:9]))


class FrozenExtraProjection:
    """One temporary continuation field; gradient/potential/Hessian agree.

    The original R1 deflection is sampled, then its cubic spline is integrated
    exactly. It is refreshed on the fitted body; it is not a physical solution
    at every intermediate optimizer iterate. Final checks use original R1.
    """
    def __init__(self, screened, *, count=180):
        self.screened=screened
        a=screened.radius_Mpc
        self.radii=np.geomspace(a*1e-3,max(20*a,.05),count)
        alpha=np.array([screened.extra_alpha(float(r)) for r in self.radii])
        self.spline=CubicSpline(self.radii,alpha,bc_type='not-a-knot',extrapolate=False)
        self.integral=self.spline.antiderivative()
        self.anchor=float(self.integral(a))

    def application(self, x, center):
        y=np.asarray(x)-np.asarray(center); r=float(np.linalg.norm(y))
        if not self.radii[0]<=r<=self.radii[-1]:
            # A rare optimizer branch far outside the interpolation interval
            # still evaluates the original frozen field rather than extrapolating.
            local=np.asarray(self.screened.centre_Mpc)+y
            return self.screened.extra_application(local)
        alpha=float(self.spline(r)); ap=float(self.spline(r,1)); unit=y/r
        potential=float(self.integral(r))-self.anchor
        hessian=alpha/r*np.eye(2)+(ap-alpha/r)*np.outer(unit,unit)
        return potential,alpha*unit,hessian


class CompositeContinuationLens:
    def __init__(self, body, extra):
        self.base=body; self.extra=extra; self.b=body.b; self.centre_Mpc=body.centre_Mpc
    def application(self,x):
        p,a,h=self.base.application(x); dp,da,dh=self.extra.application(x,self.centre_Mpc)
        return p+dp,a+da,h+dh


def macro_application(parameters, *, model='conventional', frozen_extra=None,
                      astrometry=None, require_success=True):
    p=np.asarray(parameters,float); d=frozen_distances()
    observed=load_astrometry() if astrometry is None else astrometry
    body=body_from_parameters(p)
    if model=='conventional':
        lens=body
    elif model=='r1':
        lens=(ProfileR1Lens(body,frozen_coefficients()) if frozen_extra is None
              else CompositeContinuationLens(body,frozen_extra))
    else:
        raise ValueError('model must be conventional or r1')
    geometry,U=screen_geometry(*p[4:7],distances=d)
    source=d['os']*ARCSEC*p[2:4]
    images=[]; applications=[]; failures=[]
    for seed in observed['observed_images_light_centered_arcsec']:
        def fun(q):
            return lens_application(q*ARCSEC,geometry,lens,source)['residual']/ARCSEC
        def jac(q):
            return lens_application(q*ARCSEC,geometry,lens,source)['jacobian']
        solved=root(fun,seed,jac=jac,tol=2e-10)
        defect=float(np.max(abs(fun(solved.x))))
        if defect>2e-7 or not np.isfinite(solved.x).all():
            failures.append(defect)
        images.append(solved.x)
        applications.append(lens_application(solved.x*ARCSEC,geometry,lens,source))
    images=np.asarray(images)
    separation=min(np.linalg.norm(images[i]-images[j]) for i in range(4) for j in range(i))
    if separation<.1: failures.append(float(.1-separation))
    if failures and require_success:
        raise ArithmeticError('one or more requested image branches failed or merged')
    arrivals=np.array([app['arrival_days'] for app in applications])
    relative=np.r_[(images[1:]-images[0]).ravel(),-images[0]]
    return dict(images_light_centered_arcsec=images,relative_astrometry_arcsec=relative,
                delays_days=arrivals[1:]-arrivals[0],arrival_days=arrivals,
                max_lens_residual_rad=float(max(np.max(abs(app['residual'])) for app in applications)),
                max_arrival_stationarity_Mpc=float(max(np.max(abs(app['stationarity'])) for app in applications)),
                minimum_image_separation_arcsec=float(separation),branch_failures=failures,
                geometry=geometry,source=source,lens=lens,U=U)


def _initializer(delay, whiten, observed, initial=None):
    """Training-only lens-equation initializer before solving image branches."""
    d=frozen_distances(); fixed=observed['observed_images_light_centered_arcsec']*ARCSEC
    seeds=([np.asarray(initial,float)] if initial is not None else
           [np.array([11.8,.75,0.,0.,0.,g1,g2,0.,0.]) for g1,g2 in
            ((.05,.05),(-.05,-.05),(.1,-.05),(-.1,.05))])
    def residual(p):
        body=body_from_parameters(p); geometry,U=screen_geometry(*p[4:7],distances=d)
        source=d['os']*ARCSEC*p[2:4]
        apps=[lens_application(theta,geometry,body,source) for theta in fixed]
        lensres=np.concatenate([app['residual']/ARCSEC/.003 for app in apps])
        arrivals=np.array([app['arrival_days'] for app in apps])
        return np.r_[lensres,whiten@(arrivals[1:]-arrivals[0]-delay),p[7:9]/.003]
    trials=[least_squares(residual,np.clip(seed,LOWER+1e-9,UPPER-1e-9),bounds=(LOWER,UPPER),
                          max_nfev=300,ftol=1e-10,xtol=1e-10,gtol=1e-8) for seed in seeds]
    return min(trials,key=lambda fit:np.sum(fit.fun**2)).x


def fit_macro(train_delays_days, train_delay_covariance_days2, *, model='conventional',
              initial=None, astrometry=None, max_nfev=500, refresh_steps=2,
              delay_log_likelihood=None, exact_maxiter=60):
    """Fit historical astrometry and optional exact training-light objective.

    The covariance of a local delay curvature is an initialization weight only
    when ``delay_log_likelihood`` is supplied; it is not a confidence region.
    No withheld magnitudes or nuisance refitting are consumed by this API.
    """
    delay=np.asarray(train_delays_days,float)
    covariance=np.asarray(train_delay_covariance_days2,float)
    if delay.shape!=(3,) or covariance.shape!=(3,3):
        raise ValueError('three imageB/C/D-minusA training delays and3x3 covariance required')
    if not np.allclose(covariance,covariance.T,rtol=1e-10,atol=1e-12):
        raise ValueError('training delay covariance must be symmetric')
    whiten=np.linalg.solve(np.linalg.cholesky(covariance),np.eye(3))
    observed=load_astrometry() if astrometry is None else astrometry
    p=_initializer(delay,whiten,observed,initial)
    trace=[]; lastfit=None; extra=None
    loops=1 if model=='conventional' else refresh_steps+1
    for iteration in range(loops):
        if model=='r1':
            extra=FrozenExtraProjection(ProfileR1Lens(body_from_parameters(p),frozen_coefficients()))
        def residual(parameters):
            try:
                app=macro_application(parameters,model=model,frozen_extra=extra,astrometry=observed)
                return np.r_[(app['relative_astrometry_arcsec']-observed['observed_relative_arcsec'])/observed['sigma_arcsec'],
                             whiten@(app['delays_days']-delay),parameters[7:9]/.003]
            except (ArithmeticError,ValueError):
                # A deterministic rejected branch gives finite optimizer output.
                # It is never admitted as an evaluated physical fit.
                return np.r_[np.full(11,1e6),parameters[7:9]/.003]
        lastfit=least_squares(residual,p,bounds=(LOWER,UPPER),max_nfev=max_nfev,
                              ftol=2e-11,xtol=2e-11,gtol=2e-7,
                              diff_step=2e-5,x_scale='jac')
        p=lastfit.x
        trace.append(dict(iteration=iteration,chi2=float(np.sum(lastfit.fun**2)),
                          optimizer_success=bool(lastfit.success),function_evaluations=int(lastfit.nfev),
                          optimizer_status=int(lastfit.status),message=str(lastfit.message),
                          parameters=p.tolist()))
    quadratic_parameters=p.copy()
    exact_trace=[]
    if delay_log_likelihood is not None:
        # The local delay curvature is only a numerical initializer. The final
        # refinement consumes the actual frozen-hyperparameter training likelihood.
        for iteration in range(1 if model=='conventional' else 2):
            if model=='r1':
                extra=FrozenExtraProjection(ProfileR1Lens(body_from_parameters(p),frozen_coefficients()))
            def objective(parameters):
                try:
                    app=macro_application(parameters,model=model,frozen_extra=extra,astrometry=observed)
                    astro=(app['relative_astrometry_arcsec']-observed['observed_relative_arcsec'])/observed['sigma_arcsec']
                    return float(astro@astro+np.sum((parameters[7:9]/.003)**2)-2*delay_log_likelihood(app['delays_days']))
                except (ArithmeticError,ValueError):
                    return 1e15
            optimized=minimize(objective,p,method='L-BFGS-B',jac='3-point',bounds=list(zip(LOWER,UPPER)),
                               options=dict(maxiter=exact_maxiter,ftol=2e-12,gtol=1e-6,maxls=25,
                                            finite_diff_rel_step=np.array([2e-6,2e-5,1e-4,1e-4,1e-5,1e-5,1e-5,1e-4,1e-4])))
            p=optimized.x
            exact_trace.append(dict(iteration=iteration,objective=float(optimized.fun),
                                    success=bool(optimized.success),message=str(optimized.message),
                                    evaluations=int(optimized.nfev),parameters=p.tolist()))
    final=macro_application(p,model=model,astrometry=observed)
    fitres=np.r_[(final['relative_astrometry_arcsec']-observed['observed_relative_arcsec'])/observed['sigma_arcsec'],
                  whiten@(final['delays_days']-delay),p[7:9]/.003]
    result=dict(model=model,status='CONDITIONAL_OBSERVED_MACRO_FIT',parameters=dict(zip(PARAMETER_NAMES,p.tolist())),
                parameter_vector=p.tolist(),distances=frozen_distances(),
                predicted_relative_astrometry_arcsec=final['relative_astrometry_arcsec'].tolist(),
                observed_relative_astrometry_arcsec=observed['observed_relative_arcsec'].tolist(),
                predicted_images_light_centered_arcsec=final['images_light_centered_arcsec'].tolist(),
                predicted_delays_days=final['delays_days'].tolist(),
                training_delays_days=delay.tolist(),training_delay_covariance_days2=covariance.tolist(),
                standardized_astrometry_residuals=fitres[:8].tolist(),
                whitened_training_delay_residuals=fitres[8:11].tolist(),
                mass_center_prior_residuals=fitres[11:].tolist(),
                chi2_astrometry=float(fitres[:8]@fitres[:8]),chi2_delays=float(fitres[8:11]@fitres[8:11]),
                chi2_center_prior=float(fitres[11:]@fitres[11:]),chi2_total=float(fitres@fitres),
                optimizer_trace=trace,optimizer_success=bool(lastfit.success),
                quadratic_initializer_parameters=quadratic_parameters.tolist(),
                exact_training_refinement=exact_trace,
                max_lens_residual_rad=final['max_lens_residual_rad'],
                max_arrival_stationarity_Mpc=final['max_arrival_stationarity_Mpc'],
                minimum_image_separation_arcsec=final['minimum_image_separation_arcsec'],
                numerical_branch_scope='four requested noncritical branches; not a global faint central-image census',
                astrometry_covariance=observed['assumptions'],
                physical_scope='spherical prescribed local source plus fixed-location total rear optical screen; no solved global scalar environment',
                physical_admission=False)
    result['active_parameter_bounds']={name:('lower' if abs(value-lo)<1e-5 else 'upper')
        for name,value,lo,hi in zip(PARAMETER_NAMES,p,LOWER,UPPER)
        if min(abs(value-lo),abs(value-hi))<1e-5}
    result['optimizer_success']=(bool(exact_trace[-1]['success']) if exact_trace else bool(lastfit.success))
    if delay_log_likelihood is not None:
        result['exact_training_log_likelihood']=float(delay_log_likelihood(final['delays_days']))
        result['exact_training_objective']=float(fitres[:8]@fitres[:8]+fitres[11:]@fitres[11:]-2*result['exact_training_log_likelihood'])
    # Independent full two-plane arrival action; source-only constants cancel.
    full=np.array([two_plane_arrival_length(theta*ARCSEC,final['source'],final['lens'],final['U'])[0]
                   for theta in final['images_light_centered_arcsec']])
    from strong_lensing_comparison import MPC_OVER_C_DAYS
    result['two_plane_minus_effective_delay_days']=((full[1:]-full[0])*MPC_OVER_C_DAYS-final['delays_days']).tolist()
    if model=='r1':
        screened=final['lens']; impacts=np.linalg.norm(final['images_light_centered_arcsec']*ARCSEC*frozen_distances()['ol']-np.asarray(screened.centre_Mpc),axis=1)
        result['coefficients']=asdict(frozen_coefficients())
        result['rV_kpc']=screened.rV_Mpc*1000
        result['original_source_checks']=screened.checks(impacts)
        c=screened.coefficients; length=c.hubble_length_Mpc
        sections=[]
        for r in impacts:
            mu,_,_=screened.mass_and_slope(float(r)); x=r/length
            gradient=screened.scalar_gradient_natural(float(r))
            mass_natural=screened.mass_length_Mpc/length*mu
            density=screened.base.density_Msun_per_Mpc3(float(r))
            density_natural=8*np.pi*G_MPC_KMS2_MSUN/C_KM_S**2*density*length**2
            scalar_lhs=c.C*x*x*gradient-4/c.lambda_over_MplH02*x*gradient**2
            phi_prime_natural=length*screened.mass_length_Mpc*mu/r**2*(1+screened.extra_force(float(r)))
            gauss_lhs=2*x*x*phi_prime_natural+c.b*x*x*gradient
            sections.append(dict(radius_kpc=float(r*1000),density_over_Mpl2H02=float(density_natural),
                                 consumed_reduced_source_J_eff=float(-.5*c.b*density_natural),
                                 scalar_gradient_over_MplH0=float(gradient),
                                 scalar_integrated_lhs=float(scalar_lhs),
                                 scalar_integrated_source=float(-c.b*mass_natural),
                                 metric_Gauss_lhs=float(gauss_lhs),metric_Gauss_source=float(2*mass_natural),
                                 scalar_relative_residual=float(abs(scalar_lhs+c.b*mass_natural)/(c.b*mass_natural)),
                                 metric_relative_residual=float(abs(gauss_lhs-2*mass_natural)/(2*mass_natural))))
        result['original_section_source_applications']=sections
        approximate=macro_application(p,model=model,frozen_extra=extra,astrometry=observed)
        result['continuation_projection_verification']=dict(
            max_image_difference_arcsec=float(np.max(abs(approximate['images_light_centered_arcsec']-final['images_light_centered_arcsec']))),
            max_delay_difference_seconds=float(np.max(abs(approximate['delays_days']-final['delays_days']))*86400),
            interpretation='temporary frozen radial interpolation compared against original projections after final body refresh')
        halo=DirectEffectiveHaloLens(screened)
        differences=[]
        for theta in final['images_light_centered_arcsec']:
            actual=lens_application(theta*ARCSEC,final['geometry'],screened,final['source'])
            matched=lens_application(theta*ARCSEC,final['geometry'],halo,final['source'])
            differences.append(dict(potential_Mpc=float(matched['potential']-actual['potential']),
                                    lens_equation_rad=float(np.max(abs(matched['residual']-actual['residual']))),
                                    arrival_days=float(matched['arrival_days']-actual['arrival_days'])))
        result['matched_effective_halo_control']=dict(pointwise_differences=differences,
            interpretation='independent direct positive-mass projection of same retained field; not an additional halo added on top of scalar')
    return result
