"""Positive nonlinear rear-plane geometry for the HE0435 development revision.

The local two-Plummer body and its original total-source radial R1 solution are
unchanged. A positive regular Plummer perturber and a remaining quadratic
optical potential share one rear plane. They are separate terms of that plane's
potential, not two insertions of an already-complete environmental field.
The perturber center follows its observed ray through the foreground galaxy.
No unused observation is queried or read here.
"""
from __future__ import annotations

from functools import lru_cache
from dataclasses import dataclass

import numpy as np
from scipy.integrate import quad
from scipy.optimize import root

import strong_lensing_observed_macro as old
import strong_lensing_macro_revision as local
from strong_lensing_comparison import ARCSEC, C_KM_S, MPC_OVER_C_DAYS
from strong_lensing_local_family import G_MPC_KMS2_MSUN
from strong_lensing_robustness import ProfileR1Lens


def exact_gaia_relative_geometry(ra_deg,dec_deg,source_covariance_mas2):
    """Exact gnomonic relative geometry and full common-A tangent covariance.

    Input rows are ordered A/B/C/D. Each 2x2 catalogue covariance is in its
    own local alpha-star/declination basis, in mas squared; ra_error already
    includes cos(delta). Different catalogue sources are conditionally treated
    as independent because their cross-covariances are unpublished. The exact
    6x8 Jacobian also differentiates the uncertain tangent origin at A.
    """
    ra=np.asarray(ra_deg,float);dec=np.asarray(dec_deg,float)
    source_cov=np.asarray(source_covariance_mas2,float)
    if (ra.shape!=(4,) or dec.shape!=(4,) or source_cov.shape!=(4,2,2)
            or not np.isfinite(ra).all() or not np.isfinite(dec).all()
            or not np.isfinite(source_cov).all() or np.max(abs(dec))>=90):
        raise ValueError('Four finite ordered nonpolar sky positions and 2x2 covariances required')
    for cov in source_cov:
        if not np.allclose(cov,cov.T,rtol=0.,atol=1e-12):
            raise ValueError('Symmetric tangent covariances required')
        np.linalg.cholesky(cov)
    a=np.deg2rad(ra);d=np.deg2rad(dec)
    delta=(a[1:]-a[0]+np.pi)%(2*np.pi)-np.pi
    c,s=np.cos(d[1:]),np.sin(d[1:])
    c0,s0=np.cos(d[0]),np.sin(d[0])
    cd,sd=np.cos(delta),np.sin(delta)
    u=c*sd
    v=c0*s-s0*c*cd
    w=s0*s+c0*c*cd
    if np.min(w)<=0:
        raise ValueError('Gaia source outside the reference tangent hemisphere')
    relative=np.c_[u/w,v/w]/ARCSEC
    jac=np.zeros((6,8))
    for index in range(3):
        # Derivatives of (u,v,w) wrt alpha_i and delta_i in radians.
        na=np.array([c[index]*cd[index],s0*c[index]*sd[index]])
        wa=-c0*c[index]*sd[index]
        nd=np.array([-s[index]*sd[index],c0*c[index]+s0*s[index]*cd[index]])
        wd=s0*c[index]-c0*s[index]*cd[index]
        numer=np.array([u[index],v[index]])
        da=(na*w[index]-numer*wa)/w[index]**2
        dd=(nd*w[index]-numer*wd)/w[index]**2
        d0=(np.array([0.,-w[index]])*w[index]-numer*v[index])/w[index]**2
        row=slice(2*index,2*index+2)
        # One mas of input tangent displacement gives .001 arcsec times
        # the radian/radian derivative, including each row's own cosine.
        jac[row,0]=-.001*da/c0
        jac[row,1]=.001*d0
        jac[row,2*(index+1)]=.001*da/c[index]
        jac[row,2*(index+1)+1]=.001*dd
    full=np.zeros((8,8))
    for index,cov in enumerate(source_cov):
        full[2*index:2*index+2,2*index:2*index+2]=cov
    return dict(relative_B_C_D_to_A_east_north_arcsec=relative.ravel(),
                covariance_arcsec2=jac@full@jac.T,
                jacobian_arcsec_per_tangent_mas=jac,
                scope='Exact gnomonic tangent at uncertain A; independent catalogue-source approximation with each own alpha-star basis; all shared-origin covariance retained')


@lru_cache(maxsize=16)
def distances(rear_redshift=.7821):
    """Use the already-adopted background; evaluate maps at the measured plane."""
    zp = float(rear_redshift)
    reference = old.frozen_distances()
    zl, zs = reference['z_l'], reference['z_s']
    if not zl < zp < zs:
        raise ValueError('Rear-plane redshift must lie between galaxy and source')
    if zp == reference['z_p']:
        return dict(reference)
    history = old._frozen_history()
    length = C_KM_S/reference['H0_km_s_Mpc']
    curvature = reference['curvature_per_Mpc2']
    redshifts = {'l':zl, 'p':zp, 's':zs}
    chi = {key:length*quad(lambda z:1/float(history.sol(-np.log1p(z))[2]),
                           0., z, epsabs=2e-12, epsrel=2e-12)[0]
           for key,z in redshifts.items()}
    def distance(a,b):
        separation = chi[b]-(0. if a=='o' else chi[a])
        return np.sin(np.sqrt(curvature)*separation)/np.sqrt(curvature)/(1+redshifts[b])
    d = {key:distance(key[0],key[1]) for key in ('ol','op','os','lp','ls','ps')}
    d.update(z_l=zl,z_p=zp,z_s=zs,H0_km_s_Mpc=reference['H0_km_s_Mpc'],
             curvature_per_Mpc2=curvature)
    identity = ((1+zl)*d['ol']*d['ps']+(1+zp)*d['os']*d['lp']
                -(1+zp)*d['op']*d['ls'])
    d['distance_sum_relative_residual'] = identity/((1+zp)*d['op']*d['ls'])
    return d


@dataclass(frozen=True)
class PositiveRearPlummer:
    mass_Msun: float
    scale_Mpc: float
    centre_Mpc: tuple

    def __post_init__(self):
        if (not np.isfinite([self.mass_Msun,self.scale_Mpc,*self.centre_Mpc]).all()
                or self.mass_Msun < 0 or self.scale_Mpc <= 0 or len(self.centre_Mpc)!=2):
            raise ValueError('Finite nonnegative mass, positive scale and two-vector center required')

    @property
    def mass_length_Mpc(self):
        return G_MPC_KMS2_MSUN*self.mass_Msun/C_KM_S**2

    def application(self, position_Mpc):
        x = np.asarray(position_Mpc,float)-np.asarray(self.centre_Mpc)
        denominator = x@x+self.scale_Mpc**2
        m = self.mass_length_Mpc
        potential = 2*m*np.log(denominator/(2*self.scale_Mpc**2))
        alpha = 4*m*x/denominator
        hessian = 4*m*(np.eye(2)/denominator-2*np.outer(x,x)/denominator**2)
        return float(potential),alpha,hessian

    def surface_density_Msun_per_Mpc2(self, position_Mpc):
        x = np.asarray(position_Mpc,float)-np.asarray(self.centre_Mpc)
        return self.mass_Msun*self.scale_Mpc**2/(np.pi*(x@x+self.scale_Mpc**2)**2)


def components(parameters, perturber_observed_arcsec, *, perturber_mass_Msun,
               perturber_scale_kpc=.2, rear_redshift=.7821,
               compact_mass_Msun=1e8, compact_scale_kpc=.03,
               model='conventional', frozen_extra=None):
    p = np.asarray(parameters,float)
    observed = np.asarray(perturber_observed_arcsec,float)
    if p.shape!=(9,) or observed.shape!=(2,) or not np.isfinite(p).all() or not np.isfinite(observed).all():
        raise ValueError('Nine finite local parameters and measured two-vector perturber position required')
    d = distances(rear_redshift)
    # ol is independent of the rear-plane location, so the existing body map
    # remains exactly the same physical source prescription.
    body = local.body_from_parameters(p, compact_mass_Msun, compact_scale_kpc=compact_scale_kpc)
    if model=='conventional':
        lens = body
    elif model=='r1':
        lens = (ProfileR1Lens(body,old.frozen_coefficients()) if frozen_extra is None
                else old.CompositeContinuationLens(body,frozen_extra))
    else:
        raise ValueError('model must be conventional or r1')
    gamma = np.array([[p[4]+p[5],p[6]],[p[6],p[4]-p[5]]])
    beta_lp = d['os']*d['lp']/(d['op']*d['ls'])
    if min(np.linalg.eigvalsh(np.eye(2)-gamma))<=0 or min(np.linalg.eigvalsh(np.eye(2)-beta_lp*gamma))<=0:
        raise ValueError('Quadratic remainder leaves the adopted preconjugate branch')
    U = d['os']/(d['ps']*d['op'])*gamma
    theta_g = observed*ARCSEC
    alpha_g = lens.application(d['ol']*theta_g)[1]
    center = d['op']*theta_g-d['lp']*alpha_g
    rear = PositiveRearPlummer(float(perturber_mass_Msun),perturber_scale_kpc/1000,tuple(center))
    return dict(distances=d,lens=lens,rear=rear,U=U,beta0=p[2:4]*ARCSEC,
                perturber_observed_arcsec=observed,perturber_center_Mpc=center)


def arrival_action(theta, theta_rear, system):
    """Original full two-plane action and both independent partial gradients."""
    d=system['distances']; lens=system['lens']; rear=system['rear']; U=system['U']
    theta=np.asarray(theta,float); theta_rear=np.asarray(theta_rear,float)
    beta=system['beta0']; xp=d['op']*theta_rear
    pl,al,_=lens.application(d['ol']*theta)
    pp,ap,_=rear.application(xp)
    kp=(1+d['z_p'])*d['op']*d['os']/d['ps']
    kl=(1+d['z_l'])*d['ol']*d['op']/d['lp']
    action=(.5*kl*np.sum((theta-theta_rear)**2)+.5*kp*np.sum((theta_rear-beta)**2)
            -(1+d['z_l'])*pl-(1+d['z_p'])*(.5*xp@U@xp+pp))
    grad_l=kl*(theta-theta_rear)-(1+d['z_l'])*d['ol']*al
    grad_p=-kl*(theta-theta_rear)+kp*(theta_rear-beta)-(1+d['z_p'])*d['op']*(U@xp+ap)
    return dict(ct_Mpc=float(action),arrival_days=float(action*MPC_OVER_C_DAYS),
                gradient_theta=grad_l,gradient_theta_rear=grad_p)


def ray_application(theta, system):
    """Evaluate original ordered two-plane ray equation and its exact Jacobian."""
    d=system['distances']; lens=system['lens']; rear=system['rear']; U=system['U']
    theta=np.asarray(theta,float)
    pl,al,hl=lens.application(d['ol']*theta)
    xp=d['op']*theta-d['lp']*al
    pp,ap,hp=rear.application(xp)
    residual=theta-d['ls']/d['os']*al-d['ps']/d['os']*(U@xp+ap)-system['beta0']
    d_rear=d['op']*np.eye(2)-d['lp']*d['ol']*hl
    jacobian=(np.eye(2)-d['ls']/d['os']*d['ol']*hl
              -d['ps']/d['os']*(U+hp)@d_rear)
    result=arrival_action(theta,xp/d['op'],system)
    result.update(residual=residual,jacobian=jacobian,position_rear_Mpc=xp,
                  theta_rear=xp/d['op'],alpha_local=al,alpha_rear=U@xp+ap,
                  Hessian_local=hl,Hessian_rear=U+hp)
    return result


def macro_application(parameters, perturber_observed_arcsec, *, perturber_mass_Msun,
                      perturber_scale_kpc=.2,rear_redshift=.7821,
                      compact_mass_Msun=1e8,compact_scale_kpc=.03,
                      model='conventional',frozen_extra=None,astrometry=None,
                      require_success=True):
    system=components(parameters,perturber_observed_arcsec,
                      perturber_mass_Msun=perturber_mass_Msun,
                      perturber_scale_kpc=perturber_scale_kpc,rear_redshift=rear_redshift,
                      compact_mass_Msun=compact_mass_Msun,compact_scale_kpc=compact_scale_kpc,
                      model=model,frozen_extra=frozen_extra)
    observed=old.load_astrometry() if astrometry is None else astrometry
    images=[]; applications=[]; failures=[]
    for seed in observed['observed_images_light_centered_arcsec']:
        def fun(q):
            return ray_application(q*ARCSEC,system)['residual']/ARCSEC
        def jac(q):
            return ray_application(q*ARCSEC,system)['jacobian']
        solved=root(fun,seed,jac=jac,tol=2e-11)
        defect=float(max(abs(fun(solved.x))))
        if defect>2e-7 or not np.isfinite(solved.x).all():
            failures.append(defect)
        images.append(solved.x)
        applications.append(ray_application(solved.x*ARCSEC,system))
    images=np.asarray(images)
    separation=min(np.linalg.norm(images[i]-images[j]) for i in range(4) for j in range(i))
    if separation<.1:
        failures.append(float(.1-separation))
    if failures and require_success:
        raise ArithmeticError('Requested bright branches failed or merged')
    arrivals=np.array([app['arrival_days'] for app in applications])
    return dict(images_light_centered_arcsec=images,
                relative_astrometry_arcsec=np.r_[(images[1:]-images[0]).ravel(),-images[0]],
                delays_days=arrivals[1:]-arrivals[0],arrival_days=arrivals,
                magnifications=np.array([1/np.linalg.det(app['jacobian']) for app in applications]),
                max_lens_residual_rad=float(max(np.max(abs(app['residual'])) for app in applications)),
                max_arrival_stationarity_Mpc=float(max(max(np.max(abs(app['gradient_theta'])),np.max(abs(app['gradient_theta_rear']))) for app in applications)),
                minimum_image_separation_arcsec=float(separation),branch_failures=failures,
                system=system,applications=applications,physical_admission=False)


def image_census(parameters, perturber_observed_arcsec, *, perturber_mass_Msun,
                 perturber_scale_kpc=.2,rear_redshift=.7821,
                 compact_mass_Msun=1e8,compact_scale_kpc=.03,
                 model='conventional',extra_seeds_arcsec=()):
    """Deterministic broad search; report its finite scope, never completeness.

    Search about both observed galaxy centers and continue every conventional
    root into the exact R1 solution. Unlike the old constant-screen radial
    problem, the nonlinear rear plane has no degree-nine census polynomial.
    """
    kwargs=dict(perturber_mass_Msun=perturber_mass_Msun,
                perturber_scale_kpc=perturber_scale_kpc,rear_redshift=rear_redshift,
                compact_mass_Msun=compact_mass_Msun,compact_scale_kpc=compact_scale_kpc)
    conventional_system=components(parameters,perturber_observed_arcsec,**kwargs)
    p=np.asarray(parameters,float)
    observed=old.load_astrometry()['observed_images_light_centered_arcsec']
    seeds=list(observed)+list(extra_seeds_arcsec)
    angles=np.linspace(0.,2*np.pi,20,endpoint=False)
    directions=np.c_[np.cos(angles),np.sin(angles)]
    for center in [p[7:9],np.asarray(perturber_observed_arcsec,float)]:
        seeds.append(center)
        for radius in np.geomspace(1e-5,12.,28):
            seeds.extend(center+radius*directions)
    def verified_roots(system,seeds):
        images=[]
        for seed in seeds:
            def fun(q):
                return ray_application(q*ARCSEC,system)['residual']/ARCSEC
            def jac(q):
                return ray_application(q*ARCSEC,system)['jacobian']
            try:
                solved=root(fun,seed,jac=jac,tol=2e-11)
                if (np.isfinite(solved.x).all() and max(abs(fun(solved.x)))<2e-8
                        and not any(np.linalg.norm(solved.x-known)<1e-7 for known in images)):
                    images.append(solved.x)
            except (ArithmeticError,ValueError,FloatingPointError):
                continue
        return images
    images=verified_roots(conventional_system,seeds)
    system=conventional_system
    if model=='r1':
        system=components(parameters,perturber_observed_arcsec,model='r1',**kwargs)
        # Search with a temporary original-field interpolation, then verify all
        # found roots in the original exact projection. The local source used
        # for this interpolation is the same total two-Plummer body.
        projection=old.FrozenExtraProjection(system['lens'],count=400)
        fast=components(parameters,perturber_observed_arcsec,model='r1',frozen_extra=projection,**kwargs)
        candidates=verified_roots(fast,list(images)+seeds)
        images=verified_roots(system,candidates)
    elif model!='conventional':
        raise ValueError('model must be conventional or r1')
    images=np.array(sorted(images,key=lambda q:float(np.linalg.norm(q-p[7:9]))))
    apps=[ray_application(image*ARCSEC,system) for image in images]
    if not apps:
        raise ArithmeticError('No verified image in declared search')
    determinants=np.array([np.linalg.det(app['jacobian']) for app in apps])
    return dict(images_light_centered_arcsec=images.tolist(),image_count=len(images),
                signed_magnifications=(1/determinants).tolist(),absolute_magnifications=(1/abs(determinants)).tolist(),
                arrival_days=[app['arrival_days'] for app in apps],
                max_lens_residual_rad=float(max(np.max(abs(app['residual'])) for app in apps)),
                max_arrival_stationarity_Mpc=float(max(max(np.max(abs(app['gradient_theta'])),np.max(abs(app['gradient_theta_rear']))) for app in apps)),
                min_jacobian_singular_value=float(min(np.linalg.svd(app['jacobian'],compute_uv=False)[-1] for app in apps)),
                seeds_count=len(seeds),model=model,
                scope='Finite 20-angle,28-log-radius search at both observed galaxy centers to12arcsec; R1 continues all found conventional roots and interpolated-field seeds with final original-projection verification. No theorem excludes additional remote/critical branches.',
                physical_admission=False)
