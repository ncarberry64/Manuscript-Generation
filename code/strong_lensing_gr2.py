"""Genuine 2PM GR optical control for a prescribed finite Tolman-VII body.

The independent impact coordinate b is the invariant asymptotic scattering
impact, not the closest areal radius. In the host lens calculation it is
identified with the physical thin-plane impact. This controls the ordinary-GR
part of that local optical approximation; it does not supply scalar 2PM terms
or an exact finite-distance cosmological arrival-time calculation.

For epsilon=m/R and y=b/R, delta_alpha=2 epsilon**2 I(y) and
delta_psi_hat=2 epsilon**2 R J(y), with J(1)=0 and J'=I.
The coefficients include the leading TOV pressure and areal/optical-radius
conversion. The exterior is exactly 15*pi/4*(m/b)**2 at this order.
"""
from __future__ import annotations

import math
import numpy as np
from scipy.integrate import quad


def _asin_minus_yt(y):
    """asin(y)-y*sqrt(1-y*y), evaluated without small-y cancellation."""
    if y < 0.05:
        # D'=2y^2/sqrt(1-y^2). Integrate its binomial series.
        z = y*y
        total = 0.0
        coefficient = 1.0
        power = 1.0
        for n in range(12):
            total += 2*coefficient*power/(2*n+3)
            coefficient *= (2*n+1)/(2*n+2)
            power *= z
        return y*y*y*total
    t = math.sqrt(max(0.0, (1-y)*(1+y)))
    return math.asin(y)-y*t


def coefficient_I(y):
    y = float(y)
    if y < 0:
        raise ValueError('nonnegative impact required')
    if y == 0:
        return 0.0
    if y >= 1:
        return 15*math.pi/(8*y*y)
    t = math.sqrt((1-y)*(1+y))
    return 15*_asin_minus_yt(y)/(4*y*y) + y*(7.5*t-4*t**5-48/7*t**7)


def coefficient_I_prime(y):
    y = float(y)
    if y < 0:
        raise ValueError('nonnegative impact required')
    if y == 0:
        return -6/7
    if y >= 1:
        return -15*math.pi/(4*y**3)
    t = math.sqrt((1-y)*(1+y))
    return (-15*_asin_minus_yt(y)/(2*y**3)
            +15*t+20*t**3+24*t**5-384/7*t**7)


def coefficient_J(y):
    """Primitive of I, with zero at y=1."""
    y = float(y)
    if y < 0:
        raise ValueError('nonnegative impact required')
    if y == 0:
        return 15*math.pi/8-26/3
    if y >= 1:
        return 15*math.pi/8*((y-1)/y)
    t = math.sqrt((1-y)*(1+y))
    if y > 0.8:
        # Algebraically identical form, stable as y approaches one.
        first = -15*math.pi/8*t*t/(y*(1+y))
        return (first + 15/(4*y)*_asin_minus_yt(t)
                -2.5*t**3+4/7*t**7+16/21*t**9)
    return (15*math.pi/8-15/4*(math.asin(y)/y+t)
            -2.5*t**3+4/7*t**7+16/21*t**9)


class TolmanVII2PMCorrection:
    """Pure GR 2PM correction; potentials in Mpc, angles dimensionless."""

    def __init__(self, base):
        self.mass_length_Mpc = float(base.mass_length_Mpc)
        self.radius_Mpc = float(base.radius_Mpc)
        self.centre_Mpc = np.asarray(base.centre_Mpc, dtype=float)
        self.epsilon = self.mass_length_Mpc/self.radius_Mpc

    def extra_alpha(self, impact_Mpc):
        return 2*self.epsilon**2*coefficient_I(float(impact_Mpc)/self.radius_Mpc)

    def extra_alpha_prime(self, impact_Mpc):
        return 2*self.epsilon**2/self.radius_Mpc*coefficient_I_prime(float(impact_Mpc)/self.radius_Mpc)

    def extra_potential(self, impact_Mpc):
        return 2*self.epsilon**2*self.radius_Mpc*coefficient_J(float(impact_Mpc)/self.radius_Mpc)

    def extra_application(self, x):
        displacement = np.asarray(x, dtype=float)-self.centre_Mpc
        r = float(np.linalg.norm(displacement))
        p = self.extra_potential(r)
        if r == 0:
            return p, np.zeros(2), self.extra_alpha_prime(0.)*np.eye(2)
        direction = displacement/r
        alpha = self.extra_alpha(r)
        derivative = self.extra_alpha_prime(r)
        hessian = alpha/r*np.eye(2)+(derivative-alpha/r)*np.outer(direction,direction)
        return p, alpha*direction, hessian


class TolmanVII2PMLens:
    """The original finite-body 1PM lens plus the common ordinary-GR 2PM term.

    Can be passed as the base of the existing R1ScreenedLens so that exactly
    the same GR correction enters the conventional and scalar-on models.
    """

    def __init__(self, base):
        self.base = base
        self.correction = TolmanVII2PMCorrection(base)
        for name in ('mass_Msun','radius_Mpc','centre_Mpc','mass_length_Mpc','b'):
            setattr(self, name, getattr(base,name))

    def application(self, x):
        p,a,h = self.base.application(x)
        dp,da,dh = self.correction.extra_application(x)
        return p+dp,a+da,h+dh

    def intrinsic(self, coordinate):
        original = self.base.intrinsic(coordinate)
        def partial(x):
            p,a = original(x)
            dp,da,_ = self.correction.extra_application(x)
            if coordinate == 'log_mass':
                return p+2*dp,a+2*da
            if coordinate == 'log_radius':
                r = float(np.linalg.norm(np.asarray(x)-np.asarray(self.centre_Mpc)))
                radial_extra = self.correction.extra_alpha(r)
                direction = np.zeros(2) if r == 0 else (np.asarray(x)-np.asarray(self.centre_Mpc))/r
                return (p-dp-r*radial_extra,
                        a-2*da-r*self.correction.extra_alpha_prime(r)*direction)
            return p,a
        return partial


def exact_static_deflection(body, impact_Mpc, *, epsabs=1e-19, epsrel=3e-13):
    """Flat-subtracted exact static GR scattering integral for a solved TOV body.

    body is strong_lensing_local_family.TolmanVII, with its solved areal metric.
    Use optical radius rho=r exp(-nu). The integral then has fixed endpoints
    and no turning-point square-root singularity:

      alpha = 2 int_0^(pi/2) [1/(sqrt(f)*(1-r nu'))-1] dtheta,
      rho=b/sin(theta), f=1-2 m(r)/r.

    expm1/log1p preserve the small GR deflection instead of subtracting pi.
    The returned quadrature error does not include the underlying TOV error.
    """
    b = float(impact_Mpc)
    if b < 0:
        raise ValueError('nonnegative impact required')
    if b == 0:
        return 0.0,0.0
    R = float(body.radius_Mpc)
    e = float(body.compactness)
    m = float(body.mass_length_Mpc)

    def fields(r):
        x = r/R
        if x >= 1:
            u = m/r
            logf = math.log1p(-2*u)
            return 0.5*logf, logf, u/(1-2*u)
        pbar,nubar = map(float,body.solution.sol(x))
        mu_over_x = .5*x*x*(5-3*x*x)
        u = e*mu_over_x
        logf = math.log1p(-2*u)
        rnuprime = e*(mu_over_x+4*math.pi*e*x*x*pbar)/(1-2*u)
        return e*nubar,logf,rnuprime

    def integrand(theta):
        sine = math.sin(theta)
        if sine == 0:
            return 0.0
        rho = b/sine
        r = rho
        for _ in range(12):
            nu,logf,rnuprime = fields(r)
            if rnuprime >= 1:
                raise ValueError('optical radius must be monotone on this weakly compact branch')
            residual = math.log1p((r-rho)/rho)-nu
            change = r*residual/(1-rnuprime)
            r -= change
            if abs(change) <= 3e-16*r:
                break
        else:
            raise ArithmeticError('optical-radius inversion did not converge')
        _,logf,rnuprime = fields(r)
        return math.expm1(-0.5*logf-math.log1p(-rnuprime))

    optical_surface = R/math.sqrt(1-2*e)
    points = [math.asin(b/optical_surface)] if b < optical_surface else None
    value,error = quad(integrand,0,math.pi/2,points=points,epsabs=epsabs,
                       epsrel=epsrel,limit=160)
    return 2*value,2*error
