"""Explicit source reduction and conditional complement witnesses.

The R1 equations and the supplied BHSM source laws are kept distinct. These
applications test their normalization and matching requirements; they do not
select a microscopic surface-to-R1 field map or a new physical amplitude.
Natural R1 units are Mpl=H0=c=1 unless dimensions are stated otherwise.
"""
from __future__ import annotations

from dataclasses import dataclass
import math

import numpy as np


def reduced_source(density, scalar_source, b_T, Mpl=1.0):
    """J_eff=J_T-b_T rho/(2 Mpl²), from the common metric constraint."""
    if Mpl <= 0:
        raise ValueError('positive Planck normalization required')
    return scalar_source-b_T*density/(2*Mpl*Mpl)


def spherical_scalar_flux(radius, pi_prime, C, lambda3):
    """r² F_R1·n, F_R1=C grad pi-2[(Delta pi)grad pi-Hess pi grad pi]/lambda3."""
    if radius < 0 or lambda3 <= 0:
        raise ValueError('nonnegative radius and positive cubic scale required')
    return C*radius*radius*pi_prime-4*radius*pi_prime*pi_prime/lambda3


def scalar_shell_jump(radius, surface_density, scalar_surface_source, b_T, Mpl=1.0):
    """Jump of r² F_R1·n for a proper-radial delta layer in the retained weak metric."""
    if radius <= 0:
        raise ValueError('positive shell radius required')
    return radius*radius*reduced_source(surface_density, scalar_surface_source, b_T, Mpl)


def reduced_linear_system(A, B, D, J_P, J_Q):
    """Exact Schur operator AND source on an already-owned invertible complement.

    No pseudoinverse is used. Kernel coordinates and their compatibility
    equations must be retained separately when the physical D has a kernel.
    """
    A, B, D = map(lambda x: np.asarray(x, dtype=float), (A, B, D))
    J_P, J_Q = np.asarray(J_P, dtype=float), np.asarray(J_Q, dtype=float)
    solved_B = np.linalg.solve(D, B.T)
    solved_source = np.linalg.solve(D, J_Q)
    return A-B@solved_B, J_P-B@solved_source


def neighborhood_radial_flux(radius, seam_radius, source_amplitude):
    """F_top·n for div F_top=a Theta(r-r*)/r² and F_top(r*)=0.

    This distributed exterior source has growing integrated charge. It is
    not a finite-charge delta layer or the decaying R1 exterior flux.
    """
    if seam_radius <= 0 or radius < seam_radius:
        raise ValueError('requires radius >= seam_radius > 0')
    return source_amplitude*(radius-seam_radius)/(radius*radius)


@dataclass(frozen=True)
class ComplementWitness:
    """One positive conditional two-field lift of the same leading R1 C,D.

    All momenta are q=p ell0². mass_squared=(M_c ell0)² and mixing=g ell0.
    ell0 is an unspecified common reference unit, not a chosen galaxy scale.
    This is the existing 03_reduced_response quadratic family instantiated
    as a mathematical nonuniqueness witness, not a selected BHSM solution.
    """
    C: float
    D_time: float
    mass_squared: float
    complement_gradient: float
    mixing: float

    def __post_init__(self):
        if min(self.C, self.D_time, self.mass_squared, self.complement_gradient) <= 0:
            raise ValueError('positive retained and complement coefficients required')

    @property
    def Z0(self):
        return self.C+self.mixing*self.mixing/self.mass_squared

    def complement(self, q):
        if q < 0:
            raise ValueError('nonnegative spatial eigenvalue required')
        return self.mass_squared+self.complement_gradient*q

    def hessian(self, q):
        return np.array([[self.Z0*q, self.mixing*math.sqrt(q)],
                         [self.mixing*math.sqrt(q), self.complement(q)]])

    def reduced_kernel(self, q):
        return q*(self.Z0-self.mixing*self.mixing/self.complement(q))

    @property
    def p2_coefficient_in_ell0_squared(self):
        return self.mixing*self.mixing*self.complement_gradient/self.mass_squared**2

    @property
    def normalized_B_in_ell0_squared(self):
        return self.p2_coefficient_in_ell0_squared/self.C

    @property
    def expansion_scale_q(self):
        return self.mass_squared/self.complement_gradient

    def kinetic(self):
        return np.diag([self.D_time, 1.0])
