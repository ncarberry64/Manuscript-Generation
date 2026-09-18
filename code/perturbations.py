"""Minimal harmonic-filter utilities for the BHSM/topographic bridge."""
from __future__ import annotations
import math

def nu(n: int) -> int:
    if n < 0:
        raise ValueError("n must be non-negative")
    return n*(n+2)

def p_n(n: int, a: float, R_H: float) -> float:
    return nu(n)/(a*a*R_H*R_H)

def C_n(n: int, c_s2: float, chi: float, ell: float) -> float:
    return c_s2 - chi*(1.0 - ell*ell*nu(n))

def selective_filter_condition(ell: float) -> bool:
    ell2 = ell * ell
    return (1.0 / 15.0) < ell2 < (1.0 / 8.0)

def response_per_source(n: int, a: float, R_H: float, c_s2: float, chi: float, ell: float) -> float:
    p = p_n(n, a, R_H)
    c = C_n(n, c_s2, chi, ell)
    if p == 0.0:
        raise ZeroDivisionError("n=0 has no inverse-gradient response in this massless form")
    return 1.0/(p*c)

def n2_to_n1_response_ratio(c1: float, c2: float) -> float:
    """Source-normalized S3 n=2 to n=1 susceptibility ratio.

    This is not, without observer-shell projection, a sky-quadrupole/dipole ratio.
    """
    return 3.0*c1/(8.0*c2)

def quadrupole_to_dipole_response_ratio(c1: float, c2: float) -> float:
    """Legacy compatibility alias; prefer n2_to_n1_response_ratio()."""
    return n2_to_n1_response_ratio(c1, c2)


def physical_scalar_harmonic_min() -> int:
    return 2

def n2_pure_dipole_envelope(psi):
    import numpy as np
    return np.sin(2.0 * np.asarray(psi))

def n3_to_n2_response_ratio(c2: float, c3: float) -> float:
    return 8.0 * c2 / (15.0 * c3)
