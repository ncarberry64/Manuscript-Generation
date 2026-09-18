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
    e2 = ell*ell
    return (1.0/8.0) < e2 < (1.0/3.0)

def response_per_source(n: int, a: float, R_H: float, c_s2: float, chi: float, ell: float) -> float:
    p = p_n(n, a, R_H)
    c = C_n(n, c_s2, chi, ell)
    if p == 0.0:
        raise ZeroDivisionError("n=0 has no inverse-gradient response in this massless form")
    return 1.0/(p*c)

def quadrupole_to_dipole_response_ratio(c1: float, c2: float) -> float:
    return 3.0*c1/(8.0*c2)
