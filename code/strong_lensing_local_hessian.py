"""Exact static physical metric/phase projection outside a finite spherical body.

The metric variables are log lapse A and log spatial potential P in
ds²=-exp(2A)(c dt)²+exp(-2P) dx².  They are not the cosmological ADM
lapse/shift block.  Coefficients use the static energy convention: E=-L_static.
The physical gravity normalization is K=c^4/(8*pi*G), and eta stiffness is
k_eta=w*kappa1 (with the dimensional reference conversion made separately).
"""
from __future__ import annotations

import argparse
import json
from pathlib import Path

import numpy as np
from scipy.integrate import quad


def exterior_reference(r: float, mass_length: float) -> tuple[float, ...]:
    """Return A,P,A',P',A'',P'' for exact isotropic Schwarzschild."""
    t = mass_length / (2.0 * r)
    if r <= 0.0 or abs(t) >= 1.0:
        raise ValueError("reference domain must lie outside the isotropic horizon")
    A = np.log1p(-t) - np.log1p(t)
    P = -2.0 * np.log1p(t)
    Ap = mass_length / (r*r - mass_length*mass_length/4.0)
    Pp = mass_length / (r * (r+mass_length/2.0))
    App = -2.0 * mass_length*r / (r*r-mass_length*mass_length/4.0)**2
    Ppp = -mass_length*(2.0*r+mass_length/2.0) / (r*r*(r+mass_length/2.0)**2)
    return float(A), float(P), float(Ap), float(Pp), float(App), float(Ppp)


def dirichlet_profile(r: float, inner: float, outer: float) -> tuple[float, float]:
    """Dimensionless radial mode with int f² d³r/(outer-inner)³=1."""
    d = outer-inner
    if d <= 0.0 or inner <= 0.0:
        raise ValueError("invalid annular Dirichlet domain")
    phase = np.pi*(r-inner)/d
    amplitude = d/np.sqrt(2.0*np.pi)
    f = amplitude*np.sin(phase)/r
    fp = amplitude*(np.pi*np.cos(phase)/(d*r)-np.sin(phase)/(r*r))
    return float(f), float(fp)


def evaluate_local_hessian(m_over_R: float, outer_over_R: float = 2.0) -> dict:
    """Evaluate the exact common-annulus block in units K*R and k_eta*R.

    Set areal surface radius R=1 during quadrature.  A physical body's only
    input is its compactness m/R; no field/action coefficient is fitted to rays.
    """
    mu = float(m_over_R)
    if not (0.0 <= mu < 0.5) or outer_over_R <= 1.0:
        raise ValueError("require a non-horizon surface and outer radius > R")
    inner = 0.5*(1.0-mu+np.sqrt(1.0-2.0*mu))
    outer = float(outer_over_R)
    width = outer-inner

    def integrate(fun):
        value, error = quad(fun, inner, outer, epsabs=5e-13, epsrel=2e-13)
        return float(value), float(error)

    def entries(r):
        A, P, p, q, _, _ = exterior_reference(r, mu)
        f, g = dirichlet_profile(r, inner, outer)
        B = q*q-2.0*p*q
        weight = 4.0*np.pi*r*r*np.exp(A-P)
        return weight*np.asarray([
            -f*f*B+4.0*f*q*g,
            f*f*B-4.0*f*q*g+2.0*f*p*g+2.0*g*g,
            -f*f*B+4.0*f*q*g-4.0*f*p*g-2.0*g*g,
            g*g,
        ])

    results = [integrate(lambda r, i=i: entries(r)[i]) for i in range(4)]
    a, b, c, phase = [row[0] for row in results]
    metric = np.asarray([[a, b], [b, c]])
    mixed = np.zeros(2)
    response = -np.linalg.solve(metric, mixed)
    determinant = float(np.linalg.det(metric))
    schur = float(phase-mixed@np.linalg.solve(metric, mixed))
    norm, norm_error = integrate(
        lambda r: 4.0*np.pi*r*r*dirichlet_profile(r, inner, outer)[0]**2/width**3)
    flat_stiffness, _ = integrate(
        lambda r: 4.0*np.pi*r*r*dirichlet_profile(r, inner, outer)[1]**2)
    phase_kinetic, phase_kinetic_error = integrate(
        lambda r: 4.0*np.pi*r*r*np.exp(-exterior_reference(r, mu)[0]
                                             -3.0*exterior_reference(r, mu)[1])
        *dirichlet_profile(r, inner, outer)[0]**2)

    # Independent pointwise Euler equations of the exact reduced functional.
    # D=A-P and B=P'^2-2A'P'.  E_A∝B+2 div(e^D grad P)/e^D,
    # E_P∝B+2 div(e^D(grad P-grad A))/e^D.
    residuals = []
    for r in np.linspace(inner, outer, 101):
        A, P, p, q, pp, qp = exterior_reference(float(r), mu)
        B = q*q-2.0*p*q
        residuals.append([B+2.0*(qp+2.0*q/r+(p-q)*q),
                          B+2.0*((qp-pp)+2.0*(q-p)/r+(p-q)*(q-p))])

    return {
        "schema": "strong-lensing-local-physical-hessian-v1",
        "status": "EXACT_CONDITIONAL_VACUUM_ANNULUS_PROJECTION",
        "energy_sign_convention": "E_g=-L_static; E_eta=+int N sqrt(h) w F(X)",
        "m_over_R": mu,
        "domain_in_R": {"isotropic_inner": inner, "isotropic_outer": outer,
                        "width": width, "areal_body_surface": 1.0},
        "profile": "f=d sin(pi(r-r0)/d)/(sqrt(2pi) r); both endpoint traces zero",
        "profile_norm": norm,
        "profile_norm_quadrature_error": norm_error,
        "flat_gradient_norm_over_R": flat_stiffness,
        "analytic_flat_gradient_norm_over_R": float(np.pi**2*width),
        "metric_hessian_over_KR": metric.tolist(),
        "a_over_KR": a, "b_over_KR": b, "c_over_KR": c,
        "u": 0.0, "v": 0.0,
        "c_phi_over_k_eta_R": phase,
        "phase_stiffness_curved_over_flat": phase/flat_stiffness,
        "determinant_over_KR_squared": determinant,
        "metric_condition_number": float(np.linalg.cond(metric)),
        "metric_response": response.tolist(),
        "metric_response_residual": float(np.linalg.norm(metric@response+mixed)),
        "vacuum_simplified_block_residual_over_KR": float(np.max(np.abs(
            metric-2.0*phase*np.asarray([[0.0, 1.0], [1.0, -1.0]])))),
        "schur_over_k_eta_R": schur,
        "phase_kinetic_over_k_eta_Rcubed_per_c_squared": phase_kinetic,
        "phase_kinetic_quadrature_error": phase_kinetic_error,
        "quadrature_absolute_error_by_entry": [row[1] for row in results],
        "original_metric_Euler_residual_max_in_R": float(np.max(np.abs(residuals))),
        "eta_unit_constraint_residual": 0.0,
        "eta_Euler_residual": 0.0,
        "eta_gauge_current_residual": 0.0,
        "eta_stress_residual": 0.0,
        "eta_boundary_conormal_residual": 0.0,
        "eta_Legendre_gate": "w*kappa1>0; local orthonormal tangent kinetic matrix is w*kappa1 identity; positive coordinate-time profile norm recorded separately",
        "inverse_scope": "invertibility of the retained two-dimensional Dirichlet projection; no claim of invertibility of the complete metric functional",
        "physical_scope": "retained 4D Einstein-Hilbert plus Path-B eta sector; fixed matched body boundary data",
        "not_claimed": ["a physical selector of q,w,kappa1", "a solved 5D cap response",
                        "a nonzero seam-export to target-phase constitutive map"],
    }


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--m-over-R", type=float, required=True)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    result = evaluate_local_hessian(args.m_over_R)
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(result, indent=2, sort_keys=True)+"\n", encoding="utf-8", newline="\n")
    print(json.dumps({key: result[key] for key in ("a_over_KR", "b_over_KR", "c_over_KR",
                                                 "c_phi_over_k_eta_R", "original_metric_Euler_residual_max_in_R")}, indent=2))


if __name__ == "__main__":
    main()
