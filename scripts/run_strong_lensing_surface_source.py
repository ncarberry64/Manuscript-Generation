"""Apply the derived source and junction identities to retained sections.

No background shooting, light preparation, image solve, amplitude fit or
historical artifact is repeated or changed. Complement models are explicit
mathematical witnesses of nonuniqueness, not physical BHSM admissions.
"""
from __future__ import annotations

import argparse
from dataclasses import asdict
import hashlib
import json
from pathlib import Path
import sys

import numpy as np
from scipy.integrate import quad

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT/'code'))
from strong_lensing_local_family import TolmanVII, TolmanVIILens
from strong_lensing_r1_screened import R1Coefficients, R1ScreenedLens
from strong_lensing_surface_source import ComplementWitness, reduced_source, spherical_scalar_flux, scalar_shell_jump


def evaluate():
    source = ROOT/'artifacts/strong_lensing_r1_screened_20261020/results.json'
    saved = json.loads(source.read_text(encoding='utf-8'))
    coefficients = R1Coefficients(**saved['coefficients'])
    base = TolmanVIILens(mass_Msun=saved['body']['mass_Msun'], radius_Mpc=saved['body']['radius_Mpc'])
    lens = R1ScreenedLens(base, coefficients)
    b, C, lam = coefficients.b, coefficients.C, coefficients.lambda_over_MplH02
    L_H = coefficients.hubble_length_Mpc
    R = base.radius_Mpc/L_H
    m = base.mass_length_Mpc/L_H
    D = 1+12*coefficients.H_over_H0*coefficients.Tdot_over_MplH0/lam+6*coefficients.Tdot_over_MplH0**4/lam**2

    def density(r):
        # rho/(Mpl² H0²), G=1/(8 pi), r in c/H0 units.
        return 15*m/R**3*max(1-(r/R)**2, 0.)

    rows = []
    for radius_over_R in (.1, .5, 1., 2., 10.):
        r = radius_over_R*R
        pi_prime = lens.scalar_gradient_natural(r*L_H)
        flux = spherical_scalar_flux(r, pi_prime, C, lam)
        rhs = quad(lambda s: s*s*reduced_source(density(s), 0., b),
                   0., min(r, R), epsabs=1e-27, epsrel=2e-13)[0]
        mu = float(TolmanVII.mass_fraction(radius_over_R))
        phi_prime = m*mu/(r*r)-b*pi_prime/2
        metric_flux = r*r*(phi_prime+b*pi_prime/2)
        rows.append(dict(radius_over_R=radius_over_R, radius_Mpc=r*L_H,
                         pi_prime_in_MplH0=pi_prime, metric_phi_prime_in_H0=phi_prime,
                         computed_nonlinear_scalar_flux=flux,
                         independent_integrated_effective_source=rhs,
                         scalar_flux_relative_residual=abs(flux-rhs)/abs(rhs),
                         original_metric_Gauss_relative_residual=abs(metric_flux-m*mu)/(m*mu)))

    p_R = lens.scalar_gradient_natural(base.radius_Mpc)
    # Independent direct quadratic root, compared with the computed rational root.
    p_R_interior = lam*C*R/8*(1-np.sqrt(1+16*b*m/(lam*C*C*R**3)))
    p_R_exterior = p_R
    flux_jump = spherical_scalar_flux(R, p_R_exterior, C, lam)-spherical_scalar_flux(R, p_R_interior, C, lam)
    expected_jump = scalar_shell_jump(R, 0., 0., b)
    compactness = base.mass_length_Mpc/base.radius_Mpc
    f_surface = 1-2*compactness
    nu_prime_inside = base.mass_length_Mpc/(base.radius_Mpc**2*f_surface)  # p(R)=0 radial Einstein row.
    nu_prime_outside = compactness/(base.radius_Mpc*f_surface)  # Schwarzschild row.
    K_uu_inside = -nu_prime_inside*np.sqrt(f_surface)
    K_uu_outside = -nu_prime_outside*np.sqrt(f_surface)
    leading_R1_K_inside = -(m/R**2-b*p_R_interior/2)/L_H
    leading_R1_K_outside = -(m/R**2-b*p_R_exterior/2)/L_H
    shell = dict(surface='material timelike worldtube r=R; distinct from the cosmological cap equator',
                 orientation='same +dr unit normal on both material limits',
                 scalar_boundary_data='continuous pi_T and nonlinear canonical flux; no explicit scalar layer in this R1 action',
                 density_surface_in_Mpl2H02=density(R), surface_pressure=0.,
                 pi_prime_interior_quadratic_root=p_R_interior,
                 pi_prime_exterior_computed_root=p_R_exterior,
                 nonlinear_flux_jump=flux_jump, required_scalar_shell_jump=expected_jump,
                 relative_flux_jump=abs(flux_jump)/abs(rows[2]['computed_nonlinear_scalar_flux']),
                 K_uu_inside_per_Mpc=K_uu_inside, K_uu_outside_per_Mpc=K_uu_outside,
                 delta_K_uu_per_Mpc=K_uu_outside-K_uu_inside,
                 K_uu_evaluation_scope='ordinary-GR reference: exact TOV surface p(R)=0 matched to Schwarzschild; not the rolling-scalar metric',
                 leading_R1_K_uu_inside_per_Mpc=leading_R1_K_inside,
                 leading_R1_K_uu_outside_per_Mpc=leading_R1_K_outside,
                 leading_R1_delta_K_uu_per_Mpc=leading_R1_K_outside-leading_R1_K_inside,
                 leading_R1_K_uu_scope='leading weak-metric K_uu=-Phi_prime, evaluated from the R1 metric Gauss row on both material limits; no exact modified-TOV junction claimed',
                 material_surface_Q_Sigma=0.,
                 qualification='The GR reference junction and leading R1 junction are separate applications at their stated orders. Their zero material seam functional follows from smooth matching; it does not set the volume matter-curvature source or all BHSM cap sources to zero.')

    witnesses = []
    for name, M2, g in [('A', 1., 1.), ('B', 4., 2.)]:
        witness = ComplementWitness(C, D, M2, 1., g)
        q = (np.arange(1, 21)/10)**2
        eigenvalues = np.array([np.linalg.eigvalsh(witness.hessian(float(p))) for p in q])
        witnesses.append(dict(name=name, parameters=asdict(witness), Z0=witness.Z0,
                              leading_spatial_C=witness.Z0-g*g/M2,
                              leading_temporal_D=D,
                              p2_coefficient_in_ell0_squared=witness.p2_coefficient_in_ell0_squared,
                              normalized_B_in_ell0_squared=witness.normalized_B_in_ell0_squared,
                              derivative_expansion_scale_p_ell0_squared=witness.expansion_scale_q,
                              complement_spectrum=[witness.complement(float(p)) for p in q],
                              q_spectrum=q.tolist(),
                              full_spatial_hessian_eigenvalues=eigenvalues.tolist(),
                              minimum_spatial_eigenvalue=float(eigenvalues.min()),
                              kinetic_eigenvalues=np.linalg.eigvalsh(witness.kinetic()).tolist(),
                              reduced_kernels=[witness.reduced_kernel(float(p)) for p in q]))
    mass_normalization = 4*np.pi*quad(lambda r: r*r*density(r), 0., R, epsabs=1e-25)[0]
    max_scalar = max(r['scalar_flux_relative_residual'] for r in rows)
    max_metric = max(r['original_metric_Gauss_relative_residual'] for r in rows)
    passed = (max_scalar < 1e-12 and max_metric < 1e-12 and shell['relative_flux_jump'] < 1e-12
              and abs(shell['delta_K_uu_per_Mpc']) < 1e-12
              and all(w['minimum_spatial_eigenvalue'] > 0 and min(w['kinetic_eigenvalues']) > 0 for w in witnesses))
    paths = [source, ROOT/'code/strong_lensing_surface_source.py', Path(__file__),
             ROOT/'code/strong_lensing_r1_screened.py', ROOT/'code/strong_lensing_local_family.py',
             ROOT/'manuscript/sections/03_reduced_response.tex']
    return dict(status='PASS' if passed else 'FAIL',
                classification='DERIVED_EFFECTIVE_SOURCE_APPLICATION_AND_CONDITIONAL_NONUNIQUENESS_PROOF_NOT_MICROSCOPIC_CLOSURE',
                source_input='frozen supplied finite-body R1 coefficient receipt; no background refit',
                equations=dict(metric='2 Mpl² Delta Phi=rho+rho_Sigma-b_T Delta pi_T',
                               scalar_before_elimination='A Delta pi_T-2[(Delta pi_T)²-|Hess pi_T|²]/lambda+b_T Delta Phi=J_T',
                               consumed_source='J_eff=J_T-b_T(rho+rho_Sigma)/(2 Mpl²)',
                               scalar_flux='r² F_R1·n=C r² pi_T_prime-4 r pi_T_prime²/lambda',
                               shell='[r² F_R1·n]=R²(J_T,Sigma-b_T sigma_Sigma/(2 Mpl²))'),
                normalization=dict(Mpl_H0_c=1, length_unit_Mpc=L_H, radius_natural=R,
                                   mass_length_natural=m, density_integral=mass_normalization,
                                   required_density_integral=8*np.pi*m,
                                   integrated_scalar_source=-b*m),
                computed_section_and_original_residuals=rows,
                material_matching=shell,
                effective_mass_application=dict(extra_enclosed_mass_at_R_Msun=base.mass_Msun*lens.extra_force(base.radius_Mpc),
                                               extra_total_mass_Msun=base.mass_Msun*coefficients.unscreened_extra,
                                               positive_density_proof='dln(f_M F)/dlnr=[(1+1/S)n_M+3(1-1/S)]/2>=0'),
                complement_nonuniqueness=dict(
                    status='TWO_EXPLICIT_POSITIVE_QUADRATIC_LIFTS_SAME_LEADING_C_D_DIFFERENT_B_AND_CUTOFF',
                    domain='same ball radius 10 pi ell0; both real fields regular at origin and Dirichlet at boundary; radial q_n=(n/10)²',
                    normalized_radial_basis='f_n(r)=sin(n pi r/L)/(r sqrt(2 pi L)); integral 4 pi r² f_n² dr=1',
                    source_data='same retained source J_P, explicitly unforced complement J_Q=0 for this mathematical counterexample; no quantum vacuum is chosen',
                    witnesses=witnesses,
                    proof='D_QQ=M²+c q>0 and K_red=C q+g² c q²/[M²(M²+c q)]>0 for q>0; diag(D_time,1)>0. q=0 is a formal limit excluded by this Dirichlet domain. A kernel on another physical domain would have to remain explicit.',
                    qualification='ell0 and the embedding of either complement into the BHSM domain are not supplied. Neither witness is selected for the galaxy or asserted as the microscopic cutoff.'),
                old_seam_direct_identification=dict(
                    inherited_source='S_top=a Theta(r-r_star)/r², a=epsilon_Sigma Q_Sigma/Q',
                    integrated_source='4 pi a(r-r_star)',
                    radial_flux='F_top·n=a(1/r-r_star/r²), with F_top(r_star)=0',
                    direct_rescaling_obstruction='For nonzero a the exterior flux falls as 1/r, while the isolated R1 metric correction and canonical scalar flux fall as 1/r². No constant field normalization identifies them on an unbounded exterior.',
                    scope='Rejects that direct identification only; a different owned finite-support export/operator/nonlocal map remains possible.'),
                remaining_defining_relation='An action-owned physical embedding e_T and source pullback must identify the retained cap/surface sector with canonical R1 T and its cubic vertex; the supplied Schur family does not fix its complement spectrum or physical ell0.',
                physical_admission=False, fitted_new_amplitude=False,
                sha256={p.relative_to(ROOT).as_posix():hashlib.sha256(p.read_bytes()).hexdigest() for p in paths})


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--output', type=Path, default=ROOT/'artifacts/strong_lensing_bhsm_integration_20261020/source_application.json')
    args = parser.parse_args()
    result = evaluate()
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(result, indent=2, allow_nan=False)+'\n', encoding='utf-8', newline='\n')
    print(json.dumps({key:result[key] for key in ('status', 'normalization', 'material_matching', 'effective_mass_application')}, indent=2))
    raise SystemExit(0 if result['status'] == 'PASS' else 1)


if __name__ == '__main__':
    main()
