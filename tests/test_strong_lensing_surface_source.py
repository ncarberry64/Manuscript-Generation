from pathlib import Path
import sys

import numpy as np
from scipy.integrate import quad

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT/'code'))
sys.path.insert(0, str(ROOT/'scripts'))
from strong_lensing_surface_source import (
    ComplementWitness, neighborhood_radial_flux, reduced_linear_system,
    reduced_source, scalar_shell_jump, spherical_scalar_flux,
)
from run_strong_lensing_surface_source import evaluate


def test_consumed_source_matches_the_coupled_original_metric_scalar_equations():
    A, b, Mpl, rho, J_T = 10., 3., 2., 7., 11.
    # Original rows: A Delta pi+b Delta Phi=J_T;
    # b Delta pi+2 Mpl² Delta Phi=rho.
    pi_laplacian, phi_laplacian = np.linalg.solve(
        [[A, b], [b, 2*Mpl*Mpl]], [J_T, rho])
    C = A-b*b/(2*Mpl*Mpl)
    np.testing.assert_allclose(C*pi_laplacian, reduced_source(rho, J_T, b, Mpl), rtol=2e-15)
    np.testing.assert_allclose(2*Mpl*Mpl*phi_laplacian, rho-b*pi_laplacian, rtol=2e-15)


def test_schur_source_preserves_a_nonzero_complement_force():
    A = np.array([[4., .2], [.2, 3.]])
    B = np.array([[.4], [.7]])
    D = np.array([[2.]])
    J_P, J_Q = np.array([1., -2.]), np.array([.75])
    reduced, source = reduced_linear_system(A, B, D, J_P, J_Q)
    exact = np.linalg.solve(np.block([[A, B], [B.T, D]]), np.r_[J_P, J_Q])
    np.testing.assert_allclose(np.linalg.solve(reduced, source), exact[:2], rtol=2e-15)
    assert np.linalg.norm(np.linalg.solve(reduced, J_P)-exact[:2]) > .01


def test_nonlinear_flux_reproduces_a_regular_uniform_density_source():
    C, lam, b, rho = 4.2, 1., .113, 5.
    # pi'=k r; div F=3 C k-12 k²/lambda=-b rho/2.
    k = lam/8*(C-np.sqrt(C*C+8*b*rho/(3*lam)))
    for radius in (.1, 1., 3.):
        volume_source = quad(lambda r: r*r*reduced_source(rho, 0., b), 0., radius)[0]
        np.testing.assert_allclose(spherical_scalar_flux(radius, k*radius, C, lam), volume_source, rtol=2e-13)


def test_surface_jump_has_the_correct_areal_normalization():
    R, sigma, J_T_surface, b, width = 2., 3., 5., .113, 1e-4
    integrated_layer = 4*np.pi*quad(
        lambda r: r*r*reduced_source(sigma/(2*width), J_T_surface/(2*width), b),
        R-width, R+width)[0]
    zero_width_limit = 4*np.pi*scalar_shell_jump(R, sigma, J_T_surface, b)
    np.testing.assert_allclose(integrated_layer/zero_width_limit, 1+width*width/(3*R*R), rtol=3e-12)


def test_distributed_old_seam_source_is_not_a_finite_charge_layer():
    r_star, amplitude = 2., 3.
    assert neighborhood_radial_flux(r_star, r_star, amplitude) == 0.
    for radius in (3., 10., 100.):
        integrated_source = quad(lambda r: amplitude, r_star, radius)[0]
        np.testing.assert_allclose(radius*radius*neighborhood_radial_flux(radius, r_star, amplitude), integrated_source, rtol=2e-15)
    # The asymptotic charge keeps growing, unlike the R1 isolated scalar charge.
    f1 = neighborhood_radial_flux(1e5, r_star, amplitude)
    f2 = neighborhood_radial_flux(2e5, r_star, amplitude)
    assert 1.999 < (2e5)**2*f2/((1e5)**2*f1) < 2.001


def test_two_same_domain_positive_lifts_do_not_determine_B_or_cutoff():
    C, D_time = 4.238858559433196, 4.396221028708471
    first = ComplementWitness(C, D_time, 1., 1., 1.)
    second = ComplementWitness(C, D_time, 4., 1., 2.)
    for witness in (first, second):
        assert witness.Z0-witness.mixing**2/witness.mass_squared == C
        assert np.linalg.eigvalsh(witness.kinetic()).min() > 0
        for q in (np.arange(1, 21)/10)**2:
            assert np.linalg.eigvalsh(witness.hessian(float(q))).min() > 0
            assert witness.reduced_kernel(float(q)) > C*q
        # The explicit normalized radial basis owns the same L2 domain.
        L = 10*np.pi
        for n in (1, 3):
            norm = quad(lambda r: 2/L*np.sin(n*np.pi*r/L)**2, 0., L)[0]
            np.testing.assert_allclose(norm, 1., rtol=2e-14)
    assert first.normalized_B_in_ell0_squared == 4*second.normalized_B_in_ell0_squared
    assert second.expansion_scale_q == 4*first.expansion_scale_q
    assert first.reduced_kernel(.5) != second.reduced_kernel(.5)


def test_actual_computed_sections_satisfy_volume_and_material_matching_equations():
    result = evaluate()
    assert result['status'] == 'PASS'
    assert not result['physical_admission'] and not result['fitted_new_amplitude']
    rows = result['computed_section_and_original_residuals']
    assert all(row['pi_prime_in_MplH0'] < 0 for row in rows)
    assert all(row['computed_nonlinear_scalar_flux'] < 0 for row in rows)
    assert max(row['scalar_flux_relative_residual'] for row in rows) < 1e-12
    shell = result['material_matching']
    assert shell['material_surface_Q_Sigma'] == 0.
    assert shell['relative_flux_jump'] < 1e-12
    assert shell['leading_R1_delta_K_uu_per_Mpc'] == 0.
    assert shell['leading_R1_K_uu_inside_per_Mpc'] != shell['K_uu_inside_per_Mpc']
    assert result['effective_mass_application']['extra_enclosed_mass_at_R_Msun'] > 0
    assert result['effective_mass_application']['extra_total_mass_Msun'] > result['effective_mass_application']['extra_enclosed_mass_at_R_Msun']
