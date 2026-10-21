from pathlib import Path
import sys
import numpy as np
from scipy.integrate import quad
import pytest

ROOT=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT/'code'))
from strong_lensing_local_family import (
    TolmanVII,TolmanVIILens,phase_family_receipt,saved_geometry,compute_local_comparison,
)


@pytest.fixture(scope='module')
def calculation():
    saved,g,s=saved_geometry(ROOT/'artifacts/strong_lensing_first_comparison_20261020/results.json')
    return compute_local_comparison(saved,g,s)


def test_original_tov_and_einstein_equations_and_causal_eos():
    checks=TolmanVII().checks()
    assert checks['mass_normalization_residual']<1e-13
    assert checks['max_Hamiltonian_residual_scaled_by_compactness']<1e-11
    assert checks['max_radial_Einstein_residual_scaled_by_compactness']<1e-8
    assert checks['max_angular_Einstein_residual_scaled_by_compactness_squared']<1e-7
    assert checks['max_independent_TOV_relative_scale_residual']<1e-8
    assert checks['max_independent_lapse_Einstein_scaled_residual']<1e-8
    assert checks['sound_speed_condition'] and checks['dominant_energy_condition']
    assert 0<checks['max_sound_speed_squared_over_c_squared']<1e-4
    assert checks['surface_pressure_residual']<1e-13
    assert checks['surface_lapse_log_residual']<1e-15
    assert checks['surface_nu_derivative_residual']<1e-15


def test_projected_density_and_mass_are_the_same_finite_profile():
    lens=TolmanVIILens()
    # Direct line-of-sight projection rho(r), then projected disc mass.
    for y in (.05,.4,.8,.99):
        sigma=quad(lambda z:15/(8*np.pi)*(1-y*y-z*z),
            -np.sqrt(1-y*y),np.sqrt(1-y*y),epsabs=1e-12)[0]
        analytic_sigma=5/(2*np.pi)*(1-y*y)**1.5
        assert abs(sigma-analytic_sigma)<1e-13
        fraction=quad(lambda t:5*t*(1-t*t)**1.5,0,y,epsabs=1e-13,epsrel=1e-13)[0]
        assert abs(fraction-lens.projected_fraction(y))<1e-13
    assert lens.projected_fraction(1.)==1.


def test_potential_gradient_and_hessian_include_the_surface_and_centre():
    lens=TolmanVIILens()
    for v in (np.zeros(2),np.array([.2,.3])*lens.radius_Mpc,
              np.array([.8,.6])*lens.radius_Mpc,np.array([1.3,.2])*lens.radius_Mpc):
        step=lens.radius_Mpc*1e-5
        potential,alpha,hess=lens.application(v)
        gradient=[];jac=[]
        for i in range(2):
            dx=np.eye(2)[i]*step
            p=lens.application(v+dx);m=lens.application(v-dx)
            gradient.append((p[0]-m[0])/(2*step))
            jac.append((p[1]-m[1])/(2*step))
        np.testing.assert_allclose(gradient,alpha,rtol=3e-8,atol=1e-15)
        np.testing.assert_allclose(np.array(jac).T,hess,rtol=3e-8,atol=1e-10)


def test_actual_constant_s6_sections_satisfy_every_declared_eta_row():
    receipt=phase_family_receipt()
    for section in receipt['sections']:
        assert section['unit_constraint_residual']<1e-14
        assert section['eta_Euler_residual']==0
        assert section['color_Gauss_residual']==0
        assert section['eta_energy_and_stress']==0
        assert section['derivative_target_norm_squared']==1
    assert 'not identified' in receipt['topographic_matching']


def test_local_encounter_is_nonzero_but_bhsm_reference_equivalence_is_derived(calculation):
    sol=calculation['solution']
    assert len(sol['images_rad'])==3
    assert sol['max_lens_residual_rad']<1e-13
    assert np.max(abs(sol['relative_delays_days']))>1
    assert np.max(abs(calculation['matched_optical_difference']))==0
    assert np.max(abs(calculation['intrinsic_responses']['q']['images_arcsec']))==0
    assert np.max(abs(calculation['intrinsic_responses']['q']['relative_delays_days']))==0


def test_intrinsic_encounter_derivatives_reproduce_resolved_images_and_delays(calculation):
    for check in calculation['intrinsic_response_checks'].values():
        assert check['relative_whitened_error']<2e-7
    for key in ('log_mass','log_radius'):
        terms=calculation['intrinsic_responses'][key]['arrival_response_terms_days']
        assert np.max(abs(terms[:,:3]))==0
        assert np.max(abs(terms[:,3]))>0


def test_finite_profile_has_its_own_nuisance_projection(calculation):
    info=calculation['identifiability']
    assert info['q_only']['signal_rank']==0
    assert info['cosmological_fixed_nuisance']['signal_rank']==3
    assert info['cosmological_unrestricted']['signal_rank']==0
    assert calculation['finite_profile_factorization']['relative_whitened_discrepancy']<1e-10
