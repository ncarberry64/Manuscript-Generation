"""Focused checks of the new physical-boundary finite-window control."""
from dataclasses import replace
import json
from pathlib import Path
import sys

import numpy as np
from scipy.integrate import quad

ROOT=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT/'code'))
from strong_lensing_encounter_accuracy import (
    PhysicalBoundaryGeometry,RetainedEncounterField,distributed_ray,
)


def inputs():
    report=json.loads((ROOT/'artifacts/strong_lensing_r1_screened_20261020/results.json').read_text())
    manifest=json.loads((ROOT/'preregistration/prediction_manifest.json').read_text())
    geometry=PhysicalBoundaryGeometry.from_saved(report,manifest['background_R1']['Omega_k0'])
    return report,geometry,RetainedEncounterField.from_saved(report)


def test_individual_wavefronts_follow_physical_closed_flrw_boundary_data():
    report,g,_=inputs()
    ko,ks,j=g.exterior(0.)
    # Individual curvatures differ from a freely chosen flat-space split.
    assert abs(ko-1/g.D_ol_Mpc)>1e-8
    np.testing.assert_allclose(ko+ks,g.D_os_Mpc/(g.D_ol_Mpc*g.D_ls_Mpc),rtol=3e-15)
    np.testing.assert_allclose(j,report['experiment']['source_transverse_Mpc'][0]/g.D_ls_Mpc,rtol=3e-15)
    np.testing.assert_allclose(g.S(g.observer_length_Mpc),g.D_ol_Mpc,rtol=2e-15)


def test_zero_field_distributed_ray_preserves_the_same_unlensed_preparation():
    report,g,f=inputs()
    zero=replace(f,mass_length_Mpc=0.,rV_Mpc=0.,epsilon_infinity=0.)
    theta=report['experiment']['source_transverse_Mpc'][0]/g.D_os_Mpc
    ray=distributed_ray(g,zero,2.,False,theta)
    np.testing.assert_allclose(ray['images_rad'],theta,rtol=2e-12,atol=1e-20)
    assert abs(ray['integrated_momentum_residual_rad'])<1e-18


def test_potential_differences_agree_with_an_independent_radial_integral():
    _,_,f=inputs()
    for x1,x0,z in [(0.,.002,0.),(.004,.001,.0002),(.003,.004,20.)]:
        r0=np.hypot(x0,z);r1=np.hypot(x1,z)
        dr=(x1-x0)*(x1+x0)/(r1+r0)
        for on in [False,True]:
            independent=quad(lambda q:float(f.force(r0+q*dr,on))*dr,0,1,
                             epsabs=1e-32,epsrel=2e-12,limit=160,
                             points=[(f.radius_Mpc-r0)/dr] if dr and 0<(f.radius_Mpc-r0)/dr<1 else None)[0]
            np.testing.assert_allclose(f.potential_difference(x1,x0,z,on),independent,rtol=4e-12,atol=1e-29)


def test_finite_window_alpha_has_the_original_line_of_sight_normalization():
    _,_,f=inputs();b=.002;L=.7
    for on in [False,True]:
        points=np.geomspace(b,L,22)[:-1]
        independent=4*b*quad(lambda z:float(f.force(np.hypot(b,z),on))/np.hypot(b,z),
                             0,L,points=points,epsabs=1e-22,epsrel=3e-13,limit=200)[0]
        np.testing.assert_allclose(f.alpha(b,L,on),independent,rtol=3e-12)
    assert f.alpha(0.,L,True)==0.


def test_extra_potential_tail_is_bounded_by_the_unscreened_envelope():
    _,_,f=inputs();b1=.00038;b0=.00383;L=20.
    def integrand(q):
        if q==0:return 0.
        z=L/q;r0=np.hypot(b0,z);r1=np.hypot(b1,z)
        dr=(b1-b0)*(b1+b0)/(r1+r0)
        difference=f._extra_radial_integral(np.array([r0]),np.array([dr]))[0]
        return difference*L/q**2
    actual=4*abs(quad(integrand,0,1,epsabs=1e-32,epsrel=2e-12)[0])
    bound=f.tail_bounds(b1,b0,L)
    assert actual<=bound['potential_difference_Mpc']*(1+2e-12)
    assert bound['potential_difference_Mpc']<=bound['quadratic_potential_bound_Mpc']*(1+2e-12)
    # At this large radius the bound is close, while retaining its defined scope.
    assert actual>.99999*bound['potential_difference_Mpc']


def test_distributed_section_satisfies_boundary_and_integrated_force_rows():
    report,g,f=inputs()
    theta=report['conventional']['images_rad'][0][0]
    ray=distributed_ray(g,f,f.rV_Mpc,True,theta)
    assert ray['max_scaled_boundary_residual']<1e-12
    assert abs(ray['integrated_momentum_residual_rad'])<1e-15
    assert ray['max_collocation_relative_residual']<=2.1e-11


def test_owned_force_reproduces_saved_finite_body_scalar_ratios():
    report,_,f=inputs()
    for row in report['radial_profile']:
        r=row['radius_kpc']/1000
        ratio=float(f.extra_force(r)/f.force(r,False))
        np.testing.assert_allclose(ratio,row['local_force_extra_fraction'],rtol=3e-15)
