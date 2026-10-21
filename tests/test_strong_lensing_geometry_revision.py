"""Focused equation, source and geometry checks for the nonlinear rear plane."""
import sys
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parents[1]/'code'))

import strong_lensing_geometry_revision as revision
import strong_lensing_macro_revision as old
from strong_lensing_comparison import ARCSEC
from strong_lensing_local_family import G_MPC_KMS2_MSUN
from strong_lensing_comparison import C_KM_S


P=np.array([11.689817182237041,.7541853768397996,-.007952211194155661,
            .04112658412493257,.25625719126715585,.12482978852293553,
            .06540984613892988,-.0011291968484621834,.0010014758182217533])
G1=np.array([-2.585,-3.637])


def test_zero_mass_reproduces_existing_exact_images_and_relative_delays():
    before=old.macro_application(P,compact_mass_Msun=1e8)
    after=revision.macro_application(P,G1,perturber_mass_Msun=0.,rear_redshift=.8)
    np.testing.assert_allclose(after['images_light_centered_arcsec'],before['images_light_centered_arcsec'],atol=2e-11,rtol=0.)
    np.testing.assert_allclose(after['delays_days'],before['delays_days'],atol=2e-10,rtol=0.)


def test_perturber_center_is_its_foreground_deflected_observed_ray():
    s=revision.components(P,G1,perturber_mass_Msun=3e11)
    d=s['distances']; alpha=s['lens'].application(d['ol']*G1*ARCSEC)[1]
    expected=d['op']*G1*ARCSEC-d['lp']*alpha
    np.testing.assert_allclose(s['perturber_center_Mpc'],expected,atol=1e-17,rtol=0.)
    assert np.linalg.norm(expected-d['op']*G1*ARCSEC)>1e-5


def test_two_plane_ordered_jacobian_matches_finite_difference():
    s=revision.components(P,G1,perturber_mass_Msun=8e10)
    theta=np.array([1.2,.55])*ARCSEC
    step=1e-10
    finite=np.column_stack([(revision.ray_application(theta+step*np.eye(2)[i],s)['residual']-
                             revision.ray_application(theta-step*np.eye(2)[i],s)['residual'])/(2*step)
                            for i in range(2)])
    app=revision.ray_application(theta,s)
    np.testing.assert_allclose(app['jacobian'],finite,rtol=2e-8,atol=3e-9)
    # A nonlinear rear tide and foreground Hessian need not commute.
    assert abs(app['jacobian'][0,1]-app['jacobian'][1,0])>1e-5


def test_full_action_partial_gradients_match_independent_variations():
    s=revision.components(P,G1,perturber_mass_Msun=8e10)
    theta=np.array([1.2,.55])*ARCSEC
    rear=np.array([1.1,.48])*ARCSEC
    app=revision.arrival_action(theta,rear,s); step=2e-10
    for variable,name in [(0,'gradient_theta'),(1,'gradient_theta_rear')]:
        finite=[]
        for i in range(2):
            args=[theta.copy(),rear.copy()]
            args[variable][i]+=step
            plus=revision.arrival_action(*args,s)['ct_Mpc']
            args[variable][i]-=2*step
            minus=revision.arrival_action(*args,s)['ct_Mpc']
            finite.append((plus-minus)/(2*step))
        np.testing.assert_allclose(app[name],finite,atol=5e-11,rtol=2e-8)


def test_positive_projected_source_satisfies_original_poisson_trace():
    rear=revision.PositiveRearPlummer(2e11,.0002,(-.03,.02))
    for point in [(-.03,.02),(0.,0.),(.01,-.01)]:
        potential,alpha,hessian=rear.application(point)
        density=rear.surface_density_Msun_per_Mpc2(point)
        assert density>0
        np.testing.assert_allclose(np.trace(hessian),8*np.pi*G_MPC_KMS2_MSUN/C_KM_S**2*density,rtol=2e-10,atol=1e-18)


def test_measured_redshift_distance_identity_and_stationary_bright_rays():
    d=revision.distances(.7821)
    assert abs(d['distance_sum_relative_residual'])<1e-13
    app=revision.macro_application(P,G1,perturber_mass_Msun=8e10)
    assert app['max_lens_residual_rad']<1e-17
    assert app['max_arrival_stationarity_Mpc']<1e-13
    assert not app['physical_admission']


def test_zero_perturber_original_r1_projection_is_retained():
    before=old.macro_application(P,compact_mass_Msun=1e8,model='r1')
    after=revision.macro_application(P,G1,perturber_mass_Msun=0.,rear_redshift=.8,model='r1')
    np.testing.assert_allclose(after['images_light_centered_arcsec'],before['images_light_centered_arcsec'],atol=3e-10,rtol=0.)
    np.testing.assert_allclose(after['delays_days'],before['delays_days'],atol=3e-9,rtol=0.)


def test_exact_gnomonic_gaia_jacobian_including_uncertain_origin_matches_complex_step():
    ra=np.array([69.5623046,69.5618849,69.5616035,69.562038])
    dec=np.array([-12.2873141,-12.2871606,-12.2874817,-12.2877625])
    cov=np.tile(np.diag([.1**2,.08**2]),(4,1,1))
    result=revision.exact_gaia_relative_geometry(ra,dec,cov)
    sky=np.c_[ra,dec]*np.pi/180
    def independently_project(s):
        da=s[1:,0]-s[0,0];di=s[1:,1];d0=s[0,1]
        den=np.sin(d0)*np.sin(di)+np.cos(d0)*np.cos(di)*np.cos(da)
        east=np.cos(di)*np.sin(da)/den
        north=(np.cos(d0)*np.sin(di)-np.sin(d0)*np.cos(di)*np.cos(da))/den
        return np.c_[east,north].ravel()/ARCSEC
    np.testing.assert_allclose(result['relative_B_C_D_to_A_east_north_arcsec'],independently_project(sky),atol=1e-10,rtol=0.)
    cols=[]
    for source in range(4):
        for coordinate in range(2):
            h=1e-10
            z=sky.astype(complex)
            radians_per_mas=ARCSEC/1000
            if coordinate==0:radians_per_mas/=np.cos(sky[source,1])
            z[source,coordinate]+=1j*h*radians_per_mas
            cols.append(np.imag(independently_project(z))/h)
    np.testing.assert_allclose(result['jacobian_arcsec_per_tangent_mas'],np.array(cols).T,atol=2e-13,rtol=2e-10)
    assert np.linalg.eigvalsh(result['covariance_arcsec2']).min()>0
    assert abs(result['covariance_arcsec2'][0,2])>1e-9
