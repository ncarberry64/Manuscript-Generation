"""Independent propagation/action and injected-observable checks."""
import sys
from pathlib import Path

import numpy as np
import pytest

sys.path.insert(0,str(Path(__file__).resolve().parents[1]/'code'))

import strong_lensing_observed_macro as macro
from strong_lensing_comparison import ARCSEC, lens_application
from strong_lensing_robustness import ProfileR1Lens


P=np.array([11.72,.74,-.01,.052,.09,.173,.091,0.,0.])


def test_remote_plane_forward_map_matches_explicit_two_plane_ray():
    d=macro.frozen_distances(); geometry,U=macro.screen_geometry(.07,.12,-.04)
    assert abs(d['distance_sum_relative_residual'])<1e-14
    theta=np.array([1.1,-.7])*ARCSEC
    lens=macro.body_from_parameters(P)
    alpha=lens.application(d['ol']*theta)[1]
    rear=d['op']*theta-d['lp']*alpha
    explicit=d['os']*theta-d['ls']*alpha-d['ps']*U@rear
    compact=geometry.D_os@theta-geometry.D_ls@alpha
    np.testing.assert_allclose(explicit,compact,rtol=3e-15,atol=1e-18)
    np.testing.assert_allclose(geometry.T,geometry.T.T,rtol=0,atol=1e-12)


def test_schur_arrival_action_matches_independent_full_action_differences():
    d=macro.frozen_distances(); geometry,U=macro.screen_geometry(.08,.09,.045)
    source=d['os']*np.array([.03,-.04])*ARCSEC
    lens=macro.body_from_parameters(P)
    points=np.array([[1.1,.6],[-1.3,-.2],[.2,-1.]])*ARCSEC
    full=[]; reduced=[]
    for theta in points:
        value,rear,residual=macro.two_plane_arrival_length(theta,source,lens,U)
        assert max(abs(residual))<2e-17
        full.append(value)
        reduced.append(lens_application(theta,geometry,lens,source)['ct_Mpc'])
    np.testing.assert_allclose(np.array(full)[1:]-full[0],np.array(reduced)[1:]-reduced[0],
                               rtol=2e-13,atol=1e-22)


def test_screen_rejects_conjugate_branch_instead_of_using_pseudoinverse():
    with pytest.raises(ValueError,match='pre-conjugate'):
        macro.screen_geometry(.7,.4,0.)


def test_frozen_projection_is_one_consistent_potential_gradient_hessian():
    screened=ProfileR1Lens(macro.body_from_parameters(P),macro.frozen_coefficients())
    frozen=macro.FrozenExtraProjection(screened)
    x=np.asarray(screened.centre_Mpc)+np.array([.006,.003])
    potential,alpha,hessian=frozen.application(x,screened.centre_Mpc)
    step=1e-7
    finite_gradient=[]; finite_hessian=[]
    for axis in np.eye(2):
        plus=frozen.application(x+step*axis,screened.centre_Mpc)
        minus=frozen.application(x-step*axis,screened.centre_Mpc)
        finite_gradient.append((plus[0]-minus[0])/(2*step))
        finite_hessian.append((plus[1]-minus[1])/(2*step))
    np.testing.assert_allclose(alpha,finite_gradient,rtol=3e-10,atol=1e-19)
    np.testing.assert_allclose(hessian,np.array(finite_hessian).T,rtol=3e-9,atol=1e-16)
    exact=screened.extra_application(x)
    np.testing.assert_allclose(alpha,exact[1],rtol=3e-7,atol=1e-18)


def test_injected_observables_recovered_without_claiming_unique_parameters():
    actual=macro.macro_application(P)
    observed=dict(observed_relative_arcsec=actual['relative_astrometry_arcsec'],
                  sigma_arcsec=np.full(8,.002),
                  observed_images_light_centered_arcsec=actual['images_light_centered_arcsec'],
                  assumptions='synthetic independent injection for application test')
    initial=P+np.array([.002,.001,.0001,-.0001,0.,.0002,-.0002,0.,0.])
    fit=macro.fit_macro(actual['delays_days'],np.diag([.05,.05,.05])**2,
                        astrometry=observed,initial=initial,max_nfev=300)
    assert fit['chi2_total']<1e-9
    np.testing.assert_allclose(fit['predicted_relative_astrometry_arcsec'],observed['observed_relative_arcsec'],atol=2e-8,rtol=0)
    np.testing.assert_allclose(fit['predicted_delays_days'],actual['delays_days'],atol=2e-7,rtol=0)
    assert fit['max_lens_residual_rad']<1e-17
    assert max(abs(np.array(fit['two_plane_minus_effective_delay_days'])))<1e-10


def test_invalid_delay_covariance_cannot_be_silently_regularized():
    with pytest.raises(np.linalg.LinAlgError):
        macro.fit_macro([1.,2.,3.],np.diag([1.,0.,1.]))


def test_exact_training_callback_is_consumed_and_reported_at_original_application():
    actual=macro.macro_application(P)
    observed=dict(observed_relative_arcsec=actual['relative_astrometry_arcsec'],
                  sigma_arcsec=np.full(8,.002),
                  observed_images_light_centered_arcsec=actual['images_light_centered_arcsec'],
                  assumptions='synthetic callback application')
    called=[]
    def likelihood(delays):
        called.append(np.asarray(delays).copy())
        return 137.-.5*np.sum(((delays-actual['delays_days'])/.1)**2)
    fit=macro.fit_macro(actual['delays_days'],np.eye(3)*.1**2,astrometry=observed,
                        initial=P,delay_log_likelihood=likelihood,exact_maxiter=3)
    assert len(called)>18
    expected=likelihood(np.asarray(fit['predicted_delays_days']))
    assert fit['exact_training_log_likelihood']==expected
    assert fit['exact_training_refinement']
