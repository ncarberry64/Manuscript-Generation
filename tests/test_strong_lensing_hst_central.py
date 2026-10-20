"""Scientific controls for the HST central-image inference."""
import sys
from pathlib import Path
import numpy as np
import pytest

sys.path.insert(0,str(Path(__file__).resolve().parents[1]/"code"))
from strong_lensing_hst_central import (point_template, fit_point,
    positivity_upper, nonnegative_profile, registered_light_center,
    assess_central_image_system, sky_offset_to_pixel)


def gaussian_psf():
    y,x=np.indices((21,21));p=np.exp(-((x-10)**2+(y-10)**2)/3)
    return p/p.sum(),(10,10)


def test_point_flux_and_centroid_recovered_with_local_gradient():
    psf,pc=gaussian_psf();xy=np.array([14.3,13.7]);p=point_template((31,31),xy,psf,pc)
    y,x=np.indices(p.shape);d=120*p+.03+.001*x-.002*y
    fit=fit_point(d,np.full(p.shape,.01),[14,14],psf,pc,5)
    np.testing.assert_allclose(fit["xy"],xy,atol=1e-7)
    assert abs(fit["flux"]-120)<1e-6


def test_registration_preserves_east_north_signs_and_rotation():
    sky=np.array([[1,.5],[-.3,1.1],[-1.3,0],[.2,-1.05]])
    pix=sky/.05;pix[:,0]*=-1
    t=.001;rotation=np.array([[np.cos(t),-np.sin(t)],[np.sin(t),np.cos(t)]])
    obs=pix@rotation.T+[80.2,79.7]
    center,rot,resid=registered_light_center(obs,sky,.05)
    np.testing.assert_allclose(center,[80.2,79.7],atol=1e-8)
    assert abs(rot-t)<1e-9 and resid<1e-8


def test_arbitrary_nonnegative_nucleus_can_hide_a_fainter_source():
    psf,pc=gaussian_psf();p=point_template((31,31),[15,15],psf,pc)
    d=10*p+.2
    chi,background,_=nonnegative_profile(d,np.full(p.shape,.01),p,4)
    assert chi<1e-20
    assert np.min(background)>=0
    assert nonnegative_profile(d,np.full(p.shape,.01),p,30)[0]>1000


def test_limit_respects_correlated_noise_and_nonnegative_galaxy():
    # Every pixel shares one random noise variate. Treating the pixels as
    # independent in an ordinary chi-square test would be incorrect here.
    p=np.array([.05,.2,.5,.2,.05]);s=np.full(5,.1);rng=np.random.default_rng(17)
    true=8.;covered=0
    for z in rng.normal(size=3000):
        d=true*p+z*s
        upper,_=positivity_upper(d,s,p)
        covered+=upper>=true
    assert covered/3000>.95


def test_arbitrary_positive_light_weakens_bound_monotonically():
    p=np.array([.1,.2,.4,.2,.1]);s=np.full(5,.03)
    d=3*p
    u,_=positivity_upper(d,s,p)
    v,_=positivity_upper(d+np.array([.5,.3,.2,.3,.5]),s,p)
    assert v>u


def test_bright_central_image_is_rejected_without_sersic_prior():
    p=np.array([.05,.2,.5,.2,.05]);s=np.full(5,.01)
    # Even assigning ALL measured central light to a point source cannot
    # accommodate the proposed flux. No galaxy slope or nucleus prior enters.
    d=2*p+.02
    upper,_=positivity_upper(d,s,p)
    assert upper<2.2
    assert upper<5


def test_overlapping_images_are_tested_as_a_sum():
    psf,pc=gaussian_psf();p=point_template((31,31),[15,15],psf,pc)
    d=3*p+.005;s=np.full(d.shape,.002)
    single,_=assess_central_image_system(d,s,d.shape,[15,15],.05,[[0,0]],[.02],100,psf,pc)
    pair,total=assess_central_image_system(d,s,d.shape,[15,15],.05,[[0,0],[0,0]],[.02,.02],100,psf,pc)
    assert single["upper_scale_95"]>1
    assert pair["upper_scale_95"]<1
    np.testing.assert_allclose(total,4*p,atol=1e-14)


def test_empty_support_and_invalid_noise_are_rejected():
    with pytest.raises(ValueError,match="positive template support"):
        positivity_upper(np.ones(3),np.ones(3),np.zeros(3))
    with pytest.raises(ValueError,match="positive template support"):
        positivity_upper(np.ones(3),np.zeros(3),np.ones(3))


def test_central_template_applies_fitted_registration_rotation():
    psf,pc=gaussian_psf();shape=(41,41);center=np.array([20.,20.])
    offset=[.1,.05];rot=.2
    expected=center+sky_offset_to_pixel(offset,.05,rot)
    p=point_template(shape,expected,psf,pc)
    _,model=assess_central_image_system(4*p+.1,np.full(shape,.01),shape,
        center,.05,[offset],[.04],100,psf,pc,rotation_rad=rot)
    np.testing.assert_allclose(model,4*p,atol=1e-14)
    unrotated=point_template(shape,center+sky_offset_to_pixel(offset,.05),psf,pc)
    assert np.max(np.abs(model-4*unrotated))>.01
