import sys
from pathlib import Path
import numpy as np
from numpy.testing import assert_allclose
from scipy.interpolate import RectBivariateSpline
sys.path.insert(0,str(Path(__file__).resolve().parents[1]/"code"))
from strong_lensing_read_history import (zero_flux_laplacian,group_time_moments,
    temporal_covariance,conditional_innovation,selected_history,fit_history)


def test_full_graph_charge_conservation_and_no_wrap():
    lap=zero_flux_laplacian((5,6))
    f=np.zeros(30);f[0]=4
    out=lap@f**2
    assert_allclose(out.sum(),0,atol=1e-14)
    assert out[5]==0 and out[6]==16 and out[1]==16 and out[0]==-32


def test_group_second_moment_is_not_square_of_group_time():
    m,v=group_time_moments(np.array([[1.,2.],[6.,7.]]))
    assert_allclose(v,m*m+.25)
    assert v[0]==2.5


def test_jump_starts_segment_but_is_retained_and_saturation_removed():
    idx,seg=selected_history(np.array([0,4,0,1,0,2,2],np.uint8))
    assert_allclose(idx,[0,1,2,4]);assert_allclose(seg,[0,1,1,2])


def test_shared_reads_change_innovation_and_covariance():
    t=np.array([[1.,2.],[6.,7.],[11.,12.]])
    c=temporal_covariance(t,10.,2.,4.)
    assert_allclose(c[0,0],10/2*1.25+4.)
    assert c[0,2]>0
    innovation,var=conditional_innovation(np.array([2.,3.]),np.array([5.]),np.zeros(2),np.zeros(1),c[:2,:2],c[2:,:2],c[2:,2:])
    assert abs(innovation[0]-5)>1
    assert var[0,0]<c[2,2]


def test_positive_kappa_recovery_and_held_group_predictive_fit():
    # Nonuniform compact source plus fixed padded boundary: exact nonlinear mean.
    rng=np.random.default_rng(1703)
    ny=nx=7
    yy,xx=np.indices((ny,nx));f=.2+30*np.exp(-((xx-3.1)**2+(yy-2.8)**2)/1.5)
    t=np.array([[5*j+1.,5*j+2.] for j in range(9)])
    m,v=group_time_moments(t)
    lap=zero_flux_laplacian(f.shape)
    k=8e-5
    mean=50+m[:,None,None]*f+.5*k*v[:,None,None]*(lap@(f.ravel()**2)).reshape(f.shape)
    gain=np.full(f.shape,1e8);rn=np.full(f.shape,.04)
    g=mean+rng.normal(0,.02,mean.shape)
    dq=np.zeros(g.shape,np.uint8);pdq=np.zeros(f.shape,np.uint32)
    result,arrays=fit_history(g,dq,pdq,gain,rn,t,roi=(1,6,1,6),max_nfev=150)
    assert result["optimizer_success"]
    assert abs(result["kappa"]-k)<1.5e-5
    assert result["held_raw_RMS_DN"]<.1
    assert result["positive_source_mean_verified"]
    assert not result["physical_calibration_pass"]
    assert result["profile_parameter_rank"]==26
    changed=g.copy();changed[8]+=100
    result2,arrays2=fit_history(changed,dq,pdq,gain,rn,t,roi=(1,6,1,6),max_nfev=150)
    # Deliberately corrupt the unused group: parameters, frozen noise and forecasts
    # must stay identical while its measured prediction error changes by 100 DN.
    assert result2["kappa"]==result["kappa"]
    assert_allclose(arrays2["parameter_covariance"],arrays["parameter_covariance"],rtol=0,atol=0)
    assert_allclose(arrays2["held_records"][:,2],arrays["held_records"][:,2],rtol=0,atol=0)
    assert_allclose(arrays2["held_records"][:,3]-arrays["held_records"][:,3],100,atol=1e-12)


def test_optical_profile_and_history_joint_recovery_with_full_covariance():
    rng=np.random.default_rng(942)
    shape=(7,7);yy,xx=np.indices(shape)
    os=4;gxy=(np.arange(65)-32)/os
    gy,gx=np.meshgrid(gxy,gxy,indexing="ij")
    grid=np.exp(-(gx*gx+gy*gy)/(2*.85**2))/(2*np.pi*.85**2)
    interp=RectBivariateSpline(np.arange(65),np.arange(65),grid.T)
    center=np.array([3.13,2.81]);flux=300.;background=np.array([.4,.001,-.002])
    psf=interp.ev(32+os*(xx-center[0]),32+os*(yy-center[1]))
    f=(flux*psf+background[0]+background[1]*(xx-center[0])+background[2]*(yy-center[1]))/2
    times=np.array([[5*j+1.,5*j+2.] for j in range(9)])
    m,v=group_time_moments(times);k=8e-5
    lap=zero_flux_laplacian(shape)
    source=20+m[:,None,None]*f+.5*k*v[:,None,None]*(lap@(f.ravel()**2)).reshape(shape)
    values=source+rng.normal(0,.005,source.shape)
    optical={"psf_grid":grid,"oversampling":os,"flat":np.ones(shape),"photmjsr":2.,
             "initial_center_xy":center+[.12,-.15],"initial_flux_SCI":270.,
             "initial_background_SCI":[.35,0.,0.]}
    result,arrays=fit_history(values,np.zeros(values.shape,np.uint8),np.zeros(shape,np.uint32),
         np.full(shape,1e10),np.full(shape,.01),times,roi=(0,7,0,7),optical=optical,max_nfev=100)
    assert result["optimizer_success"]
    assert_allclose(result["optical_fit"]["centroid_local_xy_pixel"],center,atol=5e-4)
    assert abs(result["kappa"]-k)<3e-6
    assert result["profile_parameter_rank"]==7
    assert arrays["parameter_covariance"].shape==(7,7)
    assert abs(arrays["parameter_covariance"][0,-1])>0
    assert result["positive_individual_read_source_verified"]
    assert result["held_raw_RMS_DN"]<.02
    assert not result["optical_fit"]["physical_optical_response_validation_pass"]
