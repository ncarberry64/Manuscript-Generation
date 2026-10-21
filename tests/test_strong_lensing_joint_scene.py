"""Synthetic operator, covariance and recovery checks; no native target reads."""
from pathlib import Path
import sys
import numpy as np
from numpy.testing import assert_allclose
from scipy.interpolate import RectBivariateSpline

sys.path.insert(0,str(Path(__file__).resolve().parents[1]/'code'))
import strong_lensing_joint_scene as js


def response(eta=(.1,-.07),size=25,over=4):
    z=(np.arange(size)-(size-1)/2)/over
    x,y=np.meshgrid(z,z)
    base=np.exp(-(x*x+y*y)/(.8**2*2))+.001
    return js.PositiveResponse(base,np.array([[-.6,.1],[.6,-.2]]),np.eye(2),.65,np.array(eta),over)


def test_delivered_positive_mean_and_analytic_derivatives():
    r=response();x=np.array([-.73,.04,.59]);y=np.array([.29,-.11,.21])
    raw=np.exp(-((x[:,None]-r.nodes[:,0])**2+(y[:,None]-r.nodes[:,1])**2)/(2*r.width**2))
    base=RectBivariateSpline(r.coord,r.coord,r.base.T).ev(x,y)
    assert_allclose(r.values(x,y),base*np.exp(raw@r.eta)/r.normalization,rtol=2e-15)
    p,dx,dy=r.values(x,y,True);h=1e-6
    assert_allclose(dx,(r.values(x+h,y)-r.values(x-h,y))/(2*h),rtol=1e-8,atol=1e-10)
    assert_allclose(dy,(r.values(x,y+h)-r.values(x,y-h))/(2*h),rtol=1e-8,atol=1e-10)
    derivative=r.cal_values(x,y);moment_derivative=r.moment_derivative.copy();eta=r.eta.copy()
    for k in range(2):
        e=np.eye(2)[k]*h
        r.update(eta+e);plus=r.values(x,y);mp=r.moment.copy()
        r.update(eta-e);minus=r.values(x,y);mm=r.moment.copy()
        assert_allclose(derivative[:,k],(plus-minus)/(2*h),rtol=2e-8,atol=1e-10)
        assert_allclose(moment_derivative[:,k],(mp-mm)/(2*h),rtol=2e-8,atol=1e-10)
    r.update(eta)
    assert_allclose(r.grid.sum()/r.over**2,1.,atol=2e-15)
    assert_allclose(r.derivative_grid.sum(axis=(1,2)),0.,atol=1e-14)
    assert_allclose(r.values(np.array([r.coord[-1]+.1]),np.array([0.])),0.)


def test_scaled_svd_preserves_unmeasured_null_directions():
    j=np.array([[1.,2.,0.]])
    fit=js.svd_fit(j,scales=np.array([2.,.5,3.]))
    assert fit['rank']==1 and not fit['full_column_rank']
    assert fit['null_vectors_dimensionless'].shape==(2,3)
    assert_allclose((j*fit['column_scales'])@fit['null_vectors_dimensionless'].T,0.,atol=1e-14)
    assert_allclose(j@fit['covariance']@j.T,[[1.]],atol=1e-14)


def test_cut_calibration_covariance_retains_shared_correlations():
    # Two independent measured centers share one uncertain response offset.
    # Their conditional errors are unit and the common calibration variance=4.
    a=np.eye(2);b=np.ones((2,1));v=np.array([[4.]])
    conditional,cut,sensitivity,fit=js.propagated_cut_covariance(a,b,v)
    assert_allclose(conditional,np.eye(2),atol=1e-14)
    assert_allclose(sensitivity,-b,atol=1e-14)
    assert_allclose(cut,[[5.,4.],[4.,5.]],atol=1e-14)
    # Profiling an additive common nuisance removes exactly its one direction.
    prof,rr=js.profile_columns(b,a,np.array([2.,3.]))
    assert_allclose(b.T@prof,0.,atol=1e-14)
    assert_allclose(b.T@rr,0.,atol=1e-14)


def test_scene_convolution_flux_and_all_calibration_batches():
    r=response();z=np.arange(41)/r.over
    yy,xx=np.meshgrid(z,z,indexing='ij')
    # A constant sky must remain constant well inside the support. The ePSF is
    # already integrated over a native pixel; exactly one 1/over^2 factor.
    rec={'response':r,'sky_fine':np.ones((9,41,41)),
         'fine_ix':np.array([20]),'fine_iy':np.array([20])}
    js.refresh_scene_projection(rec)
    assert_allclose(rec['scene_matrix'],np.ones((1,9)),atol=2e-15)
    # >24 modes tests the FFT batching regression (mode='same' first-input
    # shape must not silently discard all but the first calibration mode).
    r.derivative_grid=np.broadcast_to(r.grid,(53,*r.grid.shape)).copy()
    columns=js.scene_calibration_columns(rec,np.ones(9))
    assert columns.shape==(1,53)
    assert_allclose(columns,9.,atol=2e-14)


def test_synthetic_point_recovery_with_fixed_coherent_sky():
    r=response(size=41);yy,xx=np.mgrid[-4:5,-4:5]
    x=xx.ravel().astype(float);y=yy.ravel().astype(float)
    truth=np.array([.12,-.17,2.4]);flux=np.exp(truth[2])
    # Fit only within finite support and keep background columns independent.
    use=(abs(x)<r.coord[-1]-.5)&(abs(y)<r.coord[-1]-.5);x=x[use];y=y[use]
    bg=np.c_[np.ones(len(x)),x/10,y/10]
    data=flux*r.values(x-truth[0],y-truth[1])+bg@np.array([.2,.03,-.02])
    rec={'response':r,'x':x,'y':y,'error':np.ones(len(x))*.01,'data':data,
         'detector_background':bg,'scene_matrix':np.zeros((len(x),9)),
         'guess':np.array([.0,.0]),'flux_guess':flux*.9,'patch':0,'exposure':'synthetic','name':'A'}
    problem=js.SceneProblem([rec],scene_fixed=np.zeros(36))
    fit=problem.fit()
    assert fit.success
    assert_allclose(fit.x[:3],truth,atol=1e-8)
    assert_allclose(problem.evaluate(fit.x)[0],0.,atol=1e-8)


def test_actual_history_required_and_charge_exchange_conserves_flux():
    q=np.arange(20.).reshape(4,5)
    assert_allclose(js.conservative_charge_laplacian(q).sum(),0.,atol=1e-14)
    reads=np.stack([q,2*q]);weights=np.stack([-np.ones_like(q),np.ones_like(q)])
    assert_allclose(js.history_rate(reads,weights),q)
    import pytest
    with pytest.raises(ValueError):js.history_rate(reads,np.ones((2,1,1)))
