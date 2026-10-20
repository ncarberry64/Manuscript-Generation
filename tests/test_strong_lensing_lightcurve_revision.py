"""Independent covariance and prediction checks for the revised classical model."""
from dataclasses import replace
import sys
from pathlib import Path
import numpy as np
from numpy.testing import assert_allclose
from scipy.linalg import cho_factor, cho_solve
from scipy.special import logsumexp

sys.path.insert(0,str(Path(__file__).resolve().parents[1]/'code'))
from strong_lensing_lightcurve_revision import (
    RevisionCurves,RevisionParameters,concatenate,dense_covariance,
    log_likelihood,predict,mixture_predict,score_mixture,encode,decode,
)


def fixture():
    data=RevisionCurves(np.array([100.,103.,108.,113.,120.,130.,135.,140.]),
        np.array([0,1,0,1,0,1,0,1]),np.array([1.1,1.3,1.15,1.35,1.2,1.4,1.25,1.42]),
        np.full(8,.02),np.array([0,0,0,0,1,1,1,1]),2,('old','new'))
    p=RevisionParameters((0.,2.),extrinsic_sigma_mag=(.1,.15),
        jitter_mag=(.01,.015),brownian_mag_per_sqrt_year=(.02,.04),
        reference_time_days=90.,magnitude_zero=1.,offset_sigma_mag=2.)
    return data,p


def test_state_likelihood_matches_independent_dense_covariance_and_accelerator():
    data,p=fixture(); cov=dense_covariance(data,p)
    factor=cho_factor(cov,lower=True);r=data.magnitude-p.magnitude_zero
    dense=-.5*(len(r)*np.log(2*np.pi)+2*np.log(np.diag(factor[0])).sum()+r@cho_solve(factor,r))
    assert_allclose(log_likelihood(data,p,backend='python'),dense,atol=2e-9)
    assert_allclose(log_likelihood(data,p),dense,atol=2e-9)


def test_smoother_matches_dense_conditional_mean_and_variance():
    data,p=fixture();train=data.subset(np.arange(8)<5);test=data.subset(np.arange(8)>=5)
    pred=predict(train,test,p);cov=dense_covariance(data,p)
    cross=cov[5:,:5]; invr=np.linalg.solve(cov[:5,:5],train.magnitude-p.magnitude_zero)
    expected_mean=p.magnitude_zero+cross@invr
    expected_cov=cov[5:,5:]-cross@np.linalg.solve(cov[:5,:5],cross.T)
    assert_allclose(pred['mean'],expected_mean,atol=2e-8)
    assert_allclose(pred['variance'],np.diag(expected_cov),atol=2e-8)


def test_prediction_does_not_consume_target_values():
    data,p=fixture();train=data.subset(np.arange(8)<5);test=data.subset(np.arange(8)>=5)
    a=mixture_predict(train,test,[p]);b=mixture_predict(train,replace(test,magnitude=test.magnitude+1e5),[p])
    assert_allclose(a['mean'],b['mean'],rtol=0,atol=0)
    assert_allclose(a['variance'],b['variance'],rtol=0,atol=0)


def test_driver_updates_component_weights_by_conditional_bayes_rule():
    data,p=fixture();train=data.subset(np.arange(8)<4);driver=data.subset(np.arange(8)==4);test=data.subset(np.arange(8)>4)
    ps=[p,replace(p,source_sigma_mag=.5)]
    result=mixture_predict(train,test,ps,driver=driver)
    logw=np.array([log_likelihood(concatenate(train,driver),x)-log_likelihood(train,x) for x in ps])
    expected=np.exp(logw-logsumexp(logw))
    assert_allclose(result['weights'],expected,atol=2e-10)


def test_calibration_joint_density_has_the_required_normalization_jacobian():
    data,p=fixture();train=data.subset(np.arange(8)<5);test=data.subset(np.arange(8)>=5)
    scale=np.array([1.3,1.8]);pred=mixture_predict(train,test,[p])
    transformed=replace(test,magnitude=pred['mean']+(test.magnitude-pred['mean'])/scale[test.image])
    raw=score_mixture(train,transformed,[p]);cal=score_mixture(train,test,[p],calibration=scale)
    assert_allclose(cal['joint_conditional_log_density'],raw['joint_conditional_log_density']-
                    np.log(scale[test.image]).sum(),atol=2e-9)
    assert_allclose(cal['variance'],pred['variance']*scale[test.image]**2,atol=2e-12)


def test_brownian_and_random_slope_create_growing_future_variance():
    data,p=fixture();train=data.subset(np.arange(8)<6)
    times=np.array([141.,1000.,5000.]);future=RevisionCurves(times,np.zeros(3,int),np.zeros(3),
        np.full(3,.02),np.ones(3,int),2,data.instrument_names)
    variances=predict(train,future,p)['variance']
    assert variances[2]>variances[1]>variances[0]


def test_time_origin_mixture_offsets_are_applied_to_cross_instrument_covariance():
    data,p=fixture();shift=replace(p,instrument_time_offset_days=(.5,0.))
    shifted_data=replace(data,time_days=data.time_days+np.where(data.instrument==0,.5,0.))
    assert_allclose(dense_covariance(data,shift),dense_covariance(shifted_data,p),atol=1e-12)
    assert_allclose(log_likelihood(data,shift),log_likelihood(shifted_data,p),atol=1e-10)


def test_parameter_coordinate_round_trip():
    _,p=fixture();roundtrip=decode(encode(p),p)
    for key in ('source_sigma_mag','source_tau_days','extrinsic_tau_days'):
        assert_allclose(getattr(p,key),getattr(roundtrip,key))
    for key in ('delays_days','extrinsic_sigma_mag','jitter_mag','brownian_mag_per_sqrt_year'):
        assert_allclose(getattr(p,key),getattr(roundtrip,key))


def test_interval_score_penalizes_misses_and_includes_predictive_width():
    data,p=fixture();train=data.subset(np.arange(8)<5);test=data.subset(np.arange(8)>=5)
    prediction=mixture_predict(train,test,[p]);centered=replace(test,magnitude=prediction['mean'])
    score=score_mixture(train,centered,[p]);summary=score['overall']
    assert_allclose(summary['mean_95_interval_score_mag'],summary['mean_95_interval_width_mag'])
    outside=replace(test,magnitude=prediction['quantiles']['0.975']+1.)
    miss=score_mixture(train,outside,[p])['overall']
    assert_allclose(miss['mean_95_interval_score_mag'],miss['mean_95_interval_width_mag']+40.)
