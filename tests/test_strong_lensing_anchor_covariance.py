import numpy as np
import pytest
from pathlib import Path
import sys
from scipy import sparse
sys.path.insert(0,str(Path(__file__).resolve().parents[1]/'code'))
from strong_lensing_anchor_covariance import electronic_anchor_rate_map, anchored_noise_contraction


def test_nonlocal_anchor_matches_nonlinear_implicit_boundary_finite_difference():
    # L^-1(q)=sqrt(q); spatial mixing does not commute with local derivatives.
    h = np.array([[.9,.1],[.2,.8]]); g = np.linalg.inv(h)
    beta = np.array([20.,30.]); signal = np.array([[[1.,2.],[2.,4.]],[[15.,11.],[21.,16.]]])
    q = beta[None,None]+signal
    inv = .5/np.sqrt(q)
    adc = np.einsum('pq,jq->jp',h,np.sqrt(q).mean(axis=1))
    ld = 2*adc; w = np.array([[-.7,-.5],[.7,.5]])
    _,k = electronic_anchor_rate_map(h,g,[0,1],w,ld,inv)
    target = np.sqrt(q[0]).mean(axis=0)
    observed0 = h@target
    def rate(a):
        b=beta.copy()
        for _ in range(8):
            b-=(np.sqrt(b[None]+signal[0]).mean(axis=0)-g@a)/(.5/np.sqrt(b[None]+signal[0])).mean(axis=0)
        native = np.einsum('pq,jq->jp',h,np.sqrt(b[None,None]+signal).mean(axis=1))**2
        return np.sum(w*native,axis=0)
    delta=1e-5;fd=np.column_stack([(rate(observed0+np.eye(2)[i]*delta)-rate(observed0-np.eye(2)[i]*delta))/(2*delta) for i in range(2)])
    np.testing.assert_allclose(k.toarray(),fd,rtol=2e-7,atol=1e-8)
    assert abs(k[0,1])>1e-5


def test_same_realization_cross_covariance_is_not_independent_anchor_addition():
    rng=np.random.default_rng(13);f=rng.normal(size=(5,8));c=f@f.T
    k=np.array([[.2,.1,-.3],[.4,-.2,.1]])
    answer,t=anchored_noise_contraction(c,k,2)
    np.testing.assert_allclose(answer,(t@f)@(t@f).T,rtol=1e-14,atol=1e-14)
    assert not np.allclose(answer,c[:2,:2]+k@c[2:,2:]@k.T)
    assert np.linalg.eigvalsh(answer).min()>0


def test_native_zero_weight_requires_its_actual_chain():
    with pytest.raises(ValueError,match='ZEROFRAME'):
        electronic_anchor_rate_map(sparse.eye(2),sparse.eye(2),[0,1],np.ones((2,2)),np.ones((2,2)),np.ones((2,2,2)),native_zero_weights=np.ones(2))
