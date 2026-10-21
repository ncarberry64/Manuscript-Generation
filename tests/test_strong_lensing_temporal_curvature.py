from pathlib import Path
import sys
import numpy as np
from numpy.testing import assert_allclose
sys.path.insert(0,str(Path(__file__).resolve().parents[1]/'scripts'))
from diagnose_he0435_temporal_curvature import fit_curves


def test_actual_group_average_design_recovers_curvature_with_correlated_poisson_metric():
    times=np.array([[1+5*i,2+5*i] for i in range(9)],float)
    beta,rate,curvature=3.,2.,.002
    data=(beta+rate*times.mean(axis=1)+curvature*(times**2).mean(axis=1))[:,None,None]
    valid=np.ones_like(data,bool)
    p,cov,q,ql=fit_curves(data,valid,np.array([[rate]]),np.array([[2.]]),np.array([[4.]]),times)
    assert_allclose(p[:,0,0],[beta,rate,curvature],atol=1e-10)
    assert q[0,0]<1e-20 and ql[0,0]>0
    assert np.linalg.eigvalsh(cov[:,:,0,0]).min()>0


def test_two_groups_do_not_define_a_free_intercept_rate_curvature_triple():
    times=np.array([[1+5*i,2+5*i] for i in range(9)],float)
    data=np.ones((9,1,1));valid=np.zeros_like(data,bool);valid[:2]=True
    p,cov,q,ql=fit_curves(data,valid,np.ones((1,1)),np.ones((1,1)),np.ones((1,1)),times)
    assert np.isnan(p).all() and np.isnan(cov).all() and np.isnan(q).all()
