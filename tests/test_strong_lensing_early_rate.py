from pathlib import Path
import sys
import numpy as np
from numpy.testing import assert_allclose
sys.path.insert(0,str(Path(__file__).resolve().parents[1]/'code'))
from strong_lensing_early_rate import first_two_group_product


def test_actual_two_group_projection_and_frame_poisson_contraction():
    times=np.array([[1.,2.],[6.,7.]])
    groups=np.array([[[3.]],[[13.]]]);dq=np.zeros_like(groups,dtype=np.uint8)
    out=first_two_group_product(groups,dq,times,5.,np.array([[.5]]),3.,np.array([[2.]]),np.array([[4.]]),np.array([[2.]]))
    assert_allclose(out['rate_SCI'],[[12.]])
    w=np.array([-.5,-.5,.5,.5])/5;tt=times.ravel()
    expected=(2/2*(w@np.minimum(tt[:,None],tt[None,:])@w)+16/(2*25))*36
    assert_allclose(out['variance_SCI2'],[[expected]])
    assert out['valid'].item()


def test_common_support_rejects_second_group_jump_or_unphysical_pilot():
    times=np.array([[1.,2.],[6.,7.]])
    groups=np.ones((2,1,3));dq=np.zeros_like(groups,dtype=np.uint8);dq[1,0,0]=4
    out=first_two_group_product(groups,dq,times,5.,np.ones((1,3)),1.,np.ones((1,3)),np.ones((1,3)),np.array([[1.,-1.,1.]]))
    assert_allclose(out['valid'],[[False,False,True]])


def test_physical_read_times_remain_distinct_from_rounded_native_group_time():
    times=np.array([[1.,2.],[6.,7.]])
    out=first_two_group_product(np.ones((2,1,1)),np.zeros((2,1,1),np.uint8),times,5.1,np.ones((1,1)),1.,np.ones((1,1)),np.zeros((1,1)),np.ones((1,1)))
    w=np.array([-.5,-.5,.5,.5])/5.1;t=times.ravel()
    assert_allclose(out['photon_variance_SCI2'],[[5.1/5*(w@np.minimum(t[:,None],t[None,:])@w)]])
