from pathlib import Path
import sys
import numpy as np
from scipy import sparse
sys.path.insert(0,str(Path(__file__).resolve().parents[1]/'code'))
from strong_lensing_calibration_order import ipc_averaged_read_observation,first_group_preipc_intercept


def test_ipc_precedes_native_linearity_and_is_not_a_second_mean_convolution():
    h=sparse.csr_matrix([[.99,.01],[.02,.98]])
    c=np.array([[0.,0.],[1.,1.],[.002,.003]])
    adc=np.array([[[3.,9.],[5.,12.]]])
    q=adc+c[2]*adc**2
    got=ipc_averaged_read_observation(q,c,adc,h)
    x=(h@adc.mean(axis=1).T).T
    np.testing.assert_allclose(got['mean_DN'],x+c[2]*x*x,rtol=1e-13)
    wrong=(h@(adc+c[2]*adc**2).mean(axis=1).T).T
    assert np.max(abs(got['mean_DN']-wrong))>.001


def test_full_interpixel_read_jacobian_against_finite_difference():
    h=sparse.csr_matrix([[.97,.03],[.02,.98]])
    c=np.array([[0.,0.],[1.,1.],[.002,.003]])
    q=np.array([[[4.,11.],[6.,14.]],[[8.,17.],[10.,20.]]])
    base=ipc_averaged_read_observation(q,c,q,h,derivatives=True)
    step=1e-4
    for j in range(2):
        for r in range(2):
            for p in range(2):
                plus=q.copy();minus=q.copy();plus[j,r,p]+=step;minus[j,r,p]-=step
                a=ipc_averaged_read_observation(plus,c,q,h)['mean_DN'][j]
                b=ipc_averaged_read_observation(minus,c,q,h)['mean_DN'][j]
                np.testing.assert_allclose(base['read_to_native_group_jacobian'][j][r].toarray()[:,p],(a-b)/(2*step),rtol=1e-8,atol=1e-10)


def test_implicit_first_group_intercept_derivative_is_not_frozen():
    c=np.array([[0.,0.],[1.,1.],[.002,.003]])
    signal=np.array([[4.,11.],[6.,14.]])
    target=np.array([8.,13.]);initial=np.array([2.,3.])
    base=first_group_preipc_intercept(signal,c,target,initial)
    step=1e-4
    for r in range(2):
        for p in range(2):
            plus=signal.copy();minus=signal.copy();plus[r,p]+=step;minus[r,p]-=step
            a=first_group_preipc_intercept(plus,c,target,initial)['beta_preipc_DN']
            b=first_group_preipc_intercept(minus,c,target,initial)['beta_preipc_DN']
            want=np.zeros(2);want[p]=base['beta_signal_derivative'][r,p]
            np.testing.assert_allclose((a-b)/(2*step),want,rtol=1e-8,atol=1e-10)
