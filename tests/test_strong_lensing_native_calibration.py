from pathlib import Path
import sys
import numpy as np
from scipy import sparse
sys.path.insert(0,str(Path(__file__).resolve().parents[1]/'code'))
from strong_lensing_native_calibration import CalibrationOrderedHistory,FixedHarmonicContext


def test_harmonic_fixed_context_direction_is_linear_and_positive():
    context=FixedHarmonicContext((1,3),[True,False,True])
    np.testing.assert_allclose(context.positive([2.,-9.,8.]),[2.,5.,8.])
    np.testing.assert_allclose(context.action(np.array([[1.,2.],[999.,999.],[3.,6.]])),[[1.,2.],[2.,4.],[3.,6.]])


def model():
    shape=(1,3);n=3;times=np.array([[1.,2.],[6.,7.],[11.,12.],[16.,17.],[21.,22.],[26.,27.],[31.,32.],[36.,37.],[41.,42.]])
    g=sparse.csr_matrix([[1.01,-.01,0],[-.01,1.02,-.01],[0,-.01,1.01]])
    h=sparse.csr_matrix(np.linalg.inv(g.toarray()))
    adc=np.broadcast_to(np.array([3.,4.,5.]),(1,9,*shape)).copy()
    c=np.zeros((1,3,*shape));c[:,1]=1.;c[:,2]=.0002
    d={'gain_e_per_DN':np.ones((1,*shape))*2.,'linearity_coefficients':c,'pixeldq':np.zeros((1,*shape),np.uint32),
       'native_read_times_s':times,'frame_time_s':1.,'dark_RAPID42_DN':np.zeros((1,42,*shape)),
       'native_dark_RAPID_frame_indices':np.array([[5*j,5*j+1] for j in range(9)]),
       'native_subtracted_dark_group_DN':np.zeros((1,9,*shape)),
       'pre_linearity_group_DN':adc,'raw_group_DN':adc+100.,
       'saturation_threshold_DN':np.array([[[140.,145.,150.]]]),
       'native_corrected_zeroframe_DN':np.ones((1,*shape))}
    rate=np.array([8.,12.,10.])
    return CalibrationOrderedHistory(d,0,g,h,rate),rate


def test_full_native_direction_chain_includes_caps_and_implicit_intercept():
    obj,rate=model();kappa=1e-5
    base=obj.evaluate(rate,kappa,derivatives=True)
    assert base['physical']['clipped_read_mask'].any()
    directions=np.eye(3);dk=np.zeros(3)
    weights=np.zeros((9,3));weights[0]=-.1;weights[8]=.1
    wz=np.array([.01,0.,.02])
    actual=obj.project_directions(base,directions,dk,[0,1,2],weights,wz)
    def result(r,k):
        s=obj.evaluate(r,k)
        return np.sum(weights*s['groups_DN'].reshape(9,3),axis=0)+wz*s['zeroframe_DN'].ravel()
    eps=1e-4
    for i in range(3):
        e=np.eye(3)[i]*eps
        np.testing.assert_allclose(actual[:,i],(result(rate+e,kappa)-result(rate-e,kappa))/(2*eps),rtol=2e-6,atol=1e-8)
    actual_k=obj.project_directions(base,np.zeros((3,1)),[1.],[0,1,2],weights,wz)[:,0]
    np.testing.assert_allclose(actual_k,(result(rate,kappa+1e-7)-result(rate,kappa-1e-7))/(2e-7),rtol=2e-6,atol=1e-7)
