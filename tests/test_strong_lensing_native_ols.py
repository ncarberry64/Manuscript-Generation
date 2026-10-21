from pathlib import Path
import sys
import numpy as np
sys.path.insert(0,str(Path(__file__).resolve().parents[1]/"code"))
from strong_lensing_native_ols import piecewise_ols_pixel,piecewise_ols_arrays

OPTIONS=dict(gain=2.05,readnoise_cds_dn=12.,group_time=53.684,frame_time=10.73677)


def test_exact_two_group_contrast_and_zero_intercept_response():
    q=piecewise_ols_pixel([20.,140.],[0,0],**OPTIONS)
    expected=np.array([-1.,1.,0.])/float(np.float32(OPTIONS['group_time']))
    assert np.allclose(q['jacobian_group_and_zeroframe'],expected,rtol=0,atol=1e-17)
    assert np.isclose(q['slope_double'],120/float(np.float32(OPTIONS['group_time'])))


def test_multisegment_median_variance_derivative_matches_continuous_difference():
    y=np.array([500.,550.,612.,680.,1200.,1300.,1370.,1454.,1554.])
    dq=np.array([0,0,0,0,4,0,0,0,0])
    r=piecewise_ols_pixel(y,dq,**OPTIONS)
    assert len(r['segments'])==2
    assert abs(r['branch']['median_weight_derivative_coefficient'])>1e-5
    actual=[]
    for j in range(9):
        plus=y.copy();minus=y.copy();plus[j]+=.001;minus[j]-=.001
        actual.append((piecewise_ols_pixel(plus,dq,**OPTIONS)['slope_double']-
                       piecewise_ols_pixel(minus,dq,**OPTIONS)['slope_double'])/.002)
    assert np.allclose(actual,r['jacobian_group_and_zeroframe'][:9],rtol=1e-7,atol=2e-12)
    # A frozen segment-combination weight derivative misses a real dependency.
    assert r['branch']['locally_unique_derivative']


def test_native_saturation_single_group_and_zeroframe_branches():
    first=piecewise_ols_pixel([100.,999.,999.],[0,2,2],zeroframe=80.,**OPTIONS)
    assert first['branch']['one_good_first_group_branch']
    assert np.isclose(first['slope_double'],100./(1.5*float(np.float32(OPTIONS['frame_time']))))
    zero=piecewise_ols_pixel([999.,999.,999.],[2,2,2],zeroframe=80.,**OPTIONS)
    assert zero['branch']['zeroframe_branch']
    assert not zero['jacobian_group_and_zeroframe'][:3].any()
    assert np.isclose(zero['jacobian_group_and_zeroframe'][-1],1./float(np.float32(OPTIONS['frame_time'])))


def test_pure_jump_start_pruning_and_median_tie_visibility():
    r=piecewise_ols_pixel([0.,2.,4.,10.,12.,14.,16.,18.,20.],[0,1,0,4,0,0,0,0,0],**OPTIONS)
    assert [(x['start'],x['end']) for x in r['segments']]==[(3,9)]
    # Unique derivative of a single-segment slope survives a tied median:
    # its variance cancels out of its own normalized combination weight.
    assert r['branch']['median_tie']
    assert r['branch']['locally_unique_derivative']


def test_array_wrapper_matches_selected_pixel_and_invalid_mask():
    y=np.array([[20.,200.],[140.,400.]])
    val,jac,diag=piecewise_ols_arrays(y,np.zeros_like(y,dtype=np.uint8),
        pixeldq=np.array([0,1]),**OPTIONS)
    assert np.isfinite(val[0]) and np.isnan(val[1])
    assert jac.shape==(3,2) and not jac[:,1].any()
    assert diag[1]['status']=='invalid_pixel'
