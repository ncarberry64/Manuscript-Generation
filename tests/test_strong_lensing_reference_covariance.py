import sys
from pathlib import Path
import numpy as np
from numpy.testing import assert_allclose
sys.path.insert(0,str(Path(__file__).resolve().parents[1]/'code'))
from strong_lensing_reference_covariance import conditional_sirs


def example(seed=1):
    rng=np.random.default_rng(seed);n=32;yy,xx=np.indices((n,n));edge=(yy<4)|(yy>=n-4)|(xx<4)|(xx>=n-4)
    y,x=np.where(edge);xy=np.c_[x,y];values=rng.normal(0,3,len(x)).astype(np.float32)
    # Deliberate top and side outliers, fixed observed rejection must exclude.
    values[3]=200;dq=np.zeros(len(x),np.uint32);dq[7]=1
    out=np.array([[10,12],[11,12],[10,12],[20,18]])
    kernels=([np.array([.1,.3,.1])]*4,[np.array([.05,.15,.05])]*4)
    var=np.linspace(1,4,len(x))
    return values,xy,dq,out,kernels,var,n


def test_factored_operator_and_covariance_keep_all_cross_terms():
    a=example();o=conditional_sirs(*a[:-1],size=a[-1]);j=o['D'].toarray()+o['U']@o['T'].toarray()
    assert_allclose(o['covariance'],(j*a[5])@j.T,atol=2e-15)
    assert_allclose(o['factor']@o['factor'].T,o['covariance'],atol=2e-14)
    assert o['covariance'][0,1]!=0
    assert o['pixel_feature_index'][0]==o['pixel_feature_index'][2]
    assert o['eigenvalues'][0]>-o['eigenvalue_roundoff_tolerance']


def test_frozen_mask_reference_derivative_and_y_orientation():
    a=list(example(8));o=conditional_sirs(*a[:-1],size=a[-1]);j=o['D'].toarray()+o['U']@o['T'].toarray()
    direction=np.random.default_rng(5).normal(size=len(a[0]));h=.001
    b=a.copy();b[0]=(a[0]+h*direction).astype(np.float32);plus=conditional_sirs(*b[:-1],size=b[-1])
    b[0]=(a[0]-h*direction).astype(np.float32);minus=conditional_sirs(*b[:-1],size=b[-1])
    # Float32 native top/bottom arithmetic is explicitly retained in mean replay.
    assert_allclose((plus['native_correction']-minus['native_correction'])/(2*h),j@direction,atol=3e-4)
    assert [19,1,0] in o['features'].tolist()


def test_native_persistent_sirs_work_buffer_dtype_changes_rounding():
    a=list(example(12));a[0]=(a[0]+20000).astype(np.float32)
    first=conditional_sirs(*a[:-1],size=a[-1],native_work_dtype=np.float32)
    later=conditional_sirs(*a[:-1],size=a[-1],native_work_dtype=np.float64)
    # NIRDataset retains the float64 array returned by SIRS after group zero.
    # Pretending every group has the first group's float32 work buffer loses
    # native precision, even though the fixed-mask derivative is unchanged.
    assert first['native_work_dtype']=='float32'
    assert later['native_work_dtype']=='float64'
    first_rounding=np.max(abs(first['linear_operator_reconstruction']-first['native_correction']))
    later_rounding=np.max(abs(later['linear_operator_reconstruction']-later['native_correction']))
    assert first_rounding>1e-5
    assert later_rounding<1e-9
