from pathlib import Path
import sys
import numpy as np

ROOT=Path(__file__).resolve().parents[1]
CODE=ROOT/"code"
if str(CODE) not in sys.path:
    sys.path.insert(0,str(CODE))

import r1_2mrs_foreground_transfer as fg


def test_axis_cosine_bounded():
    x=fg.axis_cosine(np.array([0.,90.,180.]),np.array([0.,0.,0.]))
    assert np.all(np.abs(x)<=1.0+1e-15)


def test_gaussian_counts_peak_on_source():
    sn=fg.radec_to_unit([0.,90.],[0.,0.])
    gal=fg.radec_to_unit([0.],[0.])
    c=fg.gaussian_sky_counts(sn,gal,sigma_deg=10.0)
    assert c[0] > c[1]


def test_foreground_visibility_is_causal():
    # A shell beyond the SN must not contribute.
    z=np.array([0.015,0.045])
    d=np.zeros((2,5))
    d[:,4]=1.0  # 0.04-0.05 shell
    raw=fg.foreground_raw_score(z,d)
    assert abs(raw[0]) < 1e-14
    assert abs(raw[1]) > 0.0


def test_standardize_unit_variance():
    z,_,_=fg.standardize(np.array([1.,2.,3.,4.]))
    assert abs(np.mean(z)) < 1e-14
    assert abs(np.std(z,ddof=1)-1.0) < 1e-14
