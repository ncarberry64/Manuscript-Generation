from pathlib import Path
import math
import sys

ROOT = Path(__file__).resolve().parents[1]
CODE = ROOT / "code"
if str(CODE) not in sys.path:
    sys.path.insert(0, str(CODE))

import lowz_distance_hubble_corollary as lc


def test_distance_modulus_conversion():
    report = lc.diagnostics()
    expected = math.log(10.0) / 5.0 * (-0.0412)
    assert abs(report["delta_DL_over_DL"] - expected) < 1e-15


def test_distance_duality_fractional_anisotropy():
    report = lc.diagnostics()
    assert report["delta_DA_over_DA"] == report["delta_DL_over_DL"]


def test_lowz_hubble_is_opposite_fractional_distance():
    report = lc.diagnostics()
    assert report["delta_H0_over_H0"] == -report["delta_DL_over_DL"]
