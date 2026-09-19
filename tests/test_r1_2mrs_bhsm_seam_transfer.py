from pathlib import Path
import sys

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
CODE = ROOT / "code"
if str(CODE) not in sys.path:
    sys.path.insert(0, str(CODE))

import r1_2mrs_bhsm_seam_transfer as seam


def test_brighter_object_has_larger_seam_weight():
    w = seam.seam_charge_weight(
        np.array([8.0, 10.0]),
        np.array([0.02, 0.02]),
    )
    assert w[0] > w[1]


def test_weight_is_square_root_of_linear_luminosity_proxy():
    w = seam.seam_charge_weight(
        np.array([8.0, 10.0]),
        np.array([0.02, 0.02]),
    )
    ratio = w[0] / w[1]
    expected = 10.0 ** (-0.2 * (8.0 - 10.0))
    assert abs(ratio - expected) < 1e-12


def test_seam_weight_positive():
    w = seam.seam_charge_weight(
        np.array([8.0, 9.0]),
        np.array([0.01, 0.02]),
    )
    assert np.all(w > 0)
