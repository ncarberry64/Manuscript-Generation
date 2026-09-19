import math
from pathlib import Path
import sys

ROOT = Path(__file__).resolve().parents[1]
CODE = ROOT / "code"
if str(CODE) not in sys.path:
    sys.path.insert(0, str(CODE))

import bhsm_los_topographic_transfer as los


def test_reference_normalization_is_exact():
    assert math.isclose(
        los.normalized_local_transfer(los.Z_REF),
        1.0,
        rel_tol=0.0,
        abs_tol=2e-12,
    )


def test_localized_response_decorrelates_with_distance():
    f03 = abs(los.normalized_local_transfer(0.30))
    f05 = abs(los.normalized_local_transfer(0.50))
    f10 = abs(los.normalized_local_transfer(1.00))
    assert f03 < 0.2
    assert f05 < f03
    assert f10 < f05


def test_seam_candidate_reproduces_fourth_power_identity():
    GM = 2.0
    RH = 10.0
    xi = 1.3
    c = 1.0
    rs = los.seam_radius(GM, RH, xi, c)
    V2 = los.seam_charge_speed2(GM, rs)
    assert math.isclose(
        V2 * V2,
        xi * GM * c * c / RH,
        rel_tol=2e-15,
        abs_tol=2e-15,
    )


def test_exterior_profile_activates_continuously():
    V2 = 0.7
    rs = 3.0
    assert los.exterior_topographic_acceleration_magnitude(
        rs, rs, V2
    ) == 0.0
    assert (
        los.exterior_topographic_acceleration_magnitude(
            2.0 * rs, rs, V2
        )
        > 0.0
    )


def test_path_average_decorrelates():
    near = los.path_average_rms_ratio(2.0, 1.0)
    far = los.path_average_rms_ratio(20.0, 1.0)
    assert 0.0 < far < near < 1.0
