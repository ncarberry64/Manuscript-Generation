from pathlib import Path
import sys

ROOT = Path(__file__).resolve().parents[1]
CODE = ROOT / "code"
if str(CODE) not in sys.path:
    sys.path.insert(0, str(CODE))

import n2_topographic_legendre_gate as lg


def test_stationary_branch_is_inside_positive_cone():
    assert lg.parallel_eigenvalue(1.0, 1.0, 0.0) == 2.0
    assert lg.transverse_eigenvalue(1.0, 1.0) == 2.0


def test_dynamic_negative_branch_detected():
    assert lg.parallel_eigenvalue(1.0, 1.0, 0.40) < 0.0
