from pathlib import Path
import sys

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
CODE = ROOT / "code"
if str(CODE) not in sys.path:
    sys.path.insert(0, str(CODE))

import n2_metric_kkt_projection as mk


def test_metric_response_matches_explicit_inverse():
    A = mk.metric_matrix(2.0, 0.25, 3.0)
    B = mk.mixed_matrix(0.1, 0.02, -0.05, 0.03)

    R = mk.metric_response(A, B)
    assert np.allclose(A @ R + B, 0.0, atol=1e-14)


def test_schur_is_symmetric():
    A = mk.metric_matrix(2.0, 0.25, 3.0)
    B = mk.mixed_matrix(0.1, 0.02, -0.05, 0.03)
    C = np.array([[5.0, -1.0], [-1.0, 13.0]])

    K = mk.schur_scalar(C, A, B)
    assert np.allclose(K, K.T, atol=1e-14)
