from pathlib import Path
import sys

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
CODE = ROOT / "code"
if str(CODE) not in sys.path:
    sys.path.insert(0, str(CODE))

import bhsm_optical_constitutive as con


def test_static_export_has_retained_r_inverse_square_shape():
    r = np.array([2.0, 3.0, 6.0])
    s = con.static_seam_export(
        r,
        r_star=3.0,
        q_sigma=18.0,
        q_coupling=2.0,
    )
    np.testing.assert_allclose(
        s,
        np.array([0.0, 1.0, 0.25]),
        rtol=0.0,
        atol=1e-15,
    )


def test_moving_seam_distributional_coefficients():
    bulk, seam = con.dynamic_seam_export_components(
        r_star=2.0,
        q_sigma=8.0,
        q_sigma_dot=3.0,
        r_star_dot=0.5,
        q_coupling=2.0,
    )
    assert bulk == 1.5
    assert seam == -0.5


def test_metric_response_matches_direct_linear_solve():
    a, b, c = 4.0, 1.0, 3.0
    u, v = 2.0, -1.0
    response = con.local_metric_response(a, b, c, u, v)
    A = np.array([[a, b], [b, c]])
    B = np.array([u, v])
    np.testing.assert_allclose(
        A @ response + B,
        np.zeros(2),
        rtol=0.0,
        atol=1e-14,
    )


def test_weyl_and_slip_reconstruct_metric_response():
    r = con.local_metric_response(4.0, 1.0, 3.0, 2.0, -1.0)
    w = con.weyl_response_coefficient(4.0, 1.0, 3.0, 2.0, -1.0)
    s = con.slip_response_coefficient(4.0, 1.0, 3.0, 2.0, -1.0)
    assert abs((w + s) - r[0]) < 1e-15
    assert abs((w - s) - r[1]) < 1e-15


def test_covariance_pushforward_is_positive_semidefinite():
    O = np.array([[1.0, 0.3], [0.0, 1.0]])
    R = np.array([[0.5, 0.0], [0.2, 0.7]])
    G = np.eye(2)
    C = np.array([[2.0, 0.4], [0.4, 1.0]])
    Q = con.pushforward_covariance(O, R, G, C)
    eig = np.linalg.eigvalsh(Q)
    assert np.min(eig) >= -1e-14


def test_compound_poisson_covariance_scales_with_rate():
    M2 = np.diag([2.0, 3.0])
    Q = con.compound_poisson_source_covariance(4.0, M2)
    np.testing.assert_allclose(Q, np.diag([8.0, 12.0]))
