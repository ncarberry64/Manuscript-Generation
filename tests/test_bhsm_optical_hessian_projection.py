from pathlib import Path
import sys

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
CODE = ROOT / "code"
if str(CODE) not in sys.path:
    sys.path.insert(0, str(CODE))

import bhsm_optical_hessian_projection as hp


def test_projection_matches_direct_entries():
    H = np.diag([2.0, 3.0, 5.0])
    P = np.array([[1.0,0.0],[0.0,1.0],[0.0,0.0]])
    A = hp.project_symmetric_block(H, P)
    np.testing.assert_allclose(A, np.diag([2.0,3.0]))


def test_mixed_projection():
    H = np.array([[1.0,2.0],[3.0,4.0],[5.0,6.0]])
    P = np.array([[1.0,0.0],[0.0,1.0],[0.0,0.0]])
    phi = np.array([2.0,-1.0])
    B = hp.project_mixed_column(H, P, phi)
    np.testing.assert_allclose(B, np.array([0.0,2.0]))


def test_response_residual_is_zero():
    c = hp.LocalOpticalCoefficients(
        a=4.0,b=1.0,c=3.0,u=2.0,v=-1.0,c_phi=5.0
    )
    assert c.response_residual() < 1e-14


def test_basis_covariance():
    c = hp.LocalOpticalCoefficients(
        a=4.0,b=1.0,c=3.0,u=2.0,v=-1.0,c_phi=5.0
    )
    S = np.array([[1.0,0.25],[0.1,1.0]])
    assert hp.validate_basis_covariance(c,S) < 1e-13


def test_schur_matches_block_elimination():
    c = hp.LocalOpticalCoefficients(
        a=4.0,b=1.0,c=3.0,u=2.0,v=-1.0,c_phi=5.0
    )
    B = c.mixed_column
    expected = 5.0 - B.T @ np.linalg.solve(c.metric_matrix, B)
    assert abs(c.schur_curvature() - expected) < 1e-15
