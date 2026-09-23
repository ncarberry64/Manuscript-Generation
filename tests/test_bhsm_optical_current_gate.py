from pathlib import Path
import sys

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
CODE = ROOT / "code"
if str(CODE) not in sys.path:
    sys.path.insert(0, str(CODE))

import bhsm_optical_current_gate as cg


def test_zero_current_closes_mixed_source_gate():
    fp = np.ones(2)
    fpp = np.ones(2)
    j = np.zeros((2,3))
    gp = np.ones((2,3))
    A = np.tile(np.eye(3),(2,1,1))
    g = np.tile(np.eye(3),(2,1,1))
    S = cg.spatial_mixed_source_tensor(fp,fpp,j,gp,A,g)
    np.testing.assert_allclose(S,0.0)


def test_zero_texture_gradient_closes_mixed_source_gate():
    fp = np.ones(2)
    fpp = np.ones(2)
    j = np.ones((2,3))
    gp = np.zeros((2,3))
    A = np.tile(np.eye(3),(2,1,1))
    g = np.tile(np.eye(3),(2,1,1))
    S = cg.spatial_mixed_source_tensor(fp,fpp,j,gp,A,g)
    np.testing.assert_allclose(S,0.0)


def test_nonzero_spatial_current_can_source_metric_projection():
    fp = np.array([1.0])
    fpp = np.array([0.5])
    j = np.array([[1.0,0.0,0.0]])
    gp = np.array([[2.0,0.0,0.0]])
    A = np.array([np.diag([1.0,0.0,0.0])])
    g = np.array([np.eye(3)])
    S = cg.spatial_mixed_source_tensor(fp,fpp,j,gp,A,g)
    eA = np.array([np.diag([1.0,0.0,0.0])])
    eP = np.array([np.eye(3)])
    uv = cg.project_mixed_coefficients(np.array([1.0]),eA,eP,S)
    assert np.linalg.norm(uv) > 0.0


def test_rank_one_covariance_has_rank_at_most_one():
    C = cg.rank_one_source_covariance(np.array([2.0,-1.0]),3.0)
    assert np.linalg.matrix_rank(C, tol=1e-12) == 1


def test_multichannel_covariance_rank_bounded_by_channel_count():
    B = np.array([[1.0,0.0],[0.5,1.0],[0.0,2.0]])
    Cq = np.diag([2.0,3.0])
    C = cg.multichannel_source_covariance(B,Cq)
    assert np.linalg.matrix_rank(C, tol=1e-12) <= 2
