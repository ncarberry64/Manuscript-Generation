from pathlib import Path
import sys

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
CODE = ROOT / "code"
if str(CODE) not in sys.path:
    sys.path.insert(0, str(CODE))

import bhsm_optical_low_rank_covariance as lr


def test_one_channel_pushforward_is_rank_one():
    M = np.array([[1.0],[2.0],[-0.5],[3.0]])
    C = lr.pushforward_channel_covariance(M,np.array([[2.0]]))
    assert lr.covariance_rank(C) == 1


def test_rank_bound_for_two_channels():
    M = np.array([
        [1.0,0.0],
        [0.5,1.0],
        [1.0,-1.0],
        [0.0,2.0],
    ])
    Cq = np.array([[2.0,0.3],[0.3,1.0]])
    C = lr.pushforward_channel_covariance(M,Cq)
    assert lr.covariance_rank(C) <= 2


def test_rank_one_principal_minors_vanish():
    v = np.array([1.0,2.0,-3.0,0.5])
    C = np.outer(v,v)
    minors = lr.principal_minor_residuals(C)
    np.testing.assert_allclose(minors,0.0,rtol=0.0,atol=1e-14)


def test_rank_one_diagnostics_are_exact():
    v = np.array([1.0,2.0,3.0])
    C = np.outer(v,v)
    d = lr.rank_diagnostics(C)
    assert abs(d["R1"] - 1.0) < 1e-14
    assert d["rho21"] < 1e-14


def test_eigenclip_is_psd():
    C = np.array([[1.0,2.0],[2.0,1.0]])
    P = lr.nearest_psd_eigenclip(C)
    assert np.min(np.linalg.eigvalsh(P)) >= -1e-14
