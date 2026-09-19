from pathlib import Path
import sys

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
CODE = ROOT / "code"
if str(CODE) not in sys.path:
    sys.path.insert(0, str(CODE))

import bhsm_optical_transfer as opt


def test_closed_s3_background_jacobi_identity():
    for chi in (0.1, 0.5, 1.0):
        assert abs(opt.jacobi_background_residual(chi, 0.018)) < 1e-15


def test_flat_limit_SK_and_HK():
    chi = 0.7
    assert abs(opt.S_K(chi, 0.0) - chi) < 1e-15
    assert abs(opt.C_K(chi, 0.0) - 1.0) < 1e-15
    assert abs(opt.H_K(chi, 0.0) - 1.0 / chi) < 1e-15


def test_source_matrix_places_sources_in_correct_equations():
    source = np.array([1.0, 2.0, 3.0, 4.0])
    mapped = opt.B_OPT @ source
    expected = np.array([1.0, 0.0, 2.0, 0.0, 3.0, 0.0, 4.0])
    np.testing.assert_allclose(mapped, expected, rtol=0.0, atol=0.0)


def test_pure_isotropic_tidal_source_has_no_shear():
    t = np.array([[2.0, 0.0], [0.0, 2.0]])
    f, g1, g2 = opt.optical_tidal_components(t)
    assert f == -2.0
    assert g1 == 0.0
    assert g2 == 0.0


def test_trace_free_tidal_source_has_no_convergence():
    t = np.array([[2.0, 1.5], [1.5, -2.0]])
    f, g1, g2 = opt.optical_tidal_components(t)
    assert f == 0.0
    assert g1 == -2.0
    assert g2 == -1.5


def test_covariance_injection_is_positive_semidefinite():
    C = np.zeros((7, 7))
    Q = np.diag([1.0, 2.0, 3.0, 4.0])
    dC = opt.covariance_rhs(0.4, C, 0.018, Q)
    vals = np.linalg.eigvalsh(0.5 * (dC + dC.T))
    assert np.min(vals) >= -1e-14


def test_fixed_affine_distance_relation():
    assert opt.fixed_affine_fractional_distance(0.01, 0.003) == 0.017
