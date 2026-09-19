from pathlib import Path
import sys

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
CODE = ROOT / "code"
if str(CODE) not in sys.path:
    sys.path.insert(0, str(CODE))

import r1_n2_growth_lensing_transfer as gl
import r1_n2_reduced_propagator as rp


def _solutions():
    prop = rp.integrate_propagator()
    dust = [
        gl.integrate_dust_response(prop, 0),
        gl.integrate_dust_response(prop, 1),
    ]
    return prop, dust


def test_zero_dust_response_at_anchor():
    prop, dust = _solutions()
    M = gl.dust_transfer_matrix(
        prop,
        dust,
        rp.Z_START,
    )
    assert np.allclose(M, 0.0, atol=1e-12)


def test_dust_ode_dense_output_consistency():
    prop, dust = _solutions()
    assert gl.numerical_derivative_check(
        prop,
        dust,
        ngrid=10,
    ) < 5e-7


def test_weyl_kernel_is_finite():
    prop, _ = _solutions()
    for z in (2.1, 1.0, 0.0):
        k = gl.weyl_kernel(prop, z)
        assert np.all(np.isfinite(k))


def test_fsigma8_state_kernel_is_finite():
    prop, dust = _solutions()
    for z in (2.1, 1.0, 0.0):
        k = gl.fractional_fsigma8_kernel(
            prop,
            dust,
            z,
        )
        assert np.all(np.isfinite(k))
