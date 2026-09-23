from pathlib import Path
import sys

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
CODE = ROOT / "code"
if str(CODE) not in sys.path:
    sys.path.insert(0, str(CODE))

import r1_n2_luminosity_distance_kernel as dl
import r1_n2_growth_lensing_transfer as gl
import r1_n2_reduced_propagator as rp


def _setup():
    prop = rp.integrate_propagator()
    dust = [
        gl.integrate_dust_response(prop, 0),
        gl.integrate_dust_response(prop, 1),
    ]
    chi = dl.integrate_background_lightcone()
    return prop, dust, chi


def test_curved_distance_has_flat_limit():
    x = 0.4
    # S_K(x)/x = 1 + O(K x^2); R1 curvature is small.
    assert abs(dl.S_K(x)/x - 1.0) < 2e-3


def test_distance_row_is_finite_and_dual():
    prop, dust, chi = _setup()
    r = dl.lightcone_rows(prop, dust, chi, 0.02)
    assert np.all(np.isfinite(r["F_DL"]))
    assert np.allclose(r["F_DA"], r["F_DL"], atol=0.0)


def test_calibration_null_is_exact():
    report = dl.diagnostics()
    assert abs(
        report["calibration"]["null_residual"]
    ) < 1e-12


def test_single_calibration_leaves_one_state_direction():
    report = dl.diagnostics()
    audit = report["calibration_rank_audit"]
    assert audit["phase_space_dimension"] == 2
    assert audit["rank"] == 1
    assert audit["remaining_state_direction_dimension"] == 1


def test_cross_observable_rows_not_all_parallel_to_DL():
    report = dl.diagnostics()
    rows = report["calibration_rank_audit"]["cross_observable_rows"]
    assert any(
        abs(r["det_DL_fsigma8"]) > 1e-8
        for r in rows
    )
    assert any(
        abs(r["det_DL_Weyl"]) > 1e-8
        for r in rows
    )
