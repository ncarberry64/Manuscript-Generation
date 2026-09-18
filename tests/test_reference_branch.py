from __future__ import annotations

import hashlib
import json
import math
import sys
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
CODE = ROOT / "code"
if str(CODE) not in sys.path:
    sys.path.insert(0, str(CODE))

import background
import growth
import perturbations

FROZEN_MANIFEST_SHA256 = "0a995d15171f3b133edaf203869329e9eec7913a9ce8e23e8469175bc7538bfa"


def canonical_bytes(obj) -> bytes:
    return json.dumps(
        obj, sort_keys=True, separators=(",", ":"), ensure_ascii=False
    ).encode("utf-8")


def test_manifest_integrity():
    manifest_path = ROOT / "preregistration" / "prediction_manifest.json"
    obj = json.loads(manifest_path.read_text(encoding="utf-8-sig"))
    digest = hashlib.sha256(canonical_bytes(obj)).hexdigest()
    assert digest == FROZEN_MANIFEST_SHA256


def test_r1_background_reproduces_and_is_admissible():
    v0 = background.shoot_v0()
    assert math.isclose(v0, 2.3986073449708067, rel_tol=0.0, abs_tol=5e-10)

    sol = background.integrate(v0)
    assert sol.success

    zs = [10.0, 3.0, 2.1, 1.5, 1.0, 0.8, 0.5, 0.3, 0.0]
    rows = [background.diagnostics(sol, v0, z) for z in zs]

    assert all(row["D_kin"] > 0.0 for row in rows)
    assert max(row["friedmann_residual"] for row in rows) < 1.0e-9
    assert abs(rows[-1]["E"] - 1.0) < 1.0e-8


def test_harmonic_filter_selects_n1_only():
    ell = 0.40
    assert perturbations.selective_filter_condition(ell)

    # The sign of F_n = 1 - ell^2 n(n+2) determines whether chi softens
    # (positive) or stiffens (negative) the harmonic.
    f1 = 1.0 - ell * ell * perturbations.nu(1)
    assert f1 > 0.0

    for n in range(2, 25):
        fn = 1.0 - ell * ell * perturbations.nu(n)
        assert fn < 0.0


def test_r1_growth_forecast_reproduces():
    gr1 = growth.integrate_growth(growth.growth_rhs_r1)
    gl = growth.integrate_growth(growth.growth_rhs_lcdm)

    expected = {
        2.1: 0.9934646335007011,
        1.5: 0.9886863442761461,
        1.0: 0.9816359978859575,
        0.8: 0.9779545909356957,
        0.5: 0.9727064416234715,
        0.3: 0.9716770504951750,
        0.0: 0.9848217331453974,
    }

    for z, target in expected.items():
        N = math.log(1.0 / (1.0 + z))
        _, dpr = gr1.sol(N)
        _, dpl = gl.sol(N)
        ratio = float(dpr / dpl)
        assert math.isclose(ratio, target, rel_tol=0.0, abs_tol=2e-8)
