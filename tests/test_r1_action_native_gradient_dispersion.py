import json
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
BASE = ROOT / "artifacts/R1_ACTION_NATIVE_GRADIENT_DISPERSION_V1.json"
CONV = ROOT / "artifacts/R1_ACTION_NATIVE_GRADIENT_DISPERSION_CONVERGENCE_V1.json"


def load(path):
    return json.loads(path.read_text(encoding="utf-8"))


def test_gradient_dispersion_artifacts_pass_and_preserve_owner():
    base = load(BASE)
    conv = load(CONV)
    wanted = "PASS_NUMERICAL_PRINCIPAL_GRADIENT_DISPERSION_AUDIT"
    assert base["status"] == wanted
    assert conv["status"] == wanted
    for data in (base, conv):
        own = data["n2_owner_match"]
        assert own["worst_action_block_difference"] <= 5e-12
        assert own["worst_reduced_block_difference"] <= 5e-12
        assert data["observational_fits"] == 0
        assert data["parameter_retuning"] is False
        assert data["gate7_dependency"] is False


def test_extended_principal_symbol_converges_and_is_gradient_stable():
    base = load(BASE)["principal_symbol"]
    conv = load(CONV)["principal_symbol"]
    assert conv["max_growth_intercept_fit"] < base["max_growth_intercept_fit"]
    assert conv["max_growth_intercept_fit"] < 1e-3
    assert conv["min_scalar_c2_intercept"] > 0.0
    assert conv["max_radiation_c2_error"] < 5e-3
    assert conv["all_tail_growth_nonincreasing"] is True


def test_independent_gravity_scalar_diagnostics_positive():
    conv = load(CONV)
    g = conv["gravity_sector_owner"]
    assert g["min_G_S_n2"] > 0.0
    assert g["min_F_S_n2_gravity_sector_only"] > 0.0
    assert g["min_c_s2_gravity_sector_only"] > 0.0
