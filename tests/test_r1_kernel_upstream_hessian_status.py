import json
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
STATUS = ROOT / "docs" / "r1_kernel_upstream_hessian_status.json"


def test_five_kernels_share_one_upstream_hessian_gate():
    obj = json.loads(STATUS.read_text(encoding="utf-8"))
    assert obj["source_lineage"]["metric_topographic_mixing"] == "NONZERO_FORMULA_OPEN"
    assert obj["source_lineage"]["full_self_adjoint_domain"] == "OPEN"

    kernels = obj["kernel_status"]
    assert set(kernels) == {
        "F_DL",
        "F_H",
        "F_DA",
        "F_fsigma8",
        "F_lens",
    }
    assert set(kernels.values()) == {"BLOCKED_BY_COMMON_N2_RESPONSE"}

    assert (
        obj["existing_growth_py"]["classification"]
        == "REFERENCE_BACKGROUND_GROWTH_RATIO_NOT_DIRECTIONAL_N2_KERNEL"
    )
