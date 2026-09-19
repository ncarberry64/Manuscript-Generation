from pathlib import Path
import sys

ROOT = Path(__file__).resolve().parents[1]
CODE = ROOT / "code"
if str(CODE) not in sys.path:
    sys.path.insert(0, str(CODE))

import validate_r1_kernel_closure_spec as vc


def test_r1_kernel_closure_protocol_is_frozen_and_complete():
    report = vc.validate()
    assert report["status"] == "R1_KERNEL_CLOSURE_PROTOCOL_VALID"
    assert set(report["required_kernels"]) == {
        "F_DL",
        "F_H",
        "F_DA",
        "F_fsigma8",
        "F_lens",
    }
