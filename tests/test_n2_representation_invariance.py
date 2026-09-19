from pathlib import Path
import sys

ROOT = Path(__file__).resolve().parents[1]
CODE = ROOT / "code"
if str(CODE) not in sys.path:
    sys.path.insert(0, str(CODE))

import n2_representation_invariance as ri


def test_representation_invariance():
    report = ri.diagnostics()
    assert report["status"] == "N2_REPRESENTATION_INVARIANCE_VALID"
    assert report["quadratic_form_residual"] < 1e-14
    assert report["metric_response_residual"] < 1e-14
    assert report["observable_residual"] < 1e-14
