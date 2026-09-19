import math
from pathlib import Path
import sys

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
CODE = ROOT / "code"
if str(CODE) not in sys.path:
    sys.path.insert(0, str(CODE))

import n2_scalar_topographic_candidate as st


def test_canonical_n2_candidate_is_positive():
    report = st.diagnostics()
    assert report["positive"] is True
    assert report["positivity_margin"] == 64.0

    expected = [
        9.0 - math.sqrt(17.0),
        9.0 + math.sqrt(17.0),
    ]
    assert np.allclose(
        report["eigenvalues_L1"],
        expected,
        rtol=0.0,
        atol=1e-14,
    )


def test_general_positivity_margin_matches_determinant():
    M = st.candidate_matrix(lambda_rho=2.5, I2=2.0)
    assert np.isclose(
        np.linalg.det(M).real,
        st.positivity_margin(lambda_rho=2.5, I2=2.0),
    )
