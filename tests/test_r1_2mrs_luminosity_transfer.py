from pathlib import Path
import sys
import numpy as np
import pandas as pd

ROOT=Path(__file__).resolve().parents[1]
CODE=ROOT/"code"
if str(CODE) not in sys.path:
    sys.path.insert(0,str(CODE))

import r1_2mrs_luminosity_transfer as lt


def test_kmag_priority():
    df=pd.DataFrame({"kmag":[10.],"k_m_k20fe":[9.]})
    assert lt.find_kmag(df)=="k_m_k20fe"


def test_luminosity_brighter_is_heavier():
    w=lt.luminosity_weight(np.array([8.,10.]),np.array([0.02,0.02]))
    assert w[0]>w[1]


def test_weight_positive():
    w=lt.luminosity_weight(np.array([8.,9.]),np.array([0.01,0.02]))
    assert np.all(w>0)
