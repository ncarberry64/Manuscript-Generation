"""Synthetic integration checks: no numerical Gaia target rows are used."""
import sys
from pathlib import Path
import numpy as np
import pytest
sys.path.insert(0,str(Path(__file__).resolve().parents[1]/'scripts'))
from run_he0435_gaia_validation import score_catalogue_rows,start_query_lock


def synthetic():
    xy=np.array([[0.,0.],[-1.4,.5],[-2.4,-.6],[-.9,-1.6]])
    rows=[{'source_id':str(i),'ra':70.+x/3600/np.cos(np.deg2rad(-12.)),
           'dec':-12.+y/3600,'ra_error':'.2','dec_error':'.3','ra_dec_corr':'.1'}
          for i,(x,y)in enumerate(xy)]
    models={name:{'relative_B_C_D_to_A_east_north_arcsec':xy[1:].reshape(6).tolist()}
            for name in ['conventional','conditional_R1']}
    return rows,models


def test_synthetic_runner_schema_common_covariance_and_null_error():
    rows,models=synthetic();diag,observed,ordered,scored,cov=score_catalogue_rows(rows,models,[0.,2.,5.])
    assert diag['permutation_abcd']==[0,1,2,3]
    assert observed.shape==(3,2) and cov.shape==(6,6)
    assert set(scored['conventional'])=={'0.0','2.0','5.0'}
    assert scored['conventional']['0.0']['formal_chi2']<.01
    assert scored['conventional']==scored['conditional_R1']
    rows[2]['ra_dec_corr']=''
    _,obs,_,scored,cov=score_catalogue_rows(rows,models,[0.,2.,5.])
    assert cov is None and 'formal_score_undefined'in scored['conventional']
    assert np.all(np.isfinite(obs)) and len(scored['conventional']['residual_mas'])==3


def test_runner_lock_remains_after_simulated_failure(tmp_path):
    path=tmp_path/'exposure.json';start_query_lock(path,{'status':'before-network'})
    with pytest.raises(FileExistsError):start_query_lock(path,{'status':'retry'})
    assert path.exists()
