"""Scientific gates for the environmental-state manuscript integration."""
import importlib.util
from pathlib import Path
import numpy as np
import pytest

ROOT=Path(__file__).resolve().parents[1]
SPEC=importlib.util.spec_from_file_location('environment_audit',ROOT/'scripts/audit_r1_coupled_environmental_state.py')
AUDIT=importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(AUDIT)

@pytest.fixture(scope='module')
def result(): return AUDIT.build_audit(ROOT)

def test_seed_and_stored_canonical_propagator(result):
    assert result['topographic_seed_recovery_max_abs_difference']<1e-12
    assert result['stored_full_propagator_z1p5_max_abs_difference']<1e-8
    assert result['seed_inverse_residual']<1e-12

def test_environment_and_matter_alone_span_both_output_bases(result):
    assert len(result['rows'])==8
    for row in result['rows']:
        assert row['rank_rel_1e-8']==row['q2_Pi2_rank_rel_1e-8']==2
        assert row['matter_delta_v_determinant']>0
        assert row['q2_Pi2_matter_determinant']>0
        assert row['Pi2_over_zeta_dot']>0
        assert row['output_basis_conversion_error']<1e-12
    assert result['retuning'] is False and result['observational_environment_selection'] is False

def test_common_spatial_profile_needs_aligned_initial_topography(result):
    U=np.array(result['rows'][0]['physical_propagator'])
    E=np.outer([1.,0.,1.,0.],np.eye(9)[0])
    X=U[:2,2:]@E
    assert np.linalg.matrix_rank(X)==1
    Xi=np.outer([1.,0.],np.eye(9)[1])
    assert np.linalg.matrix_rank(U[:2,:2]@Xi+X)==2

def test_generic_environment_produces_two_spatial_patterns(result):
    E=np.zeros((4,9)); E[0,0]=1;E[2,1]=1
    for row in result['rows']:
        assert np.linalg.matrix_rank(np.array(row['U_XE_q2_Pi2'])@E)==2
