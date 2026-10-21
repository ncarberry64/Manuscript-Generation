from pathlib import Path
import sys
import numpy as np
from scipy import sparse
sys.path.insert(0,str(Path(__file__).resolve().parents[1]/'code'))
from strong_lensing_native_completion import right_harmonic,native_rate_adjoint,NativeCoupledScene
from strong_lensing_native_calibration import FixedHarmonicContext
from test_strong_lensing_native_calibration import model


def test_harmonic_adjoint_preserves_the_same_linear_map():
    context=FixedHarmonicContext((2,3),[True,False,True,True,False,True])
    rng=np.random.default_rng(381);a=rng.normal(size=(3,6));d=rng.normal(size=(6,4))
    np.testing.assert_allclose(right_harmonic(a,context)@d,a@context.action(d),rtol=2e-14,atol=2e-14)


def test_read_adjoint_matches_forward_chain_with_caps_and_intercept():
    obj,rate=model();state=obj.evaluate(rate,1e-5,derivatives=True)
    rng=np.random.default_rng(146);directions=rng.normal(size=(3,5));dk=rng.normal(size=5)*1e-6
    weights=rng.normal(size=(10,3))*.01
    actual,kappa=native_rate_adjoint(obj,state,np.arange(3),weights)
    forward=obj.project_directions(state,directions,dk,np.arange(3),weights[:9],weights[9])
    np.testing.assert_allclose(actual@directions+kappa[:,None]*dk,forward,rtol=2e-13,atol=2e-13)


def test_positive_small_kappa_is_not_mistaken_for_an_active_zero_bound():
    from types import SimpleNamespace
    obj=object.__new__(NativeCoupledScene)
    obj.target=SimpleNamespace(records=[]);obj.fields=[];obj.n=1;obj.eta_start=1;obj.field_start=2;obj.kappa_index=3
    obj.scales=np.array([1.,.2,1.,1e-7]);obj.evaluations=1;obj.jacobians=1
    obj.bounds=lambda:(np.array([-np.inf,-np.inf,-np.inf,0.]),np.array([np.inf,np.inf,np.inf,1e-6]))
    obj.evaluate=lambda x,jacobian:(np.ones(4),np.diag([1.,1.,1.,1e10]))
    result=obj.diagnostics(np.array([0.,0.,0.,5e-9]))
    assert result['scaled_projected_KKT_kappa']==1000.
    assert 3 not in result['active_lower_indices']
