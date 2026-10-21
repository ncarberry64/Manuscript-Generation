from pathlib import Path
import sys
import numpy as np
import pytest
sys.path.insert(0, str(Path(__file__).resolve().parents[1]/"code"))
from strong_lensing_bounded_collection import bounded_leading_collection
from strong_lensing_stochastic_response import nearest_neighbor_edges


def app(rate, caps, kind="storage_ceiling", kappa=.001, derivatives=False):
    return bounded_leading_collection(rate, [1., 4.], kappa, caps,
        nearest_neighbor_edges((1, 3)), prescription=kind,
        ceiling_source="synthetic declared trial, not detector calibration", derivatives=derivatives)


def test_original_leading_mean_and_conservation_before_ceiling():
    r=np.array([1., 5., 2.]);a=app(r,[100., 100., 100.])
    lap=np.array([r[1]**2-r[0]**2, r[0]**2+r[2]**2-2*r[1]**2, r[1]**2-r[2]**2])
    want=np.array([1.,4.])[:,None]*r+.0005*np.array([1.,16.])[:,None]*lap
    np.testing.assert_allclose(a["stored_mean_e"],want)
    np.testing.assert_allclose(a["stored_mean_e"].sum(axis=1), [8.,32.])
    np.testing.assert_allclose(a["lost_mean_e"],0.,atol=1e-14)


def test_finite_storage_and_readout_ceiling_have_different_neighbor_means():
    r=np.array([1.,8.,2.]);caps=np.array([100.,8.,100.])
    stored=app(r,caps);electronic=app(r,caps,"readout_ceiling")
    assert stored["observed_read_mean_e"][1,0] < electronic["observed_read_mean_e"][1,0]
    assert stored["lost_mean_e"][1] > 0
    assert electronic["lost_mean_e"][1] == 0
    np.testing.assert_allclose(electronic["stored_mean_e"].sum(axis=1), [11.,44.])


@pytest.mark.parametrize("kind",["storage_ceiling","readout_ceiling"])
def test_piecewise_analytic_rate_cap_and_kappa_jacobians(kind):
    rate=np.array([1.3,8.1,2.2]);caps=np.array([100.,9.7,100.]);k=.001
    base=app(rate,caps,kind,k,True);h=1e-5
    assert not base["nondifferentiable_transition_mask"].any()
    for i in range(3):
        e=np.eye(3)[i]*h
        fd=(app(rate+e,caps,kind,k)["observed_read_mean_e"]-app(rate-e,caps,kind,k)["observed_read_mean_e"])/(2*h)
        actual=np.stack([j.toarray()[:,i] for j in base["observed_rate_jacobian_s"]])
        np.testing.assert_allclose(actual,fd,rtol=2e-7,atol=1e-8)
        fd=(app(rate,caps+e,kind,k)["observed_read_mean_e"]-app(rate,caps-e,kind,k)["observed_read_mean_e"])/(2*h)
        actual=np.stack([j.toarray()[:,i] for j in base["observed_cap_jacobian"]])
        np.testing.assert_allclose(actual,fd,rtol=2e-7,atol=1e-8)
    fd=(app(rate,caps,kind,k+h)["observed_read_mean_e"]-app(rate,caps,kind,k-h)["observed_read_mean_e"])/(2*h)
    np.testing.assert_allclose(base["observed_kappa_jacobian_e2"],fd,rtol=1e-8,atol=1e-8)


def test_domain_requires_named_ceiling_and_prescription():
    with pytest.raises(ValueError,match="source"):
        bounded_leading_collection([1.],[1.],0.,[4.],[],prescription="storage_ceiling",ceiling_source="")
    with pytest.raises(ValueError,match="explicitly"):
        bounded_leading_collection([1.],[1.],0.,[4.],[],prescription="unknown",ceiling_source="trial")
