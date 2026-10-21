"""Synthetic coupled-Jacobian, active-bound and convergence tests."""
from pathlib import Path
import sys
import numpy as np
from numpy.testing import assert_allclose
sys.path.insert(0,str(Path(__file__).resolve().parents[1]/'code'))
import strong_lensing_joint_scene as js
import strong_lensing_response_completion as rc


def synthetic():
    z=np.arange(41)/4-5;xx,yy=np.meshgrid(z,z);base=np.exp(-(xx*xx+yy*yy)/2)
    nodes=np.array([[-.7,.3],[.5,-.6]])
    def response():return js.PositiveResponse(base,nodes,np.eye(2),.65,np.zeros(2))
    y,x=np.mgrid[-3:4,-3:4];x=x.ravel().astype(float);y=y.ravel().astype(float);bg=np.c_[np.ones(len(x)),x/10,y/10]
    target_response=response();field_response=response();true_target=np.array([.13,-.09,3.,.1,.02,-.03]);true_field=np.array([-.1,.17,3.3,.2,.002,-.003])
    data=np.exp(true_target[2])*target_response.values(x-true_target[0],y-true_target[1])+bg@true_target[3:]
    common={'x':x,'y':y,'error':np.full(len(x),.1),'guess':np.zeros(2),'patch':0,'exposure':'E1','name':'A','flux_guess':np.exp(3.)}
    target={**common,'response':target_response,'data':data,'detector_background':bg,
            'sky_fine':np.zeros((9,41,41)),'fine_ix':((x+5)*4).astype(int),'fine_iy':((y+5)*4).astype(int),'scene_matrix':np.zeros((len(x),9))}
    fieldbg=np.c_[np.ones(len(x)),x,y];fielddata=np.exp(true_field[2])*field_response.values(x-true_field[0],y-true_field[1])+fieldbg@true_field[3:]
    field={**common,'name':'field','response':field_response,'data':fielddata}
    problem=rc.CoupledScene(js.SceneProblem([target],scene_fixed=np.zeros(36)),[field],np.zeros(2))
    return problem,np.r_[true_target,[0.,0.],true_field]


def test_active_bound_projected_gradient_signs():
    p=np.array([0.,0.,1.,1.,.5]);g=np.array([2.,-3.,-4.,5.,6.]);lo=np.zeros(5);hi=np.ones(5)
    projected,_,_=rc.projected_gradient(p,g,lo,hi,np.ones(5))
    assert_allclose(projected,[0.,-3.,0.,5.,6.])


def test_full_coupled_jacobian_directional_derivative():
    problem,truth=synthetic();x=truth.copy();x[6:8]=[.03,-.02]
    r,j=problem.evaluate(x,True);direction=np.random.default_rng(42).normal(size=len(x))*.03;h=1e-5
    plus=problem.evaluate(x+h*direction)[0];minus=problem.evaluate(x-h*direction)[0]
    assert_allclose(j@direction,(plus-minus)/(2*h),rtol=2e-7,atol=2e-8)


def test_full_nonlinear_synthetic_recovery_and_kkt():
    problem,truth=synthetic();initial=truth.copy();initial[:2]+=[.05,-.06];initial[6:8]=[.03,-.02];initial[8:10]+=[-.04,.05]
    fit,diagnostic=problem.solve(initial,max_nfev=100)
    assert fit.success
    assert diagnostic['scaled_projected_KKT_gradient_max']<1e-6
    assert_allclose(fit.x,truth,atol=2e-7)
    assert diagnostic['target_raw_chi2']<1e-15
