import sys
from pathlib import Path
import numpy as np
import pytest

sys.path.insert(0,str(Path(__file__).resolve().parents[1]/'code'))
from strong_lensing_comparison import (
    ARCSEC,C_KM_S,MPC_OVER_C_DAYS,LensGeometry,GalaxyPotential,
    solve_images,complete_map_response,observable_vector,response_vector,
    nuisance_columns,nuisance_projected_information,
)


def example():
    geometry=LensGeometry(900*np.eye(2),1700*np.eye(2),1000*np.eye(2),.3)
    b=4*np.pi*(220/C_KM_S)**2
    galaxy=GalaxyPotential(b)
    thetaE=b*1000/1700
    beta=np.array([.2*thetaE,0.])
    return geometry,galaxy,geometry.D_os@beta,thetaE,beta


def test_sis_analytic_images_relative_arrival_and_stationarity():
    g,p,s,e,beta=example();a=solve_images(g,p,s)
    np.testing.assert_allclose(a['images_rad'],[[beta[0]+e,0],[beta[0]-e,0]],atol=1e-15)
    analytic=2*g.T[0,0]*e*beta[0]*MPC_OVER_C_DAYS
    np.testing.assert_allclose(a['relative_delays_days'],[analytic],rtol=2e-13)
    assert a['max_stationarity_Mpc']<1e-14


def test_complete_response_agrees_with_independently_resolved_images():
    g,p,s,_,_=example();a=solve_images(g,p,s)
    d=np.array([.07*g.D_ol,-.03*g.D_os,.02*g.D_ls])
    r=complete_map_response(a['images_rad'],g,p,s,d)
    differences=[]
    for h in (2e-4,1e-4):
        plus=solve_images(g.displaced(d,h),p,s)
        minus=solve_images(g.displaced(d,-h),p,s)
        differences.append((observable_vector(plus)-observable_vector(minus))/(2*h))
    np.testing.assert_allclose(response_vector(r),differences[-1],rtol=2e-8,atol=2e-9)
    assert np.linalg.norm(differences[-1]-response_vector(r))<np.linalg.norm(differences[0]-response_vector(r))*1.1
    # Holding physical s means beta changes; the source term must survive.
    assert np.linalg.norm(r['source_angle_response'])>0
    assert np.max(abs(r['arrival_response_terms_days'][:,1]))>0


def test_incompatible_distance_maps_are_not_silently_symmetrized():
    with pytest.raises(ValueError,match='not symmetric'):
        LensGeometry(900*np.eye(2),np.array([[1700,100],[0,1700]]),1000*np.eye(2),.3)


def test_mass_sheet_transformation_images_and_delays():
    g,p,s,_,_=example();a=solve_images(g,p,s);lam=.91
    sheet=(1-lam)/(1+g.z_l)*np.linalg.solve(g.D_ol.T,g.T)@np.linalg.inv(g.D_ol)
    transformed=GalaxyPotential(lam*p.b,sheet_per_Mpc=tuple(map(tuple,sheet)))
    b=solve_images(g,transformed,lam*s)
    np.testing.assert_allclose(b['images_rad'],a['images_rad'],atol=2e-15)
    np.testing.assert_allclose(b['relative_delays_days'],lam*a['relative_delays_days'],rtol=3e-13)


def test_hubble_delay_scaling_depends_on_source_convention():
    g,p,s,_,_=example();a=solve_images(g,p,s)
    d=-np.array([g.D_ol,g.D_os,g.D_ls])
    fixed_physical=complete_map_response(a['images_rad'],g,p,s,d)
    fixed_angular=complete_map_response(a['images_rad'],g,p,s,d,d_source_Mpc=-s)
    np.testing.assert_allclose(fixed_angular['relative_delays_days'],-a['relative_delays_days'],atol=2e-12)
    # SIS deflection strength has no length scale. At fixed physical source,
    # beta~H0 cancels Ddt~H0^-1 in its relative delay.
    np.testing.assert_allclose(fixed_physical['relative_delays_days'],0,atol=3e-12)


def test_free_nuisance_projection_has_unconstrained_not_zero_variance_states():
    g,p,s,_,_=example();a=solve_images(g,p,s)
    N,names=nuisance_columns(g,p,s,a)
    directions=np.array([[.04,.02,-.03],[.01,.05,.02],[-.02,.03,.01]])
    B=np.column_stack([response_vector(complete_map_response(a['images_rad'],g,p,s,
        np.array([g.D_ol*row[0],g.D_os*row[1],g.D_ls*row[2]]))) for row in directions])
    C=np.diag([.001**2]*4+[.1**2])
    result=nuisance_projected_information(B,N,C)
    assert 'independent_mass_sheet_log_scale' in names
    assert result['signal_rank']==0
    assert result['null_directions'].shape==(3,3)
    assert 'unconstrained' in result['null_direction_interpretation']


def test_smooth_core_independent_resolve_has_three_noncaustic_images():
    g,p,s,_,_=example();smooth=GalaxyPotential(p.b,core_Mpc=5e-5)
    a=solve_images(g,smooth,s)
    assert len(a['images_rad'])==3
    assert a['minimum_image_Jacobian_singular_value']>.1
    assert a['max_lens_residual_rad']<1e-15


def test_noncommuting_compatible_matrix_maps_keep_full_arrival_response():
    Dol=np.array([[900.,30.],[10.,920.]])
    Dls=np.array([[1100.,20.],[-5.,1050.]])
    T=np.array([[2200.,80.],[80.,2400.]])
    dDol=np.array([[20.,-7.],[4.,-10.]])
    dDls=np.array([[-9.,2.],[3.,15.]])
    dT=np.array([[60.,-8.],[-8.,-40.]])
    def family(t):
        o,l,clock=Dol+t*dDol,Dls+t*dDls,T+t*dT
        source_map=l@np.linalg.solve(o.T,clock)/1.3
        return LensGeometry(o,source_map,l,.3)
    g=family(0.)
    dDos=(dDls@np.linalg.solve(Dol.T,T)
        +Dls@np.linalg.solve(Dol.T,dT-dDol.T@np.linalg.solve(Dol.T,T)))/1.3
    lens=GalaxyPotential(4*np.pi*(220/C_KM_S)**2)
    source=g.D_os@np.array([.2,.05])*ARCSEC
    a=solve_images(g,lens,source)
    response=complete_map_response(a['images_rad'],g,lens,source,
        np.array([dDol,dDos,dDls]))
    np.testing.assert_allclose(response['time_delay_matrix_response'],dT,atol=1e-12)
    h=2e-5
    plus=solve_images(family(h),lens,source)
    minus=solve_images(family(-h),lens,source)
    numeric=(observable_vector(plus)-observable_vector(minus))/(2*h)
    np.testing.assert_allclose(response_vector(response),numeric,rtol=3e-7,atol=2e-8)
