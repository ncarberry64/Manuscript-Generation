from copy import deepcopy
from decimal import Decimal, localcontext
import json
from pathlib import Path
import numpy as np
import pytest
from bhtracker import fit
from bhtracker.examples import make_example
from bhtracker.forward import C,GM_SUN,PC
from bhtracker.interpretation import spherical_boundary,bhsm_interpretation
from bhtracker.report import report_html


def test_boundary_map_has_the_same_mass_at_different_enclosing_radii():
    mass=1e5
    radii=[1e-5,.01,100.]
    energies=[]
    for radius in radii:
        actual=spherical_boundary(mass,radius,2.)
        with localcontext() as ctx:
            ctx.prec=65
            m=Decimal(str(GM_SUN))*Decimal(str(mass))/(Decimal(str(C))**2*Decimal(str(PC)))
            r=Decimal(str(radius))
            expected=(1-(1-2*m/r).sqrt())/r
        assert actual['curvature_deficit_pc_inverse']==pytest.approx(float(expected),rel=2e-15)
        assert actual['reconstructed_monopole_msun']==pytest.approx(mass,rel=2e-15)
        energies.append(actual['boundary_energy_over_Msun_c2'])
    assert energies[0]>energies[1]>energies[2]


def test_boundary_errors_are_fixed_radius_pushforwards():
    mass=1e5;radius=1e-5;step=.1
    actual=spherical_boundary(mass,radius,3.)
    plus=spherical_boundary(mass+step,radius);minus=spherical_boundary(mass-step,radius)
    for value,error in [('boundary_energy_over_Msun_c2','boundary_energy_sigma_over_Msun_c2'),
                        ('curvature_deficit_pc_inverse','curvature_deficit_sigma_pc_inverse')]:
        derivative=(plus[value]-minus[value])/(2*step)
        assert actual[error]==pytest.approx(3*derivative,rel=2e-8)


@pytest.mark.parametrize('mass,radius',[(0,1),(1,0),(float('nan'),1),(1,float('inf')),(1e12,1e-10)])
def test_boundary_rejects_invalid_static_geometry(mass,radius):
    with pytest.raises(ValueError):spherical_boundary(mass,radius)


def test_gr_fit_is_preserved_and_has_its_bhsm_interpretation():
    raw,_=make_example('nearby_motion');original=deepcopy(raw)
    result=fit(raw,starts=2)
    stored=json.loads((Path(__file__).resolve().parent/'fixtures/pre_interpretation_fit.json').read_text())
    np.testing.assert_allclose(result['predictions'],stored['predictions'],rtol=2e-11,atol=1e-12)
    assert result['estimates']['mass_msun']==pytest.approx(stored['estimates']['mass_msun'],rel=2e-11)
    assert result['status']==stored['status']=='TRACKED'
    assert raw==original
    meaning=result['bhsm_interpretation']
    assert meaning['framework']=='BHSM' and meaning['spherical_boundary']['status']=='AVAILABLE'
    assert meaning['spherical_boundary']['reconstructed_monopole_msun']==pytest.approx(result['estimates']['mass_msun'])
    assert 'BHSM interpretation' in report_html(result)


def test_r1_mass_is_not_silently_given_schwarzschild_boundary_energy():
    meaning=bhsm_interpretation('r1',{'mass_msun':1e5},{'mass_msun':3.},{'r0_pc':.01},'TRACKED')
    assert meaning['spherical_boundary']['status']=='NOT_ASSUMED_FOR_R1_METRIC'
    assert 'Action metric/source' in meaning['parameter_meanings'][0]['interpretation']


def test_failed_mass_fit_does_not_export_an_interpreted_boundary_measurement():
    for status in ['FIT_REJECTED','INSUFFICIENT_INFORMATION']:
        meaning=bhsm_interpretation('gr',{'mass_msun':1e5},{'mass_msun':3.},{'r0_pc':.01},status)
        assert meaning['spherical_boundary']['status']!='AVAILABLE'
