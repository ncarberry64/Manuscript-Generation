"""Focused independent static-geodesic validation of the finite-body 2PM term."""
from pathlib import Path
import hashlib
import json
import math
import sys

import numpy as np

HERE=Path(__file__).resolve().parent
REPO=HERE.parents[1]
sys.path.insert(0,str(REPO/'code'))
from strong_lensing_local_family import TolmanVII,TolmanVIILens
from strong_lensing_gr2 import (TolmanVII2PMCorrection,coefficient_I,
                           exact_static_deflection)


def run(output=None):
    source=REPO/'artifacts/strong_lensing_r1_screened_20261020/results.json'
    original=json.loads(source.read_text())
    base=TolmanVIILens()
    body=TolmanVII()
    correction=TolmanVII2PMCorrection(base)
    rows=[]
    for b in original['image_impact_radii_Mpc']:
        p,alpha,_=base.application([b,0.])
        extra=correction.extra_alpha(b)
        exact,error=exact_static_deflection(body,b)
        rows.append(dict(impact_Mpc=b,impact_over_R=b/base.radius_Mpc,
            alpha_1PM_rad=float(alpha[0]),alpha_2PM_rad=extra,
            alpha_static_exact_rad=exact,static_quadrature_error_rad=error,
            static_exact_minus_1PM_rad=exact-alpha[0],
            static_exact_minus_1PM_minus_2PM_rad=exact-alpha[0]-extra,
            extra_potential_Mpc=correction.extra_potential(b),
            extra_alpha_prime_per_Mpc=correction.extra_alpha_prime(b),
            correction_over_1PM=extra/alpha[0]))

    # The central ray traverses the fluid. Vary compactness only to independently
    # isolate the Taylor coefficient from the exact TOV/null-geodesic calculation.
    b=original['image_impact_radii_Mpc'][1]
    central=[]
    for e in [5e-4,2.5e-4,1.25e-4]:
        mass=base.mass_Msun*e/body.compactness
        trial=TolmanVII(mass_Msun=mass,radius_Mpc=base.radius_Mpc)
        weak=TolmanVIILens(mass_Msun=mass,radius_Mpc=base.radius_Mpc)
        exact,error=exact_static_deflection(trial,b)
        alpha1=weak.application([b,0.])[1][0]
        central.append(dict(compactness=e,coefficient=(exact-alpha1)/e**2,
                            quadrature_error_in_coefficient=error/e**2))
    extrapolated=float(np.polynomial.polynomial.polyfit(
        [r['compactness'] for r in central],[r['coefficient'] for r in central],2)[0])
    target=2*coefficient_I(b/base.radius_Mpc)

    # An independent exterior identity, evaluated at both outer image impacts.
    exterior=[]
    for i in [0,2]:
        b=original['image_impact_radii_Mpc'][i]
        expected=15*math.pi/4*(base.mass_length_Mpc/b)**2
        actual=correction.extra_alpha(b)
        exterior.append(dict(image_index=i,relative_error=(actual-expected)/expected))
    result=dict(
        classification='ordinary-GR isolated-static 2PM control in the common thin-plane impact convention',
        body=dict(mass_Msun=base.mass_Msun,radius_Mpc=base.radius_Mpc,
                  mass_length_Mpc=base.mass_length_Mpc,compactness=body.compactness),
        image_rows=rows,
        central_exact_compactness_sequence=central,
        central_coefficient=dict(analytic=target,extrapolated_exact=extrapolated,
            absolute_difference=extrapolated-target,
            relative_difference=(extrapolated-target)/target),
        exterior_coefficient_checks=exterior,
        method='exact optical-radius scattering integral with expm1/log1p flat subtraction; independent analytic 2PM expansion includes Newtonian TOV pressure and impact-to-areal conversion',
        limits=['b is asymptotic invariant impact, identified with the existing physical thin-plane coordinate',
                'this does not derive a scalar 2PM response',
                'the projected primitive is an integrable thin-plane generating potential, not an exact cosmological time-transfer function',
                'quadrature error omits underlying TOV and floating-point uncertainty; compactness sequence supplies an independent coefficient check'],
        sha256={str(p.relative_to(REPO)):hashlib.sha256(p.read_bytes()).hexdigest()
                for p in [source,REPO/'code/strong_lensing_local_family.py',
                          REPO/'code/strong_lensing_gr2.py',Path(__file__)]})
    output = output or REPO/'artifacts/strong_lensing_r1_screened_20261020/review/gr_optical_control_results.json'
    output.parent.mkdir(parents=True, exist_ok=True)
    output.write_text(json.dumps(result,indent=2,allow_nan=False)+'\n', encoding='utf-8', newline='\n')
    print(json.dumps({'image_rows':rows,'central_coefficient':result['central_coefficient'],
                      'exterior_coefficient_checks':exterior},indent=2))
    print(output)


if __name__ == '__main__':
    import argparse
    parser = argparse.ArgumentParser()
    parser.add_argument('--output', type=Path, help='Write a separate replay receipt; default preserves the original command behavior.')
    args = parser.parse_args()
    run(args.output)
