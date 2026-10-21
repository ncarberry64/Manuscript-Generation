"""Re-solve the matched lens pair with the same finite-body GR 2PM term.

This is a common-conventional-term robustness control, not a complete
second-order scalar-tensor or finite-distance cosmological calculation.
"""
from __future__ import annotations

import argparse
from dataclasses import asdict
import hashlib
import json
from pathlib import Path
import sys

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / 'code'))

import numpy as np
from strong_lensing_comparison import solve_images, complete_map_response
from strong_lensing_local_family import TolmanVIILens, saved_geometry
from strong_lensing_r1_screened import R1Coefficients, R1ScreenedLens, solve_axis_images
from strong_lensing_gr2 import TolmanVII2PMLens


def encode(value):
    if isinstance(value, np.ndarray):
        return value.tolist()
    if isinstance(value, np.generic):
        return value.item()
    raise TypeError(type(value).__name__)


def difference(scalar, conventional):
    scalar_ratio = scalar['relative_delays_days'][0] / scalar['relative_delays_days'][1]
    conventional_ratio = conventional['relative_delays_days'][0] / conventional['relative_delays_days'][1]
    return dict(
        images_microarcsec=(scalar['images_arcsec'] - conventional['images_arcsec']) * 1e6,
        relative_delays_seconds=(scalar['relative_delays_days'] - conventional['relative_delays_days']) * 86400,
        conventional_delay_ratio=conventional_ratio,
        scalar_delay_ratio=scalar_ratio,
        delay_ratio_fractional_change_ppm=(scalar_ratio / conventional_ratio - 1) * 1e6,
    )


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--reference', type=Path, default=ROOT / 'artifacts/strong_lensing_first_comparison_20261020/results.json')
    parser.add_argument('--output', type=Path, default=ROOT / 'artifacts/strong_lensing_r1_screened_20261020/gr2_control.json')
    args = parser.parse_args()
    old, geometry, source = saved_geometry(args.reference)
    coefficients = R1Coefficients.frozen(geometry.z_l, old['experiment']['H0_km_s_Mpc'])
    base = TolmanVIILens()
    corrected_base = TolmanVII2PMLens(base)
    scalar1 = R1ScreenedLens(base, coefficients)
    scalar2 = R1ScreenedLens(corrected_base, coefficients)
    gr1 = solve_images(geometry, base, source)
    r1 = solve_axis_images(geometry, scalar1, source, gr1['images_rad'])
    gr2 = solve_images(geometry, corrected_base, source, seeds=gr1['images_rad'])
    r2 = solve_axis_images(geometry, scalar2, source, gr2['images_rad'])
    before = difference(r1, gr1)
    after = difference(r2, gr2)
    zero_maps = np.zeros((3, 2, 2))
    def correction_partial(x):
        p, a, _ = corrected_base.correction.extra_application(x)
        return p, a
    stationary = complete_map_response(gr1['images_rad'], geometry, base, source,
                                      zero_maps, intrinsic_encounter=correction_partial)
    conventional_change = difference(gr2, gr1)
    impacts = np.linalg.norm(gr1['images_rad'] @ geometry.D_ol.T, axis=1)
    files = [Path(__file__), ROOT/'code/strong_lensing_gr2.py',
             ROOT/'code/strong_lensing_r1_screened.py', args.reference]
    results = dict(
        status='COMMON_FINITE_BODY_GR_2PM_ROBUSTNESS_CONTROL',
        source_commit='49a1dabbe148f8ba9da4ac079a8d93fe97c5de76',
        preparation='Same prescribed density, incident light, physical source, measured redshifts, clock, compatible maps, frozen R1 coefficients and leading scalar profile.',
        coefficients=asdict(coefficients),
        body=dict(mass_Msun=base.mass_Msun, radius_Mpc=base.radius_Mpc,
                  compactness=corrected_base.correction.epsilon),
        conventional_1pm=gr1, scalar_on_1pm=r1,
        conventional_1pm_plus_2pm=gr2, scalar_on_common_gr2=r2,
        scalar_difference_before=before, scalar_difference_after=after,
        change_in_scalar_difference=dict(
            images_microarcsec=after['images_microarcsec']-before['images_microarcsec'],
            relative_delays_seconds=after['relative_delays_seconds']-before['relative_delays_seconds'],
            delay_ratio_fractional_change_ppm=after['delay_ratio_fractional_change_ppm']-before['delay_ratio_fractional_change_ppm']),
        conventional_common_correction=conventional_change,
        scalar_on_common_correction=difference(r2,r1),
        correction_at_original_images=dict(impact_Mpc=impacts,
            radial_alpha_rad=[corrected_base.correction.extra_alpha(r) for r in impacts],
            potential_Mpc=[corrected_base.correction.extra_potential(r) for r in impacts]),
        checks=dict(
            max_lens_residual_rad=max(item['max_lens_residual_rad'] for item in (gr1,r1,gr2,r2)),
            stationary_conventional_delay_change_seconds=stationary['relative_delays_days']*86400,
            finite_minus_stationary_conventional_delay_seconds=conventional_change['relative_delays_seconds']-stationary['relative_delays_days']*86400,
            scalar_profile_identical=bool(scalar1.rV_Mpc == scalar2.rV_Mpc and scalar1.coefficients == scalar2.coefficients)),
        scope=[
            'The ordinary-GR second post-Minkowskian scattering deflection is derived from the same static Tolman-VII body, including its leading TOV pressure and impact-coordinate conversion.',
            'The impact coordinate is the invariant asymptotic scattering impact, identified with the physical thin-plane impact in this control.',
            'The thin-plane potential is the radial primitive of the deflection, and the same consistent gradient/Hessian and arrival functional are used for both models.',
            'The scalar spatial profile is held at its declared leading quasi-static order. Scalar-GR second-order cross terms and scalar corrections to the 2PM fluid geometry are not calculated.',
            'This control is not a complete finite-distance 2PM cosmological arrival calculation or a bound on all physical errors.',
            'The positive effective-halo equivalence and nuisance-rank results in the manuscript are leading-order results; this common control does not claim to break them.'
        ],
        provenance_sha256={str(p.relative_to(ROOT)):hashlib.sha256(p.read_bytes()).hexdigest() for p in files})
    if not results['checks']['scalar_profile_identical']:
        raise ArithmeticError('scalar profile changed in a common-GR control')
    if results['checks']['max_lens_residual_rad'] > 1e-17:
        raise ArithmeticError('image root residual exceeds focused control bound')
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(results, indent=2, default=encode)+'\n', encoding='utf-8', newline='\n')
    print(json.dumps({k:results[k] for k in ('status','scalar_difference_after',
          'change_in_scalar_difference','conventional_common_correction','checks')}, indent=2, default=encode))
    print('output', args.output)


if __name__ == '__main__':
    main()
