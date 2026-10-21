"""Eight predeclared profile/input controls; preserve all historical results."""
from __future__ import annotations

import argparse
from dataclasses import asdict
import hashlib
import json
from pathlib import Path
import sys

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT/'code'))
from strong_lensing_local_family import saved_geometry
from strong_lensing_r1_screened import R1Coefficients
from strong_lensing_robustness import BASE_HALF_MASS_MPC, TOLMAN_HALF_X, evaluate_case


def encode(value):
    if isinstance(value, np.ndarray):
        return value.tolist()
    if isinstance(value, np.generic):
        return value.item()
    raise TypeError(type(value).__name__)


def evaluate():
    reference = ROOT/'artifacts/strong_lensing_first_comparison_20261020/results.json'
    _, geometry, source = saved_geometry(reference)
    frozen = ROOT/'artifacts/strong_lensing_r1_screened_20261020/results.json'
    prior = json.loads(frozen.read_text(encoding='utf-8'))
    coefficients = R1Coefficients(**prior['coefficients'])
    design = [('reference', {}), ('mass_plus_20_percent', {'mass_factor':1.2}),
              ('size_plus_20_percent', {'size_factor':1.2}),
              ('source_plus_20_percent', {'source_factor':1.2})]
    cases = []
    for profile in ('tolman_vii', 'plummer'):
        for name, variation in design:
            case = evaluate_case(geometry, coefficients, source, profile, detailed=name == 'reference', **variation)
            case['case'] = name
            cases.append(case)
            print(profile, name, 'images', len(case['conventional']['images_rad']),
                  'offsets_s', case['reference_relative_offsets']['delays_seconds'], flush=True)
    summaries = []
    baseline = {case['profile']:case for case in cases if case['case'] == 'reference'}
    for case in cases:
        base = baseline[case['profile']]
        ordinary = (case['conventional']['relative_delays_days']-base['conventional']['relative_delays_days'])*86400
        scalar = case['reference_relative_offsets']['delays_seconds']
        summaries.append(dict(profile=case['profile'], case=case['case'],
                              ordinary_delay_change_from_same_profile_baseline_seconds=ordinary,
                              R1_reference_relative_delay_offsets_seconds=scalar,
                              ordinary_change_to_R1_offset_max_norm_ratio=float(np.max(abs(ordinary))/np.max(abs(scalar)))))
    profile_change = (baseline['plummer']['conventional']['relative_delays_days']-baseline['tolman_vii']['conventional']['relative_delays_days'])*86400
    profile_offset_change = (baseline['plummer']['reference_relative_offsets']['delays_seconds']-baseline['tolman_vii']['reference_relative_offsets']['delays_seconds'])
    files = [reference, frozen, ROOT/'code/strong_lensing_robustness.py', Path(__file__),
             ROOT/'code/strong_lensing_r1_screened.py', ROOT/'code/strong_lensing_local_family.py',
             ROOT/'code/strong_lensing_comparison.py']
    return dict(status='BOUNDED_PROFILE_ROBUSTNESS_STUDY_EXECUTED',
                scope='same state-conditioned light transfer, eight paired applications; controlled inputs are not fitted uncertainties or empirical priors',
                action_state=dict(coefficients=asdict(coefficients), changed=False, new_scalar_amplitude=False),
                shared_preparation=dict(source_Mpc=source, source_change_only_in_declared_control=True,
                                        clocks_redshifts_and_distance_maps='identical within every response/reference pair and throughout the study'),
                profile_matching=dict(total_mass_Msun=1e11, half_mass_radius_kpc=BASE_HALF_MASS_MPC*1000,
                                      tolman_half_radius_over_support=TOLMAN_HALF_X,
                                      convention='same total 3D mass and 3D half-mass radius; Plummer a=r_half sqrt(2^(2/3)-1)',
                                      alternative='standard Plummer sphere: rho=3 M/(4 pi a³) (1+r²/a²)^(-5/2), regular center, finite total mass, no hard material edge',
                                      citation='Plummer1911, MNRAS71(5),460-470,doi:10.1093/mnras/71.5.460'),
                design=[dict(case=name, **variation) for name, variation in design],
                cases=cases, within_profile_sensitivity=summaries,
                cross_profile_baseline=dict(conventional_delay_difference_seconds=profile_change,
                                            R1_offset_difference_seconds=profile_offset_change,
                                            ordinary_profile_change_to_max_R1_offset_ratio=float(np.max(abs(profile_change))/max(np.max(abs(x['reference_relative_offsets']['delays_seconds'])) for x in baseline.values()))),
                interpretation='The positive effective-mass representation survives the profile change; ordinary profile/source freedoms alter delays much more than the action-owned offsets in this bounded comparison. This is model sensitivity, not evidence of a detection or a prior on a real lens.',
                approximation='leading weak metric, static quasi-static scalar branch, geometric optics thin plane; no GR2 or environmental equation is newly solved',
                physical_admission=False,
                sha256={p.relative_to(ROOT).as_posix():hashlib.sha256(p.read_bytes()).hexdigest() for p in files})


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--output', type=Path, default=ROOT/'artifacts/strong_lensing_acceptance_20261020/robustness.json')
    args = parser.parse_args()
    result = evaluate()
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(result, default=encode, indent=2, allow_nan=False)+'\n', encoding='utf-8', newline='\n')
    print(json.dumps(result['cross_profile_baseline'], default=encode, indent=2))


if __name__ == '__main__':
    main()
