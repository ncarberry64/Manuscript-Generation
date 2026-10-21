"""Apply nonlinear conventional nuisance fits to the saved R1 injection."""
from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path
import sys

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT/'code'))
from strong_lensing_local_family import saved_geometry
from strong_lensing_refit import fit_conventional


def encode(value):
    if isinstance(value, np.ndarray):
        return value.tolist()
    if isinstance(value, np.generic):
        return value.item()
    raise TypeError(type(value).__name__)


def evaluate():
    injection_path = ROOT/'artifacts/strong_lensing_r1_screened_20261020/results.json'
    reference_path = ROOT/'artifacts/strong_lensing_first_comparison_20261020/results.json'
    injection = json.loads(injection_path.read_text(encoding='utf-8'))
    _, geometry, source = saved_geometry(reference_path)
    image_angles = np.asarray(injection['screened']['images_rad'])[:, 0]
    delays = np.asarray(injection['screened']['relative_delays_days'])
    mass, radius = injection['body']['mass_Msun'], injection['body']['radius_Mpc']
    fits = {name: fit_conventional(geometry, mass, radius, source, image_angles, delays,
                                  free_center=free_center)
            for name, free_center in [('fixed_center', False), ('free_center', True)]}
    free = fits['free_center']
    passed = (free['optimizer_success'] and np.max(abs(free['image_errors_microarcsec'])) < 1e-5
              and np.max(abs(free['delay_errors_seconds'])) < 1e-5
              and free['forward']['max_original_lens_residual_rad'] < 1e-17)
    files = [injection_path, reference_path, ROOT/'code/strong_lensing_refit.py', Path(__file__)]
    return dict(status='PASS' if passed else 'FAIL',
                experiment='Deterministic R1 injection read from the existing receipt; finite nonlinear conventional refit, not an observational detection or likelihood.',
                preparation='Truth uses unchanged physical source, spectrum, clocks, redshifts/maps and frozen R1 state. Only inverse-model nuisance coordinates are fitted.',
                coordinates=['delta_log_M', 'delta_log_R', 'source_shift_arcsec',
                             'lens_center_arcsec', 'source_plane_convergence'],
                original_applications=['theta-beta-G alpha=0 at all three images',
                                       'relative stationary Fermat arrival times for the two image pairs'],
                numerical_row_scales='0.001 arcsec lens residual; 0.1 day delay residual. Numerical conditioning only, not observational uncertainty.',
                fits=fits,
                raw_fixed_body_offsets_seconds=injection['changes']['relative_delays_seconds'],
                mass_sheet_accounting='One physical sheet is included in the conventional encounter. Propagation maps remain the original maps and contain no duplicate sheet.',
                qualification='Local three-image axial inverse problem for this finite profile. A successful finite refit does not establish equivalence of every profile or higher-order metric; the separate positive-halo theorem supplies the complete retained-field equivalence.',
                sha256={p.relative_to(ROOT).as_posix():hashlib.sha256(p.read_bytes()).hexdigest() for p in files})


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output', type=Path, default=ROOT/'artifacts/strong_lensing_acceptance_20261020/nonlinear_refit.json')
    args = parser.parse_args()
    result = evaluate()
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(result, default=encode, indent=2, allow_nan=False)+'\n',
                           encoding='utf-8', newline='\n')
    print(json.dumps(result, default=encode, indent=2))
    raise SystemExit(0 if result['status']=='PASS' else 1)


if __name__ == '__main__':
    main()
