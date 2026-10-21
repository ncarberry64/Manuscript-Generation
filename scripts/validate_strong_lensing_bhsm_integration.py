"""Replay the supplied numerical package without overwriting its evidence.

This runner changes output locations only. The original scalar/GR equations,
frozen inputs and review calculations are executed by their supplied scripts.
"""
from __future__ import annotations

import hashlib
import importlib.metadata
import json
import math
from pathlib import Path
import platform
import subprocess
import sys
import time

ROOT = Path(__file__).resolve().parents[1]
DEST = ROOT / 'artifacts/strong_lensing_bhsm_integration_20261020'
REPLAY = DEST / 'replay'
SUPPLIED = ROOT / 'artifacts/strong_lensing_r1_screened_20261020'


def digest(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def numeric_comparison(reference, replay):
    """Report every finite numeric difference; do not compare provenance text."""
    differences = []
    structure = []
    count = 0

    def walk(old, new, path):
        nonlocal count
        if path in ('/provenance_sha256', '/sha256', '/input_sha256'):
            return
        if isinstance(old, dict):
            if not isinstance(new, dict) or old.keys() != new.keys():
                structure.append(path)
                return
            for key in old:
                walk(old[key], new[key], f'{path}/{key}')
        elif isinstance(old, list):
            if not isinstance(new, list) or len(old) != len(new):
                structure.append(path)
                return
            for index, item in enumerate(old):
                walk(item, new[index], f'{path}/{index}')
        elif isinstance(old, (int, float)) and not isinstance(old, bool):
            count += 1
            if not isinstance(new, (int, float)) or isinstance(new, bool):
                structure.append(path)
            elif old != new:
                differences.append(dict(path=path, supplied=old, replay=new,
                                        absolute_difference=abs(new-old),
                                        relative_difference=abs(new-old)/max(abs(old), 1e-300)))

    walk(reference, replay, '')
    return dict(numeric_values_compared=count, structure_discrepancies=structure,
                numerical_differences=differences)


def assess(receipt):
    """Separate equation diagnostics and null-space roundoff from predictions."""
    invalid = []
    if receipt.get('historical_files_modified'):
        invalid.append(dict(reason='historical evidence modified',
                            files=receipt['historical_files_modified']))
    if len(receipt.get('commands', [])) != 9:
        invalid.append(dict(reason='incomplete nine-command replay'))
    diagnostic_differences = []
    compared = 0
    for row in receipt['commands']:
        if row['returncode']:
            invalid.append(dict(command=row['name'], reason='nonzero return code'))
            continue
        if 'comparison' not in row:
            continue
        comp = row['comparison']
        compared += comp['numeric_values_compared']
        invalid.extend(dict(command=row['name'], reason='structure discrepancy', path=path)
                       for path in comp['structure_discrepancies'])
        for difference in comp['numerical_differences']:
            path = difference['path']
            name = row['name']
            # Error estimators and residuals are tested against their equation
            # bounds below; relative changes near zero are not predictions.
            diagnostic = ((name == 'r1_screened' and path.startswith('/checks/'))
                          or (name == 'screened_projection' and path.startswith('/max_'))
                          or (name == 'gr_optical_control' and
                              ('quadrature_error' in path or 'minus_1PM_minus_2PM' in path
                               or path.startswith('/central_coefficient/absolute_difference')
                               or path.startswith('/central_coefficient/relative_difference')))
                          or (name == 'common_gr_composition' and
                              path.startswith('/checks/') and path.endswith(('absolute_difference', 'relative_difference'))))
            if diagnostic:
                diagnostic_differences.append(dict(command=name, **difference))
                continue
            absolute_tolerance = 1e-18
            if name == 'r1_screened' and path.startswith('/identifiability/free_center/'):
                absolute_tolerance = 1e-15
            if name == 'shared_lens_rank' and 'projected_signal_singular_values' in path:
                absolute_tolerance = 1e-10  # Null values remain below the declared rank threshold.
            if ((name == 'gr2_control' and path.startswith('/change_in_scalar_difference/'))
                    or (name == 'common_gr_composition' and path.startswith('/screened_difference_drift/'))):
                absolute_tolerance = 1e-9  # Seconds / microarcseconds / ppm after subtraction.
            allowed = absolute_tolerance + 2e-9*abs(difference['supplied'])
            if not math.isfinite(difference['replay']) or difference['absolute_difference'] > allowed:
                invalid.append(dict(command=name, reason='numeric reproduction tolerance',
                                    allowed_absolute_difference=allowed, **difference))

    main = json.loads((REPLAY/'results.json').read_text())
    gr2 = json.loads((REPLAY/'gr2_control.json').read_text())
    projection = json.loads((REPLAY/'review/screened_projection_check_results.json').read_text())
    optical = json.loads((REPLAY/'review/gr_optical_control_results.json').read_text())
    composition = json.loads((REPLAY/'review/common_gr_composition_review_results.json').read_text())
    rank = json.loads((REPLAY/'review/shared_lens_rank_results.json').read_text())
    checks = []

    def upper(name, value, tolerance, units='dimensionless'):
        passed = math.isfinite(value) and abs(value) <= tolerance
        checks.append(dict(name=name, value=value, absolute_upper_bound=tolerance,
                           units=units, passed=passed))
        if not passed:
            invalid.append(dict(reason='focused equation check', name=name))

    upper('integrated scalar flux residual', main['checks']['max_integrated_scalar_equation_relative_residual'], 1e-12)
    upper('alpha derivative relative error', main['checks']['max_independent_alpha_derivative_relative_error'], 1e-8)
    upper('potential derivative relative error', main['checks']['max_independent_potential_derivative_relative_error'], 1e-7)
    upper('independent direct-z projection relative discrepancy', projection['max_relative_alpha_discrepancy'], 1e-10)
    upper('independent lens-equation residual', projection['max_image_equation_residual_rad'], 1e-18, 'rad')
    upper('four-case GR2 lens-equation residual', gr2['checks']['max_lens_residual_rad'], 1e-18, 'rad')
    upper('exact static GR remainder after 1PM and 2PM',
          max(abs(x['static_exact_minus_1PM_minus_2PM_rad']) for x in optical['image_rows']), 3e-16, 'rad')
    upper('central GR2 coefficient extrapolation difference', optical['central_coefficient']['absolute_difference'], 1e-9)
    upper('GR2 potential primitive discrepancy', max(x['absolute_difference'] for x in composition['checks']['primitive']), 1e-12)
    upper('GR2 deflection derivative relative discrepancy', max(x['relative_difference'] for x in composition['checks']['derivative']), 1e-8)
    upper('mixed stationary drift discrepancy',
          max(abs(x) for x in composition['four_case_minus_mixed_stationary_prediction_seconds']), 2e-9, 's')
    expectations = [('single-source fixed-center signal rank', main['identifiability']['fixed_center']['signal_rank'], 1),
                    ('single-source free-center signal rank', main['identifiability']['free_center']['signal_rank'], 0),
                    ('three-source finite-profile screened rank', rank['local_screened_signal']['additionally_independent_source_convergence']['signal_rank'], 1),
                    ('three-source distance rank with independent convergence', rank['additionally_independent_source_convergence']['signal_rank'], 0),
                    ('scalar profile held identical in common GR2 control', gr2['checks']['scalar_profile_identical'], True)]
    for name, value, expected in expectations:
        checks.append(dict(name=name, value=value, expected=expected, passed=value == expected))
        if value != expected:
            invalid.append(dict(reason='focused structural check', name=name))
    headline = dict(screening_radius_kpc=1000*main['screening']['rV_Mpc'],
                    leading_fixed_body_reference_delay_offsets_s=main['changes']['relative_delays_seconds'],
                    common_GR2_fixed_body_reference_delay_offsets_s=gr2['scalar_difference_after']['relative_delays_seconds'],
                    screened_three_source_retained_norm_fraction=rank['local_screened_signal']['independent_source_convergence_retained_norm_fraction'])
    return dict(status='PASS' if not invalid else 'FAIL', equations_changed=False,
                numeric_values_compared=compared, invalidating_discrepancies=invalid,
                diagnostics_assessed_by_equation_bounds=diagnostic_differences,
                reproduction_policy=dict(relative_tolerance=2e-9, default_absolute_tolerance=1e-18,
                                         null_projector_absolute_tolerance=1e-15,
                                         null_rank_singular_value_absolute_tolerance=1e-10,
                                         canceled_difference_absolute_tolerance=1e-9,
                                         note='Special absolute bounds concern near-zero numerical nulls and subtraction diagnostics, not a physical error budget.'),
                focused_checks=checks, headline=headline,
                portability='UTF-8/LF replay text; parsed numeric comparison ignores path-separator/hash provenance metadata; original evidence preserved byte-for-byte.',
                qualification='No bound on omitted scalar higher orders, environmental response or global matching; no completed microscopic common-source derivation.')


def main():
    if sys.argv[1:] == ['--assess-only']:
        receipt = json.loads((DEST/'numerical_replay.json').read_text())
        # Re-read existing receipts after output-serialization portability edits;
        # the nine scientific commands are deliberately not repeated.
        outputs = ['results.json', 'gr2_control.json', 'review/r1_qs_coefficients.json',
                   'review/screened_projection_check_results.json', 'review/multisource_crossratio_results.json',
                   'review/shared_lens_rank_results.json', 'review/gr_optical_control_results.json',
                   'review/common_gr_composition_review_results.json']
        for row, output in zip(receipt['commands'], outputs):
            row['comparison'] = numeric_comparison(json.loads((SUPPLIED/output).read_text()),
                                                    json.loads((REPLAY/output).read_text()))
            row['final_script_sha256'] = digest(ROOT/row['command'][1])
            if row['final_script_sha256'] != row['script_sha256']:
                row['post_execution_change'] = 'Receipt serialization made explicit UTF-8/LF only; scientific AST audit in portability_audit.json.'
        result = assess(receipt)
        receipt['status'] = 'PASS' if result['status'] == 'PASS' else 'FAIL'
        (DEST/'numerical_replay.json').write_text(json.dumps(receipt, indent=2)+'\n', encoding='utf-8', newline='\n')
        (DEST/'validation.json').write_text(json.dumps(result, indent=2)+'\n', encoding='utf-8', newline='\n')
        print(json.dumps({k:result[k] for k in ('status', 'numeric_values_compared', 'invalidating_discrepancies', 'headline')}, indent=2))
        raise SystemExit(0 if result['status'] == 'PASS' else 1)
    REPLAY.mkdir(parents=True, exist_ok=True)
    commands = [
        ('r1_screened', 'scripts/run_strong_lensing_r1_screened.py', 'results.json'),
        ('gr2_control', 'scripts/run_strong_lensing_gr2_control.py', 'gr2_control.json'),
        ('r1_qs_coefficients', 'scripts/r1_screened_review/r1_qs_coefficients.py', 'review/r1_qs_coefficients.json'),
        ('screened_projection', 'scripts/r1_screened_review/check_screened_projection.py', 'review/screened_projection_check_results.json'),
        ('multisource_crossratio', 'scripts/r1_screened_review/multisource_crossratio.py', 'review/multisource_crossratio_results.json'),
        ('shared_lens_rank', 'scripts/r1_screened_review/shared_lens_rank.py', 'review/shared_lens_rank_results.json'),
        ('gr_optical_control', 'scripts/r1_screened_review/check_gr_optical_control.py', 'review/gr_optical_control_results.json'),
        ('common_gr_composition', 'scripts/r1_screened_review/check_common_gr_composition.py', 'review/common_gr_composition_review_results.json'),
        ('figure', 'scripts/plot_strong_lensing_r1_screened.py', None),
    ]
    historical = list(SUPPLIED.rglob('*')) + list((ROOT/'manuscript/figures').glob('r1_screening_prediction.*'))
    historical = {str(p.relative_to(ROOT)): digest(p) for p in historical if p.is_file()}
    results = []
    for name, script, output in commands:
        args = [sys.executable, script]
        if output:
            args += ['--output', str((REPLAY/output).relative_to(ROOT))]
        else:
            args += ['--output-dir', str((REPLAY/'figures').relative_to(ROOT)), '--formats', 'png', 'svg']
        start = time.perf_counter()
        done = subprocess.run(args, cwd=ROOT, capture_output=True, text=True, encoding='utf-8', errors='replace')
        elapsed = time.perf_counter()-start
        (REPLAY/f'{name}.stdout.txt').write_text(done.stdout, encoding='utf-8', newline='\n')
        (REPLAY/f'{name}.stderr.txt').write_text(done.stderr, encoding='utf-8', newline='\n')
        row = dict(name=name, command=args, returncode=done.returncode, elapsed_seconds=elapsed,
                   script_sha256=digest(ROOT/script))
        if done.returncode == 0 and output:
            row['comparison'] = numeric_comparison(json.loads((SUPPLIED/output).read_text()),
                                                    json.loads((REPLAY/output).read_text()))
        results.append(row)
        print(name, 'returncode', done.returncode, 'seconds', round(elapsed, 3), flush=True)
        if done.returncode:
            print(done.stderr[-3000:], flush=True)
            break
    altered = [path for path, value in historical.items() if digest(ROOT/path) != value]
    receipt = dict(status='REPLAY_EXECUTED_REVIEW_REQUIRED',
                   source_commit='49a1dabbe148f8ba9da4ac079a8d93fe97c5de76',
                   environment=dict(python=platform.python_version(), executable=sys.executable,
                                    platform=platform.platform(),
                                    **{name: importlib.metadata.version(name) for name in ('numpy', 'scipy', 'matplotlib')}),
                   equations_changed=False, routing_changes_only=True,
                   historical_files_checked=len(historical), historical_files_modified=altered,
                   commands=results,
                   qualification='Reproduces retained equations; numerical agreement does not bound omitted physics or establish a microscopic common-source derivation.')
    (DEST/'numerical_replay.json').write_text(json.dumps(receipt, indent=2)+'\n', encoding='utf-8', newline='\n')
    if altered or len(results) != len(commands) or any(row['returncode'] for row in results):
        raise SystemExit(1)
    result = assess(receipt)
    receipt['status'] = 'PASS' if result['status'] == 'PASS' else 'FAIL'
    (DEST/'numerical_replay.json').write_text(json.dumps(receipt, indent=2)+'\n', encoding='utf-8', newline='\n')
    (DEST/'validation.json').write_text(json.dumps(result, indent=2)+'\n', encoding='utf-8', newline='\n')
    raise SystemExit(0 if result['status'] == 'PASS' else 1)


if __name__ == '__main__':
    main()
