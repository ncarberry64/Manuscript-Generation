"""Portable joint-research exports; numerical replay is distinct from validation.

No reported solver mass is promoted to a black-hole classification here.
An exported internal state can be replayed without repeating optimization;
profiles, noise controls and systematic checks remain separate evidence.
"""
import base64
import csv
import hashlib
import html
import json
from pathlib import Path

import numpy as np

from .provenance import REPOSITORY_URL, project_citation

BUNDLE_SCHEMA = 'bhtracker.joint_analysis_bundle.v1'


def is_joint_result(result):
    return isinstance(result, dict) and result.get('schema') in {
        'bhtracker.native_joint_mass_result.v1', 'bhtracker.native_joint_mass_result.v2'}


def validate_result(result):
    """Check saved predictions AND physical summaries against the recorded state."""
    if not is_joint_result(result) or 'dataset_input' not in result:
        raise ValueError('A completed versioned joint result with its native input is required')
    if result['schema'].endswith('.v2'):
        from .joint_mass_v2 import JointData, JointModel
        from .joint_workflow_v2 import replay
    else:
        from .joint_mass import JointData, JointModel
        from .joint_workflow import replay
    checked = replay(result)
    if not checked['passed']:
        raise ValueError('Joint state replay differs; retain the matching scientific source and unmodified result')
    dataset = JointData(result['dataset_input'])
    model = JointModel(dataset, **result['model_options'])
    actual = model.evaluate(np.asarray(result['best_internal_state']), details=True)

    def same(a, b):
        if isinstance(b, dict):
            return isinstance(a, dict) and set(a) == set(b) and all(same(a[k], b[k]) for k in b)
        if b is None:
            return a is None
        if isinstance(b, (float, int, np.number)):
            return not isinstance(a, bool) and isinstance(a, (int, float, np.number)) and bool(np.isclose(a, b, rtol=2e-10, atol=2e-8))
        return a == b

    if not same(result.get('parameters'), actual['parameters']):
        raise ValueError('Reported joint physical parameters differ from the saved state')
    for key in ('photometric_residuals_magnitude', 'whitened_photometric_residuals',
                'astrometric_residuals_NE_mas', 'whitened_astrometric_residuals'):
        if not np.allclose(np.asarray(result[key]), actual[key], rtol=2e-10, atol=2e-8):
            raise ValueError(f'Reported joint {key} differ from the native data and saved state')
    for name, source in dataset.document.get('native_files', {}).items():
        raw = base64.b64decode(source['raw_bytes_base64'], validate=True)
        if hashlib.sha256(raw).hexdigest() != source['sha256']:
            raise ValueError(f'Native source byte hash mismatch: {name}')
    checked['same_physical_parameters'] = True
    checked['same_native_residuals'] = True
    checked['method'] = 'FORWARD_STATE_REPLAY; optimization and scientific acceptance are separate'
    return dataset, model, checked


def report_html(result):
    dataset, _, checked = validate_result(result)
    esc = lambda value: html.escape(str(value))
    parameters = ''.join(f'<tr><th>{esc(k)}</th><td>{esc(v)}</td></tr>' for k, v in result['parameters'].items() if not isinstance(v, dict))
    return f'''<!doctype html><html lang="en"><meta charset="utf-8"><meta name="viewport" content="width=device-width">
<title>BHSM Blackhole Tracker · Joint astrometry research</title>
<style>body{{font:16px/1.6 system-ui;max-width:1050px;margin:3rem auto;padding:0 1.5rem;color:#172b40}}table{{border-collapse:collapse;width:100%}}td,th{{padding:.5rem;text-align:left;border-bottom:1px solid #dce3e9}}pre{{white-space:pre-wrap;overflow-wrap:anywhere;background:#f1f5f8;padding:1rem}}.status{{background:#fff3d8;padding:1rem}}</style>
<h1>BHSM Blackhole Tracker · Joint research</h1><p>Formal source and citation: <a href="{REPOSITORY_URL}">BHSM GitHub repository</a>.</p>
<p class="status">{esc(result['mass_status'])}. Solver mass is a conditional research estimate. This export does not certify a black hole, a discovery, noise adequacy or confidence-interval coverage.</p>
<p>Target: {esc(dataset.document.get('target', 'not supplied'))}. {len(dataset.document['photometry'])} native photometric rows; {len(dataset.document['astrometry'])} native angular rows. Selected rows, exclusions, native bytes, time assumptions and frame provenance are retained in dataset.json.</p>
<h2>What the calculation measures</h2><p>Under the compact-lens model, θE² = κ M πrel and πE = πrel/θE, so M = θE/(κ πE); μrel = θE/tE. Angular centroid data and annual photometric parallax must independently constrain the two scales. Duration alone is not a mass measurement. A boundary-limited or degenerate fit leaves mass unidentified.</p>
<h2>BHSM interpretation</h2><p>One ordinary compact-lens optical potential represents the BHSM Einstein/GHY spacetime-side response, counted once. This calculation supplies no separate optical term for the proposed nongeometric core, envelope or membrane.</p><p>Absolute lens motion, absolute source distance, lens distance, horizon and core state remain unavailable unless independent frame, distance and dark-lens evidence are supplied. A conditional mass above a remnant threshold does not identify a horizon.</p>
<h2>Conditional solver parameters</h2><table>{parameters}</table>
<h2>Separate acceptance checks</h2><p>Inspect all starts, fixed-mass profiles, both tails and dependence of noise residuals, timing and calibration sensitivity, lens-light assumptions and prior literature. A scalar goodness-of-fit probability or a successful replay is insufficient. The benchmark validation receipt is separate from this individual fit.</p>
<pre>{esc(json.dumps(result.get('diagnostics'), indent=2))}</pre><h2>Model options</h2><pre>{esc(json.dumps(result['model_options'], indent=2))}</pre>
<h2>Native clocks and frame</h2><pre>{esc(json.dumps(dataset.document['metadata'], indent=2))}</pre>
<h2>Numerical replay</h2><pre>{esc(json.dumps(checked, indent=2))}</pre><p>Forward replay evaluates the saved state with the exact scientific source. It does not repeat the optimizer or establish scientific reliability. See REPRODUCE.md for refitting and profile commands.</p></html>'''


def export_result(result, directory):
    dataset, model, checked = validate_result(result)
    dest = Path(directory)
    dest.mkdir(parents=True, exist_ok=True)
    for name, value in [('dataset.json', dataset.document), ('result.json', result),
                        ('forward_replay_check.json', checked), ('model_options.json', result['model_options'])]:
        (dest / name).write_text(json.dumps(value, indent=2, allow_nan=False) + '\n', encoding='utf-8')
    (dest / 'report.html').write_text(report_html(result), encoding='utf-8')
    with (dest / 'photometry.csv').open('w', newline='', encoding='utf-8') as stream:
        writer = csv.writer(stream)
        writer.writerow(['source_row_id', 'instrument', 'JD_TDB_model_epoch', 'native_magnitude', 'supplied_sigma_mag', 'model_magnitude', 'residual_mag', 'whitened_innovation'])
        for i, row in enumerate(dataset.document['photometry']):
            writer.writerow([row['source_row_id'], row['instrument'], model.tp[i], row['magnitude'], row['magnitude_error'], result['photometric_predictions_magnitude'][i], result['photometric_residuals_magnitude'][i], result['whitened_photometric_residuals'][i]])
    epochs = model.average(model.ta[:, None])[:, 0]
    with (dest / 'astrometry.csv').open('w', newline='', encoding='utf-8') as stream:
        writer = csv.writer(stream)
        writer.writerow(['source_row_id', 'filter', 'JD_TDB_model_mean_epoch', 'native_north_mas', 'native_east_mas', 'supplied_sigma_north_mas', 'supplied_sigma_east_mas', 'model_north_mas', 'model_east_mas'])
        for i, row in enumerate(dataset.document['astrometry']):
            writer.writerow([row['source_row_id'], row['filter'], epochs[i], row['north_mas'], row['east_mas'], row['north_error_mas'], row['east_error_mas'], *result['astrometric_predictions_NE_mas'][i]])
    native_index=[]
    for i, (name, source) in enumerate(sorted(dataset.document.get('native_files', {}).items()), 1):
        filename=f'native_source_{i:02}.bin'
        (dest / filename).write_bytes(base64.b64decode(source['raw_bytes_base64'], validate=True))
        native_index.append({'exported_filename':filename,'original_source_name':name,'sha256':source['sha256']})
    (dest / 'native_sources.json').write_text(json.dumps(native_index,indent=2)+'\n',encoding='utf-8')
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    fig, axes = plt.subplots(2, 2, figsize=(12, 8), layout='constrained')
    for instrument in model.names:
        indices = np.flatnonzero(np.array([r['instrument'] == instrument for r in model.pr]))
        indices = indices[np.argsort(model.tp[indices], kind='stable')]
        axes[0, 0].errorbar(model.tp[indices] - 2450000, model.yp[indices], yerr=model.sp[indices], fmt='.', ms=1.5, alpha=.4, label=instrument)
        axes[0, 0].plot(model.tp[indices] - 2450000, np.array(result['photometric_predictions_magnitude'])[indices], lw=.7)
        axes[1, 0].plot(model.tp[indices] - 2450000, np.array(result['whitened_photometric_residuals'])[indices], '.', ms=1.5, alpha=.5, label=instrument)
    axes[0, 0].invert_yaxis()
    axes[0, 0].set(ylabel='Native magnitude', xlabel='JD TDB − 2450000')
    axes[0, 0].legend(fontsize=7)
    axes[1, 0].set(ylabel='Whitened innovation', xlabel='JD TDB − 2450000')
    axes[1, 0].axhline(0, color='.5', lw=.7)
    for k, label in enumerate(['North', 'East']):
        observed = np.array([[r['north_mas'], r['east_mas']] for r in model.ar])[:, k]
        sigma = np.array([[r['north_error_mas'], r['east_error_mas']] for r in model.ar])[:, k]
        predicted = np.array(result['astrometric_predictions_NE_mas'])[:, k]
        axes[0, 1].errorbar(epochs - 2450000, observed, yerr=sigma, fmt='.', label=label)
        axes[0, 1].plot(epochs - 2450000, predicted, 'x', ms=3)
        axes[1, 1].errorbar(epochs - 2450000, observed - predicted, yerr=sigma, fmt='.', label=label)
    axes[0, 1].set(ylabel='Relative angular position (mas)', xlabel='JD TDB − 2450000')
    axes[1, 1].set(ylabel='Angular residual (mas)', xlabel='JD TDB − 2450000')
    axes[0, 1].legend(fontsize=8)
    axes[1, 1].axhline(0, color='.5', lw=.7)
    fig.suptitle('BHSM joint research · all rows · conditional estimate')
    fig.savefig(dest / 'native_joint_diagnostics.png', dpi=160)
    plt.close(fig)
    workflow='bhtracker.joint_workflow_v2' if result['schema'].endswith('.v2') else 'bhtracker.joint_workflow'
    (dest / 'REPRODUCE.md').write_text('''# Joint research replay and refit

Install the matching Blackhole Tracker wheel with its astrometry extra:

    python -m pip install "./bhsm_blackhole_tracker-0.2.0rc4-py3-none-any.whl[astrometry]"

Use the supplied local wheel when this release is not on a package index.
All commands below run in the extracted bundle directory.

    python -m bhtracker replay . --output replayed
    python -m WORKFLOW_MODULE replay result.json
    python -m WORKFLOW_MODULE fit dataset.json refitted.json --options "$(cat model_options.json)"
    python -m WORKFLOW_MODULE profile dataset.json refitted.json mass_profile.json

The refit example's shell substitution is for POSIX shells. On PowerShell use
`--options (Get-Content -Raw model_options.json)`, or call `fit_joint` from Python.
Forward replay checks the exact scientific-source hash, dataset hash, saved
state, predictions, residuals and derived parameters. It does not rerun the
optimizer. Refitting and mass profiling are separate potentially long jobs.
No confidence-interval coverage, adequate noise model, dark lens or horizon
is certified by replay. Preserve the separate validation/control receipts.
Native bytes and exclusions remain in dataset.json, with source row identities,
clock and calibration provenance. CSV epochs are model evaluation coordinates,
not a replacement for native clocks or an author's unknown averaging weights.
BHSM interpretation uses one spacetime-side optical potential counted once.
'''.replace('WORKFLOW_MODULE',workflow), encoding='utf-8')
    files = sorted(path for path in dest.iterdir() if path.is_file() and path.name != 'manifest.json')
    manifest = {'schema': BUNDLE_SCHEMA, 'dataset_sha256': dataset.sha256,
                'scientific_source_sha256': result['scientific_source_sha256'],
                'replay_method': 'FORWARD_STATE_REPLAY', 'project_citation': project_citation(),
                'files': {p.name: hashlib.sha256(p.read_bytes()).hexdigest() for p in files}}
    (dest / 'manifest.json').write_text(json.dumps(manifest, indent=2) + '\n', encoding='utf-8')
    return str(dest.resolve())


def replay_bundle(read, manifest, output):
    from .data import read_json
    required = {'dataset.json', 'result.json', 'model_options.json', 'forward_replay_check.json'}
    if not required <= set(manifest.get('files', {})):
        raise ValueError('Joint bundle manifest lacks required replay files')
    for name, digest in manifest['files'].items():
        if Path(name).name != name or '/' in name or '\\' in name:
            raise ValueError('Joint manifest filenames must be single path components')
        if hashlib.sha256(read(name)).hexdigest() != digest:
            raise ValueError(f'Joint replay hash mismatch for {name}')
    result = read_json(read('result.json'))
    dataset, _, checked = validate_result(result)
    if read_json(read('dataset.json')) != dataset.document or dataset.sha256 != manifest.get('dataset_sha256'):
        raise ValueError('Joint bundle input disagrees with the fitted dataset')
    if read_json(read('model_options.json')) != result['model_options']:
        raise ValueError('Joint bundle model options disagree with the saved result')
    if manifest.get('scientific_source_sha256') != result['scientific_source_sha256']:
        raise ValueError('Joint bundle source pin differs from the result')
    export_result(result, output)
    (Path(output) / 'replay_check.json').write_text(json.dumps(checked, indent=2) + '\n', encoding='utf-8')
    return checked
