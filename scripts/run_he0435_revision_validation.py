"""Commit a revised development freeze, then open independent LCO targets once.

No target parser is called by prepare(). A durable exclusive-create marker
locks every prediction input before the first numerical target read.
"""
from __future__ import annotations

import argparse
from dataclasses import asdict, replace
from hashlib import sha256
import json
from pathlib import Path
import subprocess
import sys

import numpy as np
from scipy.special import logsumexp

ROOT = Path(__file__).resolve().parents[1]
sys.path[:0] = [str(ROOT/'code'), str(ROOT/'scripts')]
from fit_he0435_revision import INSTRUMENTS, load_development
from strong_lensing_lightcurve_revision import RevisionCurves, RevisionParameters, log_likelihood, score_mixture

OUT = ROOT/'artifacts/strong_lensing_revision_20261020'
PROTOCOL = ROOT/'preregistration/he0435_unused_validation_20261020.json'
TARGET = ROOT/'data/strong_lensing_revision_20261020/sealed/he0435r_2019onward.dat'
FREEZE = OUT/'frozen_revision_predictions.json'
MARKER = OUT/'revision_scoring_started.json'
SCORES = OUT/'revision_unused_scores.json'


def digest(path):
    return sha256(Path(path).read_bytes()).hexdigest()


def encode(value):
    if isinstance(value, np.ndarray):
        return value.tolist()
    if isinstance(value, np.generic):
        return value.item()
    raise TypeError(type(value).__name__)


def save(path, value):
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(value, indent=2, default=encode, allow_nan=False)+'\n', encoding='utf-8', newline='\n')


def assert_unexposed():
    if MARKER.exists() or SCORES.exists():
        raise RuntimeError('Independent targets already exposed; preserve results and do not retune this experiment')


def transported_components(train, particles, weights, delays):
    """Declared finite development-only nuisance transport, not full evidence."""
    source = [RevisionParameters(**p) for p in particles]
    transformed = [replace(p, delays_days=tuple(delays)) for p in source]
    old = np.array([log_likelihood(train, p) for p in source])
    new = np.array([log_likelihood(train, p) for p in transformed])
    logweights = np.log(np.asarray(weights))+new-old
    logweights -= logsumexp(logweights)
    return {'particles': [asdict(p) for p in transformed], 'weights': np.exp(logweights),
            'development_log_likelihood': new,
            'transport_effective_sample_size': float(1/np.sum(np.exp(2*logweights))),
            'delay_vector_days': list(delays)}


def prepare():
    assert_unexposed()
    protocol = json.loads(PROTOCOL.read_text())
    if digest(TARGET) != protocol['target_sha256']:
        raise RuntimeError('Independent target file changed')
    development_path = OUT/'revision_development_fit.json'
    development = json.loads(development_path.read_text())
    conventional_path = OUT/'revised_conventional_training.json'
    r1_path = OUT/'revised_r1_training.json'
    conventional = json.loads(conventional_path.read_text())
    r1 = json.loads(r1_path.read_text())
    if conventional['compact_mass_Msun'] != r1['compact_mass_Msun']:
        raise RuntimeError('Paired models require the same selected compact source')
    train = load_development()
    particles, weights = development['particles'], development['weights']
    models = {'empirical_development_delays': {
        'particles': particles, 'weights': weights,
        'transport_effective_sample_size': float(1/np.sum(np.asarray(weights)**2))}}
    for name, fit in [('revised_conventional_macro', conventional), ('revised_conditional_R1_macro', r1)]:
        models[name] = transported_components(train, particles, weights, [0., *fit['predicted_delays_days']])
    fixed_delays = np.asarray(conventional['predicted_delays_days'])+np.asarray(
        r1['fixed_conventional_body_R1_activation']['delay_offsets_seconds'])/86400.
    models['fixed_body_R1_activation'] = transported_components(train, particles, weights, [0., *fixed_delays])
    models['matched_conventional_effective_mass'] = {
        'identical_to': 'revised_conditional_R1_macro',
        'reason': 'Source-derived positive effective conventional mass has identical retained optical metric; structural score identity, not independent confirmation'}
    models['zero_delay_control'] = transported_components(train, particles, weights, [0.]*4)
    # Freeze all relevant local applications, inputs, candidate/HST results and
    # calibration. Score output and exposure marker are intentionally excluded.
    paths = [PROTOCOL, TARGET, development_path, conventional_path, r1_path,
             ROOT/'preregistration/he0435_revision_development_20261020.json',
             ROOT/'preregistration/he0435_revision_long_horizon_20261020.json',
             ROOT/'data/strong_lensing_observed_20261020/candidates/HE0435-1223/HE0435_all_Millon2020.dat',
             ROOT/'data/strong_lensing_observed_20261020/candidates/HE0435-1223/astrometry_redshifts.json',
             ROOT/'data/strong_lensing_revision_20261020/calibration/he0435r_pre2019.dat',
             ROOT/'preregistration/prediction_manifest.json', ROOT/'code/background.py',
             ROOT/'code/action_native_matter_n2.py', Path(__file__)]
    paths += sorted((ROOT/'code').glob('strong_lensing*.py'))
    paths += [ROOT/'code/strong_lensing_revision_filter.c']
    paths += [ROOT/'scripts'/name for name in ['fit_he0435_revision.py', 'test_he0435_central_image.py', 'fit_he0435_macro_revision.py']]
    paths += sorted((ROOT/'data/strong_lensing_revision_20261020/hst').glob('*cutout*.fits'))
    paths += sorted((ROOT/'data/strong_lensing_revision_20261020').glob('*.json'))
    paths += [ROOT/'data/strong_lensing_revision_20261020/hst/provenance.json',
              ROOT/'data/strong_lensing_revision_20261020/sealed/ReadMe',
              ROOT/'preregistration/he0435_macro_revision_development_20261020.json']
    paths += [p for p in sorted(OUT.glob('*.json')) if p.name not in
              {FREEZE.name, MARKER.name, SCORES.name, 'revision_publication_manifest.json'}]
    paths = sorted(set(paths))
    missing = [str(p) for p in paths if not p.exists()]
    if missing:
        raise RuntimeError('Incomplete pre-score applications: '+str(missing))
    frozen = {'status': 'REVISED_MODELS_FROZEN_BEFORE_UNUSED_TARGET_PARSE',
              'target_photometry_consumed': False, 'models': models,
              'calibration': development['calibration'],
              'posterior_diagnostics': development['posterior_diagnostics'],
              'source_sha256': {p.relative_to(ROOT).as_posix(): digest(p) for p in paths},
              'qualification': protocol['uncertainty'],
              'physical_guard': 2e-10, 'physical_admission': False}
    save(FREEZE, frozen)
    print(json.dumps({'freeze_sha256': digest(FREEZE), 'models': {k: v.get('transport_effective_sample_size') for k, v in models.items()}}), flush=True)


def load_targets():
    if not MARKER.exists():
        raise RuntimeError('Numerical target parser requires the persistent exposure lock')
    rows = np.loadtxt(TARGET, usecols=range(9))
    if rows.shape != (54, 9) or np.min(rows[:, 0]) < 58484:
        raise RuntimeError('Unexpected independent target domain')
    return RevisionCurves(np.repeat(rows[:, 0], 4), np.tile(np.arange(4), len(rows)),
                          rows[:, 1::2].ravel(), rows[:, 2::2].ravel(),
                          np.full(len(rows)*4, INSTRUMENTS.index('LCO_R')), 4, INSTRUMENTS)


def committed_bytes(path):
    return subprocess.check_output(['git', 'show', 'HEAD:'+path.relative_to(ROOT).as_posix()], cwd=ROOT)


def verify_frozen_inputs(frozen):
    if sha256(committed_bytes(FREEZE)).hexdigest() != digest(FREEZE):
        raise RuntimeError('Commit reviewed prediction freeze before target opening')
    for relative, expected in frozen['source_sha256'].items():
        path = ROOT/relative
        if digest(path) != expected:
            raise RuntimeError('Changed frozen input: '+relative)
        if sha256(committed_bytes(path)).hexdigest() != expected:
            raise RuntimeError('Frozen input bytes must also be committed: '+relative)


def score():
    assert_unexposed()
    frozen = json.loads(FREEZE.read_text())
    verify_frozen_inputs(frozen)
    commit = subprocess.check_output(['git', 'rev-parse', 'HEAD'], cwd=ROOT, text=True).strip()
    with MARKER.open('x', encoding='utf-8', newline='\n') as stream:
        json.dump({'status': 'UNUSED_TARGET_EXPOSURE_STARTED_NO_RETUNING',
                   'freeze_sha256': digest(FREEZE), 'freeze_commit': commit}, stream, indent=2)
        stream.write('\n')
    # This is the first numerical parse of independent confirmatory photometry.
    train, later = load_development(), load_targets()
    driver, targets = later.subset(later.image == 2), later.subset(later.image != 2)
    result = {'status': 'ONE_FROZEN_INDEPENDENT_EXPOSURE_RETROSPECTIVE_TEST',
              'freeze_commit': commit, 'freeze_sha256': digest(FREEZE),
              'primary_target_count': len(targets.time_days), 'driver_count': len(driver.time_days),
              'secondary_target_count': len(later.time_days), 'primary': {}, 'secondary': {},
              'target_design': {k: getattr(later, k) for k in ['time_days', 'image', 'magnitude', 'error_mag']},
              'physical_guard': 2e-10, 'physical_admission': False,
              'qualification': 'Independent LCO exposures, retrospective shared DIA reference; qualitative publication awareness disclosed. C-driver primary and all-image secondary overlap. No target fitting. Fixed macro and approximate finite nuisance uncertainty; correlated epochs preclude iid-binomial coverage significance.'}
    for name, model in frozen['models'].items():
        if 'identical_to' in model:
            for mode in ('primary', 'secondary'):
                result[mode][name] = {'identical_to': model['identical_to'], 'reason': model['reason']}
            continue
        for mode, test, allowed, calmode in [('primary', targets, driver, 'primary_C_driver'),
                                            ('secondary', later, None, 'secondary_forecast')]:
            scales = frozen['calibration']['modes'][calmode]['per_image_scale']
            result[mode][name] = score_mixture(train, test, model['particles'], model['weights'],
                                              driver=allowed, calibration=scales)
            print(json.dumps({'mode': mode, 'model': name, 'log_score': result[mode][name]['joint_conditional_log_density'],
                              **result[mode][name]['overall']}), flush=True)
        # Persist substantive progress after every model, without clearing lock.
        save(OUT/'revision_score_progress.json', result)
    save(SCORES, result)
    print('SCORED', digest(SCORES), flush=True)


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('mode', choices=('prepare', 'score'))
    args = parser.parse_args()
    prepare() if args.mode == 'prepare' else score()
