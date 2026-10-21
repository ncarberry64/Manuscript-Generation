"""Freeze train-only physical predictions, then score reserved HE0435 light once.

Run --prepare first, review and commit the frozen artifacts, then run --score.
The latter refuses to overwrite a score or consume a changed prediction input.
The explicit later-image-B driver is allowed conditioning, never a refit.
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

ROOT = Path(__file__).resolve().parents[1]
sys.path[:0] = [str(ROOT/'code'), str(ROOT/'scripts')]
from fit_he0435_lightcurves import load_training
from strong_lensing_observed_lightcurve import (
    LightCurves, OUParameters, concatenate, log_likelihood, score_predictions,
)
from strong_lensing_observed_macro import fit_macro

OUT = ROOT/'artifacts/strong_lensing_observed_20261020'
DATA = ROOT/'data/strong_lensing_observed_20261020/candidates/HE0435-1223/HE0435_all_Millon2020.dat'
PROTOCOL = ROOT/'preregistration/he0435_heldout_20261020.json'
CUT = 57447.03918


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
    path.write_text(json.dumps(value, indent=2, default=encode, allow_nan=False)+'\n',
                    encoding='utf-8', newline='\n')


def prepare():
    if (OUT/'withheld_scores.json').exists() or (OUT/'scoring_started.json').exists():
        raise RuntimeError('Outcomes already scored: do not retune this held-out experiment')
    protocol = json.loads(PROTOCOL.read_text())
    assert digest(DATA) == protocol['data_sha256']
    train = load_training()
    receipt = json.loads((OUT/'lightcurve_training.json').read_text())
    empirical = OUParameters(**receipt['parameters'])
    curvature = json.loads((OUT/'training_curvature.json').read_text())
    local = next(x for x in curvature['step_sensitivity'] if x['step_days'] == .1)
    covariance = np.asarray(local['local_inverse_curvature_covariance_days2'])

    def training_log_density(delays):
        parameters = replace(empirical, delays_days=tuple(np.r_[0., delays]))
        return log_likelihood(train, parameters)

    conventional = fit_macro(empirical.delays_days[1:], covariance,
                             delay_log_likelihood=training_log_density)
    save(OUT/'conventional_training.json', conventional)
    print(json.dumps({'model': 'conventional', 'delays': conventional['predicted_delays_days'],
                      'astrometry_chi2': conventional['chi2_astrometry'] }), flush=True)
    screened = fit_macro(empirical.delays_days[1:], covariance, model='r1',
                         initial=conventional['parameter_vector'],
                         delay_log_likelihood=training_log_density)
    save(OUT/'r1_training.json', screened)
    print(json.dumps({'model': 'r1', 'delays': screened['predicted_delays_days'],
                      'astrometry_chi2': screened['chi2_astrometry']}), flush=True)
    models = {
        'empirical_training_transfer': empirical,
        'conventional_macro': replace(empirical, delays_days=tuple(
            np.r_[0., conventional['predicted_delays_days']])),
        'conditional_R1_macro': replace(empirical, delays_days=tuple(
            np.r_[0., screened['predicted_delays_days']])),
        # The equality follows from the source-derived effective mass, and
        # its independent original projection is checked in r1_training.
        'matched_conventional_effective_mass': replace(empirical, delays_days=tuple(
            np.r_[0., screened['predicted_delays_days']])),
        'zero_delay_control': replace(empirical, delays_days=(0., 0., 0., 0.)),
    }
    paths = [PROTOCOL, DATA, OUT/'lightcurve_training.json', OUT/'training_curvature.json',
             OUT/'conventional_training.json', OUT/'r1_training.json', Path(__file__),
             ROOT/'code/strong_lensing_observed_lightcurve.py',
             ROOT/'code/strong_lensing_observed_macro.py', ROOT/'scripts/fit_he0435_lightcurves.py',
             ROOT/'code/strong_lensing_comparison.py', ROOT/'code/strong_lensing_robustness.py',
             ROOT/'code/strong_lensing_r1_screened.py', ROOT/'code/strong_lensing_local_family.py',
             ROOT/'code/background.py', ROOT/'code/action_native_matter_n2.py',
             ROOT/'preregistration/prediction_manifest.json',
             ROOT/'data/strong_lensing_observed_20261020/candidates/HE0435-1223/astrometry_redshifts.json']
    frozen = {
        'status': 'FROZEN_BEFORE_TARGET_SCORING',
        'protocol_commit': 'e55725d19135c3a434b9e12621845eae122e1b8f',
        'training_epoch_count': len(train.time_days)//4,
        'training_observation_count': len(train.time_days),
        'target_values_consumed': False,
        'empirical_nuisance_fit': 'Two-stage conditional empirical Bayes: fit common variability/delays on training, then hold variability hyperparameters fixed while each macro model fits astrometry plus the exact training light-curve likelihood.',
        'macro_delay_curvature': '0.1day fixed-hyperparameter curvature is a numerical initializer only; final macro objective uses original training light-curve log density, not a claimed observed Gaussian delay posterior.',
        'parameter_uncertainty': 'All predictive intervals/scores condition on fitted nuisance and macro parameters; no full hyperparameter or physical posterior marginalization.',
        'models': {name: asdict(value) for name, value in models.items()},
        'training_log_likelihood_by_model': {name: log_likelihood(train, value)
                                            for name, value in models.items()},
        'source_sha256': {p.relative_to(ROOT).as_posix(): digest(p) for p in paths},
        'structural_control': 'The matched conventional effective-mass field has the same transfer as conditional R1. Its shared score is a structural identity, not a separate fitted confirmation.',
        'physical_admission': False,
    }
    save(OUT/'frozen_predictions.json', frozen)
    print('FROZEN', digest(OUT/'frozen_predictions.json'), flush=True)


def score():
    output = OUT/'withheld_scores.json'
    marker = OUT/'scoring_started.json'
    if output.exists() or marker.exists():
        raise RuntimeError('Withheld outcomes have already been scored; preserve that result')
    frozen_path = OUT/'frozen_predictions.json'
    frozen = json.loads(frozen_path.read_text())
    for path, expected in frozen['source_sha256'].items():
        if digest(ROOT/path) != expected:
            raise RuntimeError('Prediction input changed after freeze: '+path)
    # Require that the reviewed freeze itself was committed before this call.
    git_path = frozen_path.relative_to(ROOT).as_posix()
    committed = subprocess.check_output(['git', 'show', 'HEAD:'+git_path], cwd=ROOT)
    if sha256(committed).hexdigest() != digest(frozen_path):
        raise RuntimeError('Commit the reviewed prediction freeze before opening the targets')
    # Persist first exposure even if a later prediction/score call aborts.
    # An aborted evaluation must never silently reopen model preparation.
    with marker.open('x', encoding='utf-8', newline='\n') as stream:
        json.dump({'status': 'TARGET_OPENING_STARTED_PARAMETERS_LOCKED',
                   'frozen_prediction_sha256': digest(frozen_path),
                   'freeze_commit': subprocess.check_output(
                       ['git', 'rev-parse', 'HEAD'], cwd=ROOT, text=True).strip()}, stream, indent=2)
        stream.write('\n')
    train = load_training()
    rows = []
    for line in DATA.read_text().splitlines():
        tokens = line.split()
        if tokens and float(tokens[0]) > CUT:
            rows.append([float(x) for x in tokens[:9]])
    rows = np.asarray(rows)
    assert rows.shape == (100, 9)
    later = LightCurves(np.repeat(rows[:, 0], 4), np.tile(np.arange(4), len(rows)),
                        rows[:, 1::2].ravel(), rows[:, 2::2].ravel(), 4)
    driver = later.subset(later.image == 1)
    targets = later.subset(later.image != 1)
    conditioning = concatenate(train, driver)
    result = {'status': 'ONE_FROZEN_RETROSPECTIVE_WITHHELD_APPLICATION',
              'lens': 'HE0435-1223', 'frozen_prediction_sha256': digest(frozen_path),
              'publication_freeze_commit': subprocess.check_output(
                  ['git', 'rev-parse', 'HEAD'], cwd=ROOT, text=True).strip(),
              'primary': {}, 'secondary': {},
              'primary_target_count': len(targets.time_days),
              'allowed_driver_count': len(driver.time_days),
              'secondary_target_count': len(later.time_days),
              'qualification': 'Retrospective batch cross-image prediction conditions on all later B points, including driver points later than some target dates. The secondary all-image forecast uses training only. No refitting, per-image magnitude shift or new microlensing fit consumes target values. Intervals are point-parameter conditional; measurement-noise covariance is approximated as independent while latent brightness covariance is correlated.',
              'native_epoch_origin': 'Use native day differences only; published half-day origin discrepancy does not affect delays or split.',
              'target_design': {'time_days': targets.time_days, 'image': targets.image,
                                'magnitude': targets.magnitude, 'error_mag': targets.error_mag},
              'driver_design': {'time_days': driver.time_days, 'image': driver.image,
                                'magnitude': driver.magnitude, 'error_mag': driver.error_mag},
              'physical_admission': False}
    for name, values in frozen['models'].items():
        parameters = OUParameters(**values)
        if name == 'matched_conventional_effective_mass':
            # Equality was derived and projection-checked before the test.
            # Reuse the exact shared transfer without pretending independent data.
            result['primary'][name] = result['primary']['conditional_R1_macro']
            result['secondary'][name] = result['secondary']['conditional_R1_macro']
        else:
            result['primary'][name] = score_predictions(conditioning, targets, parameters)
            result['secondary'][name] = score_predictions(train, later, parameters)
        print(json.dumps({'model': name,
                          'primary_joint_log_density': result['primary'][name]['joint_conditional_test_log_density'],
                          'primary_summary': result['primary'][name]['overall']}), flush=True)
    save(output, result)


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    mode = parser.add_mutually_exclusive_group(required=True)
    mode.add_argument('--prepare', action='store_true')
    mode.add_argument('--score', action='store_true')
    args = parser.parse_args()
    prepare() if args.prepare else score()
