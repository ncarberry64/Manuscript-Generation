"""Training-only finite-step and conditional delay-slice diagnostics, no refit."""
from __future__ import annotations

from concurrent.futures import ProcessPoolExecutor
from dataclasses import replace
from hashlib import sha256
import json
from pathlib import Path
import sys

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT/'code'))
from strong_lensing_observed_lightcurve import OUParameters, _decode, _encode, log_likelihood
from fit_he0435_lightcurves import load_training


OUTPUT = ROOT/'artifacts/strong_lensing_observed_20261020'


def delay_slice(image):
    training = load_training()
    receipt = json.loads((OUTPUT/'lightcurve_training.json').read_text())
    parameters = OUParameters(**receipt['parameters'])
    base = log_likelihood(training, parameters)
    points = []
    for offset in np.linspace(-5., 5., 101):
        delays = np.asarray(parameters.delays_days).copy()
        delays[image] += offset
        value = log_likelihood(training, replace(parameters, delays_days=tuple(delays)))
        points.append({'delay_days': float(delays[image]), 'offset_days': float(offset),
                       'delta_log_likelihood_from_frozen_best': float(value-base)})
    return {'image_index': image, 'all_other_parameters_fixed': True, 'points': points}


def main():
    training = load_training()
    training_receipt = json.loads((OUTPUT/'lightcurve_training.json').read_text())
    parameters = OUParameters(**training_receipt['parameters'])
    coordinates = _encode(parameters)
    bounds = training_receipt['coordinate_bounds']
    gradients = []
    for delay_step, log_step in ((.0001, .0001), (.001, .0001), (.01, .001)):
        steps = np.r_[np.full(3, delay_step), np.full(5, log_step)]
        gradient = []
        for index, step in enumerate(steps):
            delta = np.zeros_like(coordinates)
            delta[index] = step
            plus = -log_likelihood(training, _decode(coordinates+delta, parameters))
            minus = -log_likelihood(training, _decode(coordinates-delta, parameters))
            gradient.append(float((plus-minus)/(2*step)))
        # KKT-compatible projection at a lower/upper active bound. Symmetric
        # diagnostics may evaluate just outside bounds but do not refit there.
        projected = np.asarray(gradient).copy()
        for i, (low, high) in enumerate(bounds):
            if coordinates[i] <= low+1e-6 and projected[i] > 0:
                projected[i] = 0.
            if coordinates[i] >= high-1e-6 and projected[i] < 0:
                projected[i] = 0.
        gradients.append({'delay_step_days': delay_step, 'log_hyperparameter_step': log_step,
                          'negative_log_likelihood_central_gradient': gradient,
                          'projected_gradient': projected.tolist(),
                          'maximum_absolute_projected_gradient': float(np.max(abs(projected)))})
    with ProcessPoolExecutor(max_workers=3) as pool:
        slices = list(pool.map(delay_slice, (1, 2, 3)))
    result = {'training_only': True, 'parameters_unchanged': True,
              'gradient_step_sensitivity': gradients,
              'conditional_delay_slices': slices,
              'interpretation': 'These are fixed-other-parameter slices, not profiled likelihoods or calibrated confidence intervals; OU delay likelihood has sample-coincidence cusps.',
              'training_receipt_sha256': sha256((OUTPUT/'lightcurve_training.json').read_bytes()).hexdigest(),
              'code_sha256': sha256((ROOT/'code/strong_lensing_observed_lightcurve.py').read_bytes()).hexdigest(),
              'runner_sha256': sha256(Path(__file__).read_bytes()).hexdigest()}
    (OUTPUT/'training_diagnostics.json').write_text(json.dumps(result, indent=2)+'\n')
    print(json.dumps({'gradients': gradients,
                      'slice_maxima_offsets_days': [max(row['points'], key=lambda p:p['delta_log_likelihood_from_frozen_best'])['offset_days'] for row in slices]}))


if __name__ == '__main__':
    main()
