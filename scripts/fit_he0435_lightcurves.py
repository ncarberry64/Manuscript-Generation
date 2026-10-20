"""Frozen-protocol, training-only HE0435 OU fits; never scores later targets."""
from __future__ import annotations

from concurrent.futures import ProcessPoolExecutor, as_completed
from dataclasses import asdict
from hashlib import sha256
import json
from pathlib import Path
import sys
import time

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT/'code'))
from strong_lensing_observed_lightcurve import (
    LightCurves, fit_training, training_delay_curvature,
)


DATA = ROOT/'data/strong_lensing_observed_20261020/candidates/HE0435-1223/HE0435_all_Millon2020.dat'
OUTPUT = ROOT/'artifacts/strong_lensing_observed_20261020'
CUT = 57447.03918


def load_training(path=DATA):
    # Only timestamp metadata is examined for later rows. Their magnitude and
    # uncertainty tokens are not converted, summarized, or returned.
    rows = []
    for line in path.read_text().splitlines():
        tokens = line.split()
        if not tokens or float(tokens[0]) > CUT:
            continue
        rows.append([float(token) for token in tokens[:9]])
    rows = np.asarray(rows)
    return LightCurves(np.repeat(rows[:, 0], 4), np.tile(np.arange(4), len(rows)),
                       rows[:, 1::2].ravel(), rows[:, 2::2].ravel(), 4)


def one_start(index, delays, magnitude_zero):
    training = load_training()
    start = time.perf_counter()
    parameters, receipt = fit_training(training, [delays], maxiter=200,
                                        magnitude_zero=magnitude_zero)
    return index, parameters, receipt, time.perf_counter()-start


def main():
    OUTPUT.mkdir(exist_ok=True, parents=True)
    training = load_training()
    magnitude_zero = float(np.median(training.magnitude))
    rng = np.random.default_rng(20261020)
    starts = [np.zeros(4).tolist()]+[np.r_[0., rng.uniform(-20., 20., 3)].tolist()
                                     for _ in range(3)]
    results = []
    started = time.perf_counter()
    with ProcessPoolExecutor(max_workers=4) as pool:
        jobs = [pool.submit(one_start, index, delays, magnitude_zero)
                for index, delays in enumerate(starts)]
        for job in as_completed(jobs):
            index, parameters, receipt, elapsed = job.result()
            result = {'protocol_start_index': index, 'fit': receipt, 'wall_seconds': elapsed}
            (OUTPUT/f'training_start_{index}.json').write_text(json.dumps(result, indent=2)+'\n')
            results.append((index, parameters, receipt, elapsed))
            print(json.dumps({'completed_start': index, 'log_likelihood': receipt['training_log_likelihood'],
                              'optimizer': receipt['starts'][0]['optimizer_message'],
                              'wall_seconds': elapsed}), flush=True)
    results.sort(key=lambda result: result[0])
    best = max(results, key=lambda result: result[2]['training_log_likelihood'])
    parameters = best[1]
    receipt = dict(best[2])
    receipt.update(starts=[result[2]['starts'][0] for result in results],
                   selected_start_index=best[0],
                   input_path=DATA.relative_to(ROOT).as_posix(),
                   input_sha256=sha256(DATA.read_bytes()).hexdigest(),
                   consumed_epoch_cut_inclusive=CUT,
                   consumed_epoch_count=len(training.time_days)//4,
                   magnitude_zero_rule='Median of pooled training magnitudes only',
                   wall_seconds=time.perf_counter()-started,
                   withheld_values_consumed=False,
                   code_sha256=sha256((ROOT/'code/strong_lensing_observed_lightcurve.py').read_bytes()).hexdigest(),
                   runner_sha256=sha256(Path(__file__).read_bytes()).hexdigest())
    (OUTPUT/'lightcurve_training.json').write_text(json.dumps(receipt, indent=2)+'\n')
    curvatures = [training_delay_curvature(training, parameters, step_days=step)
                  for step in (.01, .02, .05, .1)]
    curvature_receipt = {'parameters': asdict(parameters), 'training_only': True,
                         'step_sensitivity': curvatures,
                         'code_sha256': receipt['code_sha256']}
    (OUTPUT/'training_curvature.json').write_text(json.dumps(curvature_receipt, indent=2)+'\n')
    print(json.dumps({'selected_start': best[0], 'parameters': asdict(parameters),
                      'training_log_likelihood': receipt['training_log_likelihood'],
                      'wall_seconds': receipt['wall_seconds']}), flush=True)


if __name__ == '__main__':
    main()
