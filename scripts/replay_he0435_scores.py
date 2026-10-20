"""Reproduce a known held-out score without fitting or claiming a fresh test.

The original exposure marker and result are preserved. This explicit replay
writes separately and cannot reopen preparation or overwrite original data.
"""
import argparse
from hashlib import sha256
import json
from pathlib import Path
import sys

import numpy as np

ROOT=Path(__file__).resolve().parents[1]
sys.path[:0]=[str(ROOT/'code'),str(ROOT/'scripts')]
from fit_he0435_lightcurves import DATA, load_training, CUT
from strong_lensing_observed_lightcurve import LightCurves, OUParameters, concatenate, score_predictions


def main(output):
    directory=ROOT/'artifacts/strong_lensing_observed_20261020'
    original=directory/'withheld_scores.json'
    if not original.is_file() or not (directory/'scoring_started.json').is_file():
        raise RuntimeError('Replay is only for the already-exposed archived experiment')
    if output.resolve()==original.resolve() or output.resolve().is_relative_to((ROOT/'data').resolve()):
        raise ValueError('Replay output must not overwrite the original score or observations')
    freeze=json.loads((directory/'frozen_predictions.json').read_text())
    for path,expected in freeze['source_sha256'].items():
        assert sha256((ROOT/path).read_bytes()).hexdigest()==expected,path
    rows=np.loadtxt(DATA,usecols=range(9)); rows=rows[rows[:,0]>CUT]
    later=LightCurves(np.repeat(rows[:,0],4),np.tile(np.arange(4),len(rows)),
                     rows[:,1::2].ravel(),rows[:,2::2].ravel(),4)
    train=load_training(); driver=later.subset(later.image==1); target=later.subset(later.image!=1)
    old=json.loads(original.read_text()); result={'status':'KNOWN_OUTCOME_NUMERICAL_REPLAY_NOT_AN_INDEPENDENT_TEST','scores':{}}
    for name,values in freeze['models'].items():
        if name=='matched_conventional_effective_mass':
            result['scores'][name]=result['scores']['conditional_R1_macro']; continue
        parameters=OUParameters(**values)
        primary=score_predictions(concatenate(train,driver),target,parameters)
        secondary=score_predictions(train,later,parameters)
        result['scores'][name]={'primary_joint_log_density':primary['joint_conditional_test_log_density'],
                               'secondary_joint_log_density':secondary['joint_conditional_test_log_density']}
        for mode,score in [('primary',primary),('secondary',secondary)]:
            np.testing.assert_allclose(score['predicted_magnitude'],old[mode][name]['predicted_magnitude'],rtol=1e-9,atol=1e-9)
            np.testing.assert_allclose(score['joint_conditional_test_log_density'],old[mode][name]['joint_conditional_test_log_density'],rtol=0,atol=1e-5)
    output.parent.mkdir(parents=True,exist_ok=True)
    output.write_text(json.dumps(result,indent=2)+'\n',encoding='utf-8',newline='\n')
    print('PASS known-outcome replay:',output)


if __name__=='__main__':
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output',type=Path,default=ROOT/'manuscript/build/he0435_known_outcome_replay.json')
    main(parser.parse_args().output)
