"""Exposure-lock and conditional nuisance-weight regressions; no real targets."""
from dataclasses import asdict
import importlib.util
import json
from pathlib import Path
import sys

import numpy as np
import pytest

ROOT = Path(__file__).resolve().parents[1]
sys.path[:0] = [str(ROOT/'code'), str(ROOT/'scripts')]
import run_he0435_revision_validation as run
from strong_lensing_lightcurve_revision import RevisionParameters


def test_unopened_target_parser_does_not_touch_photometry(tmp_path, monkeypatch):
    monkeypatch.setattr(run, 'MARKER', tmp_path/'unexposed.json')
    monkeypatch.setattr(run.np, 'loadtxt', lambda *a, **k: pytest.fail('Parsed target before lock'))
    with pytest.raises(RuntimeError, match='persistent exposure lock'):
        run.load_targets()


def test_failed_evaluation_keeps_exposure_lock(tmp_path, monkeypatch):
    freeze = tmp_path/'freeze.json'
    freeze.write_text('{}')
    marker = tmp_path/'exposed.json'
    monkeypatch.setattr(run, 'FREEZE', freeze)
    monkeypatch.setattr(run, 'MARKER', marker)
    monkeypatch.setattr(run, 'SCORES', tmp_path/'scores.json')
    monkeypatch.setattr(run, 'verify_frozen_inputs', lambda x: None)
    monkeypatch.setattr(run.subprocess, 'check_output', lambda *a, **k: 'frozen-test-commit\n')
    monkeypatch.setattr(run, 'load_development', lambda: (_ for _ in ()).throw(ValueError('synthetic abort')))
    with pytest.raises(ValueError, match='synthetic abort'):
        run.score()
    assert json.loads(marker.read_text())['freeze_commit'] == 'frozen-test-commit'
    with pytest.raises(RuntimeError, match='already exposed'):
        run.prepare()


def test_transported_delay_weights_include_development_likelihood_ratio(monkeypatch):
    a = RevisionParameters((0., 1., 2., 3.), source_sigma_mag=.1)
    b = RevisionParameters((0., 4., 5., 6.), source_sigma_mag=.3)
    def likelihood(data, p):
        return -p.source_sigma_mag*p.delays_days[1]**2
    monkeypatch.setattr(run, 'log_likelihood', likelihood)
    result = run.transported_components(None, [asdict(a), asdict(b)], [.25, .75], [0., 2., 3., 4.])
    ratio = .25/.75*np.exp((-.1*4+.1*1)-(-.3*4+.3*16))
    np.testing.assert_allclose(result['weights'][0]/result['weights'][1], ratio)
    assert result['particles'][1]['delays_days'] == (0., 2., 3., 4.)


def test_changed_input_aborts_before_exposure(tmp_path, monkeypatch):
    freeze = tmp_path/'freeze.json'
    freeze.write_text('{}')
    dependency = tmp_path/'dependency.dat'
    dependency.write_bytes(b'changed')
    monkeypatch.setattr(run, 'ROOT', tmp_path)
    monkeypatch.setattr(run, 'FREEZE', freeze)
    monkeypatch.setattr(run, 'committed_bytes', lambda p: p.read_bytes())
    with pytest.raises(RuntimeError, match='Changed frozen input'):
        run.verify_frozen_inputs({'source_sha256': {'dependency.dat': 'incorrect'}})
