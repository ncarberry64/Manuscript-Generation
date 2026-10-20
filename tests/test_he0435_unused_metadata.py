"""Check the metadata-only acquisition boundary without parsing real targets."""

import importlib.util
from pathlib import Path


SCRIPT = Path(__file__).resolve().parents[1] / "scripts/acquire_he0435_unused_metadata.py"
SPEC = importlib.util.spec_from_file_location("he0435_unused_acquirer", SCRIPT)
ACQUIRE = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(ACQUIRE)


def isolated_paths(monkeypatch, tmp_path):
    out = tmp_path / "data"
    sealed = out / "sealed"
    sealed.mkdir(parents=True)
    monkeypatch.setattr(ACQUIRE, "ROOT", tmp_path)
    monkeypatch.setattr(ACQUIRE, "OUT", out)
    monkeypatch.setattr(ACQUIRE, "SEALED", sealed)
    return out, sealed


def test_acquisition_never_parses_non_time_columns(monkeypatch, tmp_path):
    """Nonnumeric target sentinels would fail any unintended photometry parse."""
    _, sealed = isolated_paths(monkeypatch, tmp_path)
    raw = b"58300.0000 DO_NOT_PARSE_TARGET_VALUES\n58488.2085 STILL_NOT_NUMERIC\n60326.2011 SEALED_OUTCOME\n"
    (sealed / "he0435r.dat").write_bytes(raw)
    monkeypatch.setattr(ACQUIRE, "FILES", {"he0435r.dat": 3})
    monkeypatch.setattr(ACQUIRE.requests, "get", lambda *a, **k: (_ for _ in ()).throw(AssertionError("Unexpected network access")))
    receipt = ACQUIRE.acquire()
    item = receipt["files"]["he0435r.dat"]
    assert item["rows"] == 3
    assert item["post_cut_rows"] == 2
    assert item["photometry_parsed"] is False
    assert (sealed / "he0435r.dat").read_bytes() == raw


def test_split_preserves_rows_without_parsing_targets(monkeypatch, tmp_path):
    out, sealed = isolated_paths(monkeypatch, tmp_path)
    early = b"58300.0000 PRECALIBRATION_NOT_NUMERIC\r\n"
    late = b"58484.0000 BOUNDARY_TARGET_SENTINEL\r\n60326.2011 LATER_TARGET_SENTINEL\r\n"
    (sealed / "he0435r.dat").write_bytes(early + late)
    split = ACQUIRE.split_r_without_disclosing_targets()
    assert (out / "calibration/he0435r_pre2019.dat").read_bytes() == early
    assert (sealed / "he0435r_2019onward.dat").read_bytes() == late
    assert split["development"]["rows"] == 1
    assert split["confirmatory"]["rows"] == 2
    assert split["primary_test"].startswith("Condition on new image C")


def test_existing_acquisition_metadata_remains_stable(monkeypatch, tmp_path):
    _, sealed = isolated_paths(monkeypatch, tmp_path)
    (sealed / "he0435r.dat").write_bytes(b"58300.0000 NOT_A_PHOTOMETRIC_NUMBER\n60326.2011 SEALED\n")
    monkeypatch.setattr(ACQUIRE, "FILES", {"he0435r.dat": 2})
    first = ACQUIRE.acquire()
    second = ACQUIRE.acquire()
    assert second["acquired_utc"] == first["acquired_utc"]
    assert second["files"] == first["files"]
