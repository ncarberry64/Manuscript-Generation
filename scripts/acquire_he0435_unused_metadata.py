"""Acquire independent LCO files without disclosing photometric outcomes.

Only timestamps, row/byte counts, hashes and descriptive metadata are emitted.
The original files remain sealed until a separate, frozen prediction protocol
authorizes their photometric use. This is a retrospective public-data seal, not
a claim that the publishing authors or a blind reduction never saw the data.
"""

from __future__ import annotations

import hashlib
import json
from datetime import UTC, datetime, timedelta
from pathlib import Path

import requests


ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / "data/strong_lensing_revision_20261020"
SEALED = OUT / "sealed"
MIRRORS = (
    "https://cdsarc.cds.unistra.fr/ftp/J/A+A/703/A250/",
    "https://cdsarc.cds.unistra.fr/ftp/J/A%2BA/703/A250/",
    "https://cdsarc.u-strasbg.fr/ftp/J/A+A/703/A250/",
    "https://cdsarc.cfa.harvard.edu/ftp/J/A+A/703/A250/",
)
FILES = {"ReadMe": None, "he0435r.dat": 79, "he0435v.dat": 80}


def acquire() -> dict:
    SEALED.mkdir(parents=True, exist_ok=True)
    receipt_path = OUT / "unused_observation_metadata.json"
    previous = json.loads(receipt_path.read_text(encoding="utf-8")) if receipt_path.exists() else {}
    rows = {}
    errors = []
    for name, expected in FILES.items():
        path = SEALED / name
        url = previous.get("files", {}).get(name, {}).get("url") or MIRRORS[0] + name
        if not path.exists():
            for base in MIRRORS:
                candidate = base + name
                try:
                    response = requests.get(candidate, timeout=40)
                    response.raise_for_status()
                    data = response.content
                    if b"<html" in data.lower()[:500]:
                        raise ValueError("HTML response")
                    if expected is not None:
                        lines = [line for line in data.splitlines() if line.strip()]
                        if len(lines) != expected:
                            raise ValueError("Unexpected row count")
                        # Validate the time field only. Do not parse magnitudes,
                        # uncertainty columns or any predictive target.
                        for line in lines:
                            float(line[:10])
                    path.write_bytes(data)
                    url = candidate
                    break
                except Exception as exc:
                    errors.append({"url": candidate, "error_type": type(exc).__name__})
            else:
                raise RuntimeError(f"No successful metadata-only acquisition: {name}")
        data = path.read_bytes()
        item = {
            "relative_path": path.relative_to(ROOT).as_posix(),
            "url": url,
            "bytes": len(data),
            "sha256": hashlib.sha256(data).hexdigest(),
        }
        if expected is not None:
            lines = [line for line in data.splitlines() if line.strip()]
            times = [float(line[:10]) for line in lines]
            cut = 58484.0  # MJD 2019-01-01; deliberately beyond old 2018 epochs.
            later = [t for t in times if t >= cut]
            item.update({
                "rows": len(times),
                "time_column": "Bytes 1-10, MJD, unshifted observation epoch",
                "minimum_mjd": min(times),
                "maximum_mjd": max(times),
                "candidate_conservative_cut_mjd": cut,
                "pre_cut_rows": sum(t < cut for t in times),
                "post_cut_rows": len(later),
                "post_cut_minimum_mjd": min(later) if later else None,
                "post_cut_maximum_mjd": max(later) if later else None,
                "all_observation_mjd": times,
                "first_utc_approximation": (datetime(1858, 11, 17, tzinfo=UTC) + timedelta(days=min(times))).isoformat(),
                "last_utc_approximation": (datetime(1858, 11, 17, tzinfo=UTC) + timedelta(days=max(times))).isoformat(),
                "photometry_parsed": False,
            })
        rows[name] = item
    receipt = {
        "schema": "he0435-independent-observation-metadata-v1",
        "acquired_utc": previous.get("acquired_utc") or datetime.now(UTC).isoformat(),
        "publication": "Sorgenfrei, Schmidt & Wambsganss (2025), A&A 703, A250",
        "publication_doi": "10.1051/0004-6361/202555933",
        "bibcode": "2025A&A...703A.250S",
        "catalogue_url": "https://cdsarc.cds.unistra.fr/viz-bin/cat/J/A+A/703/A250",
        "primary_methods_url": "https://arxiv.org/html/2509.09341v1",
        "facility": "Las Cumbres Observatory, 1 m CTIO telescopes; independent of prior COSMOGRAIL exposure set",
        "source_observation_dates": "2014-07-24 through 2024-01-17; 677 R and 697 V exposures combined into nightly epochs",
        "time_uncertainty": "CDS labels MJD without a declared time scale or heliocentric correction; do not silently equate to COSMOGRAIL MHJD for sub-day timing",
        "source_reduction": "Difference imaging, Gaia proper-motion alignment, GALFIT PSF extraction, common high-S/N reference, Gaia photometric zero point",
        "shared_reference_caution": "Published reduction uses a common reference image with substantial 2022 contribution; a chronological split is not an independently blinded reduction and can carry correlated reference-flux errors",
        "scientific_reuse": "Public CDS A&A catalogue for scientific reuse with original-author/publication/publisher and CDS/VizieR citation; no blanket CC license asserted",
        "usage_policy_url": "https://cds.unistra.fr/vizier-org/licences_vizier.html",
        "old_consumed_release_max_native_day": 58354.41787,
        "proposed_cut_status": "Metadata-only candidate. Parent researcher must freeze the final protocol before photometric use.",
        "target_values_disclosed_or_parsed": False,
        "published_summary_awareness": "The metadata search displayed the primary paper's abstract/method text and qualitative published variability findings. No catalogue photometry, plots, our predictions or confirmatory scores were inspected. This is public retrospective withholding, not complete scientific blinding.",
        "files": rows,
        "failed_fetch_metadata": errors,
    }
    receipt_path.write_text(json.dumps(receipt, indent=2) + "\n", encoding="utf-8", newline="\n")
    return receipt


def split_r_without_disclosing_targets() -> dict:
    """Create the parent-selected split using timestamps, preserving row bytes."""
    original = SEALED / "he0435r.dat"
    lines = original.read_bytes().splitlines(keepends=True)
    cut = 58484.0
    calibration = [line for line in lines if line.strip() and float(line[:10]) < cut]
    target = [line for line in lines if line.strip() and float(line[:10]) >= cut]
    calibration_dir = OUT / "calibration"
    calibration_dir.mkdir(exist_ok=True)
    calibration_path = calibration_dir / "he0435r_pre2019.dat"
    target_path = SEALED / "he0435r_2019onward.dat"
    calibration_path.write_bytes(b"".join(calibration))
    target_path.write_bytes(b"".join(target))
    assert b"".join(calibration + target) == original.read_bytes(), "Chronological row concatenation must retain original bytes"
    result = {
        "cut_mjd": cut,
        "cut_utc": "2019-01-01T00:00:00",
        "selection_basis": "Parent instruction; timestamp-only conservative cut beyond all consumed COSMOGRAIL observations",
        "primary_test": "Condition on new image C as the declared driver; predict held-out images A, B and D",
        "secondary_test": "Forecast all four images without conditioning on new observations",
        "driver_selection_basis": "Old consumed development-set cross-validation, before LCO target inspection",
        "development": {
            "relative_path": calibration_path.relative_to(ROOT).as_posix(),
            "rows": len(calibration),
            "sha256": hashlib.sha256(calibration_path.read_bytes()).hexdigest(),
            "photometric_use": "Authorized for instrument calibration and model development",
        },
        "confirmatory": {
            "relative_path": target_path.relative_to(ROOT).as_posix(),
            "rows": len(target),
            "sha256": hashlib.sha256(target_path.read_bytes()).hexdigest(),
            "photometric_use": "Sealed; no photometric fields parsed or inspected during acquisition/splitting",
        },
        "whitespace_split_columns": ["MJD", "mag_A", "sigma_mag_A", "mag_B", "sigma_mag_B", "mag_C", "sigma_mag_C", "mag_D", "sigma_mag_D", "Band"],
        "fixed_width_one_based_columns": {"MJD": [1, 10], "mag_A": [12, 18], "sigma_mag_A": [20, 25], "mag_B": [27, 33], "sigma_mag_B": [35, 40], "mag_C": [42, 48], "sigma_mag_C": [50, 55], "mag_D": [57, 63], "sigma_mag_D": [65, 70], "Band": [72, 72]},
        "photometric_convention": "Apparent magnitudes, Gaia-star transformed R zero point; quoted marginal 1-sigma magnitude errors. Relative to arbitrary COSMOGRAIL magnitudes a per-image instrument/reduction offset is needed, not an assumed common numeric zero point.",
        "time_convention": "MJD as declared by CDS; no time scale or heliocentric convention declared. Do not silently identify with COSMOGRAIL MHJD; common half-day-origin ambiguity retained in time-map nuisance/model scope.",
        "v_status": "Original V file retained sealed, all 80 rows; unopened and unscored",
        "published_summary_awareness": "Primary paper abstract/method text exposed qualitative variability findings and aggregate published modeling results, but no catalogue targets, plots or confirmatory scores were inspected. The study is retrospective withholding with publication-level awareness.",
        "shared_reference": "Common LCO DIA reference image includes substantial 2022 contribution. Common reference/reduction errors may correlate calibration and test rows. Independent photons/new post-2018 epochs are established; an independently blinded reduction is not.",
    }
    (OUT / "r_confirmatory_split.json").write_text(json.dumps(result, indent=2) + "\n", encoding="utf-8", newline="\n")
    return result


if __name__ == "__main__":
    receipt = acquire()
    print(json.dumps({name: {key: value for key, value in item.items() if key != "all_observation_mjd"} for name, item in receipt["files"].items()}, indent=2))
    print(json.dumps(split_r_without_disclosing_targets(), indent=2))
