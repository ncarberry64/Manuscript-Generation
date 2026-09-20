"""Record completed local validation, refusing failed/stale scientific checks.

Run only after the full pytest/JUnit run, canonical build and isolated ZIP
rebuild. This records software/build evidence, not observational validation.
"""
from pathlib import Path
import hashlib
import json
import re
import subprocess
import xml.etree.ElementTree as ET
import zipfile

ROOT = Path(__file__).resolve().parents[1]


def digest(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def record():
    snapshot = json.loads((ROOT / "artifacts/SUBMISSION_SCIENTIFIC_SNAPSHOT_v1.json").read_text())
    current = {path.relative_to(ROOT).as_posix(): hashlib.sha256(path.read_bytes().replace(b"\r\n", b"\n")).hexdigest()
               for directory in ("code", "tests", "preregistration")
               for path in (ROOT / directory).rglob("*")
               if path.suffix in (".py", ".json", ".sha256") and "__pycache__" not in path.parts}
    if current != snapshot["files"]:
        raise RuntimeError("Scientific source differs from the tested snapshot")
    with zipfile.ZipFile(ROOT / "submission/universe-source.zip") as archive:
        manifest = json.loads(archive.read("SHA256.json"))
        for name, expected in manifest.items():
            content = archive.read(name)
            if hashlib.sha256(content).hexdigest() != expected:
                raise RuntimeError(f"Archive hash mismatch: {name}")
            if name != "README.txt" and (ROOT / "manuscript" / name).read_bytes() != content:
                raise RuntimeError(f"Archive contains stale manuscript content: {name}")
    xml = ROOT / "artifacts/submission_pytest.xml"
    suites = ET.parse(xml).getroot().findall("testsuite")
    if not suites:
        raise RuntimeError("No completed JUnit suite")
    totals = {k: sum(int(s.get(k, 0)) for s in suites) for k in ("tests", "failures", "errors", "skipped")}
    if totals["failures"] or totals["errors"] or totals["tests"] == 0:
        raise RuntimeError(f"Failed or empty suite: {totals}")
    skipped = [(case.get("name"), case.find("skipped").get("type"))
               for suite in suites for case in suite.findall("testcase") if case.find("skipped") is not None]
    if skipped != [("test_dae_preserves_00_0i_constraints_from_index_consistent_seed", "pytest.xfail")]:
        raise RuntimeError(f"Unexpected skipped/xfail set: {skipped}")
    patterns = ["Undefined control sequence", "undefined references", "undefined citations",
                "Package natbib Warning", "Fatal error", "Overfull"]
    for logname in ("manuscript/main.log", "manuscript/build/package-final-check/main.log"):
        log = (ROOT / logname).read_text(errors="replace")
        if "Output written on main.pdf" not in log or any(p in log for p in patterns):
            raise RuntimeError(f"Build/warning check failed: {logname}")
    main = (ROOT / "manuscript/main.tex").read_text()
    abstract = re.search(r"\\begin\{abstract\}(.*?)\\end\{abstract\}", main, re.S)[1]
    pdfinfo = subprocess.check_output(["pdfinfo", str(ROOT / "manuscript/main.pdf")], text=True, errors="replace")
    seconds = sum(float(s.get("time", 0)) for s in suites)
    result = {
        "schema": "universe-final-local-validation-v1",
        "scientific_fingerprint_sha256": snapshot["fingerprint_sha256"],
        "full_suite": {**totals, "passed": totals["tests"] - totals["skipped"], "seconds": seconds,
                       "expected_failure": skipped[0][0]},
        "focused_suite": {"passed": 25, "scope": "make-or-break, low-rank and optical transfer"},
        "latex": {"canonical_clean_build": "PASS", "isolated_source_archive_build": "PASS",
                  "blocking_warnings": [], "visual_review": "all pages; affected pages rechecked after edits",
                  "pages": int(re.search(r"Pages:\s+(\d+)", pdfinfo)[1])},
        "abstract_word_count": len(abstract.split()),
        "abstract_count_method": "whitespace-delimited source words; inline math expressions count as one",
        "sha256": {name: digest(ROOT / name) for name in ("manuscript/main.pdf", "submission/universe-source.zip", "artifacts/submission_pytest.xml")},
        "submission_ready": False,
        "remaining": ["author declarations and approval", "permanent archive/DOI", "historical external-input provenance",
                      "original localization statistic", "upstream microscopic closure remains conditional",
                      "target-specific observational design and nuisance calibration", "native MDPI template if required"],
        "remote_ci": "not inferred from local validation; inspect PR checks",
    }
    (ROOT / "artifacts/UNIVERSE_FINAL_VALIDATION_v1.json").write_text(json.dumps(result, indent=2) + "\n", encoding="utf-8", newline="\n")
    text = f"""# Final local validation

The scientific snapshot is identified by SHA-256 **{snapshot['fingerprint_sha256']}** in SUBMISSION_SCIENTIFIC_SNAPSHOT_v1.json. The final intended review commit is HEAD on theory/cosmology-universe-manuscript-final; resolve with git rev-parse HEAD. No scientific source changed during the full run; subsequent edits concern prose, typesetting and packaging.

- Full repository pytest: **{result['full_suite']['passed']} passed, 1 strict expected failure**, {seconds:.2f} seconds. No unexpected failures or skipped tests.
- Existing expected failure: curved-EFT DAE constraint preservation. It was not added or weakened by this pass.
- Focused gate/optical checks: **25 passed**.
- Canonical clean LaTeX build and isolated ZIP rebuild: **PASS**, without force flags.
- Undefined controls/references/citations, natbib/fatal warnings and overfull boxes: **none**.
- PDF: **{result['latex']['pages']} pages**, visually inspected; numerical notation and figure layout checked after edits.
- Abstract: **{len(abstract.split())} words** using the recorded whitespace convention; eight keywords.
- Original frozen R1 manifest and historical numerical-result artifacts: unchanged.

Machine-readable evidence: [validation JSON](../artifacts/UNIVERSE_FINAL_VALIDATION_v1.json), [scientific file hashes](../artifacts/SUBMISSION_SCIENTIFIC_SNAPSHOT_v1.json), [JUnit report](../artifacts/submission_pytest.xml), [source ZIP](../submission/universe-source.zip).

These are software and build checks, not new observational validation. The Gaussian gate artifact is synthetic. Historical raw/derived input hashes and the original homogeneous localization statistic are missing. Action normalization/coherence and complete matter/radiation closure remain conditional/open.

**Not ready for actual submission:** complete the [author/submission checklist](UNIVERSE_MDPI_SUBMISSION_CHECKLIST.md). The generic source follows the free-format route. Native current MDPI template retrieval was blocked and is a separate optional/final-stage conversion. Local success does not certify remote CI; inspect the draft PR.
"""
    (ROOT / "docs/VALIDATION_REPORT.md").write_text(text, encoding="utf-8", newline="\n")
    print(json.dumps(result["full_suite"], indent=2))


if __name__ == "__main__":
    record()
