"""Record scrutiny evidence without running, weakening, or bypassing tests.

The long suite runs in the preceding integration worktree. Reuse is permitted
only when every recorded scientific/test/protocol input and collected test ID
matches. The new plotting code is executed separately, not covered by that run.
"""
from pathlib import Path
import hashlib
import importlib.metadata
import json
import re
import shutil
import subprocess
import xml.etree.ElementTree as ET
from pypdf import PdfReader

ROOT = Path(__file__).resolve().parents[1]
BASE = "083e14654d1b81b3853484b48e132bd5bc41567c"
TITLE = "Geometry Before Fields: Constraint-Reduced Cosmology and Optical Predictions in the Berger--Hopf Framework"

def digest(p): return hashlib.sha256(p.read_bytes()).hexdigest()
def load(p): return json.loads((ROOT / p).read_text(encoding="utf-8"))
def text(p): return (ROOT / p).read_text(encoding="utf-8", errors="replace")
def write(p, s): (ROOT / p).write_text(s, encoding="utf-8", newline="\n")

def main():
    eq = load("artifacts/provenance/FINAL_SCRUTINY_INPUT_EQUIVALENCE.json")
    other = Path(eq["reference_worktree"])
    assert all(digest(ROOT/p) == h == digest(other/p) for p,h in eq["input_sha256"].items())
    collection = {s.strip() for s in text("artifacts/final_scrutiny_collection.txt").splitlines() if "::" in s}
    prior_collection = {s.strip() for s in (other/"artifacts/coupled_environment_collection.txt").read_text().splitlines() if "::" in s}
    assert collection == prior_collection and len(collection) == 102
    suite = {"status": "RUNNING_OR_PENDING", "expected_count": 102,
             "source_worktree": str(other), "source_commit": eq["reference_commit"],
             "reuse_basis": "Identical 188 recorded scientific/test/protocol/sightline inputs and identical 102-test collection; new plotter exercised separately."}
    source_xml = other/"artifacts/coupled_environment_full_pytest.xml"
    if source_xml.exists():
        tree = ET.parse(source_xml)
        cases = tree.findall(".//testcase")
        actual = {c.get("classname").replace(".","/")+".py::"+c.get("name") for c in cases}
        failures = [c.get("name") for c in cases if c.find("failure") is not None or c.find("error") is not None]
        skipped = [c.get("name") for c in cases if c.find("skipped") is not None]
        allowed = ["test_dae_preserves_00_0i_constraints_from_index_consistent_seed"]
        xfails = [c.get("name") for c in cases if c.find("skipped") is not None and c.find("skipped").get("type") == "pytest.xfail"]
        ok = actual == collection and not failures and skipped == xfails == allowed
        suite.update(status="PASS" if ok else "FAIL", count=len(cases), passed=len(cases)-len(skipped)-len(failures),
                     failures=failures, expected_failures=xfails, unexpected_skips=[s for s in skipped if s not in allowed],
                     elapsed_seconds=tree.find("testsuite").get("time"), source_xml_sha256=digest(source_xml))
        shutil.copyfile(source_xml, ROOT/"artifacts/final_scrutiny_full_pytest.xml")
        shutil.copyfile(other/"artifacts/coupled_environment_full_pytest.log", ROOT/"artifacts/final_scrutiny_full_pytest.txt")
    focused = ET.parse(ROOT/"artifacts/final_scrutiny_focused_pytest.xml").findall(".//testcase")
    assert len(focused) == 14 and all(len(c) == 0 for c in focused)
    visual = load("artifacts/provenance/FINAL_SCRUTINY_VISUAL_AUDIT.json")
    assert visual["all_pages_visually_inspected"] and visual["pdf_sha256"] == digest(ROOT/"manuscript/main.pdf")
    pdf = PdfReader(ROOT/"manuscript/main.pdf")
    rebuilt = PdfReader(ROOT/"manuscript/build/final-source-rebuild/main.pdf")
    assert [p.extract_text() for p in pdf.pages] == [p.extract_text() for p in rebuilt.pages]
    for src, dest in [("manuscript/main.log", "artifacts/final_scrutiny_latex_build.txt"),
                      ("manuscript/build/final-source-rebuild/main.log", "artifacts/final_scrutiny_isolated_build.txt")]:
        assert not re.search(r"Overfull|undefined|Warning|! LaTeX Error", text(src))
        # Preserve the diagnostic content while keeping repository text clean.
        write(dest, "\n".join(line.rstrip() for line in text(src).splitlines()).rstrip()+"\n")
    main_tex = text("manuscript/main.tex")
    abstract = main_tex.split(r"\begin{abstract}")[1].split(r"\end{abstract}")[0]
    title_files = ["manuscript/main.tex", "README.md", "submission/UNIVERSE_SPECIAL_ISSUE_COVER_LETTER.md",
                   "submission/UNIVERSE_PLANNED_PAPER_ABSTRACT_250_WORDS.md"]
    assert all(TITLE in text(p) for p in title_files)
    inputs = re.findall(r"\\input\{([^}]+)\}", main_tex.split(r"\appendix")[0])
    figure_count = sum(text("manuscript/"+p+".tex").count(r"\begin{figure}") for p in inputs)
    assert figure_count == 6
    a = load("artifacts/provenance/FINAL_SCRUTINY_COUPLED_REPLAY.json")
    assert all(r["rank_rel_1e-8"] == r["q2_Pi2_rank_rel_1e-8"] == 2 for r in a["rows"])
    assert not a["retuning"] and not a["observational_environment_selection"]
    archived = load("artifacts/provenance/BHSM_R1_COUPLED_ENVIRONMENTAL_STATE_AUDIT_V1.json")["temporal_result"]["rows"]
    differences = []
    relative = []
    for old, new in zip(archived, a["rows"], strict=True):
        assert old["z"] == new["z"]
        for u,v in zip([old["s1"],old["s2"],old["matter_det"]],new["singular_values"]+[new["matter_delta_v_determinant"]]):
            differences.append(abs(u-v)); relative.append(abs(u-v)/abs(u))
    # Record the exact baseline-relative review list, without auto-staging it.
    names = subprocess.check_output(["git","diff","--name-only","850b1ab"], cwd=ROOT, text=True).splitlines()
    new = subprocess.check_output(["git","ls-files","--others","--exclude-standard"], cwd=ROOT, text=True).splitlines()
    outputs = sorted(set(names+new))
    statuses = dict(PUBLIC_BASELINE=BASE, TITLE_ALIGNED=True, STALE_MATTER_RADIATION_CLAIMS="CLEARED",
                    COUPLED_ENVIRONMENT_SECTION="INTEGRATED", MAIN_PUBLICATION_FIGURES=6,
                    FULL_TESTS=suite["status"], FOCUSED_TESTS="14 PASSED", LATEX_BUILD="PASS",
                    VISUAL_PAGE_AUDIT="PASS", ISOLATED_SOURCE_REBUILD="PASS", RETUNING=False,
                    READY_FOR_AUTHOR_FINAL_CONFIRMATION=suite["status"]=="PASS")
    receipt = dict(statuses=statuses, branch="review/cosmology-final-scrutiny-20260924", full_suite=suite,
                   page_count=len(pdf.pages), abstract_whitespace_words=len(abstract.split()),
                   references=len(re.findall(r"\\bibitem",text("manuscript/main.bbl"))),
                   scientific_input_count=eq["input_count"], coupled_replay=a,
                   archived_plot_values_vs_current_replay={"bitwise_identical":False,"max_absolute_difference":max(differences),"max_relative_difference":max(relative),"note":"The figure retains the original committed artifact; the current-code replay is separately preserved. Both give the same sampled rank conclusions. A trial 1e-12 absolute identity check failed; no frozen value or scientific gate was altered."},
                   raw_latex_log_sha256={p:digest(ROOT/p) for p in ["manuscript/main.log","manuscript/build/final-source-rebuild/main.log"]},
                   versions={n:importlib.metadata.version(n) for n in ["numpy","scipy","matplotlib","pytest","pypdf"]},
                   reviewed_changed_files=outputs,
                   hashes={p:digest(ROOT/p) for p in ["manuscript/main.pdf","submission/universe-source.zip",
                           "preregistration/prediction_manifest.json","code/final_scientific_figures.py",
                           "artifacts/final_scrutiny_focused_pytest.xml"]},
                   merge_performed=False, author_declarations_confirmed=False)
    write("artifacts/provenance/FINAL_SCRUTINY_VALIDATION_20260924.json",json.dumps(receipt,indent=2)+"\n")
    result = f"{suite.get('passed', 0)} passed, one documented strict expected failure" if suite["status"]=="PASS" else suite["status"]
    report = f"""# Final scrutiny validation — 2026-09-24

Branch: `review/cosmology-final-scrutiny-20260924`. Public baseline: `{BASE}`.
The coupled provenance and corrected integration were carried forward as `e0794fa` and `850b1ab`. The final review commit is obtained with `git rev-parse HEAD` on this branch; this report does not embed its own self-referential hash.

Title: {TITLE}

The review corrects stale matter/radiation claims to the audited frozen-R1 linear/quadratic result, separates temporal rank from spatial selection, identifies the retained GP equations as historical compatibility audits, and confines harmonic-filter extrapolations to the illustrative realization. Reference kernels and every frozen scientific input remain unchanged. No failed identity was patched and no observational fit was run.

The coupled table and figure use `(zeta,zeta_dot)` and a hatted transfer. Rank two is reproduced in both output bases at eight post-anchor epochs. Singular-value magnitudes are basis dependent; the anchor block is zero. `Pi2` is reference normalized and is not generally the coupled canonical momentum. Environmental alignment alone does not align an independent initial topographic state.

| Gate | Result |
|---|---|
| Title alignment / stale matter-radiation claims | PASS / CLEARED |
| Frozen scientific/test/protocol/sightline files | {eq['input_count']} byte-identical to the full-suite worktree |
| Focused action-native, environment and gradient checks | 14 passed |
| Full 102-test suite | {result} |
| Main figures | Six, plus two appendix figures |
| PDF / abstract / bibliography | {len(pdf.pages)} pages / {len(abstract.split())} whitespace-delimited words / {receipt['references']} references |
| LaTeX / isolated ZIP rebuild | PASS; no warnings, unresolved references/citations or overfull boxes |
| Visual page audit | PASS; every page rendered and visually checked |
| Retuning | FALSE |
| Ready for final author confirmation | {str(statuses['READY_FOR_AUTHOR_FINAL_CONFIRMATION']).upper()} |

The full suite is the run launched in `{other.name}`, not a second fresh run on this review branch. All 102 test identities and 188 recorded scientific inputs match exactly. Its evidence is imported only after completion and checked for unexpected failures/skips. The new plotter was run separately. The focused tests were rerun here; the gradient tests validate the committed audit artifacts, not a newly rerun principal-symbol integration. Both figure generators ran successfully; the four regenerated legacy PNGs are unchanged.

The environmental replay recovers both topographic seeds to {a['topographic_seed_recovery_max_abs_difference']:.6g}, inverts the seed basis to {a['seed_inverse_residual']:.6g}, and reproduces the saved z=1.5 propagator with maximum difference {a['stored_full_propagator_z1p5_max_abs_difference']:.6g}. This replays committed propagators through current constraint/basis code; it is not a new independent ODE integration or reconstruction of physical environmental initial data.

The original plotted singular values/determinants and the current replay are not bitwise identical: maximum absolute difference {max(differences):.6g}, maximum relative difference {max(relative):.6g}. A trial 1e-12 absolute identity check failed. Both artifacts are preserved without modifying their values, the figures use the declared original input, and both audits retain the same sampled rank-two conclusion. This is not reported as exact numerical reproduction of every archived scalar.

`FINAL_SCIENTIFIC_FIGURES_RECEIPT.json` hashes the analytic/artifact inputs and all new figures. `FINAL_SCRUTINY_VISUAL_AUDIT.json` records PDF and per-page render hashes. All pages were inspected in 100-dpi contact sheets, with individual detail views for revised geometry/text. Crowded analytic longitude labels were corrected and rechecked. The isolated rebuild has identical extracted text and page count. The source ZIP verifies every internal hash.

Existing Python 3.14 warnings about legacy non-raw TeX string escapes were observed during collection. No scientific owner was changed to suppress them. The LaTeX builds have zero warnings.

Open science remains explicit: independently reconstructed environmental anchor data, spatial profile/axis selection, microscopic and absolute local-source normalization, nonlinear-fluid/UV completion, joint observational covariance, exact radial BAO and survey-ready RSD. No observational signal or prospective gate success is claimed.

Author-only CRediT, funding and conflict declarations remain unfilled. No archive DOI, journal submission, external message or merge was performed. Review and merge approval remain with the author, as explicitly required by the final scrutiny handoff.

Reproduce: run `python code/figures.py`, `python code/final_scientific_figures.py`, the focused pytest files listed in the JUnit receipt, `latexmk -pdf -interaction=nonstopmode -halt-on-error main.tex` in `manuscript`, and `python scripts/package_submission.py`. With the sibling suite complete, run `python scripts/record_final_scrutiny_validation.py` to refresh the checked completion receipt. Full-suite reproduction: set `PYTHONPATH=code` and run `python -m pytest -q` (long integrations).
"""
    write("docs/FINAL_SCRUTINY_VALIDATION_20260924.md",report)
    readiness = f"""# Submission readiness

Title: {TITLE}

Current review branch: `review/cosmology-final-scrutiny-20260924`, based on public `{BASE}`. {len(pdf.pages)} pages; {len(abstract.split())} abstract words; six main and two appendix figures; {receipt['references']} references.

The corrected action-native quadratic stability and coupled environmental results are integrated at their audited scope. All 188 recorded scientific/test/protocol/sightline inputs match the integration worktree. Direction, amplitude, redshift law, R1 cosmology, sample cuts and p<0.01 gates are unchanged. RETUNING=false.

Validation: 14 focused tests pass; full suite: **{result}**. The full-suite evidence is reused only for identical scientific inputs and test identities; the new plotter was executed separately. Both figure generators ran. Clean LaTeX and isolated source-ZIP builds pass with no warnings or overfull boxes; all {len(pdf.pages)} pages were rendered and visually inspected. Full scope, provenance and limitations are in [the final scrutiny report](../docs/FINAL_SCRUTINY_VALIDATION_20260924.md).

All frozen SN results remain unchanged: Pantheon released +2.985, source/environment −0.883, fixed external-PV sensitivity −4.909; DES released +0.240, transfer +1.321. No coherent signal is established; fixed external PV has incomplete joint uncertainty, and aggregate tomography does not localize residual correlations. No observational likelihood was rerun in this manuscript pass.

Open science: environmental anchor reconstruction and spatial selection, microscopic/local-source normalization, nonlinear-fluid/UV control, complete joint covariance, radial BAO and survey-ready RSD, and independent prospective survey design. Historical missing localization statistics and original foreground hashes remain disclosed. The source ZIP is a manuscript build package, not the raw observational archive.

Author actions remain: confirm CRediT, funding and conflicts in `AUTHOR_CONFIRMATION_REQUIRED.md`; review retained claims and AI disclosure; verify special-issue routing, preprint relationship and submission declarations; approve any merge, deposit or submission. No permanent DOI or external submission has been created.

READY_FOR_AUTHOR_FINAL_CONFIRMATION = {str(statuses['READY_FOR_AUTHOR_FINAL_CONFIRMATION']).upper()}
"""
    write("submission/SUBMISSION_READINESS.md",readiness)
    print(json.dumps(statuses))

if __name__ == "__main__": main()
