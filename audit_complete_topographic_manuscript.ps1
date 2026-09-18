$ErrorActionPreference = 'Stop'

function Invoke-Checked {
    param(
        [Parameter(Mandatory=$true)][string]$Command,
        [Parameter(Mandatory=$true)][string[]]$Arguments
    )
    & $Command @Arguments
    if ($LASTEXITCODE -ne 0) {
        throw "$Command failed with exit code $LASTEXITCODE"
    }
}

if (-not (Test-Path ".git")) {
    throw "Run this script from the root of C:\Users\carbe\Manuscript-Generation."
}

$branch = (git branch --show-current).Trim()
if ($branch -ne "cosmology/topographic-universe-paper") {
    throw "Expected branch cosmology/topographic-universe-paper, found $branch"
}

Write-Host "== Repository hygiene =="

@"
__pycache__/
*.py[cod]
.pytest_cache/
.venv/
venv/

# LaTeX temporary/build products
*.aux
*.bbl
*.blg
*.fdb_latexmk
*.fls
*.log
*.out
*.synctex.gz
*.toc
manuscript/build/
"@ | Set-Content .gitignore -Encoding UTF8

@"
* text=auto eol=lf
*.py text eol=lf
*.tex text eol=lf
*.bib text eol=lf
*.md text eol=lf
*.json text eol=lf
*.yml text eol=lf
*.yaml text eol=lf
*.ps1 text eol=crlf
"@ | Set-Content .gitattributes -Encoding UTF8

git rm -r --cached --ignore-unmatch code/__pycache__ | Out-Host
Remove-Item -Recurse -Force code\__pycache__ -ErrorAction SilentlyContinue
Remove-Item -Recurse -Force .pytest_cache -ErrorAction SilentlyContinue

Write-Host "== Add validation tests =="

New-Item -ItemType Directory -Force -Path tests | Out-Null

@'
from __future__ import annotations

import hashlib
import json
import math
import sys
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
CODE = ROOT / "code"
if str(CODE) not in sys.path:
    sys.path.insert(0, str(CODE))

import background
import growth
import perturbations

FROZEN_MANIFEST_SHA256 = "0a995d15171f3b133edaf203869329e9eec7913a9ce8e23e8469175bc7538bfa"


def canonical_bytes(obj) -> bytes:
    return json.dumps(
        obj, sort_keys=True, separators=(",", ":"), ensure_ascii=False
    ).encode("utf-8")


def test_manifest_integrity():
    manifest_path = ROOT / "preregistration" / "prediction_manifest.json"
    obj = json.loads(manifest_path.read_text(encoding="utf-8-sig"))
    digest = hashlib.sha256(canonical_bytes(obj)).hexdigest()
    assert digest == FROZEN_MANIFEST_SHA256


def test_r1_background_reproduces_and_is_admissible():
    v0 = background.shoot_v0()
    assert math.isclose(v0, 2.3986073449708067, rel_tol=0.0, abs_tol=5e-10)

    sol = background.integrate(v0)
    assert sol.success

    zs = [10.0, 3.0, 2.1, 1.5, 1.0, 0.8, 0.5, 0.3, 0.0]
    rows = [background.diagnostics(sol, v0, z) for z in zs]

    assert all(row["D_kin"] > 0.0 for row in rows)
    assert max(row["friedmann_residual"] for row in rows) < 1.0e-9
    assert abs(rows[-1]["E"] - 1.0) < 1.0e-8


def test_harmonic_filter_selects_n1_only():
    ell = 0.40
    assert perturbations.selective_filter_condition(ell)

    # The sign of F_n = 1 - ell^2 n(n+2) determines whether chi softens
    # (positive) or stiffens (negative) the harmonic.
    f1 = 1.0 - ell * ell * perturbations.nu(1)
    assert f1 > 0.0

    for n in range(2, 25):
        fn = 1.0 - ell * ell * perturbations.nu(n)
        assert fn < 0.0


def test_r1_growth_forecast_reproduces():
    gr1 = growth.integrate_growth(growth.growth_rhs_r1)
    gl = growth.integrate_growth(growth.growth_rhs_lcdm)

    expected = {
        2.1: 0.9934646335007011,
        1.5: 0.9886863442761461,
        1.0: 0.9816359978859575,
        0.8: 0.9779545909356957,
        0.5: 0.9727064416234715,
        0.3: 0.9716770504951750,
        0.0: 0.9848217331453974,
    }

    for z, target in expected.items():
        N = math.log(1.0 / (1.0 + z))
        _, dpr = gr1.sol(N)
        _, dpl = gl.sol(N)
        ratio = float(dpr / dpl)
        assert math.isclose(ratio, target, rel_tol=0.0, abs_tol=2e-8)
'@ | Set-Content tests\test_reference_branch.py -Encoding UTF8

Write-Host "== Add reproducible figures =="

@'
"""Generate manuscript figures from the canonical R1 implementation.

No observational comparison data are read by this script.
"""
from __future__ import annotations

import math
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np

import background
import growth
import perturbations

ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / "manuscript" / "figures"
OUT.mkdir(parents=True, exist_ok=True)


def save_alpha_b():
    v0 = background.shoot_v0()
    sol = background.integrate(v0)
    zs = np.linspace(0.0, 10.0, 300)
    vals = [background.diagnostics(sol, v0, float(z))["alpha_B"] for z in zs]

    fig, ax = plt.subplots(figsize=(7.2, 4.6))
    ax.plot(zs, vals)
    ax.set_xlabel("Redshift z")
    ax.set_ylabel(r"$\alpha_B$")
    ax.set_title("Reference Branch R1: late-time kinetic braiding")
    ax.set_yscale("log")
    fig.tight_layout()
    fig.savefig(OUT / "r1_alphaB.png", dpi=220)
    plt.close(fig)


def save_omega_t():
    v0 = background.shoot_v0()
    sol = background.integrate(v0)
    zs = np.linspace(0.0, 10.0, 300)
    vals = [background.diagnostics(sol, v0, float(z))["Omega_T"] for z in zs]

    fig, ax = plt.subplots(figsize=(7.2, 4.6))
    ax.plot(zs, vals)
    ax.set_xlabel("Redshift z")
    ax.set_ylabel(r"$\Omega_T$")
    ax.set_title("Reference Branch R1: scalar energy fraction")
    fig.tight_layout()
    fig.savefig(OUT / "r1_omegaT.png", dpi=220)
    plt.close(fig)


def save_growth():
    gr1 = growth.integrate_growth(growth.growth_rhs_r1)
    gl = growth.integrate_growth(growth.growth_rhs_lcdm)
    zs = np.linspace(0.0, 2.1, 250)
    vals = []
    for z in zs:
        N = math.log(1.0 / (1.0 + float(z)))
        _, dpr = gr1.sol(N)
        _, dpl = gl.sol(N)
        vals.append(float(dpr / dpl))

    fig, ax = plt.subplots(figsize=(7.2, 4.6))
    ax.plot(zs, vals)
    ax.axhline(1.0, linewidth=1.0)
    ax.set_xlabel("Redshift z")
    ax.set_ylabel(r"$f\sigma_{8,\mathrm{R1}}/f\sigma_{8,\Lambda\mathrm{CDM}}$")
    ax.set_title("R1 prospective growth ratio")
    fig.tight_layout()
    fig.savefig(OUT / "r1_growth_ratio.png", dpi=220)
    plt.close(fig)


def save_filter():
    ell = 0.40
    ns = np.arange(1, 11)
    f = np.array([1.0 - ell * ell * perturbations.nu(int(n)) for n in ns])

    fig, ax = plt.subplots(figsize=(7.2, 4.6))
    ax.plot(ns, f, marker="o")
    ax.axhline(0.0, linewidth=1.0)
    ax.set_xlabel("S3 harmonic index n")
    ax.set_ylabel(r"$1-\ell^2 n(n+2)$")
    ax.set_title(r"Harmonic filter sign for $\ell=0.40$")
    ax.set_xticks(ns)
    fig.tight_layout()
    fig.savefig(OUT / "harmonic_filter.png", dpi=220)
    plt.close(fig)


def main():
    save_alpha_b()
    save_omega_t()
    save_growth()
    save_filter()
    print(f"Figures written to {OUT}")


if __name__ == "__main__":
    main()
'@ | Set-Content code\figures.py -Encoding UTF8

@'
\section{Reference-Branch Figures}

\begin{figure}[ht]
\centering
\includegraphics[width=0.78\linewidth]{figures/r1_alphaB.png}
\caption{Kinetic-braiding function \(\alpha_B(z)\) for Reference Branch R1. The logarithmic vertical scale emphasizes the suppression of the braiding sector toward high redshift.}
\label{fig:r1alphab}
\end{figure}

\begin{figure}[ht]
\centering
\includegraphics[width=0.78\linewidth]{figures/r1_omegaT.png}
\caption{Scalar energy fraction \(\Omega_T(z)\) for R1. The scalar component is subdominant at high redshift and becomes dynamically important only at late times.}
\label{fig:r1omegat}
\end{figure}

\begin{figure}[ht]
\centering
\includegraphics[width=0.78\linewidth]{figures/harmonic_filter.png}
\caption{Sign of the Schur-response factor \(1-\ell^2 n(n+2)\) for the illustrative filtering value \(\ell=0.40\). The \(n=1\) branch lies on the softening side while every displayed \(n\ge2\) mode lies on the stiffening side.}
\label{fig:harmonicfilter}
\end{figure}

\begin{figure}[ht]
\centering
\includegraphics[width=0.78\linewidth]{figures/r1_growth_ratio.png}
\caption{Prospective R1 growth ratio relative to a matched curved-\(\Lambda\)CDM background with the same early normalization.}
\label{fig:r1growth}
\end{figure}
'@ | Set-Content manuscript\sections\09a_reference_figures.tex -Encoding UTF8

$main = Get-Content manuscript\main.tex -Raw
if ($main -notmatch '09a_reference_figures') {
    $main = $main.Replace(
        '\input{sections/09_results_R1}',
        "\input{sections/09_results_R1}`n\input{sections/09a_reference_figures}"
    )
    Set-Content manuscript\main.tex -Value $main -Encoding UTF8
}

Write-Host "== Add build/audit documentation =="

@'
# Build and validation

## Scientific reproducibility

The current R1 reference state is anchored by:

- Git parent commit containing the first populated manuscript: `3923570`
- preregistration manifest SHA-256:
  `0a995d15171f3b133edaf203869329e9eec7913a9ce8e23e8469175bc7538bfa`

The manifest test intentionally fails if the frozen JSON changes without an explicit protocol/version update.

## Local validation

From the repository root:

```powershell
python -m pip install -r requirements.txt
python -m pip install pytest matplotlib
python code\background.py
python code\growth.py
python -m pytest -q
python code\figures.py
```

## LaTeX build

If `latexmk` is installed:

```powershell
cd manuscript
latexmk -pdf -interaction=nonstopmode -halt-on-error main.tex
```

Otherwise use `pdflatex` + `bibtex`:

```powershell
cd manuscript
pdflatex -interaction=nonstopmode -halt-on-error main.tex
bibtex main
pdflatex -interaction=nonstopmode -halt-on-error main.tex
pdflatex -interaction=nonstopmode -halt-on-error main.tex
```

## Claim discipline

See:

- `docs/claim_status.md`
- `docs/derivation_ledger.md`

Reference Branch R1 is exploratory and must not be described as a unique BHSM prediction unless its currently free bridge parameters are independently derived.
'@ | Set-Content docs\build_and_validation.md -Encoding UTF8

Write-Host "== Add CI validation =="

New-Item -ItemType Directory -Force -Path .github\workflows | Out-Null

@'
name: validate-manuscript

on:
  push:
  pull_request:

jobs:
  numerical-validation:
    runs-on: ubuntu-latest
    steps:
      - uses: actions/checkout@v4
      - uses: actions/setup-python@v5
        with:
          python-version: "3.12"
      - run: python -m pip install --upgrade pip
      - run: python -m pip install -r requirements.txt pytest matplotlib
      - run: python -m pytest -q
      - run: python code/figures.py

  latex-build:
    runs-on: ubuntu-latest
    steps:
      - uses: actions/checkout@v4
      - run: sudo apt-get update
      - run: sudo apt-get install -y latexmk texlive-latex-extra texlive-fonts-recommended
      - run: |
          cd manuscript
          latexmk -pdf -interaction=nonstopmode -halt-on-error main.tex
'@ | Set-Content .github\workflows\validate.yml -Encoding UTF8

Write-Host "== Ensure local dependencies =="

python -c "import numpy, scipy, matplotlib, pytest" 2>$null
if ($LASTEXITCODE -ne 0) {
    Invoke-Checked "python" @("-m", "pip", "install", "-r", "requirements.txt")
    Invoke-Checked "python" @("-m", "pip", "install", "pytest", "matplotlib")
}

Write-Host "== Run numerical tests =="
Invoke-Checked "python" @("-m", "pytest", "-q")

Write-Host "== Generate figures =="
Invoke-Checked "python" @("code\figures.py")

Write-Host "== Check frozen manifest digest =="
$digest = python -c "import json,hashlib; p='preregistration/prediction_manifest.json'; o=json.load(open(p,encoding='utf-8-sig')); b=json.dumps(o,sort_keys=True,separators=(',',':'),ensure_ascii=False).encode(); print(hashlib.sha256(b).hexdigest())"
$digest = $digest.Trim()
Write-Host "manifest sha256=$digest"
if ($digest -ne "0a995d15171f3b133edaf203869329e9eec7913a9ce8e23e8469175bc7538bfa") {
    throw "Frozen manifest hash changed unexpectedly."
}

Write-Host "== Compile manuscript if TeX is installed =="

$latexmk = Get-Command latexmk -ErrorAction SilentlyContinue
$pdflatex = Get-Command pdflatex -ErrorAction SilentlyContinue

if ($latexmk) {
    Push-Location manuscript
    try {
        Invoke-Checked "latexmk" @("-pdf", "-interaction=nonstopmode", "-halt-on-error", "main.tex")
    }
    finally {
        Pop-Location
    }
}
elseif ($pdflatex) {
    Push-Location manuscript
    try {
        Invoke-Checked "pdflatex" @("-interaction=nonstopmode", "-halt-on-error", "main.tex")
        if (Get-Command bibtex -ErrorAction SilentlyContinue) {
            Invoke-Checked "bibtex" @("main")
            Invoke-Checked "pdflatex" @("-interaction=nonstopmode", "-halt-on-error", "main.tex")
            Invoke-Checked "pdflatex" @("-interaction=nonstopmode", "-halt-on-error", "main.tex")
        }
    }
    finally {
        Pop-Location
    }
}
else {
    Write-Host "TeX not found locally; compilation skipped. GitHub Actions will perform the LaTeX build after push."
}

Write-Host "== Git diff check =="
Invoke-Checked "git" @("diff", "--check")

Write-Host "== Commit and push =="

git add .
$pending = git status --porcelain
if ($pending) {
    Invoke-Checked "git" @("commit", "-m", "Add reproducibility tests figures and manuscript validation")
    Invoke-Checked "git" @("push")
} else {
    Write-Host "No changes to commit."
}

Write-Host "== Final status =="
git status
git log -1 --oneline
