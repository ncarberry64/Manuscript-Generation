# Historical build notes

For current instructions see [REPRODUCIBILITY](REPRODUCIBILITY.md). Original notes follow for provenance.

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

- `docs/CLAIM_STATUS.md`
- `docs/derivation_ledger.md`

Reference Branch R1 is exploratory and must not be described as a unique BHSM prediction unless its currently free bridge parameters are independently derived.
