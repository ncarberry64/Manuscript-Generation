# R1 coupled environmental-state manuscript integration

The manuscript now distinguishes environmental realization of the two-component
temporal state from selection of its spatial profile or axis. The new section
follows the reference propagator and precedes observable-specific kernels;
the introduction, claim ledger and conclusion refer to the same distinction.

Worktree: `C:\Users\carbe\Manuscript-Generation\.worktrees\r1-coupled-environment-integration`

Branch: `codex/r1-coupled-environment-integration`

Base: `e32b13a`, the committed coupled-environment provenance branch, which
includes the submission integration. The original working directory's dirty
manuscript files were not copied or modified. No remote branch was changed.

## Corrections to the supplied draft

1. The supplied table describes `(zeta,zeta_dot)`, not `(q2,Pi2)`. The caption
   now says so explicitly. The audit also evaluates the full physical
   six-state propagator and verifies the positive invertible conversion
   `U_XE(q2,Pi2)=diag(1,2*a^3*G_S2) U_XE(zeta,zeta_dot)`.
2. `Pi2` is a reference-normalized coordinate. It is not identified with the
   full coupled canonical `p_zeta`.
3. Aligned environmental coefficients preserve one common spatial profile
   when the initial topographic state is zero or aligned with that profile.
   An independent initial profile can give rank two. A rank-one snapshot does
   not prove that one fixed profile applies throughout the history.
4. Temporal rank is asserted at the eight audited post-anchor epochs, not
   over an unproved continuous interval or at the anchor itself.
5. The text retains the distinction between existing reference kernels and
   a general environmental realization; no forecast kernel is silently replaced.

The original insertion draft and JSON handoff are preserved under
`docs/research_packets/2026-09-24/`. They supplied the proposed integration;
their numerical and scientific statements were checked rather than assumed.

## Scientific reproduction

The coordinate map reproduces the retained two-column seed with maximum
absolute difference `6.01005e-14` and inverse residual `2.22045e-16`.
The current-code reintegration reproduces the stored canonical propagator at
`z=1.5` exactly at the retained binary64 output precision (maximum difference
zero). This is a replay of the frozen solver, not a new independent integrator
or an improvement to its stated accuracy.

The environmental and matter-only transfers have rank two at
`z=1.5,1,0.8,0.5,0.3,0.1,0.02,0`, in both declared output bases. The audit
evaluates basis responses and does not choose a physical environmental state
from observations. It does not derive that state from earlier cosmological
history, select a global axis, or close the microscopic normalization bridge.

The frozen R1 manifest has raw SHA256
`afaf697100a70d43eab2a6086931bad5efe1ced02f9a56d422b709896be51113`.
The corresponding canonical-JSON hash remains governed by the existing frozen
manifest test. No owner code, observational results, reference kernel, or
preregistration parameter was changed.

## Validation and reproduction commands

```powershell
Set-Location 'C:\Users\carbe\Manuscript-Generation\.worktrees\r1-coupled-environment-integration'
$env:PYTHONPATH = (Join-Path (Get-Location) 'code')
& 'C:\Python314\python.exe' scripts/audit_r1_coupled_environmental_state.py --repo . --output artifacts/provenance/R1_COUPLED_ENVIRONMENT_INTEGRATION_REPLAY.json
& 'C:\Python314\python.exe' -m pytest -q tests/test_coupled_environmental_state.py tests/test_action_native_matter_n2.py --junitxml=artifacts/coupled_environment_focused_pytest.xml
& 'C:\Python314\python.exe' -m pytest -q --junitxml=artifacts/coupled_environment_full_pytest.xml -o cache_dir=manuscript/build/pytest-cache
Set-Location manuscript
latexmk -pdf -interaction=nonstopmode -halt-on-error main.tex
```

The first focused collection attempt lacked `PYTHONPATH=code` and failed to
import `background`; the corrected command passes all 11 focused tests.
The full-suite result, final log checks, source hashes and status lines are
recorded in `artifacts/provenance/R1_COUPLED_ENVIRONMENT_INTEGRATION_RECEIPT.json`.
The new tests check both temporal bases and demonstrate the spatial-profile
counterexample; they do not alter any existing expected failure.

The 50-page PDF builds without undefined references/citations, LaTeX warnings
or overfull boxes. Pages 5, 11--13 and 30 were rendered and visually inspected:
the claim table, new equations/table, transition to existing kernels and
conclusion are readable without clipping or overlap.
