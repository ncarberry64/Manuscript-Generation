# Integration checkpoint: full validation running

Worktree: `C:\Users\carbe\Manuscript-Generation\.worktrees\r1-coupled-environment-integration`

Branch: `codex/r1-coupled-environment-integration`

Base: `e32b13ad86d1b66889a4a0483fed9ef20ae84561`

The corrected section, introduction, claim ledger and conclusion are integrated.
Both temporal bases and the matter-only rank are reproduced. All 11 focused
checks pass. The 50-page PDF builds with no LaTeX warnings, unresolved citations/
references or overfull boxes; changed pages were visually reviewed. Frozen
owner code, kernels, R1 parameters, observations and preregistration are unchanged.

The **102-test full suite is still running**. Its historical JUnit record took
6696.88 seconds, dominated by three unchanged distance-integral tests. Do not
restart this run merely because the progress line stays at 70% for a long time.
The current receipt says `MANUSCRIPT_INTEGRATION=PENDING_FULL_SUITE`, not PASS.

At checkpoint, the full-suite Python process is PID 13600. A hidden, one-shot
PowerShell helper waits for that process and reruns the evidence recorder when
it exits. The helper only updates the local receipt and completion log. It does
not change scientific files, commit, merge or push. A completed receipt counts
as PASS only if all collected identities are present, no tests fail, and the
only skip is the already registered curved-EFT strict expected failure.

Read:

```powershell
Set-Location 'C:\Users\carbe\Manuscript-Generation\.worktrees\r1-coupled-environment-integration'
Get-Content artifacts/coupled_environment_full_pytest.log -Tail 12
Get-Content artifacts/coupled_environment_completion.log -ErrorAction SilentlyContinue
Get-Content artifacts/provenance/R1_COUPLED_ENVIRONMENT_INTEGRATION_RECEIPT.json
```

If the suite finished but the helper did not update the receipt, run:

```powershell
& 'C:\Python314\python.exe' scripts/record_coupled_environment_integration.py
```

After PASS, review and commit the new full JUnit record and updated receipt on
this integration branch. Remote main remains untouched. If the suite fails,
preserve the failures and investigate without changing any frozen science or
turning unexpected failures into expected ones. If the machine shut down and
no JUnit exists, rerun the full command documented in
`docs/R1_COUPLED_ENVIRONMENT_INTEGRATION.md`.
