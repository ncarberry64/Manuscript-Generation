param([Parameter(Mandatory=$true)][int]$TestProcessId)
$ErrorActionPreference = 'Stop'
$repoPath = Split-Path -Parent $PSScriptRoot
Set-Location -LiteralPath $repoPath
$logPath = Join-Path $repoPath 'artifacts/coupled_environment_completion.log'
try {
    $testProcess = Get-Process -Id $TestProcessId -ErrorAction SilentlyContinue
    if ($null -ne $testProcess) { $testProcess.WaitForExit() }
    & 'C:\Python314\python.exe' scripts/record_coupled_environment_integration.py *> $logPath
    if ($LASTEXITCODE -ne 0) { throw 'Receipt generation failed.' }
    $receipt = Get-Content -Raw -LiteralPath 'artifacts/provenance/R1_COUPLED_ENVIRONMENT_INTEGRATION_RECEIPT.json' | ConvertFrom-Json
    if ($receipt.statuses.MANUSCRIPT_INTEGRATION -eq 'PENDING_FULL_SUITE') {
        'Test process exited without a complete JUnit record; validation remains incomplete.' | Add-Content -LiteralPath $logPath
    }
    # Record evidence only. No Git commit, merge, push or manuscript edits occur here.
    $receipt.statuses.MANUSCRIPT_INTEGRATION | Add-Content -LiteralPath $logPath
} catch {
    $_ | Out-String | Add-Content -LiteralPath $logPath
    exit 1
}
