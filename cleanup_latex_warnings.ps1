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
    throw "Run this script from the Manuscript-Generation repository root."
}

$branch = (git branch --show-current).Trim()
if ($branch -ne "cosmology/topographic-universe-paper") {
    throw "Expected branch cosmology/topographic-universe-paper, found $branch"
}

Write-Host "== Fix PDF-string warnings in section titles =="

$sec2 = Get-Content "manuscript\sections\02_hyperspherical_parent.tex" -Raw
$sec2 = $sec2.Replace(
    '\subsection{Exact light-cone form of the \(n=1\) mode}',
    '\subsection{Exact light-cone form of the n=1 mode}'
)
Set-Content "manuscript\sections\02_hyperspherical_parent.tex" -Value ($sec2.TrimEnd("`r","`n") + "`n") -Encoding UTF8

$sec7 = Get-Content "manuscript\sections\07_predictions_falsification.tex" -Raw
$sec7 = $sec7.Replace(
    '\subsection{High-\(n\) recovery}',
    '\subsection{High-n recovery}'
)
Set-Content "manuscript\sections\07_predictions_falsification.tex" -Value ($sec7.TrimEnd("`r","`n") + "`n") -Encoding UTF8

Write-Host "== Fix BibTeX warning for the preprint entry =="

$bib = Get-Content "manuscript\references.bib" -Raw
$pattern = '(?s)@article\{Carberry2026Topographic,.*?\n\}'
$replacement = @'
@misc{Carberry2026Topographic,
  author = {Norman P. Carberry},
  title = {A Hyperspherical Scalar--Topographic Framework for Late-Time Cosmological Anomalies},
  year = {2026},
  howpublished = {Preprints.org},
  note = {Posted 20 January 2026, doi:10.20944/preprints202601.1427.v1}
}
'@
$bib2 = [regex]::Replace($bib, $pattern, $replacement, 1)
if ($bib2 -eq $bib) {
    Write-Warning "Carberry2026Topographic entry was not replaced; inspect references.bib manually."
} else {
    Set-Content "manuscript\references.bib" -Value ($bib2.TrimEnd("`r","`n") + "`n") -Encoding UTF8
}

Write-Host "== Fix overfull boxed result in R1 section =="

$r1 = Get-Content "manuscript\sections\09_results_R1.tex" -Raw
$old = @'
\[
\boxed{
\text{enhanced horizon-scale }n=1\text{ susceptibility}
+
\text{rapid high-}n\text{ recovery}
+
\text{few-percent late-time growth shift}.
}
\]
'@
$new = @'
\[
\boxed{
\begin{gathered}
\text{enhanced horizon-scale }n=1\text{ susceptibility}\\
+\ \text{rapid high-}n\text{ recovery}\\
+\ \text{few-percent late-time growth shift}
\end{gathered}
}
\]
'@
if ($r1.Contains($old)) {
    $r1 = $r1.Replace($old, $new)
} else {
    Write-Warning "Long boxed R1 result was not found exactly; inspect 09_results_R1.tex if the hbox warning remains."
}
Set-Content "manuscript\sections\09_results_R1.tex" -Value ($r1.TrimEnd("`r","`n") + "`n") -Encoding UTF8

Write-Host "== Re-run scientific validation =="

Invoke-Checked "python" @("-m", "pytest", "-q")

$digest = python -c "import json,hashlib; p='preregistration/prediction_manifest.json'; o=json.load(open(p,encoding='utf-8-sig')); b=json.dumps(o,sort_keys=True,separators=(',',':'),ensure_ascii=False).encode(); print(hashlib.sha256(b).hexdigest())"
$digest = $digest.Trim()
Write-Host "manifest sha256=$digest"
if ($digest -ne "0a995d15171f3b133edaf203869329e9eec7913a9ce8e23e8469175bc7538bfa") {
    throw "Frozen manifest hash changed unexpectedly."
}

Write-Host "== Clean and rebuild manuscript =="

Push-Location manuscript
try {
    Invoke-Checked "latexmk" @("-C")
    Invoke-Checked "latexmk" @("-pdf", "-interaction=nonstopmode", "-halt-on-error", "main.tex")
}
finally {
    Pop-Location
}

Write-Host "== Inspect remaining warnings =="

$log = Get-Content "manuscript\main.log" -Raw
$warningLines = $log -split "`r?`n" | Where-Object {
    $_ -match "Warning" -or $_ -match "Overfull"
}
if ($warningLines) {
    $warningLines | ForEach-Object { Write-Host $_ }
} else {
    Write-Host "No LaTeX warnings or overfull boxes detected in main.log."
}

Write-Host "== Git validation =="

Invoke-Checked "git" @("diff", "--check")

Write-Host "== Commit and push =="

git add manuscript\sections\02_hyperspherical_parent.tex `
        manuscript\sections\07_predictions_falsification.tex `
        manuscript\sections\09_results_R1.tex `
        manuscript\references.bib `
        manuscript\main.pdf

$pending = git status --porcelain
if ($pending) {
    Invoke-Checked "git" @("commit", "-m", "Clean LaTeX warnings and rebuild manuscript")
    Invoke-Checked "git" @("push")
} else {
    Write-Host "No changes to commit."
}

Write-Host "== Final status =="
git status
git log -1 --oneline
