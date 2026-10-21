# Codex task: integrate the proven repair and complete the measurement response

We are repairing the HE 0435 calculation before claiming precision validation.
BHSM remains the governing framework. The reference and R1-inclusive columns
are controls; their equivalence to a matched effective-mass representation at
retained order is not a competition between two fundamental theories.

## Start from the delivered work

Repository: ncarberry64/Manuscript-Generation, existing draft PR #23.
Reviewed immutable base: 7ec9cea7b4d1c642e04407a6ef2f0a4eaef1baf9.

Read current repository instructions and inspect status/head. Preserve newer
work. Read HE0435_RESPONSE_REPAIR_REPORT.md, the experiment receipts, and
mask_patch/docs/CODEX_HANDOFF.md. Apply the supplied ten-file additive patch in
an isolated worktree. Run its four regression groups and canonical native
replay. Check the original input hashes and actual pipeline versions.

Expected mask-only G1 reference RMS is 0.448457087207 mas; 303 valid slopes are
restored. All sixteen absolute profile checks still fail. The original frozen
experiment remains immutable. Commit the verified correction and its provenance
through the already authorized Git workflow, without merging main or submitting
the manuscript.

## Continue from the new diagnosis

Do not restart with arbitrary PSF flexibility. Executed trials show:

- A positive field-star response partly transfers and gives raw RMS ≈ 0.362 mas,
  but its held star and target profiles fail.
- A free transverse spatial response is not stably constrained by the nearly collinear
  calibration-star positions.
- Finer field-only response overfits; twelve local quasar training stamps remove
  the rank defect but their held-exposure response still fails.
- Dated STPSF with fitted common color/diffusion remains inadequate.
- Eight bounded common pupil-phase modes improve some profiles but still fail
  exact nonlinear held-field checks.
- Remaining field residuals are predominantly at 2–6 native pixels from the star,
  not in the restored core.

The detector-history check has now been executed. Read the definitive
`ramp_diagnostic/native_pipeline/common_time_field_stcal118_results.json` and
its adjacent `NATIVE_READ_HISTORY_FINDINGS.md` reproduction record. Do not repeat the preliminary
approximate group calculation as though full processing were still pending.

The original star core retains only two usable groups while the 3–6 pixel
wings retain all nine. Reconstructing and calibrating the first two individual
reads gives every retained pixel the same 10.73677 s interval. On 289 identical
pixels with a fixed dated STPSF, raw normalized profile residual falls from
7.6653% to 5.1787% (32.44%). But early p = 0.00984653 remains below the declared
0.01 threshold; formal centroid errors are about three times larger, and several
wing residuals worsen. This is not an adequate precision response yet.

Numerical provenance is resolved: stcal 1.18.0 ramp fitting reproduces the
archived field rate/ERR to 8.97e-8 / 1.65e-7 maximum relative difference. The
initial 1.19.1 mismatch is retained and occurs at jump-segment pixels. The
archive does not record its stcal version, so numerical matching does not prove
an unrecorded installation history. Preserve both replays and their gates.

Continue from the measured read-history difference and the profile's sensitivity
to the extraction interval. Use the already consumed F115W exposures to assess a
response with the actual pixel-by-pixel accumulation
and group-selection history. Test the smallest physically supported extension,
with flux conservation where charge redistribution is proposed. A common-time
image is a diagnostic, not a license to discard precision until a test passes.
The remaining wing mismatch and the shared photons/reference covariance must
be accounted for. Do not assert a unique brighter-fatter or charge-migration
coefficient from this one-star comparison.

## Required functional calculation

Implement the native-pixel observational operator explicitly:

    expected pixels = detector-and-optical response(
        BHSM-propagated sky brightness,
        source spectrum and extended source/galaxy light,
        exposure, pixel and calibrated ramp response)
      + detector/sky background.

Use a coherent extended-light scene across exposures instead of relying only
on independent quadratic backgrounds near every quasar. The physical lens
coordinates must not be the objective used to calibrate the detector or PSF.
During measurement development, point positions may be independent measurement
nuisances; keep their relationship to the physical prediction explicit.

Identify the smallest response extension supported by calibration data. Compute
any proposed time/illumination contribution from reads and instrument
information, and assess it jointly with the local optical response and
extended scene. The current comparison supplies a concrete hypothesis, not a
unique identification of its mechanism. Do not claim an unmeasured
detector coefficient, invoke an extra gravitational term to fit a field-star
profile, or double-count a retained mass/scalar contribution.

Learned response parameters need their own uncertainty and correlations.
Propagate them through the joint likelihood into centroid covariance, using a
consistent rank/SVD convention. Do not report conditional fixed-PSF information
as fully calibrated precision. Native ERR scaling by a chosen multiplier is
not a replacement for this calculation.

## Avoid already settled implementation traps

- The rate-mask policy applies to audited native CAL_VER 2.0.1 rate products,
  not raw GROUPDQ or mosaics. Extend the release guard only with evidence.
- The four existing pipelines already use OLS_C and suppress_one_group=False.
- SAT does not identify which individual slopes use ZEROFRAME.
- STPSF 2.2.0 exploratory grids must use equal FFT/detector oversampling unless
  its IPC/grid metadata behavior is explicitly repaired and verified.
- Pixel-integrated ePSFs must not be integrated a second time.
- Preserve the original GWCS and fixed angular scale; do not hide a geometry
  mismatch with arbitrary per-image offsets or a free shear in registration.
- A learned response cache must verify version, OPD, position, wavelength,
  sampling and code/configuration provenance before reuse.

## Acceptance and reporting

Demonstrate measurement-model adequacy on calibration or deliberately held
development observations before another physical validation claim. A training
fit is insufficient. Retain raw residual maps, optimizer/rank diagnostics,
calibration uncertainty and exposure consistency. Synthetic recovery verifies
implementation only.

The already consumed F115W photons are development data. The remaining reserved
bands/observations must stay unopened until the physical predictions and the
complete measurement procedure are frozen. A re-reduction of the same photons
does not create another independent test.

Keep the 2e-10 physical guard and physical_admission=false unchanged unless the
actual separate admission conditions are satisfied. A classical optical fit is
not a microscopic U,V preparation solve. Keep the manuscript for Galaxies and
do not introduce a submission date.

Return the published repair commit/PR, the identified response mechanism (or
the precise remaining dependency if not identified), executed quantitative
results, and the actual validation state. Report tests/builds actually rerun.
Do not present incomplete calibration as a successful observational result.
