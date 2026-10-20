# HE 0435: HST constraint, revised uncertainty and unused observations

This continues the original [observed-light checkpoint](STRONG_LENSING_OBSERVED_20261020.md). Its source, predictions, 85.33% coverage failure and all numerical artifacts remain unchanged. The old 1002 COSMOGRAIL epochs are now consumed development data. They are never called a new validation sample.

The regenerated real-MDPI Galaxies draft has **32 pages**, **559415 bytes**, PDF SHA256 `fafe605fa9cd9b74c87d2fc760aeb1ed6585c3609fd63ff747390542963675c6`. All pages were rendered and reviewed; the final LaTeX log has no warnings, undefined references or box warnings. The [publication manifest](../artifacts/strong_lensing_revision_20261020/revision_publication_manifest.json) separates these release checks from the physical-model failures. The immutable science archive is revision **6b9a0f06e9850e8c796762e86334f7b9e970f0fb**. The final internal scientific/coherence review passes; this is not external human peer review or physical admission.

## Actual HST application and source revision

The original unattenuated central/A macro ratio, **0.063197**, was applied to real ACS F814W and WFC3 F160W science pixels. Four empirical PSFs, registration rotation, free central position, arbitrary nonnegative galaxy/host/nuclear light, and noise/normalization controls are explicit. The F814W numerical conditional envelope is **0.0293492**, or **0.0412312** with fivefold marginal noise and a 20% lower A normalization. Conservatively rounded bounds are 0.030 and 0.042; F160W is inconclusive. These are conditional on empirical PSFs, estimated marginal errors, sky subtraction and the unattenuated macro-to-observed flux relation. The finite PSF/continuous centroid search is numerical, not a globally certified arbitrary-PSF bound. Stress tests are not calibrated confidence probabilities. See [HST methods and pixels](strong_lensing_revision_20261020/HST_CENTRAL_IMAGE.md).

Five positive outer-plus-compact Plummer sources were fitted on consumed development data and assessed in the **sum** of all extra-image templates. The representative selected source adds **10^8 solar masses at a 30 pc scale**. Its conventional/R1 central ratios are **0.00264414 / 0.00264389**, both allowed by the HST application. All allowed nuclear masses have nearly tied development objectives; this is no nuclear-mass or black-hole measurement. All five candidates and the earlier R1 optimizer drift remain in the archive.

The actual source is

    rho(r) = sum_j 3 M_j a_j^2 / [4 pi (r^2+a_j^2)^(5/2)]
    M(<r) = sum_j M_j r^3 / (r^2+a_j^2)^(3/2).

The original nonlinear scalar consumes this **total enclosed mass once**. In the original reduced units, `x=r H0`, `p=pi_prime/(Mpl H0)`, `m=G M(<r) H0`, `bhat=b_T/Mpl`, `lambdahat=lambda/(Mpl H0^2)` and `c=1`, the executed applications are

    C x^2 p - 4 x p^2/lambdahat = -bhat m(x)
    2 x^2 dPhi/dx + bhat x^2 p = 2 m(x)
    J_eff/(Mpl H0^2) = -bhat rho/(2 Mpl^2 H0^2).

At **r=0.030 kpc**, the evaluated R1 section gives `rho/(Mpl^2 H0^2)=3.7344283856e9`, `J_eff/(Mpl H0^2)=-1.7120543202e8` and `p=-3.5987997269e-5`. Scalar left/source sides both equal `-3.49480848337e-17`; metric Gauss sides equal `7.62307121243e-16`. Maximum section residuals are **3.42e-16 scalar** and **2.17e-16 Gauss**, relative. The independent positive effective-mass deflection agrees within **3.38e-16 relative**. Numerical derivative checks and matching-tail diagnostics remain separately visible; none is a complete physical-admission test.

R1 activation at the **same conventional body, source, maps and clock** changes B/C/D-minus-A delays by **20.60284, 1.00422, 31.21607 s**. Separately nuisance-refitted differences are **-6.02919, -1.06520, -14.18106 s**. The latter are not the fixed-body physical signal. Exact projection endpoint selection prevents temporary continuation interpolation or optimizer drift from deciding the accepted development objective.

## Revised classical brightness application and freeze

The new model has one shared source OU, per-image extrinsic OUs, Gaussian random slopes, Brownian increments, jitter and proper image/instrument offsets. Its exact normalized state-space likelihood conditions these latent quantities jointly. The original modules remain unchanged. This is a classical empirical brightness law, not the BHSM quantum statistical-state kernel.

All 1002 old epochs plus 25 pre-2019 LCO R epochs establish development and instrument calibration. Six rolling-origin windows, including 2000-day horizons, use the declared higher 97.5th-percentile calibration rule. Short-window results are retained; the longer windows expose additional B undercoverage. Final primary A/B/C/D scales are **[1, 3.6132095, 1, 1]**; all-image scales are **[1.6101234, 1.9072587, 1.8775404, 1.9247902]**.

A finite 72-component nuisance ensemble includes development-only catalogue epoch-origin uncertainty. It contains correlated draws; basic between-chain Rhat reaches **5.19**. It is not a converged posterior, nor 72 independent samples. Macro-delay transport and the allowed C driver update weights using their declared development/driver likelihoods; no target is used to refit a curve, offset or parameter. The affine calibration includes its joint-density Jacobian.

Pre-score freeze commit: **b4cf9a04a7962c03441356ee8f7e18e3da157165**. Freeze SHA256: `44ba7a298f9c29b1fe3fe5d20b1432ceebcd58032727eca08c4780fe94d013a9`. All **79 committed dependency hashes** were checked before exposure. The persistent marker records opening only after the commit. The full [protocol](../preregistration/he0435_unused_validation_20261020.json), development amendments, candidate results and freeze remain immutable.

## One unused-catalogue application

The separate LCO campaign supplies 54 post-2018 R epochs: **162 A/B/D targets**, with 54 C observations as an explicitly allowed batch driver. The secondary forecast reserves all 216 values. The full 80-epoch V photometry remains unscored.

These numerical catalogue target values were unused by our fitting and selection. However, the earlier calibration rows and later targets share a retrospective difference-imaging reference mostly from 2022. Its exact constituent exposure list and full reference-noise covariance are unavailable. **Photon-disjoint calibration and full reduction independence are not established.** The frozen workflow identifier containing “independent exposure” must be read with that qualification; it does not establish a blinded independent image reduction. Qualitative published microlensing findings were also known before numerical target opening. These limits were disclosed before scoring and are strengthened in the final interpretation. Gaia geometry is assessed separately under an additional pre-target protocol.

| Frozen transfer | Joint log density | RMSE (mag) | 95% coverage | Mean 95% width (mag) | 95% interval score (mag) |
|---|---:|---:|---:|---:|---:|
| Empirical development delays | 292.466597 | 0.078029 | 98.765% | 0.656568 | 0.663991 |
| Conventional macro | 291.070684 | 0.077244 | 98.765% | 0.669945 | 0.677856 |
| Conditional R1 macro | 291.070650 | 0.077244 | 98.765% | 0.669945 | 0.677856 |
| Fixed-body R1 activation | 291.070835 | 0.077244 | 98.765% | 0.669941 | 0.677852 |
| Zero-delay control | 289.015811 | 0.079471 | 98.765% | 0.624534 | 0.633425 |

Conventional intervals contain **160/162** targets. Image A covers **52/54**, B **54/54**, D **54/54**; their mean widths are **0.36672, 1.30140, 0.34171 mag**. Nominal 68.27% intervals cover **94.44%** overall: uncertainty is conservative, particularly in B. The original undercoverage is not reproduced on this test, but precise nominal calibration and sharpness are not established. The changed data, driver and model prevent treating old/new coverage as a controlled single-variable comparison.

The primary R1-minus-conventional joint score is **-0.0000345961**; fixed-body activation changes it by **+0.0001507423**. These are conditional predictive densities, not Bayes factors, detections or evidence for a uniquely identified state. The complete matched conventional effective mass retains an identical transfer by construction. The conventional model exceeds the zero-delay control by only **2.05487** joint log units; the empirical-delay ensemble scores **1.39591** higher than the restricted macro family. The interval-score ordering differs, so coverage alone does not demonstrate superior predictions.

Without C conditioning, conventional/R1 coverage is **100%**, RMSE **0.21697 mag**, mean 95% width **1.25336 mag**, with a conventional signed bias **+0.17612 mag**. The all-image R1-minus-conventional score is only **+0.0000024967**. These are overlapping primary/secondary observations, not independent confirmations. Correlated epochs preclude an independent-binomial coverage significance claim.

## Separate previously unused Gaia geometry

The same physical predictions were frozen again for Gaia before any numerical quasar position row was queried, at commit **3032a85**. This protocol was added after the LCO outcome was known; no overarching prospective registration is claimed. A preceding metadata-only query found four entries. All 24 assignments were tested against the fixed conventional pattern, with no continuous frame/rotation/scale fit, and the same labels were retained for R1. The runner-up squared-pattern residual is 4.6197 arcsec squared versus 0.000063594 for the chosen match.

The first numerical ESA Gaia DR3 query and six-coordinate application were executed once. The four sources have RUWE 1.00–1.11, no duplicate flags and 0–1% multi-peak fractions. All remain in the result. Source RA uncertainties already include the declination cosine; the relative covariance includes common-A off-diagonal blocks. Unpublished between-source covariance is not manufactured.

**Both frozen point predictions fail Gaia's catalogue-error-only precision test.** Conventional relative-coordinate residuals in east/north mas are B=(+0.623,+1.701), C=(+1.915,-5.026), D=(-5.572,+0.580). Coordinate RMS is **3.25560 mas** (R1 **3.25557 mas**), with maximum absolute coordinate residual about 5.572 mas.

| Added relative-coordinate floor | Conventional formal Q | R1 formal Q | Conventional log density | R1 log density |
|---|---:|---:|---:|---:|
| 0 mas | 3465.9053 | 3465.6671 | -1727.1767 | -1727.0576 |
| 2 mas | 15.7966 | 15.7963 | -17.5910 | -17.5908 |
| 5 mas | 2.54112 | 2.54107 | -16.4441 | -16.4440 |

The 0/2/5 mas controls were specified before position exposure. Floors add independent variance per relative coordinate after propagating Gaia's shared-A covariance. Their better formal scores are sensitivity results, not a data-fitted cure, a derived physical posterior or proof of precision validity. No uniformly successful precise physical model is established. R1-minus-conventional log differences are only +0.1191, +0.0001531 and +0.0000243 under the three assumptions; a small difference between inadequate point predictions is not a state detection.

Gaia collected different photons in 2014–2017 from HST 2003/2012. The NIC2 astrometry used for the point fit was published in 2011. New HST central mosaics nevertheless have a Gaia-EDR3 field calibration; complete independence of every calibration cannot be certified without its match list. The **numerical Gaia quasar positions were unused by all model fits and selection**. The new test diagnoses the fixed geometry while exposing source/frame/photocentre and macro-uncertainty limitations. Gaia values are now consumed and cannot be reused after a refit as a second confirmation.

## Physical scope and preservation

The remote optical screen does not determine the same-domain scalar environmental Hessian. The total scalar equation contains its explicit cross-Hessian contraction; the projected two-dimensional optical tide does not close it. Prescribed density alone does not supply the material EOS/support. Global rolling-background Cauchy/matching conditions and the microscopic cap/source embedding/complement-state relation remain separately derived physical dependencies. The exact source trail is recorded in [postscore_source_audit.json](../artifacts/strong_lensing_revision_20261020/postscore_source_audit.json).

The empirical Gaussian brightness law does not execute a U/V preparation current, factorize a dressed quantum state, or replace the exclusive two-point, inclusive four-point or complete matched-occupation statistics. No later four-point calculation is imposed as a prerequisite for these classical source/optical applications.

The **2e-10 physical guard is unchanged**; `physical_admission=false`. The inherited main manuscript equation blocks, original observed modules/results, calibration and old source controls remain intact. **73 original/revision application tests, eight Gaia tests and 12 symbolic identities pass.** Reviews are internal agent/source checks, not external human peer review. The journal submission and cover letter remain unsent.
