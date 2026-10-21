# Executed HE 0435 response repair

The user authorized execution of the attached repair package. Its task file was treated as task material within that request, not as authority to replace newer work, merge main, submit the manuscript, open reserved observations, or admit a physical state. Work continues in the existing isolated PR #23 worktree. The delivery, original predictions, old failures and unrelated work are preserved.

## Verified integration

The ZIP SHA-256 is `28ee0e7d343a428fc81a0ee91ab907111d02170372caa2e9db358ef70a07dad3`. All 109 manifest entries are preserved byte-for-byte under `data/strong_lensing_response_repair_20261021/delivery/`. The supplied ten-file mask patch and its reviewed integration were published in commit `f77af4f6a192f519d411866569bdef5853d795cb` on [draft PR #23](https://github.com/ncarberry64/Manuscript-Generation/pull/23).

The delivered legacy-fitter digest did not match the actual immutable source at the revision it cited. The default guard failed before applying data. Its failed receipt is preserved. A separate Git-bound integration manifest verifies the actual original source object and the same ten JSON inputs; it does not modify the fitter or disable its default guard. See `HE0435_RATE_MASK_INTEGRATION_20261021.md`.

The current-runtime canonical replay restores **303** valid slopes, reproduces all 16 original centroids exactly before the mask change, and obtains reference/R1-inclusive raw geometry RMS **0.448457055/0.448562655 mas** after the change. Exposure scatter is **Q=253.027695/15**; all 16 absolute profiles still fail. The reference formal geometry statistic rises to **1033.114/5**, because restoring information tightens errors despite the smaller raw discrepancy. No ERR multiplier or physical refit is used. The small difference from the supplied Linux RMS is recorded as a runtime difference, not concealed as exact numerical identity.

The delivered native processing and earlier optical/response trials were read and audited, not executed again. Numerical `stcal 1.18.0` agreement is not proof of the unrecorded archived installation version. The supplied common-time normalized residual improves by 32.4393%, but its original p=0.00984653 remains below its 0.01 gate.

## New direct-group calculation

`code/strong_lensing_read_history.py` evaluates actual supplied two-read supports, native corrected groups and group selections on the delivered 66×66 field domain. No target science values enter this application.

For each pixel p and group j it evaluates

    E[G_DN,pj] = a_p,segment(j) + mean(t)_j f_DN,p
                + kappa_e mean(t²)_j/(2 g_p) L((g f_DN)²)_p.

`L` is symmetric nearest-neighbor exchange with a zero-flux external boundary; `g` is the recorded electron/DN gain. The implied individual-read mean transport conserves electrons, before different pixel histories are selected. It does not require a conserved sum of the final mixed-history slopes. Positivity is checked at individual read times. DNU/SAT exclude a group; jumps or gaps create separate intercept segments. The actual read times are `(5j+r)*10.73677 s`, r=1,2. Exact archived OLS_C slope weights are not invented: the executed likelihood applies the group means directly.

The training-only working covariance is

    C_p,jk = max(f_train,p,0)/g_p mean_rs min(t_jr,t_ks)
             + delta_jk READNOISE_CDS_DN,p²/4.

Segment intercepts are profiled by GLS. The first eight groups train the rate and history mean; the ninth is excluded from both that fit and the covariance pilot. The held prediction includes `C_ht C_tt^-1(y_t-mu_t)`, the conditional innovation covariance, the profiled-intercept uncertainty and the full joint fitted-parameter covariance. Training and held photons are never treated as independent. This is an internal temporal holdout of historically consumed data.

| Rate mean | History term | Training Q / dof | Held Q / 257 | Held raw RMS, DN |
|---|---|---:|---:|---:|
| 289 free pixel rates | zero | 1621.402 / 1648 | 272.942 | 10.20884 |
| 289 free pixel rates | electron exchange | 1566.524 / 1647 | 260.747 | 9.67245 |
| Fixed dated 4000 K STPSF, fitted position/flux/plane | zero | 11909.177 / 1931 | 1542.413 | 25.23791 |
| Same optical mean jointly fitted with exchange | electron exchange | 11569.456 / 1930 | 1573.909 | 25.58811 |

The free-rate fit gives `kappa_e=(2.3627144 ± 0.3351756)e-7 electron^-1`. Its 5.25% held RMS reduction is a temporal-curvature diagnostic, not PSF calibration. The seven-coordinate coupled optical fit gives `kappa_e=(2.3665262 ± 0.1281599)e-7 electron^-1`, but worsens the original held optical residual. Its held working p≈1.02e-187. Neither curvature error identifies a physical detector coefficient under an inadequate mean. Final-group full-domain charge residuals are −3.69e-11 and −2.05e-11 electrons for these two exchange fits; finite-ROI boundary exchange is separately retained.

An initial DN-conserving implementation and all its sources/tests/results are preserved in `read_history/superseded_dn_surrogate/`. The recorded gain happens to be uniformly 2.05 electron/DN, so the corrected predictions and scores coincide after coefficient rescaling. The general operator now conserves electrons with spatially varying gain, verified by a nonuniform-gain regression. No actual electron leak is inferred for the uniform-gain delivered example.

This is a declared leading mean hypothesis, not a unique derivation of NIRCam charge migration. It lacks a stochastic exchange law and its induced interpixel/time covariance. Exact native linearity Jacobians, correction-before/after-read-average effects, SIRS spatial covariance, dark/reference/gain uncertainty, fixed context-rate uncertainty and selection effects remain explicit. Those defining inputs are not supplied by a small formal coefficient error.

## New coherent sky and shared calibration calculation

`code/strong_lensing_joint_scene.py` implements a single observed-sky scene, transported through each exact archived GWCS. Four disjoint image neighborhoods each have nine positive Gaussian sky coefficients shared across exposures (width 0.10 arcsec, spacing 0.18 arcsec, declared coefficient prior scale 10 SCI). One detector background plane per exposure replaces independent local quadratic backgrounds. Point positions and fluxes are measurement nuisances; the application never loads or scores a physical lens prediction.

The response retains the delivered coarse **318-coordinate** positive field family, `spline(positive base)(u)*exp(B(u) eta)`, its gauge removal and the declared 0.2 coefficient prior. Its already pixel-integrated response is convolved once with sky quadrature at one quarter native pixel. Finite 25-pixel response support is zero outside; padding and canonical 17×17 likelihood support are explicit. No arbitrary PSF enlargement, additional pixel integration, PAM multiplication or chosen ERR scaling is introduced.

The evaluated objective is the sum of native `(SCI-mean)/ERR` squared residuals, the four field-star residual sums and the two declared prior penalties. Field centroid/flux/background directions are profiled. The baseline field-only cut covariance is `C_p=F_p^-1+S C_eta,field S^T`, `S=-F_p^-1 A^T B`. A distinct **actually executed joint mean step** updates sky, point and response parameters and re-evaluates the original nonlinear field-plus-target objective with line search. All Jacobians use one dimensionless-column SVD convention, relative cutoff 1e-10, complete nullspace reporting and separate data/prior-supported ranks.

On all four target exposures, that one coupled step changes target raw **Q=50121.532→34712.081** on 4620 native pixels. The original complete objective changes **57432.279→45065.902**. The full profiled joint problem has rank **414/414**. The largest unconstrained scaled residual gradient is **14477**, including active-bound directions, so it is not a projected KKT statistic. A separate no-fit saved-mean derivative evaluation finds free detector-plane scaled half-objective gradient **12254** and projected positive-sky violation **7605**. Stationarity therefore genuinely fails; **13** positive sky coefficients are at zero. Its covariance is explicitly **unconstrained local curvature**, not a bound-truncated or calibrated posterior. The full 32×32 centroid and 20×20 intrinsic-geometry covariances retain sky–response and inter-image/inter-exposure correlations. Maximum absolute off-diagonal centroid correlation is **0.333**.

A separate three-target-exposure training application leaves target exposure four out of the shared scene/response fit, while all four field stars remain calibration observations. The held target is fitted only for its independent point/flux/one-plane nuisances. Its original native **Q=10906.483/1141** nominal conditional dof fails absolute adequacy. The trained sky/response covariance and their cross block are propagated into the held centroid covariance; they are not used to alter ERR or certify a failed mean. This is a consumed-data development check with a provisional training iterate, not new physical validation or a conclusive rejection of all converged models in this family.

The first product-grid-interpolation implementation is preserved in `joint_scene/development_product_grid_interpolation/`; the final evaluator matches the exact delivered interpolate-base-then-exponentiate family. Raw-center, full finite-grid moment and baseline-relative moment outputs are saved separately. None by itself establishes physical centroid truth, and no unmatched convention is scored against the frozen gravitational prediction.

## What remains determined, and what remains open

The rate-selection defect is corrected and the new equations have actual numerical applications, residuals, covariance arrays and independent review. A temporal effect is a plausible diagnostic lead, but the optical-coupled held residual worsens. A coherent scene lowers training residuals and provides the required shared uncertainty structure, but the present nonconverged held-exposure application remains inadequate.

An adequate instrumental mean and stochastic/reference covariance remain needed before precision inference. The compact delivered group array covers one field star in one exposure, not the quasar histories/selections in all four exposures. A fully coupled target history/sky calculation consumes those actual histories and a transferable calibrated law; absence from a compact array does not show that the histories or a state relation do not exist. The failed field optical transfer prevents silently exporting its coefficient. Reserved observations remain closed until the physical predictions and complete measurement procedure are frozen.

BHSM remains the interpretive framework. Reference and R1-inclusive results are controls; the effective-mass representation is a retained-order degeneracy, not a contest between fundamental theories. Instrumental residuals are not new gravitational amplitudes. A classical optical fit neither computes the microscopic U,V sections nor evaluates their original source current. The physical cap/canonical source pairing and same-domain state/Cauchy prescription remain separate. Exclusive chronological two-point, inclusive contour four-point and complete matched occupations retain their distinction. No vacuum, Wick or dressed-state factorization is introduced. The original once-only local/native accounting, frozen preparation and predictions, **2e-10** physical guard and **physical_admission=false** are unchanged.

## Commands and verification actually executed

On the existing local `.venv-geometry` runtime:

```powershell
.venv-geometry/Scripts/python.exe -X utf8 scripts/audit_he0435_repair_delivery.py --delivery data/strong_lensing_response_repair_20261021/delivery --output artifacts/strong_lensing_response_repair_20261021/delivery_compact_audit.json
.venv-geometry/Scripts/python.exe -X utf8 scripts/replay_he0435_verified_rate_mask.py --root . --baseline-manifest artifacts/strong_lensing_response_repair_20261021/mask_verified_baseline_manifest.json --native-dir data/strong_lensing_geometry_20261021/native --output artifacts/strong_lensing_response_repair_20261021/mask_verified_native_replay.json
.venv-geometry/Scripts/python.exe -X utf8 scripts/run_he0435_read_history.py
.venv-geometry/Scripts/python.exe -X utf8 scripts/run_he0435_joint_scene.py --resume-development
.venv-geometry/Scripts/python.exe -X utf8 scripts/export_he0435_joint_scene_residuals.py
.venv-geometry/Scripts/python.exe -X utf8 scripts/check_he0435_joint_scene_stationarity.py
.venv-geometry/Scripts/python.exe -X utf8 scripts/r1_screened_review/check_remaining_field_algebra.py --output artifacts/strong_lensing_response_repair_20261021/remaining_field_algebra_check.json
```

The first joint-scene runtime failures (unsigned-DQ memmap access, an FFT derivative batch indexing error) and the superseded evaluator are retained, not counted as successful model applications. The final residual arrays are evaluations of the saved final coefficients; they are not refits. All 21 focused test files listed in `.github/workflows/strong-lensing-manuscript.yml` were actually run together: **146 passed in 26.79 s**. This includes the four delivered mask groups, six Git-bound guard groups, seven history groups and six joint-scene groups. The 12 exact continuation algebra identities also pass. Synthetic recovery checks implementation, not empirical calibration.

The first published delivery exposed a new whole-repository pytest import collision: the immutable delivered test copy and its canonical test counterpart have the same module name. A root `conftest.py` excludes only the new delivered and result-receipt directories from collection, retaining their bytes and all canonical/unrelated test paths. This is a package-integration fix, not a relaxation of science checks or a replay of earlier CI repairs.
Whole-repository local collection then reaches 281 tests with that collision resolved, but stops on three unrelated pandas-import errors in the geometry-only local environment. No whole-repository pass or installation of unrelated dependencies is claimed; the focused 146-test result above is the executed application-suite result.

`independent_response_review.json` independently checks the new source equations, held-group statistics, electron conservation, positivity, ranks, exact response family, active-bound limitations and saved residuals. `joint_scene/independent_covariance_review.json` independently checks the saved full covariance, cross correlations, moment transforms and analytic derivatives without refitting. Original failed adequacy flags are preserved. The manuscript build and rendering audit are recorded separately after regeneration.

## Publication and regenerated paper

The executed numerical applications and reviews are published at immutable revision `34694083c22d32261c6f0ac27211d8e3f984803d`; the focused Linux source/algebra and real-MDPI build jobs passed in [run 39634766236](https://github.com/ncarberry64/Manuscript-Generation/actions/runs/39634766236). The subsequent manuscript pins that archive and integrates the corrected current/superseded units, explicit native residual equations, actual held failures and precise remaining state/history/covariance dependencies.

The final local Galaxies draft is **44 pages**, with a 187-word abstract, zero LaTeX warnings, no overfull/underfull boxes, no unresolved references/citations and all 42 original main math environments unchanged. All pages were rendered and reviewed. The old diagnostic figure was replayed from saved JSON with `scripts/plot_he0435_jwst_geometry.py --figures-only`, stacking its panels and enlarging labels without changing the original scientific receipt. The new residual figure uses full journal width. `paper_review/publication_manifest.json` binds the final sources, PDF and rendered contact sheets. These are internal computational/layout reviews, not external human peer review. The draft is not merged or submitted to a journal, and physical admission remains false.
