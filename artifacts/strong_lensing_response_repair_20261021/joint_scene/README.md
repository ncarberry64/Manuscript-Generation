# Consumed F115W joint-scene application

This directory records an **executed conditional optical application**, followed by one actual coupled Gauss–Newton continuation. It is development on already-consumed F115W data. It does not validate either lens model or admit a physical envelope preparation.

The mean has independent point positions and positive fluxes in each exposure, one shared observed-sky scene, one shared response coefficient vector, and a detector background plane per exposure. The scene comprises four finite sky patches with nine nonnegative Gaussian coefficients each, shared across exposures through the archived nonlinear GWCS. It does not uniquely separate foreground galaxy light from an extended source.

The supplied coarse positive field family is evaluated as

`P(u; eta) = spline(positive_base)(u) exp(B(u) eta) / Z(eta)`.

The positive-base pedestal, 318 supplied response modes, first-order flux/moment gauge, and coefficient prior SD 0.2 are retained. Finite-grid normalization is absorbed by free point flux. The sky basis is sampled every quarter native pixel and convolved once with the pixel-integrated response. The 17 × 17 native likelihood uses the corrected SAT policy and original ERR, with no error rescaling or extra pixel integration.

At fixed supplied `eta0`, the target positions, fluxes, sky coefficients, and backgrounds are fitted. Four actual field-star stamps supply response information with their independent point/background nuisance parameters profiled out. The subsequent joint step changes response, sky, and point parameters together, then re-evaluates the original nonlinear target and field means. Raw field likelihood, both declared priors, and raw target likelihood determine acceptance. No lens prediction enters the fit.

| Application | Target chi-square | Nominal conditional dof |
|---|---:|---:|
| All data, supplied response | 50121.5321 | 4524 |
| All data, one coupled step | 34712.0814 | 4524 |
| Exposures 1–3, supplied response | 40336.1522 | 3383 |
| Exposures 1–3, one coupled step | 24106.1452 | 3383 |
| Exposure 4, trained response and sky fixed | 10906.4833 | 1141 |

These conditional reference dof count the target mean Jacobian with response fixed; they are not a calibrated goodness-of-fit distribution for the joint constrained nonlinear estimator. Every application remains far from native-pixel adequacy. Lower raw chi-square does not establish a correct instrument response or lens geometry.

The all-data step has accepted factor 1. Its raw field + target + prior objective changes from 57432.2788 to 45065.9022. The response update has `norm(delta_eta / 0.2) = 10.8097`. This is one provisional step, with convergence and coverage unestablished. The reported unconstrained scaled gradient includes active sky directions and is not a projected KKT convergence test.

A separate no-fit replay of the original target residuals and sky prior gives scaled gradients of the half-objective: free point coordinates 368.511, log flux 646.504, detector plane 12254.320, and projected nonnegative-sky KKT violation 7605.238. These nonzero free/projected directions establish that the saved one-step mean has not reached target-nuisance stationarity. This check fixes response at saved `eta`; it does not infer a complete joint response gradient without field nuisance application.

Two covariance estimands remain separate:

- The baseline **field-only/cut** propagation adds shared response uncertainty to the conditional measurement covariance without allowing target pixels to condition the response mean or precision.
- The updated covariance is the **unconstrained local curvature** of the specified linearized joint field + target + prior problem at the actual provisional iterate. It is not an exact nonlinear, truncated, or active-bound posterior. All-data rank is 414/414 and training rank is 399/399 under the common dimensionless SVD cutoff `1e-10`; these ranks include declared priors. Thirteen all-data and fourteen training sky coefficients sit at zero. Their bound directions are recorded rather than assigned zero variance.

The full 32-coordinate covariance keeps correlations between images and exposures. Raw nuisance centers, centers plus full finite-grid response moment, and centers plus the delivered baseline-relative moment are saved separately. These are declared finite-grid conventions, not established physical centroid accuracy.

Exposure 4 is held from the **target fit only**. All four separate field-star exposures constrain calibration, and inherited target regions/initial guesses were already available. This is an internal development split, not a genuinely unused validation set.

The compact inputs do not cover the actual per-read target electron histories and pixel-specific slope weights required by a gain-aware history application. This is a statement about delivered input coverage, not proof that the defining physical state or an upstream history does not exist. `history_rate` requires actual matching arrays; the unused charge-exchange helper assumes electrons or common-gain units. No target history or transferred charge coefficient is fabricated.

Run from the worktree root using the pinned environment:

```powershell
.venv-geometry/Scripts/python.exe scripts/run_he0435_joint_scene.py --resume-development
.venv-geometry/Scripts/python.exe -m pytest tests/test_strong_lensing_joint_scene.py -q
.venv-geometry/Scripts/python.exe scripts/export_he0435_joint_scene_residuals.py
.venv-geometry/Scripts/python.exe scripts/check_he0435_joint_scene_stationarity.py
```

The application command explicitly re-executes this consumed-data development calculation; it was run once after the delivered-family point evaluator was corrected. Six focused synthetic tests pass. The earlier product-grid interpolation application is preserved only in `development_product_grid_interpolation/` with its superseded scope; it supplies no final numbers.

The exporter command only replays saved target means. `raw_residuals/` retains full 17 × 17 native SCI/ERR/DQ/mask and coordinate arrays, plus valid-row mean, raw residual, and standardized residual maps for the baseline, updated, training, and held applications. Excluded-pixel model/residual values are explicitly NaN. Every raw chi-square replays within `1.5e-11`, and the original result files remain byte-identical. The display panel clips colors at ±10; the saved arrays remain unclipped. Field residual maps cannot be reconstructed from the saved field centroids alone because the original field flux/background nuisance values were not retained; this no-refit exporter does not manufacture them.

`equations_and_scope.json` states the actual applications and remaining dependencies. `execution_summary.json` binds the source, supplied response, native-file receipts, PSFs, and GWCS. `execution_checks.json` records saved-array arithmetic checks. Native arrays, earlier lens predictions, and inherited physical preparation records are unchanged. Physical validation and admission remain false; the physical guard remains `2e-10`.
