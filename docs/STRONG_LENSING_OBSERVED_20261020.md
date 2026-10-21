# HE 0435-1223: frozen observed-light baseline and withheld comparison

This checkpoint adds an actual observed-lens application to the strong-lensing work. It trains a conventional finite-body optical model and the conditional fixed-state R1/BHSM response on historical HE 0435-1223 monitoring and published HST point astrometry, freezes their predictions, and evaluates reserved later observations once. It does **not** establish a BHSM preference, a unique envelope state or an adequate complete physical model of this lens.

The frozen protocol is [he0435_heldout_20261020.json](../preregistration/he0435_heldout_20261020.json). The numerical predictions were committed before target scoring at [`e6e62bd`](https://github.com/ncarberry64/Manuscript-Generation/commit/e6e62bddd793a6b5545040b2fc37b775876384da); the [freeze receipt](../artifacts/strong_lensing_observed_20261020/frozen_predictions.json) and [one withheld application](../artifacts/strong_lensing_observed_20261020/withheld_scores.json) identify their exact hashes. The earlier [acceptance checkpoint](STRONG_LENSING_ACCEPTANCE_20261020.md) and its 19-page PDF audit remain historical synthetic/control evidence. They do not describe the new observed-data result.

## Release and reproducibility record

The manuscript title is **Strong-Lensing Time Delays in the Berger–Hopf Framework: An R1 Calculation and an Observed-Light Test**. Immutable science checkpoint [`14ad777`](https://github.com/ncarberry64/Manuscript-Generation/tree/14ad777eaf8befda387b9c3cb4504cf499e76cfc) includes the observed code, original data and inherited source ancestors, frozen training/scoring records, and comparison figure. The protocol is separately fixed at [`e55725d`](https://github.com/ncarberry64/Manuscript-Generation/commit/e55725d19135c3a434b9e12621845eae122e1b8f); the pre-exposure prediction freeze remains [`e6e62bd`](https://github.com/ncarberry64/Manuscript-Generation/commit/e6e62bddd793a6b5545040b2fc37b775876384da).

The regenerated real-MDPI Galaxies build has **24 pages**, **427799 bytes** and no LaTeX warnings; PDF SHA-256 is `5a031e394a1b9626f9eec4e5ba0cec6c20f97cc5320419771c7f8136e78c43f6`. The [publication manifest](../artifacts/strong_lensing_observed_20261020/publication_manifest.json) records release verification and preservation checks; the [24-page contact sheet](../artifacts/strong_lensing_observed_20261020/pdf_contact_sheet.png) accompanies the rendered review. Build completion and rendered-review status are separate checks; the manifest identifies the final review status.

An explicitly labeled known-outcome reproduction is available through `python scripts/replay_he0435_scores.py`, which writes its default output under ignored `manuscript/build/`. That replay was executed once and reproduced the saved predictions/scores within its declared tolerances. It does not refit, replace the original exposure/score record, or count as a second independent test. The fitted module, parameters and original result remain unchanged. No additional numerical solve or test replay was required for this documentation/release update.

## Data and independent chronological boundary

HE 0435-1223 was chosen for its long public COSMOGRAIL monitoring and four resolved images before fitting or inspecting reserved magnitudes. The original [Bonvin et al. (2017)](https://doi.org/10.1093/mnras/stw3006) 13-year release contains 884 epochs and ends at native numerical day **57447.03918**. That independently selected endpoint partitions the common [Millon et al. (2020)](https://doi.org/10.1051/0004-6361/202037740) release into **902 training epochs** (3608 magnitudes) and **100 later epochs** (400 magnitudes). The later release was re-reduced, so old and new photometry are not concatenated. This is retrospective chronological validation inside one common reduction, not a prospectively blinded acquisition or reduction.

The [data record](../data/strong_lensing_observed_20261020/README.md) gives original-byte hashes, source URLs, licensing/attribution, coordinate signs and the half-day disagreement in published epoch-origin metadata. All splitting and delays use native numerical day differences; no absolute Julian-date conversion is made. Original catalogue files are unchanged. Per-image magnitude errors are available; full same-exposure/reduction covariance is not. The noise likelihood therefore approximates measurement errors as diagonal quoted variance plus fitted jitter. The complete brightness covariance remains correlated through the shared source, image processes and offsets.

HST/NIC2 point astrometry is from [Courbin et al. (2011)](https://doi.org/10.1051/0004-6361/201015709), Table 2, with its required additional 2 mas systematic uncertainty. Relative-position covariance is approximated as diagonal because a full covariance is not published. The mass-light center prior is separately declared. The application adopts lens/source redshifts **0.4546 / 1.693**; the latter follows Sluse et al. (2012), A&A 544, A62, [doi:10.1051/0004-6361/201219125](https://doi.org/10.1051/0004-6361/201219125), as explicitly adopted by Bonvin et al. (2017). No published delay posterior, cosmological posterior, HST pixel reconstruction or central-image detection limit is treated as a new observation here.

## Equations, training and freeze

The empirical classical brightness likelihood is

`m_i(t) = m0 + b_i + q(t-d_i) + e_i(t) + noise_i(t)`.

Here `q` is one shared stationary OU source process, `e_i` independent slower OU image processes and `b_i` proper Gaussian offsets. Quoted independent errors and training-fitted jitter enter the measurement noise. Exact Kalman likelihood evaluation and train-only RTS smoothing marginalize these latent classical processes. Positive delay means arrival later than image A. The mean, four deterministic starts, bounds, hyperparameters and delay fits consume training observations only. This likelihood is an empirical variability assumption; it is not a BHSM quantum statistical contraction or an authorized Wick replacement.

A positive spherical Plummer body supplies the conventional local source. A remote symmetric quadratic optical screen at fixed redshift 0.8 supplies asymmetric two-plane Jacobi maps. The R1 comparison recomputes the original spherical scalar source and coefficients on the same declared frozen R1 background at the actual lens redshift. Both use H0=67.4 and the same source/channel maps and observer clocks. The screen is a conditional total optical propagation input, not a derived scalar environmental state. This is not an arbitrary spherical scalar inserted into an elliptical local density.

Common variability hyperparameters are fitted first; each macro model then fits the historical astrometry and **original training light-curve likelihood**, keeping those hyperparameters fixed. Local delay curvature is an initializer, not an observed Gaussian delay posterior. The empirical training delays (B,C,D relative to A) are `(6.41984062,-0.97473042,12.48799483)` days; conventional macro predictions are `(6.02001440,0.67269613,13.61142341)` days, and conditional R1 predictions `(6.02265914,0.67277986,13.61649789)` days. The macro restriction loses 9.576/9.572 units of training log likelihood relative to the free-delay empirical transfer, so optimizer termination is not a claim of complete lens adequacy. Both source and extrinsic OU time scales reach their upper bounds; the retained [training diagnostics](../artifacts/strong_lensing_observed_20261020/training_diagnostics.json) and local-curvature step sensitivity remain visible. No global optimum or fully marginalized parameter confidence is claimed.

The freeze stores physical/nuisance predictions and source/data/protocol hashes. Before target parsing, scoring checks these hashes and exact committed freeze bytes, then exclusively creates a persistent `scoring_started.json` marker. Preparation and a second automatic scoring run refuse either that marker or the final score. A score interruption cannot silently reopen model tuning. No parameters are fitted to reserved target values.

## One withheld application

The primary test conditions on historical training data **plus all 100 later B measurements** as an explicitly allowed observed source driver, then predicts the **300 later A/C/D magnitudes**. Smoothing may use B observations later than a particular target date, so this is retrospective batch cross-image prediction, not a causal online forecast. The secondary test uses historical training only and predicts all **400 later magnitudes**, including uncertainty in new intrinsic quasar fluctuations.

The joint score is the exact normalized correlated conditional density `log p(conditioning,targets)-log p(conditioning)`, with all fitted parameters frozen. It is not the sum of marginal scores. Log densities depend on the chosen magnitude unit; comparisons below share that unit. No model evidence, Bayes factor or significance is inferred from a point-fit score difference.

| Test/model | Joint conditional log density | RMSE (mag) | Nominal 95% marginal coverage |
|---|---:|---:|---:|
| Primary: empirical training transfer | 795.198466 | 0.101520 | 85.00% |
| Primary: conventional macro | 791.623845 | 0.101333 | 85.33% |
| Primary: conditional R1 macro | 791.628828 | 0.101333 | 85.33% |
| Primary: frozen zero-delay control | 755.975529 | 0.101436 | 86.67% |
| Secondary: empirical training transfer | 1032.875197 | 0.196156 | 75.25% |
| Secondary: conventional macro | 1029.291488 | 0.194945 | 75.50% |
| Secondary: conditional R1 macro | 1029.296502 | 0.194945 | 75.50% |
| Secondary: frozen zero-delay control | 993.390943 | 0.194942 | 75.50% |

The R1-minus-conventional primary score difference is only **+0.00498269**; the secondary difference is **+0.00501436**. This is no persuasive separation or evidence for a BHSM state response. Both primary macro predictions cover only **256/300 (85.33%)** of targets within nominal 95% marginal intervals and **142/300 (47.33%)** within nominal 68.27% intervals. The shared empirical transfer also undercovers. The conditional OU prediction is therefore not demonstrated to be adequately calibrated. Correlation between targets means these fractions are diagnostics, not independent-binomial significance tests.

The conventional-model image breakdown locates the failure:

| Test/image | Count | RMSE (mag) | Mean residual observed minus predicted (mag) | Nominal 95% coverage |
|---|---:|---:|---:|---:|
| Primary A | 100 | 0.133518 | +0.112380 | 62% |
| Primary C | 100 | 0.074072 | +0.048092 | 99% |
| Primary D | 100 | 0.086553 | +0.063765 | 95% |
| Secondary A | 100 | 0.119641 | -0.098605 | 99% |
| Secondary B | 100 | 0.276789 | -0.246190 | 34% |
| Secondary C | 100 | 0.183883 | -0.162410 | 83% |
| Secondary D | 100 | 0.165151 | -0.135894 | 86% |

R1 has the same reported coverage in every row and nearly identical residuals; full precision for every model/image is retained in the score JSON. Primary calibration is especially poor in A; future-only prediction is especially poor in B. These data do not identify whether the shortfall comes from the OU/microlensing model, missing observational covariance, conditional parameter uncertainty, calibration or mass/environment restrictions. No one cause is selected by this experiment.

## Original field application, numerical residuals and additional image

The [macro forward audit](../artifacts/strong_lensing_observed_20261020/macro_forward_audit.json) turns on the original R1 response at **identical** fitted conventional body, source, screen and clock. Its fixed-body delay changes are **19.29034420, 0.93856058 and 29.03669445 seconds**, with screening radius 313.50843142 kpc. The separately nuisance-refitted R1-minus-conventional changes are **228.50609665, 7.23433396 and 438.43556314 seconds**. The latter include parameter reallocation and finite optimizer termination; they must not be represented as the fixed-body R1 signal. The fitted R1 body's screening radius is separately 313.49958785 kpc. The earlier synthetic radius 212.02922498 kpc and its two fixed-body delay offsets remain unchanged historical controls.

An actual original-section/source application at the fitted R1 A-ray impact `r=7.741611730637014 kpc` gives `rho/(Mpl^2 H0^2)=885519.9402922966`, `J_eff=-40596.79508804083` and `pi_prime/(Mpl H0)=-0.00018351972574832728`. The original integrated scalar left/source sides are `-2.368678682663907e-13` / `-2.3686786826639066e-13` (relative residual `2.13e-16`); the original metric Gauss left/source sides are both `5.1666940730059555e-12` at printed precision. These are achieved reduced-source applications, not inserted targets or a completed microscopic preparation.

At the final fitted R1 body the maximum integrated scalar-equation relative residual is `1.10e-15`; pointwise metric residuals are at most `3.80e-16`. Across requested four-image branches, original lens residuals are at most `3.29e-18 rad` and arrival stationarity `6.68e-15 Mpc`. Temporary continuation interpolation differs from the final original projection by at most `7.68e-13 arcsec` and `4.83e-7 s`. Independent deflection/potential derivative comparisons give relative errors `1.19e-12` / `3.57e-9`; these individual numerical diagnostics are retained without treating them as a complete physical error bound or global admission test. Direct positive effective-mass projection agrees with the retained R1 transfer on these channels; its shared withheld score is a structural equality, **not an independent empirical confirmation**, and its mass is counted once.

A separate dense 1060-seed image search inside an analytic outer radial bound finds **five stationary images** in both fitted models. The additional central image has macro flux about **6.32% of A** (conventional) or **6.31%** (R1), at relative-to-A positions approximately `(-1.16864,-0.59294)` / `(-1.16865,-0.59290)` arcsec. No central-image upper limit, galaxy-light subtraction or non-detection likelihood was used, so this is an important untested prediction/adequacy constraint. The radial bound is analytic; the finite multiseed search is not a proof of every interior root basin. The earlier four-branch fit scope remains explicitly distinguished from this subsequent forward audit.

## Scientific outcome and next independent test

The experiment builds an executable observed-light baseline and a reproducible test of the conditional transfer. Its reported outcome is **no conventional/BHSM discrimination and inadequate conditional predictive calibration**, with a potentially consequential extra image left unconstrained by the supplied observations. It does not map a unique physical envelope state or validate the adopted background cosmology. The original BHSM embedding, material support, global environmental state and full preparation kernel are not supplied by this empirical likelihood. The exclusive two-point, inclusive four-point and matched-occupation distinctions remain intact.

The physical-admission flag remains false and the original **2e-10 physical guard** and calibration are unchanged. Small application residuals do not convert this provisional conditional fit into a physically admitted state.

These targets are now consumed. Any improved variability, covariance, profile, screen or central-image treatment may be developed transparently, but these same observations cannot then be presented as a new untouched test. A subsequent confirmatory comparison needs a new declared protocol and genuinely unused observations or independently reserved channels. Preserve this first result and all of its failures.

The [artifact index](../artifacts/strong_lensing_observed_20261020/README.md) separates data, training, freeze, source/census audit and one score. The supporting source/likelihood and macro reviews were completed before target exposure. Fifteen focused new tests (eight likelihood and seven macro tests) were executed separately; no fitting or tests were repeated to write this checkpoint.
