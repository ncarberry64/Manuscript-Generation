# Actual histories, coupled response convergence and covariance source trail

This continues the existing PR #23 worktree at `06cf7c1b85d6605bfee0cc0d07fb2af86a66cae9`. The delivered repair package, original failed applications and unrelated newer work are preserved. All lens observations used below are the already consumed four F115W exposures. No reserved band is opened. A smaller residual or a changed conditional noise model is not physical validation.

## Actual domain and native applications

`scripts/acquire_he0435_target_histories.py` and `scripts/recover_he0435_target_histories.py` recover A/B/C/D and a field star in each of the four actual exposures, on matching native 65 by 65 domains. The downlinked SCI values are **nine two-frame group averages**, with the first individual read retained separately as ZEROFRAME. The complete later individual-read sequence was not downlinked; its times, not invented constituent values, are reconstructed from the native readout pattern:

    t_jr = (5 j + r) TFRAME, j=0,...,8, r=1,2;
    TFRAME = 10.73677 s.

The rounded `TGROUP=53.684 s` is retained where the actual native estimator consumes it; it is not substituted for physical read times. The actual observation date is **12 October 2022**, established by the science FITS header. Full-detector SIRS processing precedes extraction. All 20 rate/ERR replay checks pass: worst relative differences are 1.8394e-7 and 1.7881e-7, respectively, against the unchanged 1e-4/1e-3 gates. Compact official OLS_C applications reproduce full-detector native slopes exactly. These new target/exposure applications do not rerun the old delivered field experiment as a milestone.

The native runtime is jwst 2.0.1 with `stcal 1.18.0+winportable1`. Its explicitly archived Windows portability patch changes platform/debug definitions and constant array dimensions, not slope/variance equations. It does not establish the unrecorded archived installation version. Environment and patch receipts are under `artifacts/strong_lensing_response_completion_20261021/runtime/`.

Raw reference DQ encodings have each reference's `DQ_DEF`. Separate mapped-DQ sidecars use the official mapping; science GROUPDQ/PIXELDQ already use the native global flags. The four original recovered history files remain unchanged.

## Original mean and actual convergence

The original optical application remains

    I(theta,lambda) = sum_i F_i(lambda) delta(theta-q_i) + sum_k a_k(lambda) B_k(theta),
    f_ep = integral T_e(lambda) P_ep(theta,lambda;eta) I(theta,lambda) dtheta dlambda.

`P` already includes the pixel integral. The same 318-coordinate positive response, 36 nonnegative coherent-sky coefficients, independent point measurements, one detector plane per exposure, and all 24 field-star nuisances are solved together. No physical lens positions enter calibration. The native-rate continuation minimizes the original residual and declared priors:

    r_ep = (mean_ep - SCI_ep)/ERR_ep,
    Q = r_target^T r_target + r_field^T r_field + ||a/10||^2 + ||eta/0.2||^2.

Convergence requires the maximum projected gradient of **Q/2**, in declared dimensionless parameter scales, to be at most 1e-3 **and** last relative objective change at most 1e-10. Bound signs are included, rather than interpreting an active-bound raw gradient as stationarity failure. `code/strong_lensing_response_completion.py` implements the application and records all field nuisances explicitly.

| Application | Full Q | Target original-ERR Q | Projected KKT maximum | Status |
|---|---:|---:|---:|---|
| All four targets plus fields | 44420.297987 | 33639.007883 | 4.896e-5 | Converged conditional mean |
| E1-E3 targets plus all fields | 34631.543686 | 23295.018884 | 3.535e-4 | Converged conditional mean |
| Held E4 independent point/flux/plane fit | - | 11018.981762 / 1141 nominal dof | 1.661e-4 | Converged, inadequate |

The held residual worsens from the old one-step 10906.483. Broad structured residuals remain; only 15-28% of each image's Q is within three pixels. Convergence therefore removes an optimization ambiguity and **does not repair the held mean**. Full centroid covariance, including train/held cross blocks, is preserved as unconstrained local curvature, not calibrated posterior coverage.

An actual matched first-two-GROUP application also converges: all-data target Q=7085.723398, field Q=1875.832461, total Q=9409.145267 and projected KKT=5.636e-6. Its E1-E3 training and E4 nuisance applications converge as well. The early E4 statistic is 2277.353794 under its explicitly frozen conditional count/CDS metric. This is a different temporal projection of the same photons, with larger native sampling uncertainty; it remains inadequate. Its raw full-native comparator is retained separately. The two likelihoods are never summed as independent observations.

## Group map and bounded mean applications

The pipeline linearizes an ADC **group average**, so the appropriate conditional map from linearized individual-read means is

    H_j(q) = L( mean_r L^-1(q_jr) ) - D_native,j,
    q_jr = beta + mu_e(t_jr)/gain + D_read,jr.

The electronic intercept is applied before inverse linearity and averaging. `code/strong_lensing_group_calibration.py` evaluates this map, its exact read Jacobian and Hessian on an actual anchored monotone polynomial branch. The Hessian contraction is a stated second-order noise correction, not an exact transformed probability law. In actual data-anchored diagnostics, nonlinear group averaging reaches 87.624 DN in a bright core but at most 0.004558 DN/group in the 3-6-pixel wings; it cannot explain the broad wing failure alone. This per-pixel map is the earlier restricted application; the source-derived IPC composition below replaces it in the final native fit.

The original leading electron mean is

    mu(t) = lambda t + kappa_e t^2/2 L(lambda^2).

It is not valid to extrapolate a saturated core as measured stored charge lambda*t. `code/strong_lensing_bounded_collection.py` implements two explicit conditional continuations. For a finite storage ceiling Q, define

    H_p(t) = integral_0^t min(lambda_p s, Q_p) ds,
    D_p(t) = sum_q (lambda_p+lambda_q) [H_q(t)-H_p(t)],
    mu_p(t) = min(lambda_p t + kappa_e D_p(t), Q_p).

The antisymmetric internal D conserves charge before clipping; lost overflow is individually reported. The electronic-ceiling alternative keeps H=lambda*t^2/2 and unbounded stored charge, clipping only the observed read. Both recover the original leading mean before ceilings, and provide analytic rate, ceiling and kappa derivatives on each active branch. These are first-order hypotheses, not identified physical full-well or stochastic overflow laws.

Actual raw cores plateau below 65534 DN, above their SAT calibration cutoff. Neither the digital rail nor that cutoff identifies a physical well. The executed native applications use explicit initial-intercept-translated cutoff trials and a 10% larger trial; missing positive context ceilings receive a recorded provisional harmonic continuation. Positive point-plus-coherent-sky light is assigned to photon arrivals; the fitted detector plane is separately treated as an electronic rate offset. This partition remains conditional. Group0 supplies a per-pixel electronic-intercept boundary datum; it is not independently predicted. Actual observed GROUPDQ/JUMP selection is fixed, so deposited cosmic-ray charge after jumps is not reconstructed by this application.

`scripts/apply_he0435_bounded_native_history.py` forwards converged early optical means through the actual inverse-L/read-average map and official adaptive OLS_C, retaining native ZEROFRAME and unchanged full-rate ERR. It evaluates 240 channel/ceiling/kappa applications. This is an actual full-history mean application, not an assumed transfer of the old field coefficient. The canonical saved residual arrays and independent source review are under `bounded_native/`.

## Stochastic relation and calibration source

The leading mean alone does not determine a counting law. One explicit **conditional collection-boundary** law fixes A_pq=kappa_e(lambda_p+lambda_q) for neighboring p,q and A_pp=-sum_q A_pq. Its local moments obey

    mu' = lambda + A mu,
    C' = A C + C A^T + diag(lambda + A mu),
    K(t,s) = exp(A(t-s)) C(s), t >= s.

`code/strong_lensing_stochastic_response.py` projects the full interpixel and unequal-time covariance through caller-supplied actual read maps. Its kappa=0 limit is Poisson min(t,s); source births outside selected pixels are retained on the declared finite context. It is a local affine/linear-noise closure, not the exact nonlinear global counting process or a BHSM state contraction. Nonnegative means and small-signal/charge domains are checked without covariance floors. A parameter-dependent Gaussian instrument likelihood requires both logdet(C) and covariance derivatives; its full score is implemented and tested. Frozen-C continuation is explicitly distinguished from that optimum.

For SIRS, `code/strong_lensing_reference_covariance.py` constructs the observed-mask conditional reference derivative J=D+UT, preserving odd/even top/bottom offsets, calibrated side kernels, rejected controls and cross-stamp columns. All **36** native group corrections reproduce the actual pre-to-post pixels bitwise. This is an actual operator, not an assumed independent reference-noise array.

The READNOISE0282 HISTORY states that its CDS marginal map was made **after** ncdhas reference, bias and IPC correction from dark programs 1453 and 4442. Consequently it constrains only

    diag( D_CDS P C_raw P^T D_CDS^T ),

not the full raw covariance. In particular,

    C_post = C_science + J C_reference J^T
             - C_science,reference J^T - J C_reference,science.

Assigning CDS^2/2 to independent raw ADC noise and adding J C_reference J^T would omit the cross terms and may double count noise. Actual source darks, rather than missing numerical arrays or lens residuals, were followed to measure this relation.

Two public full-frame NRCB3 RAPID exposures from program 1453 sample two epochs within the reference pedigree (31 August 2022 and 18 October 2023), with two integrations apiece. They do not bracket the full pedigree beginning 14 April 2022. Exact 18-read supports matching all nine SHALLOW2 groups were acquired, totaling 603979776 raw bytes. The same bias/SIRS correction is applied without IPC deconvolution. Integration differences remove fixed detector means under an explicit independence assumption. Per-amplifier spectral covariance uses declared frequency-bin pooling; earliest data train and the latest epoch tests temporal transfer. Measured science/reference cross identities close to 1.09e-10 DN^2. Near-neighbor covariance of 20-41 DN^2 and horizontal lag-64 covariance of 5-12 DN^2 remain after correction: diagonal-only covariance is demonstrably incomplete. Epoch variance ratios of 0.698-0.819 remain explicit instead of being tuned on the lens.

These measured covariance components do not identify an unrestricted per-pixel, all-epoch stochastic law. Their stationarity, finite sampling, native operator stage and independent-integration assumptions remain part of the model. Unfiltered later-dark particle/outlier scores are retained as well as the declared screened comparison; the screened epoch ratios are not universal noise validation. The collection-boundary stochastic completion, physical/electronic clipping alternative, source/background photon partition and any source-dependent electronic IPC response are separate defining relations; none follows solely from a small coefficient error or existence of an old retarded array.

## Derived IPC/linearity composition

The archived CRDS context gives `jwst_nircam_ipc_0022.fits`, a ground-derived **deconvolution** G, not the forward carrier-sharing kernel. Its actual stored sum is 0.9999999913852662; no renormalization is hidden. The positive Neumann inverse H has omitted infinite-grid kernel mass at most 2.07573e-15 at order eight and reproduces G H=I to 2.22e-16 on the tested interior. This calibration normalization is not substituted for, or admitted through, the physical 2e-10 guard.

The linearity0054 header records IPC, reference-pixel and bias corrections completed **before fitting its coefficients**, while the consumed science channel skips IPC. The source-derived nominal composition is therefore

    preIPC photon rate = gain * G (effective photon rate / gain),
    native group mean = L_native[ H mean_r L_cal^-1(beta_pre + mu_pre/gain + D_pre,read) ]
                        - D_native,group.

H acts before group averaging/native L. Applying H to a final linearized mean would generally be a different operation. The empirical effective ePSF already contains its low-charge mean IPC; the first line prevents a second optical convolution. The electronic intercept is conditioned by solving `mean_r L^-1(beta+signal_r) = G actual_ADC_group0`. Its implicit derivative is included in the coupled directional chain, together with the cap defining relation, fixed harmonic context, charge response, inverse L, H, native L and ZEROFRAME dark.

`code/strong_lensing_calibration_order.py` and `code/strong_lensing_native_calibration.py` implement the original applications and their analytic directions. Five focused map/chain tests verify the full interpixel Jacobian, implicit-intercept derivative and active-ceiling rate/kappa derivatives against numerical changes, rather than merely mirroring formula constants.

All twenty central preIPC source domains are positive. The finite field-template exterior has negative deconvolved values and therefore does not define a physical incident counting rate there. A declared fixed-mask harmonic continuation fixes that boundary initializer, leaving it free as a physical input. Its central low-charge perturbation is at most 4.81955e-12 e/s in the initial saved-pilot application; that number is not asserted for the subsequent fitted boundary. Dark0337 does not record its IPC stage: nominal reference values are zero over 19 complete inner stamps, with one nonzero pixel in E1D. Nonzero values also occur in some wider IPC-support halos. Assigning nonzero readwise dark to the preIPC linearized stage remains an explicit conditional convention, rather than an inferred calibration fact.

The actual kappa=0 native application completes all forty patch/ceiling trials. The native reported-rate amplitudes are converted to physical arrival rates by TGROUP/(5 TFRAME); the initial unconverted three-ppm clock comparison is retained separately. At the nominal cutoff, target original-ERR Q=47398.057117, held E4 Q=14058.572106 and field Q=15412.923119. The 10% larger ceiling gives held Q=14058.277430. The old pixelwise map's held Q was 14034.117525. Thus the noncommuting composition is a real consumed calibration relation but **does not by itself repair the held mean**; this source trail is tested, not asserted as a successful explanation.

The actual nine-group temporal contrast independently profiles static intercept and rate before measuring curvature. A shared quadratic collection basis, common pointing/diffusion derivatives, and their joint 13-rank model all leave most held curvature unexplained. The joint model worsens held raw-curvature RMS by 2.05%; its reduced conditional quadratic score is not accepted as a repair. In an affine read model, a static slope/PSF error alone cannot generate those profiled temporal curves. Static illumination propagated through nonlinear collection/IPC/linearity can generate curvature; these tests do not uniquely identify collection, signal-dependent IPC, coefficient calibration drift or time-dependent optics.

The former finite h=1/2 DN native slope Jacobian crosses adaptive estimator branches on some pixels. `code/strong_lensing_native_ols.py` now implements the source-derived continuous piecewise OLS_C equations, including model-dependent SNR powers, median order, multisegment inverse-variance derivatives and native timing branches. All 11552 finite slopes in forty actual compact comparison applications match official C bitwise after float32 conversion. Median ties and branch boundaries remain individually reported. The native value application is preserved; differentiating quantization is not pretended to be a smooth operation.

## Converged coupled native response

`code/strong_lensing_native_completion.py` now solves the complete conditional native mean jointly. The body retains the 318 response coordinates, 36 coherent sky coefficients, independent point/flux coordinates, exposure planes, all field nuisances and a nonnegative shared collection coefficient. The actual applications are the bounded electron history, implicitly solved group0 electronic intercept, calibrated inverse linearity, interpixel forward map, read averaging, native linearity/dark subtraction, native weighted OLS, and flat/PHOTMJSR conversion above. Their full derivative includes the ceiling relation and both harmonic context maps; it does not freeze the consumed electronic intercept.

The fixed incident-rate boundary trusts the inner 33-square target and 17-square FIELD source domains, with harmonic continuation over the full 65-square context. The initial attempt trusting the compact-template edge stalled and is preserved under `native_coupled/all_data/attempt_1_compact_edge_mask/`. The final prescription gives zero domain rejections. A native-only bound test uses dimensionless scales, preventing small positive kappa from being misclassified as zero. The original static results are unchanged.

| Application | Original-ERR target Q, continuous mean | Official native float32 Q | Kappa, electron^-1 | Projected KKT maximum |
|---|---:|---:|---:|---:|
| All targets/fields, rank 439/439 | 33330.090404 | 33330.089424 | 1.594530921e-7 | 1.32149e-4 |
| E1-E3 targets, all four fields, rank 424/424 | 23005.301878 | 23005.302019 | 1.654003367e-7 | 1.64766e-4 |
| Held E4 targets, trained shared body fixed | 10974.065189 | 10974.065420 | Fixed from training | 1.91219e-4 |

Both full fits pass the declared gradient and objective-change convergence gates, with five active sky bounds. All-data total Q=44040.612639 and training total Q=34291.526487 include field and prior terms. The held target nuisance fit is also converged. E4 FIELD remains in calibration; observed E4 target group0 conditions its electronic intercept. This is an internal conditional held-target test, not a historically unused observation or an unconditional prediction of every group.

The original-error held statistic improves only **0.41%** from the converged static 11018.981762. About 79% of held Q lies outside three pixels, and un-JUMPed pixels alone contribute Q=9867.836391. Missing deposited cosmic-ray steps cannot explain the whole failure. The formal conventional geometry statistic remains 114.49/15 and fails. The full 32-coordinate train/held covariance retains sky/response/kappa crossblocks (minimum eigenvalue 2.11590e-6 pixel^2); it remains unconstrained local curvature rather than calibrated coverage.

Actual target and FIELD full-chain derivative checks agree to relative L2 below 7e-10. Native source/replay, saved-map and independent result reviews are under `native_coupled/`. The paper's revised figure compares saved static and native residuals at the same original ERR; displayed color clipping does not clip their numerical scores.

### Source covariance applications and a direct noise-mean test

The exact source-derived native derivative replaces finite h=1/2 DN differences, which crossed adaptive estimator branches. All 11552 actual finite slopes in forty compact applications match official C bitwise after float32 conversion. The final cross-stamp/amplifier dark covariance projections are in `dark_calibration/analytic_native_projection/`. They replace the old read-variance component, without adding native ERR or a second read-noise copy. Under this source-fixed conditional covariance, the saved static held target gives **Q=3901.2783/1156 pixels**, raw RMS **1.525146 SCI units unchanged**. It still fails; changing the diagnostic covariance is not validation.

An actual E4A application closes the direct nonlocal noisy-group0 intercept path at fixed conditional charge history/caps. Write B=d(predicted native rate)/d(beta_pre). Then

    d beta/d ADCobserved0 = diag(1/mean_read invLprime_0) G,
    J_residual,0 = W_observed,0 diag(Lprime_observed,0) - B diag(1/mean_read invLprime_0) G,
    C_residual = [I,-K] C_same_realization(observed_rate,ADC0_context) [I,-K]^T.

The actual conversion factors and ZEROFRAME branches are retained in the implementation. The contraction has 289 outputs and 1225 G0 context pixels; its component identity closes to 1.11e-16 SCI^2 and its minimum eigenvalue is positive. READ+DARK-only Q=3422.54025 is a component diagnostic, not a full likelihood. Source photons and the translated-cap noise path remain separate. In particular, for a=G(SAT-raw0+postSIRS0), b=G(postSIRS0), the cutoff consumes

    dCap = cap_multiplier * gain * [diag Lprime(a) G(du-dr) - diag Lprime(b) G du],
    du=d(postSIRS0), dr=d(raw0), at fixed source.

Its unconditional contraction needs the joint raw/reference/post-SIRS covariance, not just a delivered read-noise marginal. The physical full-well boundary is additionally not identified by this calibration cutoff.

To test an actual correction to the mean, a protocol was frozen before extracting sixteen fixed perimeter controls outside E4A's 17-square ROI. Seven orthonormal temporal contrasts at each control remove constant and actual mean-read-time linear terms, leaving 112 controls. The earliest independent dark spectrum predicts noise in the target via

    m_noise = C_target,control C_control,control^-1 z_control.

Target ROI residuals do not enter this predictor and its covariance is not retuned on the lens. The READ+DARK Gaussian source law gives predicted-noise RMS **0.0105645 SCI**, maximum 0.0261238. Raw target RMS changes only **1.29927705 -> 1.29906960**; original-ERR **Q worsens 3385.7513 -> 3425.0588**. A separately evaluated positive-source Poisson control contamination variant gives RMS 1.29907141, Q=3424.48285. That mixed non-Gaussian variant is a covariance-based best **linear** predictor, not generally an exact conditional expectation. Neither result repairs the mean. The frozen protocol, actual vectors, predictions and independent review are preserved under `preregistration/he0435_E4A_temporal_control_20261021.json` and `dark_calibration/E4A_temporal_control_mean/`.

### Background counterfactual and remaining defining law

On the converged native body, a direct source-domain counterfactual transfers the entire fitted late electronic plane into photon arrivals before IPC and charge response, removes the late plane, and re-solves beta from actual G0. All 24 patch applications stay in the positive trusted domain. Baseline full65 group replay is exact. Held Q changes **10974.065420 -> 10972.215765**, only **0.017%**; near-wing curvature changes approximately 0.02-0.05%. The native adaptive estimator is reapplied, rather than holding its weights fixed. This bounded alternative is insufficient, so no fitted background fraction is claimed as a repair. `native_coupled/background_partition/` contains the source application and independent review.

A final full-support dark-stage sensitivity removes the nominal preIPC RAPID-dark mean while retaining actual native dark subtraction and the same G0 boundary. The entire nonlocal support is included. Held target continuous changes are at most 7.53e-14 SCI, FIELD 2.456e-9 SCI; official float32 slopes are identical and held Q is unchanged. This specific stage convention is negligible on this saved body, while the general dark-noise and illuminated-state laws remain conditional. The actual application and arrays are under `native_coupled/dark_stage_support/`.

The remaining relation is the illuminated detector's conditional joint law

    P(ADC_0,...,ADC_8, ZEROFRAME, raw/reference controls, overflow
      | incident photon state, accumulated carriers, physical well/interface state),

with its signal-dependent collection/IPC and physical boundary prescription. The exact mean applications and source-calibrated READ+DARK component above fix parts of this relation. The nominal ground IPC fixes a low-charge linear map; dark data constrain a near-zero-illumination law. Neither determines an arbitrary nonlinear illuminated response or stochastic overflow. The trial ceiling, harmonic boundary and background partition specify a conditional iterate, not a measured detector state.

This is a derivable nonuniqueness, not an inference from an absent array. On any neighbor edge v=e_p-e_q, equal-rate conservative exchanges +v and -v contribute zero drift but **2 h v v^T** to the covariance rate for any admissible h>=0. Mean drift and low-charge IPC alone therefore do not fix that state-dependent exchange rate. This construction proves what is free; it is not silently added to the physical model. The adopted local affine moment law fixes one conditional small-signal closure only. A complete illuminated likelihood would also need the translated-cap joint-noise contraction above and consistent source/background state. No unique missing physical mechanism or successful held-exposure repair is established by these calculations.

The completed finite applications show that reaching conditional stationarity, group averaging alone, the tested ground-IPC composition alone, jump deposits alone, this extreme background repartition and this measured correlated-noise mean do not suffice to repair the held exposure. They establish local stationarity, not global uniqueness in a nonlinear fit. They do not authorize an arbitrary PSF extension fitted to gravitational residuals. Reserved science observations remain closed until a source-justified mean passes the already-consumed held-exposure check without variance inflation.

## Unchanged physical interpretation

The conventional potential and BHSM/R1 response remain declared matched controls. Instrumental calibration residuals are not new gravitational parameters. This work does not compute microscopic U,V sections or the original current j_c. The fixed incoming preparation, exact pair statistics, exclusive chronological two-point, inclusive contour four-point and complete matched-occupation distinctions remain intact. No vacuum, Wick/Gaussian dressed-state factorization, replacement statistical kernel, repeated local/native physical accounting or physical admission is introduced. The guard remains **2e-10** and `physical_admission=false`.
