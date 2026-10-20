# HE0435-1223: frozen classical light-curve likelihood

This note defines the empirical light-curve layer used to test the prepared lens-transfer predictions. It does not define a BHSM quantum state, a statistical state kernel, a Wick factorization, or a microscopic envelope preparation. Conventional and conditional R1 macro models use the same classical photometric likelihood, source variability, extrinsic variability, error treatment, and permitted conditioning data.

## Historical split and permitted conditioning

The input is the public COSMOGRAIL HE0435-1223 light-curve release associated with [Millon et al. (2020), COSMOGRAIL XIX](https://doi.org/10.1051/0004-6361/202037740): `data/strong_lensing_observed_20261020/candidates/HE0435-1223/HE0435_all_Millon2020.dat`. Its SHA256 is `880ed1cbc23d884ffec9a0d59411c946f3cb4256d962fc506464b080fb172ce3`. The columns supply native epoch and relative magnitude and quoted error for images A, B, C, and D.

The inclusive historical cut is native epoch **57447.03918**, the endpoint of the earlier thirteen-year release. Fit the **902** pre-cut epochs in the 2020 release, giving **3608** image measurements. The **100** later epochs are reserved. The older release is used to establish the endpoint, not to supply differently reduced training magnitudes. The 2020 release re-reduces historical observations, so this is a retrospective chronological validation within a common reduction. It is not a prospectively blind experiment, and it does not establish complete independence of photometric reduction choices between periods. All calculations use native day differences; the source records contain a half-day ambiguity in their absolute epoch-origin description that does not affect this split or relative delays.

The primary application conditions on the historical training data plus the **100 later image-B measurements**. It predicts the **300 later A, C, and D measurements** without fitting or shifting them. Image B is the prescribed driver before target scoring. Conditioning uses all later B measurements, so some permitted driver measurements occur after some target timestamps. This is a batch conditional prediction, not a causal real-time forecast. The secondary application predicts all **400 later measurements** using historical training only. Its source variability after the training period is inherently uncertain under the declared stochastic model.

The training loader examines later timestamps only to exclude their rows; it does not convert or summarize their magnitude/error tokens. Target magnitudes enter only the final scoring application. Known timestamps, image identities, and measurement uncertainties define the prediction design. The frozen-protocol runner checks file hashes and requires a committed prediction freeze before target scoring.

## Model and delay convention

For image index i and observer time t in days, the relative magnitude is

\[
m_i(t)=m_0+b_i+q(t-d_i)+e_i(t)+\epsilon_i(t).
\]

Here d_A=0 fixes the delay reference. Positive d_i means image i arrives later than A for the same intrinsic source feature. The processes and offsets have the proper Gaussian priors

\[
\begin{aligned}
\operatorname{Cov}[q(u),q(u')]&=\sigma_q^2e^{-|u-u'|/\tau_q},\\
\operatorname{Cov}[e_i(t),e_j(t')]&=\delta_{ij}\sigma_e^2e^{-|t-t'|/\tau_e},\\
\operatorname{Cov}[b_i,b_j]&=\delta_{ij}\sigma_b^2,\\
\operatorname{Var}[\epsilon_{ik}]&=s_{ik}^2+\sigma_j^2.
\end{aligned}
\]

All q, e_i, b_i, and measurement-noise terms are mutually independent apart from the explicitly shared q. The scale sigma_q is stationary RMS magnitude, not the asymptotic structure-function amplitude. A common extrinsic amplitude and timescale are used for the independent image processes. The quoted photometric errors s are treated as independent conditional measurement errors; the release does not supply the full inter-image/deblending noise covariance. Shared latent source covariance is retained. Neither a measured full noise covariance nor a physically complete microlensing model is claimed.

The proper broad offset prior has sigma_b=3 mag. The numerical mean level m_0 is the pooled median of training magnitudes, **-10.3459**; later observations do not set this level. The offsets are analytically marginalized with their uncertainty, rather than removed using target data. The hyperparameters are fitted empirical Bayes quantities, not externally known physical constants.

## Exact finite-state application

For each observation define source time u=t-d_i and sort by u. Stationarity and independence imply that e'_i(u)=e_i(u+d_i) has the same within-image covariance as e_i(t). The exact nine-component state is

\[
x=(q,e_A,e_B,e_C,e_D,b_A,b_B,b_C,b_D)^T.
\]

Set

\[
v=(\sigma_q^2,\sigma_e^2,\sigma_e^2,\sigma_e^2,\sigma_e^2,
\sigma_b^2,\sigma_b^2,\sigma_b^2,\sigma_b^2),
\quad
r=(\tau_q^{-1},\tau_e^{-1},\tau_e^{-1},\tau_e^{-1},\tau_e^{-1},0,0,0,0).
\]

The initial state has mean zero and covariance diag(v). Between successive sorted source times separated by nonnegative Delta u,

\[
A=\operatorname{diag}(e^{-r_a\Delta u}),\qquad
Q=\operatorname{diag}\{v_a[1-e^{-2r_a\Delta u}]\}.
\]

The observation row H_i has ones at q, e_i, and b_i and zeros elsewhere. With predicted state mean x^- and covariance P^-, the original scalar observation applications are

\[
\begin{aligned}
R_k&=m_{i_k}(t_k)-m_0-H_{i_k}x_k^-,\\
F_k&=H_{i_k}P_k^-H_{i_k}^T+s_k^2+\sigma_j^2,\\
K_k&=P_k^-H_{i_k}^T/F_k,\\
x_k&=x_k^-+K_kR_k,\\
P_k&=P_k^--P_k^-H_{i_k}^TH_{i_k}P_k^-/F_k.
\end{aligned}
\]

The exact normalized marginalized training log likelihood is

\[
\log p(D_{\rm train}\mid\theta)
=-\frac12\sum_{k\in D_{\rm train}}
\left[\log(2\pi)+\log F_k+R_k^2/F_k\right].
\]

There is no large-observation covariance approximation, spline binning, or diagonal replacement of the shared source covariance. The filter costs observation count times a fixed nine-state matrix update. Test timestamps can be inserted without assimilation; the OU transition semigroup leaves the training likelihood unchanged. This also handles target observations whose shifted source times interleave with training times.

For conditioning set C (training alone, or training plus permitted later B), run a filter assimilating C only and then a Rauch--Tung--Striebel smoother. The smoother gain and covariance application are

\[
J_k=P_{k|k}A_{k+1}^T(P_{k+1|k})^{-1},\qquad
P_{k|C}=P_{k|k}+J_k(P_{k+1|C}-P_{k+1|k})J_k^T.
\]

The corresponding conditional means and variances for each target are

\[
\hat m_i=m_0+H_i x_{k|C},\qquad
V_i=H_iP_{k|C}H_i^T+s_i^2+\sigma_j^2.
\]

Target magnitude values are not read by this prediction operation. The full correlated held-out score uses

\[
\log p(D_{\rm target}\mid C,\theta)
=\log p(C,D_{\rm target}\mid\theta)-\log p(C\mid\theta).
\]

The product of the individual marginal predictive densities is **not** substituted for this joint density. Marginal negative log densities, RMSE, mean residual, predictive sigma, and empirical 68.27%/95% interval coverage are reported separately, overall and per target image. The intervals and scores condition on fitted parameter points; they do not integrate uncertainty in the OU hyperparameters, macro parameters, mass profile, background, or envelope state.

## Training bounds, starts, and actual application

All bounds and initializers were chosen before target scoring. Delay bounds for B/C/D relative to A are [-30,30] days. The common initial hyperparameters are (sigma_q,tau_q,sigma_e,tau_e,sigma_j)=(0.15,200,0.05,1500,0.005), in magnitudes and days as appropriate. Logarithmic hyperparameter bounds are the logarithms of

| Parameter | Frozen interval |
|---|---:|
| sigma_q | 0.02 to 0.5 mag |
| tau_q | 10 to 1500 days |
| sigma_e | 0.005 to 0.5 mag |
| tau_e | 200 to 10000 days |
| sigma_j | 0.0001 to 0.03 mag |

Use four prescribed deterministic starts: all delays zero, followed by three non-reference delay vectors drawn uniformly from [-20,20] with NumPy generator seed **20261020**. The training-only L-BFGS-B objective is minus the normalized log likelihood, with explicit absolute finite-difference step 1e-5 and at most 200 iterations. Every start and outcome is retained in `lightcurve_training.json`; no target score selects a start.

The selected training application has

| Quantity | Fitted value |
|---|---:|
| d_B | +6.419840617785979 days |
| d_C | -0.9747304186323459 days |
| d_D | +12.487994831847848 days |
| sigma_q | 0.17385702623679347 mag |
| tau_q | 1500 days, upper bound |
| sigma_e | 0.1605004908096214 mag |
| tau_e | 10000 days, upper bound |
| sigma_j | 0.01003829730923635 mag |
| log p(training) | 8791.55342889895 |

These are estimates within this declared empirical variability model; they are not an independent precision-delay measurement with calibrated uncertainty. The four start log likelihoods are 8791.55342232394, 8727.397030216867, 8790.026037670279, and 8791.55342889895. Two starts reach nearly the same best objective, while two remain in poorer local modes. The source and extrinsic timescales reach their prescribed upper bounds. Those boundary results and alternative modes are retained, and the bounds are not enlarged after looking at held-out scores.

`training_diagnostics.json` records central finite-difference gradients at delay steps 0.0001, 0.001, and 0.01 days and fixed-other-parameter delay slices on a 101-point grid over plus/minus five days. All three slices have their sampled maximum at the frozen delay. The raw largest optimizer gradient includes an outward derivative at the active source-timescale bound; projected central-gradient maxima are approximately 0.142, 0.128, and 0.289 at the three diagnostic steps. The stationary OU covariance has derivative cusps when shifted sample times coincide, so small-step smooth-optimizer diagnostics need this qualification.

`training_curvature.json` records fixed-hyperparameter delay Hessians at 0.01, 0.02, 0.05, and 0.1 day steps. They are positive definite in this application, but their eigenvalues and inverse-curvature widths depend materially on the step. They are local numerical diagnostics, **not calibrated delay confidence limits or a measured Gaussian delay posterior**. The 0.1-day matrix is used only to initialize the physical macro fit. Its final objective uses the original training light-curve likelihood.

## Two-stage macro comparison and verification

First estimate the empirical OU parameters and delays using historical magnitudes only. Then hold the shared variability hyperparameters, offset prior and mean level fixed. For each physical macro model, fit its astrometry and original predicted delays against the same historical normalized light-curve likelihood and declared center prior. The local delay curvature initializes this optimization; it does not replace the final likelihood. Freeze each model's predicted delay vector and all input hashes before evaluating any reserved targets. This is a two-stage conditional empirical Bayes comparison, not a joint posterior marginalization over all physical and stochastic parameters.

The implementation is `code/strong_lensing_observed_lightcurve.py`, the repeatable training runner is `scripts/fit_he0435_lightcurves.py`, and the training-only diagnostic runner is `scripts/diagnose_he0435_training.py`. The eight focused tests check:

1. normalized Kalman likelihood against an independent dense Gaussian covariance calculation;
2. driver-conditioned RTS means and variances against exact dense Gaussian conditioning at interleaved shifted times;
3. invariance of predictions to arbitrary changes in hidden target magnitudes;
4. exact correlated conditional test density against the independent dense formula;
5. synthetic recovery and sign of a known +6.5-day delay, without target observations;
6. invariance under a common native time translation;
7. copying and immutability of observation inputs;
8. copying a JSON delay list so later external mutation cannot alter frozen parameters.

All eight pass. The last input-normalization fix was made after the numerical training application; `code_normalization_audit.json` retains the executed and final code hashes and verifies that removing only those declared tuple-copy lines restores the executed source exactly. Every executed initializer and optimizer decode already supplied an immutable tuple, so no numerical training result changed and no duplicate fit was run.

This note records the method and historical training application only. It contains no held-out outcome or result-driven choice.
