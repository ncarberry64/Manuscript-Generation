# HE 0435 predictive uncertainty revision

This note defines the revised **classical empirical brightness likelihood**.
It is a development record, not an independent validation result and not a
microscopic BHSM state contraction. The original model, predictions, score and
85.3% primary-target coverage remain in their original experiment unchanged.

## Observational domains and development protocol

The entire 1002-epoch Millon2020 COSMOGRAIL table is now consumed development
data, including the 100 epochs scored in the first experiment. The additional
25 LCO R epochs strictly before native MJD 58484 establish cross-instrument
calibration. None of those 1027 epochs is an independent confirmation. The
post-2019 LCO R target magnitudes are not opened by the development runner.

The released COSMOGRAIL table includes a tenth telescope field. It is retained:
Euler, EulerCAM, Maidanak, Mercator and SMARTS are separate instrument levels,
as is LCO R. The original first-experiment loader used only the first nine
numeric fields. Instrument labels exist; a full measurement-noise covariance
does not accompany these marginal error columns.

`preregistration/he0435_revision_development_20261020.json` was written before
revised model fitting and calibrates no future target. Its SHA256 is
`6b9b668c4570a35cc4520e8d3227d3ed3707a3f99d1d901c755ae9a76eff5feb`.
The sealed validation domain and macro-model freeze are separately controlled
by the parent experiment.

## Original applications executed

For image i, instrument k and observer epoch t, the revised application is

    u = t + delta_k - d_i
    m_i,k(t) = m0 + b_i,k + q(u) + e_i(u)
               + a_i (u-t0)/365.25 + w_i(u) + epsilon_i(t).

The A image fixes d_A=0. The common source is stationary OU:

    Cov[q(u),q(v)] = sigma_q^2 exp(-|u-v|/tau_q).

Independent per-image stationary extrinsic OUs obey

    Cov[e_i(u),e_j(v)] = delta_ij sigma_e,i^2 exp(-|u-v|/tau_e).

Independent random slopes and Brownian increments contribute

    a_i ~ N(0, (0.03 mag/year)^2)
    Cov[w_i(u),w_j(v)]
       = delta_ij D_i^2 min(u-t0,v-t0)/365.25.

Here D_i is measured in mag/sqrt(year), and w_i(t0)=0. A proper broad static
instrument/image offset prior, b_i,k ~ N(0,(30 mag)^2), allows the arbitrary
COSMOGRAIL magnitude zero and calibrated LCO apparent-magnitude zero to differ.
It also propagates uncertainty in those offsets. That broad convention is not
a measurement of macro magnification or color. Noise variance is

    Var[epsilon_i(t)] = quoted_error_i(t)^2 + jitter_i^2.

The model does not invent off-diagonal photometric error terms. A common
difference-imaging reference can induce such terms; that reduction limitation
must remain explicit in interpreting the new observations.

All latent variables are jointly conditioned, rather than fitting a new
microlensing curve to withheld targets. For time increment h, the source and
extrinsic transitions are exp(-h/tau), with innovations sigma^2(1-exp(-2h/tau)).
The slopes and offsets have identity transitions. The Brownian innovations
have variance D_i^2 h/365.25. The observation row selects the common source,
image extrinsic OU, image random slope with coefficient (u-t0)/365.25, image
Brownian process and image/instrument offset. The normalized innovations
likelihood and RTS smoothing equations are applied to this state, including
nonstationary Brownian initialization at the first source epoch.

The independent dense covariance implementation tests all of these
applications, including cross-instrument offsets, Brownian and slope terms.
An optional small C accelerator evaluates the same innovations likelihood;
the Python filter and dense covariance remain independently tested references.
No generated DLL or shared library is tracked.

## Priors and uncertainty integration

Delay support is [-30,30] days relative to A. Independent log-uniform supports
are source amplitude [.02,.8] mag and time scale [30,30000] days; per-image
extrinsic amplitude [.002,.8] mag and shared extrinsic time scale [200,30000]
days; per-image jitter [.0003,.08] mag; per-image Brownian amplitude [.001,.2]
mag/sqrt(year). These are declared statistical prior supports, not physical
constraints or measured upper limits. The reference epoch is floor(earliest
training native epoch minus 31 days). The median training magnitude defines
m0, with the offsets subsequently marginalized.

Finite retained random-walk MCMC samples of delays and all fitted variability
scales define a predictive mixture. The between-component mean spread and
within-component latent/noise variances both contribute. Two deterministic
seeds and all diagnostics are retained. The finite, correlated, short-chain
approximation does not establish posterior convergence. Basic between-chain
Rhat reaches 5.19 for an extrinsic/Brownian decomposition coordinate and about
2.16 for source amplitude/time scale, while delay Rhat is approximately
1.001--1.096 and jitter Rhat approximately 1.0005--1.011. These diagnostics
exclude a claim of a converged joint posterior. The frozen finite mixture is
used as an explicitly limited nuisance ensemble and independently calibrated
predictive algorithm, rather than as established complete Bayesian uncertainty.
Long source/extrinsic OU time scales and amplitudes are weakly identifiable;
the short-increment diffusion combinations can be much better constrained.
No MAP long-time-scale value is promoted to a measured physical state.

The COSMOGRAIL native time-origin convention differs by a possible half day
between available descriptions; the new LCO table labels MJD without a
specified heliocentric correction. The old-release common shifts
delta_old in {-0.5,0,+0.5} day have equal prior probabilities, while LCO is the
reference convention. Their weights are updated using development data only.
For a retained zero-shift nuisance sample theta_k, the shifted component has
weight proportional to

    w_k exp[log p(dev | theta_k,delta_old)
            - log p(dev | theta_k,0)].

This is a finite importance approximation, not three independently fitted
confirmations. Its reported importance ESS measures discrete weight
concentration, not MCMC effective sample size. Unknown heliocentric versus
topocentric corrections at roughly 0.006 day remain below this sparse cadence
and are a timing limitation; no shift is tuned to withheld magnitudes.

## Rolling-origin development calibration

The three fixed training cuts are native days 56000, 56800 and 57600; each
predicts its subsequent 500-day window. Every fold re-estimates its own delay
and variability parameters using only the preceding observations. An
instrument must have at least 12 earlier epochs to supply fold targets, so
an unconstrained new zero point cannot masquerade as failed time extrapolation.

The primary revised mode conditions on the C image in the later window and
withholds A/B/D. This is a batch cross-image prediction: some allowed C points
are later than the target epochs. It is not a real-time forecast. Choosing C
uses already-consumed first-experiment information; the former primary mode
used B. Those are different protocols, and their test coverages cannot be
treated as a controlled one-variable comparison. A secondary mode withholds
all four later images and receives no driver.

For each mode and each target image, let z=(observed-mixture_mean)/mixture_sigma
over the pooled rolling-origin windows. The prespecified conservative scale
is max(1, higher_empirical_97.5th_percentile(|z|)/1.95996398454). The calibration
uses no new validation photometry. Because the errors are temporally
correlated, this rule carries no exchangeability coverage guarantee.

The initial short-window calibration has 222 residuals per eligible image. Primary unscaled
Gaussian-moment 95% coverages are A=100%, B=92.79%, D=100%; the initial primary
scales for A/B/C/D are [1,1.3483272298,1,1]. C is the allowed driver and has no
primary targets. Secondary all-image forecast scales are
[1.2674779627,1.4458321609,1.5790635871,1.4795949118]. These original development
results are retained in `development_calibration_short_windows.json`.

Before opening validation values, inspection of the known target timestamps
identified an additional issue: new observations extend about five years
beyond calibration, whereas the initial development windows were 500 days.
The explicit pre-target amendment
`preregistration/he0435_revision_long_horizon_20261020.json` therefore adds
2000-day windows following cuts 54800, 55600 and 56400. The model, priors,
eligibility and quantile rule remain unchanged. All six windows are pooled;
they overlap in time and do **not** supply independent repeated validations.

The resulting 901 development residuals per target image have unscaled
primary Gaussian-moment coverages A=100%, B=67.15%, D=98.45%. This reveals a
substantial long-horizon B problem that a short-window-only calibration would
miss. The final frozen primary A/B/C/D scales are
[1,3.6132095008,1,1]; secondary all-image scales are
[1.6101234032,1.9072587333,1.8775403905,1.9247902264]. These larger intervals are
a cost of the empirical correction. Their usefulness must be assessed with
the actual independent coverage, average interval width, 95% proper interval
score and joint log score; development calibration does not certify coverage.

The predictive calibration transforms the complete conditional mixture about
its common mixture mean mbar:

    Y_cal = mbar + S (Y_raw-mbar),

where S is the positive image-dependent diagonal scale. Every component mean
and covariance is transformed consistently. The joint calibrated log density
equals the raw density at mbar+S^-1(y-mbar), minus log det(S). Marginal mixture
quantiles are solved from the full mixture CDF, rather than assuming that a
mixture is Gaussian. The proper joint score therefore includes both temporal
correlation and the calibration Jacobian.

Allowed later C observations update latent states and component probabilities
by their conditional likelihood. They do not refit hyperparameters, offsets
using other images' future values, delays or calibration. Macro-model delays
are imposed with development-only nuisance-likelihood reweighting; conventional
and conditional BHSM alternatives share the same brightness model and
calibration. All parameters and scoring rules must be frozen before opening
new A/B/D targets.

## Interpretation

The added extrinsic components and parameter mixture address concrete reasons
for the original conditional intervals being too narrow: a common extrinsic
amplitude, stationary extrapolation and point-hyperparameter conditioning.
Brownian and random-slope uncertainty can grow beyond the final training epoch.
This is an empirical noise/variability model, not a derived BHSM interface law.
Its success is decided by the new scored observations, with interval width,
residual bias and proper scores reported alongside actual coverage. Any
remaining undercoverage is a model failure to retain, not a reason to retune
the same held-out sample.
