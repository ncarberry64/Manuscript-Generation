# Gaia DR3: unused relative-image astrometry rejects the restricted point model

**Both frozen restricted macro models fail the catalogue-error-only astrometric
test.** Conventional relative-coordinate RMS is **3.2556013 mas**, with a
maximum absolute coordinate residual of **5.5716340 mas** and a conditional
catalogue-covariance quadratic residual **3465.9053 for six coordinates**.
The R1 result is essentially the same: RMS **3.2555702 mas**, maximum coordinate
residual **5.5718281 mas**, quadratic residual **3465.6671**. A code/test audit
passing is not a pass of this physical-model adequacy test.

The test evaluates previously unused Gaia DR3 numerical positions and errors.
No macro parameter, frame orientation, angular scale, image-specific offset,
source quality cut or uncertainty floor was fitted after opening these rows.
The improved HST central-light application and conservative LCO brightness
coverage therefore do **not** validate a complete HE0435 galaxy/envelope map.

## Freeze, provenance and exposure chronology

The conventional and R1 macro parameters and image predictions are copied
byte-for-byte from their previously frozen source calculation at
`b4cf9a04a7962c03441356ee8f7e18e3da157165`. Protocol, predictions, pure helper,
query runner and synthetic tests were committed in
**`3032a852fdb497fcb7711991d236897637d55caf` before querying numerical Gaia rows**.
Earlier metadata inquiries returned only a source count and column descriptions.
Three explicitly preserved pre-exposure implementation protocol versions do not
contain numerical Gaia targets. This physical validation was added after the
LCO primary result was known; no prospective registration of the whole program
is claimed.

The persistent exposure marker was created before the network request on
20 October 2026, 22:34:24 UTC (17:34:24 America/Chicago). The primary
ESA Gaia Archive `gaiadr3.gaia_source` cone query returned four sources. The
original CSV is retained at
`data/strong_lensing_revision_20261020/gaia/GaiaDR3_target_rows.csv`, SHA256
`f16654b688a122ebf7733970e9f6dca39a2f18aa4b5a2c8b16c445689c4cec3a`.
The protocol and prediction SHA256s are respectively
`748b3b40369b220588de22c61d18448f8a94c99885a16f09d2c257ea18c1da6e`
and `70f61abfa8068d48cccb56f3b7990c3c203621d8116788132898e0cf3a68d0e6`.
The saved score and marker record the same hashes and pre-query commit. Neither
the target query nor the scoring application was replayed for this review.

The [primary Gaia DR3 release publication](https://doi.org/10.1051/0004-6361/202243940)
is Gaia Collaboration; Vallenari, A.; Brown, A.G.A.; Prusti, T.; de Bruijne,
J.H.J.; et al., 2023, *Astronomy & Astrophysics* **674**, A1, *Gaia Data Release 3:
Summary of the content and survey properties*. The
[ESA release documentation](https://www.cosmos.esa.int/web/gaia/dr3)
specifies a 2014–2017 observing span and J2016.0 reference epoch. Numerical Gaia
quasar positions were not fitted in the earlier 2011 NIC2 source astrometry.

The four labels were assigned by the predeclared exhaustive 24 permutations,
minimizing unweighted residuals to the conventional relative-position pattern.
The R1 comparison uses this same assignment. The best association has residual
sum **0.0000635936392 arcsec²**; the runner-up has **4.61970294 arcsec²**, a gap
of **4.61963934 arcsec²**. No continuous transformation was fitted. Source IDs
in A/B/C/D order are `3178020716638059264`, `3178020716640423680`,
`3178020716638059136` and `3178020716638059392`.

## Actual relative-position residuals

Positive coordinates mean east and north. Residuals are observed minus frozen
prediction in mas, relative to A; A is an observed common origin, not an
independently fitted translation.

| Image | Conventional east | Conventional north | Conditional R1 east | Conditional R1 north |
|---|---:|---:|---:|---:|
| B | +0.622956 | +1.701362 | +0.622921 | +1.701178 |
| C | +1.914697 | −5.026486 | +1.914766 | −5.026188 |
| D | −5.571634 | +0.579826 | −5.571828 | +0.579850 |

The observed B/C/D relative coordinates in arcseconds are
`(−1.4762800332,+0.5527895269)`,
`(−2.4664004733,−0.6030558599)` and
`(−0.9381491117,−1.6141093406)`.

## Covariance, scores and prespecified sensitivity controls

Each Gaia position uses its published `ra_error`, `dec_error` and
`ra_dec_corr`. `ra_error` is already the standard error in right ascension
times cos(declination); there is no second cosine factor. The small-field
propagation gives relative covariance blocks `C_i+C_A` on the diagonal and
`C_A` between B/C/D. The common A error is not silently dropped. Unknown
between-source covariance is not supplied by Gaia and is not fabricated here.

The primary statistic uses those six-coordinate catalogue errors only, with
zero fitted parameters and frozen point predictions. The normalized Gaussian
joint densities have units mas⁻⁶. They are conditional comparison densities,
not a fully marginalized physical posterior or calibrated global significance.
The fixed 2 mas and 5 mas controls add variance independently to the **relative
coordinates after** propagating the common-A covariance, exactly once. They
were specified before exposure and are not estimated from these residuals.

| Added relative-coordinate floor | Conventional Q | Conditional R1 Q | Conventional log density | Conditional R1 log density |
|---|---:|---:|---:|---:|
| 0 mas | 3465.905265 | 3465.667061 | −1727.176710 | −1727.057607 |
| 2 mas | 15.796645 | 15.796339 | −17.590989 | −17.590836 |
| 5 mas | 2.541121 | 2.541072 | −16.444055 | −16.444031 |

The R1-minus-conventional log-density differences are **0.1191024**,
**0.000153105** and **0.0000243144**, respectively. The catalogue-only score
cannot be presented as evidence for R1 when both models badly miss the
catalogue-error scale. Nor can the 5 mas control be selected after the fact as
a cure or new validation. The source/screen family remains too restricted for
the measured point geometry at the published precision, or needs independently
justified additional astrometric uncertainty; this test does not identify which
physical/nuisance extension accounts for the residuals.

## All source diagnostics retained

| Image | Gaia RA*/Dec errors (mas) | RUWE | Duplicated | IPD multiple peaks (%) | IPD odd/truncated windows (%) | Excess-noise significance |
|---|---|---:|---|---:|---:|---:|
| A | 0.09822 / 0.08139 | 1.0664 | false | 0 | 1 | 1.2553 |
| B | 0.12586 / 0.10194 | 1.0685 | false | 1 | 64 | 0.8465 |
| C | 0.13382 / 0.10549 | 1.1071 | false | 0 | 19 | 2.5714 |
| D | 0.18918 / 0.14853 | 1.0038 | false | 0 | 69 | 0 |

All four have six-parameter solutions (`astrometric_params_solved=95`) and
17–23 visibility periods. No source is removed. The ordinary RUWE values and
absence of duplicated-source flags do not justify waving away the mismatch.
Conversely, truncated-window fractions, nonzero excess-noise diagnostics and
possible band-dependent photocentres show why these rows do not establish an
ideal, systematics-free six-coordinate likelihood. Proper-motion entries are
retained in the raw rows but were not used to shift the target coordinates or
retune the stationary macro model after seeing the discrepancy.

Gaia observations are distinct photons from the HST 2003/2012 observations.
The new HST mosaics' global WCS was nevertheless aligned to a Gaia EDR3 field
frame. Their central-light test used relative pixel registration and did not
fit numerical Gaia quasar positions. Independence of every global calibration
is not certified. This caveat and the unprovided between-source covariance
limit formal significance claims; they do not erase the raw 3.26 mas residual.

## Review and preservation

`artifacts/strong_lensing_revision_20261020/gaia/post_exposure_review.json`
records read-only checks of the saved raw rows, score, source hashes, marker and
pre-query Git objects. All checks pass. Six pure-helper and two synthetic
runner tests passed before exposure, including common-A covariance, angle-unit
handling, all-permutation association, normalized Gaussian scoring, null-error
raw-residual reporting and a persistent repeat-query lock. No target-dependent
code correction, source fit, new Gaia query or scoring replay was performed in
this review. Original macro and GP inputs, HST constraints, earlier scored
results, exact pair-statistics distinctions, calibration and the unchanged
`2e-10` physical guard remain preserved. Physical admission remains false.
