# HE 0435−1223: executed HST central-image pixel test

The frozen soft-body conventional model predicts a positive-parity stationary
image at `(east,north)=(+0.0037359034,−0.0271602072)` arcsec relative to lens light,
with an absolute macro flux of **0.0631971756 times image A**. The conditional R1
fit predicts `(+0.0037117881,−0.0270829239)` arcsec and **0.0631314768 times A**.
These predictions were read unchanged from the already published
`artifacts/strong_lensing_observed_20261020/macro_forward_audit.json` before this
pixel application. No light-curve validation targets were read by this analysis.

**Result:** the original nonattenuated central-image prediction exceeds the
total central F814W light allowed by the tested empirical PSFs. A conditional
Gaussian 95% upper bound, allowing arbitrary nonnegative galaxy, host, quasar
wings and unresolved nuclear light, is numerically **central/A < 0.02935** after profiling
the tested PSF family, a ±50 mas square of light-centroid freedom and a 5%
normalization nuisance. Multiplying every marginal noise by five and lowering
the A normalization by 20% raises the numerical bound to **0.04124**, still below
both frozen macro predictions. Conservatively round these numerical envelopes
to **0.030** and **0.042**, respectively. **F160W is inconclusive:** its corresponding
bound is **0.55377**, because the infrared lens galaxy dominates the nucleus.

This is an imaging test of the macro model **conditional on its central/A macro
ratio being the observed ratio**. It is not an unconditional exclusion of a
mass profile: differential extinction, source variability across the central
image's delay, and microlensing have not been derived. Retaining the old macro
model requires at least **0.833 mag** additional relative dimming under the
nominal F814W test (about 0.46 mag under the deliberately inflated-noise
stress test), or a revised central mass distribution. A revised body must be
inserted into the original total-mass scalar source, rather than adding
independently screened scalar solutions. The central image is not a BHSM-only
prediction: both tested models share this defect.

## Actual archival pixels and calibration

The data were queried and downloaded directly from the primary
[MAST CAOM API](https://mast.stsci.edu/api/v0/invoke):

| Product | Observation | Exposure | Native output scale |
|---|---|---:|---:|
| `hst_9744_12_acs_wfc_f814w_j8oi12_drc.fits` | HE0435−1223; HST GO-9744, ACS/WFC, MJD 52869.95151–52869.98285 | 1445 s | 0.05 arcsec |
| `hst_12889_01_wfc3_ir_f160w_ibyg01_drz.fits` | HE0435−1223; HST GO-12889, WFC3/IR, MJD 56228.55191–56228.77802 | 9939.589728 s | 0.128250003 arcsec |

The ACS product combines ten chip images from five exposures; WFC3 combines 24
exposures. Both are native archive science mosaics with `SCI`, `WHT`, `CTX` and
input-exposure `HDRTAB`. The archive reprocessing dates are 2026-06-18 and
2026-07-01; these are not observation dates. `SCI` is electrons/s. No display
JPEG, synthetic lens image or four-point catalogue was substituted for pixels.

Full-product URLs, SHA256s, header calibration fields and exact native-pixel
cutout locations are in `data/strong_lensing_revision_20261020/hst/provenance.json`.
The full products remain local and are excluded from Git because the ACS mosaic
is 326.7 MB. The four exact SCI/WHT/CTX cutouts per band are saved as FITS with
shifted WCS: a lens cutout and three isolated field point-source PSF candidates.
There is no resampling, deconvolution or change of pixel values in these cutouts.
The HST proposal metadata and `PHOTFLAM`, `PHOTFNU` where applicable, `PHOTPLAM`,
and `PHOTZPT` are retained in their copied headers. Absolute magnitude
calibration cancels in the same-filter central/A comparison.

For a download, use the exact URL in the provenance JSON; for example
`https://mast.stsci.edu/api/v0.1/Download/file?uri=mast:HST/product/hst_9744_12_acs_wfc_f814w_j8oi12_drc.fits`.
Verify the complete SHA256 before applying the recorded zero-based cutout
centres and square sizes to each SCI/WHT/CTX plane. This recovers the native
pixel arrays, even if FITS CHECKSUM comments contain a different creation time.

## Original applications and nuisance treatment

The pixel equation actually evaluated is

\[
 d_i=f_A r_c P_i(x_c,y_c)+B_i+e_i,\qquad B_i\geq0.
\]

Here `P` is an empirical same-filter PSF, `f_A` is obtained by fitting image A
with that PSF and a local background plane, and `r_c` is the frozen macro ratio.
The four quasar point positions are fitted independently, then registered to
the published east/north astrometry by translation and a small rotation. The
same fitted rotation is applied to every central-image offset and template;
it is not dropped after registering the outer images. The
largest registration residuals are 1.35–1.48 mas in ACS for the three field PSFs
and 2.35–4.39 mas in WFC3. The centroid search of ±50 mas in each direction is
substantially wider. It is a nuisance envelope, not an astrometric measurement.

The PSF family contains three isolated field point-source candidates and an
empirical image-A template. All are background-subtracted, clipped to positive
intensity, truncated consistently and aperture normalized. This explicitly
checks stellar PSF colour and the quasar's own detector sampling. Interpolation
and subpixel registration are applied to the actual templates. For the field
templates, the fitted A flux ranges from 1069.5–1106.7 electrons/s (ACS) and
357.4–373.9 electrons/s (WFC3). Point-source fits have residual PSF mismatch;
their photon-only chi squares are **not** acceptable global image models. The
5% normalization envelope is a declared nuisance choice, checked against a
20% alternative. The finite empirical PSF family does not certify every
possible instrumental PSF.

The maximization over centroid uses an 11×11 grid followed by continuous Powell
refinement, with the same bounds for the noise stress test. Thus the stated
upper envelope is a numerical profile over this domain, not a mathematical
certificate of a global continuum maximum. The margin between 0.04124 and
0.06313 remains much larger than the changes among the tested PSFs and grid
refinements.

### Arbitrarily flexible stellar, host and nuclear light

For a fixed central flux, every nuisance pixel is profiled independently:

\[
 \widehat B_i(f)=\max(d_i-fP_i,0),\qquad
 Q(f)=\sum_i[\min(d_i-fP_i,0)/\sigma_i]^2.
\]

This allows much greater nuclear and galaxy freedom than a Sérsic model. A
positive unresolved nucleus and any positive host emission can absorb a faint
central source; they cannot cancel a source brighter than all the measured
central light. `Q` is recorded as a diagnostic, **not** converted to a
chi-square detection significance at the nonnegative boundary.

The confidence construction uses marginal Gaussian upper limits. For `N`
pixels with `P_i > 0.03 max(P)`, set

\[
 U=\min_i\frac{d_i+\Phi^{-1}(1-0.05/N)\sigma_i}{P_i}.
\]

If the template and marginal errors are specified, the union bound gives at
least 95% coverage for a nonnegative arbitrary `B`, without assuming independent
drizzled pixels. Maximizing `U` over the tested source-position and PSF family
retains a conservative nuisance envelope. The bound is expressly conditional
on that template/noise family and the stated normalization allowance.

### Weight maps and drizzle noise

The archive `WHT` values are exposure-like: `DnnnWTSC=1`, ACS median weight near
1445 s, WFC3 near 10 ks. They are **not inverse variances**. The marginal
variance applied here is

\[
 \sigma_i^2=s_{\rm blank}^2
 \frac{\operatorname{median}(W_{\rm blank})}{W_i}
 +\frac{\max(d_i-d_{\rm sky},0)}{W_i}.
\]

The blank level and robust scatter are estimated in the same lens cutout's
2.5–3.7 arcsec annulus, with adequate weight. Shot noise is added explicitly.
The ACS blank scatter is 0.01030 electrons/s and the WFC3 scatter is 0.07048
electrons/s. Crowding can increase these estimates. This is an approximate
marginal noise model for mosaics, not a reconstructed full drizzle covariance.
The confidence construction handles arbitrary correlations conditional on
correct marginal errors; a fivefold marginal-error stress test assesses their
scale sensitivity. This follows the primary
[STScI explanation of weight maps and correlated drizzle noise](https://hst-docs.stsci.edu/drizzpac/chapter-3-description-of-the-drizzle-algorithm/3-3-weight-maps-and-correlated-noise).

## Smooth-light and unresolved-nucleus comparison

Separate actual least-squares applications fit pixels within 0.65 arcsec to a
positive shared-ellipticity Gaussian-mixture galaxy, with free centroid,
ellipticity and position angle, plus a free local background plane. An extended
family allows an intrinsic 0.02 arcsec cusp and a separate nonnegative nuclear
point source. These fits use the second field PSF. The central template remains
at its frozen registered sky location while the galaxy centroid varies; the
positivity analysis above separately profiles its positional uncertainty.

| Band and stellar-light family | No central source: Q / pixels | Frozen central source: Q / pixels | Residual RMS, no / frozen (electrons/s) |
|---|---:|---:|---:|
| ACS F814W, smooth positive MGE | 1027.42 / 532 | 46943.64 / 532 | 0.03889 / 0.39258 |
| ACS F814W, compact cusp + unresolved nucleus | 617.58 / 532 | 46943.64 / 532 | 0.02007 / 0.39258 |
| WFC3 F160W, smooth positive MGE | 509.90 / 80 | 465.46 / 80 | 0.20689 / 0.19577 |
| WFC3 F160W, compact cusp + unresolved nucleus | 441.74 / 80 | 435.29 / 80 | 0.19039 / 0.18876 |

The smooth families are descriptive central-light fits, not a lensed host-source
reconstruction or complete accepted galaxy image model. Their residuals are
retained rather than hidden by arbitrary noise rescaling. In the infrared,
galaxy/nucleus freedom readily trades against central light, so the small
chi-square improvement does not identify a central quasar image. The ACS result
uses the separate positivity inequality and therefore does not rest on the
chosen smooth stellar profile.

## Injection and focused validation

`python scripts/test_he0435_central_image.py` applies all fits and writes
`artifacts/strong_lensing_revision_20261020/hst_central_image.json` plus a
two-band pixel/model/residual/prediction figure. Positive PSF injections of
1%, 6.3197% and 10% of A into the actual images shift the fixed-position
upper bound by exactly their injected flux (floating-point precision).
This recovers the **limit's response**, not a unique source amplitude in the
presence of nuclear-light degeneracy. A separate 3000-realization perfectly
correlated Gaussian boundary simulation achieves 99.87% coverage for the
nominal 95% construction; it validates the conservative correlation treatment,
not the astronomical noise model or the whole PSF uncertainty.

Nine focused tests verify point-flux/centroid recovery with a gradient, sky-sign
and rotation registration, arbitrary positive nuclear degeneracy, correlated
noise coverage, monotonicity under positive galaxy light and rejection of an
overbright source without a Sérsic prior, and the necessity of summing
overlapping central images before applying a constraint, plus empty-support and
invalid-noise guards and consistent rotation of central templates. All pass. Original HST pixels and
original macro predictions are preserved; this test is no physical admission
of BHSM and does not alter the `2e-10` guard.

## Revised multi-image systems

For subsequent development candidates, the same runner accepts
`--prediction-json INPUT --output OUTPUT`. The input is
`{"models":{"name":{"images":[{"position_light_centered_arcsec":[east,north],"absolute_macro_flux_relative_A":ratio},...]}}}`.
Every additional image must be supplied. The original pixel equation is then
evaluated with the **sum** `f_A Σ_j r_j P_i(x_j,y_j)`, before fitting the
nonnegative nuisance light or computing the bound. Output `upper_scale_95`
multiplies this whole predicted image system; a value below one is a conditional
rejection. Applying a single-image flux bound separately to overlapping images
would be incorrect. The module's `assess_central_image_system` exposes the same
array-level calculation.

Because the HST result can inform the revised central mass profile, these later
candidate applications are explicitly labelled **development HST constraints**,
not unused validation. Independent light-curve targets must remain sealed while
these physical and instrumental choices are revised.
