# HE 0435: unused geometry observations identified before new target exposure

Status: candidate inventory and proposed protocol, not a completed validation. On
2026-10-21 only MAST observation and product metadata were requested. No new
science pixels, target centroids, numerical target fluxes, image previews or
published numerical astrometry from these candidate reductions were opened.
The full query, replies, product filenames and hashes are saved under
`data/strong_lensing_geometry_20261021/metadata/`. The protocol and model must be
frozen separately before any numerical target application. This metadata
document does not claim prospective registration of the whole research program.

## What is already consumed

The 2011 Courbin et al. Table 2 training astrometry is explicitly HST/NIC2 F160W
simultaneous deconvolution, as recorded in the original observational receipt.
The later ACS F814W (program 9744) and WFC3/IR F160W (12889) science mosaics were
used in the central-light and outer-image registration application. Their
individual component exposures and alternative reductions are consequently not
fresh independent photons. Gaia DR3 quasar positions and all old COSMOGRAIL R
and released LCO R target photometry are consumed. Gaia DR2/EDR3, Gaia focused
product copies of the same astrometry, catalogue reformatting and a reprocessing
of these HST images cannot supply the requested new observational test.

The 80 LCO V epochs remain numerically unused, but share a reference reduction
with the consumed R campaign. They test a cross-band light model only after its
chromatic source/extrinsic response and zero-point relation are specified.
They do not independently test milliarcsecond lens geometry.

## Accessible, distinct HST photons

The primary [MAST API](https://mast.stsci.edu/api/v0/MastApiTutorial.html)
metadata query finds public ACS/WFC F555W program 9744 observation
`hst_9744_12_acs_wfc_f555w_j8oi12`, obsid **26356063**, MJD
**52869.88645176–52869.91671834** (2003-08-18). The current single-visit mosaic
reports **2030 s**. The often-cited 760 s is one subset, `j8oi12010`, not all of
the current visit: the other archived subsets are 700 s and 570 s. Choose the
complete declared input set before target use; do not select a subset by its
agreement with a prediction.

Five original calibrated charge-transfer-corrected exposures are available:

| Dataset | Product | Approximate full file size |
|---|---|---:|
| `j8oi12l9q` | `j8oi12l9q_flc.fits` | 168 MB |
| `j8oi12laq` | `j8oi12laq_flc.fits` | 168 MB |
| `j8oi12lcq` | `j8oi12lcq_flc.fits` | 168 MB |
| `j8oi12leq` | `j8oi12leq_flc.fits` | 168 MB |
| `j8oi12lgq` | `j8oi12lgq_flc.fits` | 168 MB |

Resolve each filename through
`https://mast.stsci.edu/api/v0.1/Download/file?uri=mast:HST/product/FILENAME`.
The product-metadata JSON records the actual sizes. The corresponding public
mosaic is `hst_9744_12_acs_wfc_f555w_j8oi12_drc.fits` (326727360 bytes). Single
FLC exposures avoid drizzle covariance and permit an explicit distortion/PSF
analysis, although DQ, CTE residuals and PSF/host errors still require modelling.
They are different exposures from consumed ACS F814W and NIC2 F160W, while
sharing HST calibration infrastructure and the 2003 ACS visit with F814W.

A second independent-epoch HST option is WFC3/UVIS F275W, program 11732 visit
04, obsid **26506487**, MJD **55072.07934855–55072.14592278** (2009-08-29),
2524 s. Original FLC files are `ib1304j8q_flc.fits`, `ib1304j9q_flc.fits`,
`ib1304jfq_flc.fits`, `ib1304jhq_flc.fits`, each 168681600 bytes. The public
visit mosaic is `hst_11732_04_wfc3_uvis_f275w_ib1304_drc.fits`. Other metadata
show F275W visits under 11732/12324/13113 through 2013; those targets remain
unopened. The primary observational methods identify these filters/programs
in [Zahedy et al. (2016)](https://academic.oup.com/mnras/article/458/3/2423/2589311).
UV point-source centroids need an explicitly tested PSF model; a wavelength
dependent host/microlensing centroid effect must not be silently declared absent.

## Strongest cross-instrument candidate: JWST NIRCam

Public JWST GTO 1198 observation `o007`, target `t005`, contains four untouched
NIRCam bands, each with an archive exposure time **1803.776 s**:

| Band | MAST obsid | MJD range | Stage-3 product |
|---|---:|---|---|
| F115W | 97278424 | 59864.18393802–59864.20841948 | `jw01198-o007_t005_nircam_clear-f115w_i2d.fits` |
| F356W | 97278410 | 59864.18393802–59864.20841948 | `jw01198-o007_t005_nircam_clear-f356w_i2d.fits` |
| F150W | 97278411 | 59864.21189873–59864.23575872 | `jw01198-o007_t005_nircam_clear-f150w_i2d.fits` |
| F277W | 97278409 | 59864.21189799–59864.23575797 | `jw01198-o007_t005_nircam_clear-f277w_i2d.fits` |

These are 2022-10-12 photons from a different telescope/detector and are the
strongest available independence candidate. Product metadata for F115W and
F150W list 32 original `_cal.fits` products each (four integrations/exposures
across eight short-wave detectors), plus the resampled products. This count
does not identify which detector contains the lens. Detector/WCS header-only
footprints can identify that without opening science values. Do not average
these 32 files as 32 lens exposures: most contain other detector fields.
MAST download uses `uri=mast:JWST/product/FILENAME` at the same endpoint.

Header-only follow-up identifies **NRCB3** as the detector containing all four
images in the F115W exposure set. HTTP Range responses were streamed in 2880-byte
blocks through the end of the PRIMARY and SCI headers, then closed before any
SCI data bytes. The target detector's four files are
`jw01198007001_02101_0000{1,2,3,4}_nrcb3_cal.fits`; each header request read
37440 bytes. The approximate already-known field centre projects to native
pixels (796.775,1292.985), (988.934,1104.088), (1180.889,914.948) and
(1372.553,725.622), respectively. These are **pointing/footprint metadata**,
not measured quasar positions. The other seven detector footprints miss this
field. The large roughly 270-pixel steps preclude claiming that the four files
sample only a single unchanging subpixel PSF location. Full header receipts,
including HTTP content ranges and exact bytes read, are saved as
`*_headers_only.json`. The science array targets remain unopened by this
metadata operation.

The short-wave bands give the more direct point-centroid test; a new lensed
host reconstruction in long-wave bands additionally requires a specified host
source and PSF. The [STScI NIRCam calibration status](https://jwst-docs.stsci.edu/jwst-calibration-status/nircam-calibration-status/nircam-imaging-calibration-status)
reports filter/detector distortion calibration residuals. Those are calibration
evidence to propagate under a declared model, not permission to choose a floor
from the lens residuals. The target pixels and numerical literature models of
these NIRCam observations have not been inspected in this metadata search.

## Proposed geometry-only confirmatory application

1. Use the now-consumed Gaia and old HST data only for diagnosis/development.
   Freeze both revised physical models, every predicted stationary image,
   total-mass scalar source, all optical response terms, source positions,
   nuisance priors and numerical hashes before opening the new targets.
   Freeze a primary candidate/band and any secondary bands; do not replace a
   failed primary with an easier one after exposure.
2. Freeze a pixel measurement model separately. Fit same-frame field PSFs,
   smooth host/galaxy nuisance and point fluxes with free point centroids, not
   lens-model-constrained centroids. Use native SCI/ERR/DQ (including masked
   pixels), full instrument distortion and pixel-integrated PSFs. Test the
   extractor on synthetic injections and field stars before the lens target.
   Reject rather than reweight after seeing an instrumental failure under
   a previously frozen quality rule. Record all exclusions and all bands.
3. To prevent a Gaia frame alignment from masquerading as fresh geometry,
   restore an available distortion-corrected **unaligned `IDC_*`** HST WCS, or
   reconstruct it using its reference distortion files. Do **not** use raw
   `OPUS` as though distortion corrected. Inspect and save `WCSNAME`, headerlet,
   IDCTAB/non-polynomial distortion, velocity-aberration scale and reference
   versions. [STScI's WCS conventions](https://hst-docs.stsci.edu/drizzpac/chapter-4-astrometric-information-in-the-header/4-5-absolute-astrometry)
   distinguish distorted OPUS, unaligned IDC and Gaia-aligned FIT solutions.
   At JWST, use detector distortion calibration and save WCS/alignment
   references; the complete independence of every shared sky calibration is
   not established merely by having different photons.
4. A frame-independent primary statistic can retain calibrated angular scale
   using the intrinsic coordinates of the four-image configuration. With
   `u=(B-A)/|B-A|`, `v=J u` in a fixed handed tangent plane, use
   `q=(|B-A|, (C-A).u, (C-A).v, (D-A).u, (D-A).v)`.
   This removes only global translation and rotation (five geometric
   dimensions remain). No target-derived scale, shear or affine correction is
   allowed. A separately declared dimensionless `q/|B-A|` shape test has only
   four dimensions and cannot validate the absolute lens scale. Preserve
   chirality and the predeclared image-assignment rule.
5. Propagate the full measured joint centroid covariance through the exact
   Jacobian or frozen Monte Carlo extractor. The common-A and common-AB
   reference introduces off-diagonal covariance. Report the actual residual
   vector, RMS, normalized likelihood and predictive score with this
   covariance. Model uncertainty must come from the frozen development fit;
   detector uncertainty must come from ERR/PSF/calibration evidence established
   independently of target agreement. No 2/5 mas floor, covariance inflation
   tuned to residuals, or wider-interval validation claim is permitted.
6. Report every model failure at the resulting precision. Agreement could
   establish only this positional geometry conditional on the measurement
   and chromatic centroid model; it cannot alone identify an envelope state
   or a BHSM-specific interface response. Target values are consumed after
   the first numerical application, including a negative result.

F555W is the simplest accessible continuation of the current HST pipeline.
NIRCam offers stronger instrumental/epoch independence but needs a verified
native-detector PSF/WCS extraction. Whichever is actually executable must be
chosen before target exposure. None of these metadata records supplies the
missing physical response relation, if the revised model still lacks one.
