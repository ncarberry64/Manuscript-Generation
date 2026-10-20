# Observed HE 0435-1223 monitoring data

This directory preserves original public COSMOGRAIL light-curve bytes for a retrospective chronological prediction test. HE 0435-1223 was selected for its long monitoring record and four separately resolved images before training or withheld-outcome scoring. No published delay estimate, mass-model posterior or cosmological posterior is used as a fit observation.

## Original releases and attribution

The raw files are in `candidates/HE0435-1223/`; `provenance.json` records URLs, original SHA-256 values, byte counts and the selection/split record. Preserve their byte representation; do not normalize line endings.

- `HE0435_Bonvin2016.rdb.txt`: 884 epochs from the official [COSMOGRAIL released-data page](https://www.epfl.ch/labs/lastro/scientific-activities/cosmograil/data/), associated with Bonvin et al. (2017), *H0LiCOW V. New COSMOGRAIL time delays of HE 0435-1223*, MNRAS 465, 4914-4930, [doi:10.1093/mnras/stw3006](https://doi.org/10.1093/mnras/stw3006).
- `HE0435_all_Millon2020.dat`: 1002 epochs from [CDS catalogue J/A+A/640/A105](https://cdsarc.cds.unistra.fr/ftp/J/A+A/640/A105/lcabcd/HE0435_all.dat), associated with Millon et al. (2020), *COSMOGRAIL XIX. Time delays in 18 strongly lensed quasars from 15 years of optical monitoring*, A&A 640, A105, [doi:10.1051/0004-6361/202037740](https://doi.org/10.1051/0004-6361/202037740).
- The two `ReadMe_*.txt` files are the corresponding original CDS catalogue metadata.

The extended file contains native time, relative magnitude and quoted magnitude uncertainty for each of images A, B, C and D, followed by telescope identification. Magnitudes have arbitrary relative zero points. Telescope labels are preserved in the raw data; the present empirical likelihood does not fit separate instrument calibration parameters, so remaining instrument/reduction errors are a limitation.

## Chronological split and common-reduction limitation

The independently chosen boundary is native numerical day **57447.03918**, the final epoch of the earlier 13-year release. Using the extended 2020 release consistently gives **902 training epochs** at or before the boundary and **100 later epochs**. The original release supplies the boundary, not photometry concatenated into the extended release.

The releases differ beyond constant magnitude zero points: 835 epochs are shared, 49 old epochs are absent from the extension, and 67 pre-boundary epochs are added. Some common historical photometry and uncertainties have changed. These differences were assessed using pre-boundary observations only. The fit therefore uses only the common 2020 reduction, split chronologically. This is retrospective withholding within a common published reduction; it is not a prospectively blinded observation or reduction pipeline.

The primary frozen-parameter test conditions on all historical training measurements plus later image B as the explicitly allowed observed driver, then predicts the 300 withheld A/C/D magnitudes. Because it conditions on the complete later B curve and uses smoothing, this is a retrospective cross-image conditional prediction, not a strictly causal online forecast. The secondary test conditions on training data only and predicts all 400 later magnitudes, including uncertainty in new intrinsic quasar variability. The frozen protocol is in `preregistration/he0435_heldout_20261020.json`.

## Time coordinate

The EPFL release page defines MHJD as HJD minus 2400000.5, while the 2020 CDS ReadMe labels its epoch column HJD minus 2400000. The numerical epochs agree between releases. This analysis retains native numerical days for splitting and differences and performs no absolute Julian-date or calendar conversion. A common half-day origin shift does not change the inferred delays.

## Measurement covariance and empirical likelihood

The release provides per-image marginal uncertainties but no full covariance between images from a common exposure or between reduction epochs. The measurement-noise likelihood approximates these errors as independent, with a training-fitted common jitter added in quadrature. This is an explicit approximation; a shared photometric systematic need not be absorbed correctly by independent jitter.

The total brightness covariance is correlated: one shared stationary OU source process, independent slower OU image processes, and proper Gaussian image-offset priors are analytically marginalized. It is therefore not a diagonal total likelihood. This is a classical empirical variability model, independent of the quantum state contractions in BHSM preparation equations. Predictive intervals condition on fitted delays/hyperparameters and do not include their full posterior uncertainty or all physical mass-model uncertainty.

## Astrometry and redshifts

`astrometry_redshifts.json` transcribes HST/NIC2 F160W relative point astrometry from Courbin et al. (2011), A&A 536, A53, Table 2, [doi:10.1051/0004-6361/201015709](https://doi.org/10.1051/0004-6361/201015709). Coordinates preserve the published delta-alpha/delta-delta signs in arcseconds relative to image A. Their quoted internal errors require an additional **2 mas systematic uncertainty**. A full astrometric covariance was not published; any diagonal relative-position covariance used in the macro fit must be described as an approximation. No HST pixel-level Einstein-ring reconstruction is supplied here.

The new application adopts lens redshift 0.4546 and source redshift 1.693, following the spectroscopy summarized explicitly in Bonvin et al. (2017). The updated source-redshift reference is Sluse et al. (2012), *Microlensing of the broad line region in 17 lensed quasars*, A&A 544, A62, [doi:10.1051/0004-6361/201219125](https://doi.org/10.1051/0004-6361/201219125). The older 1.689 value is recorded only to explain historical provenance; it is not silently mixed into the new distance calculation.

## Scientific reuse

The [CDS/VizieR rules of usage](https://cds.unistra.fr/vizier-org/licences_vizier.html) permit scientific use with citation of the original authors, publication and publisher; A&A catalogue data are free for scientific usage. No blanket Creative Commons licence is asserted for these observations. These are third-party data, not observations owned or newly acquired by this repository. Retain COSMOGRAIL, Bonvin et al., Millon et al., Courbin et al., their journals/publishers, and CDS/VizieR attribution.

Recommended service acknowledgment: This research has made use of the VizieR catalogue access tool, CDS, Strasbourg, France ([doi:10.26093/cds/vizier](https://doi.org/10.26093/cds/vizier)).
