# Gaia DR3 astrometry: availability checked before position exposure

Status: metadata and a COUNT-only query have been executed. No numerical Gaia
position, uncertainty, proper motion, magnitude or source-quality target row
has been retrieved in this step. The already fitted conventional and R1 models
must be frozen before any such row is queried. Existing HST and light-curve
model inputs are not changed by this note.

The primary ESA Gaia Archive TAP endpoint is
`https://gea.esac.esa.int/tap-server/tap/sync`. The saved ADQL COUNT-only query
finds **four** `gaiadr3.gaia_source` records within 3.6 arcsec of the existing
HE0435 field coordinate. This is a feasibility count, not a verified image
assignment or use of numerical astrometry. The published primary
[Ducourant et al. 2018 abstract](https://doi.org/10.1051/0004-6361/201833480)
also lists HE0435 among the quadruples with all four images detected in Gaia
DR2. A later DR3 query is therefore feasible without using DR2 coordinates.

[ESA's DR3 release summary](https://www.cosmos.esa.int/web/gaia/dr3)
specifies Gaia observations from 25 July 2014 through 28 May 2017 and astrometric
reference epoch J2016.0 in ICRS. These are different photons from the HST
2003/2012 observations. The source astrometry used to fit the lens was published
in 2011, before Gaia observations. However, the newly retrieved 2026-reprocessed
HST mosaics have Gaia-EDR3-aligned WCS headers. Their central-image test used
relative pixel registration, not a numerical Gaia quasar-coordinate fit.
Without the complete HST alignment-match list, independence of every
astrometric calibration cannot be certified. This caveat must accompany any
claim of new physical validation.

## Relevant catalogue fields and covariance

`data/strong_lensing_revision_20261020/gaia/GaiaDR3_schema_only.json` saves a
TAP_SCHEMA query for column descriptions only. A post-freeze numerical query
can request `source_id`, `ra`, `dec`, `ref_epoch`, `ra_error`, `dec_error`, and
`ra_dec_corr`. Preserve quality fields such as `astrometric_params_solved`,
`astrometric_excess_noise`, `astrometric_excess_noise_sig`, `ruwe`,
`visibility_periods_used`, `ipd_frac_multi_peak`, `ipd_frac_odd_win`,
`duplicated_source` and `astrometric_n_good_obs_al`. Do not remove disagreeing
images after seeing their residuals. Null quality fields can occur for
position-only solutions and are not zero noise.

The primary [Gaia source data model](https://gea.esac.esa.int/archive/documentation/GDR3/Gaia_archive/chap_datamodel/sec_dm_main_source_catalogue/ssec_dm_gaia_source.html)
defines `ra_error` as the tangent-plane standard deviation
`sigma_alpha_star = sigma_alpha cos(dec)`, in mas, and `dec_error` in mas.
The cosine must not be applied to `ra_error` twice. The two position errors and
`ra_dec_corr` give each source's 2x2 tangent-plane covariance. Positions in
degrees must first be converted to a consistently oriented east/north tangent
plane; the frozen lens predictions use positive east and positive north.

For independent source catalogue errors, forming B/C/D relative to A yields
diagonal covariance blocks `C_i+C_A` and off-diagonal blocks `C_A`. Thus the
six relative coordinates cannot be treated as independent simply because A is
chosen as their origin. The primary
[ESA correlation documentation](https://gea.esac.esa.int/archive/documentation/GDR3/Data_processing/chap_cu3ast/sec_cu3ast_quality/ssec_cu3ast_quality_corr.html)
states that between-source correlations are not provided by the processing.
A Gaussian residual statistic using the published within-source errors must
therefore be labelled an approximation; no complete joint catalogue likelihood
is supplied by these rows.

## Freeze requirements for the subsequent application

The physical prediction is the already fitted B/C/D image offsets relative to
A, for both macro models, without a fitted transformation against Gaia targets.
Freeze the cone, assignment rule, covariance convention, numerical predictions
and their hashes before the position query. Any astrometric-floor control must
also be specified before opening rows, using the existing 2 mas HST training
systematic as a declared example rather than estimating a floor from Gaia
residuals. Such a control is not a newly derived macro posterior. Report the
actual residuals and source-quality flags even if the model fails. A match is
an additional point-geometry test, not an Einstein-ring reconstruction, a
unique envelope-state measurement or an empirical BHSM discrimination.

If this test is added after the LCO result is known, that chronology must remain
visible. Unused Gaia numerical target values can still be held out, but a
prospectively registered overall scientific program is not claimed.
