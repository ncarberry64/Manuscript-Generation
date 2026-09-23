"""Write the final, result-driven report without selecting or optimizing a model."""
from pathlib import Path
import json,hashlib,platform
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
ROOT=Path(__file__).resolve().parent; O=ROOT/'interpretation'
def table(df):
 def fmt(x):
  if isinstance(x,(float,np.floating)): return 'unavailable' if not np.isfinite(x) else f'{x:.5g}'
  return str(x)
 return '| '+' | '.join(df.columns)+' |\n|'+'|'.join(['---']*len(df.columns))+'|\n'+'\n'.join('| '+' | '.join(fmt(v) for v in row)+' |' for row in df.itertuples(index=False,name=None))
def read(name): return json.loads((O/name).read_text())
def main():
 l=pd.read_csv(O/'LIKELIHOODS.csv'); ix=l.set_index('test'); sight=pd.read_csv(ROOT/'sightlines/SN_SIGHTLINES_COMPLETE.csv.gz'); lc=pd.read_csv(ROOT/'lightcurves/LIGHTCURVE_RECONSTRUCTION.csv'); sys=pd.read_csv(ROOT/'lightcurves/LIGHTCURVE_SYSTEMATIC_REALIZATIONS.csv'); jk=pd.read_csv(ROOT/'lightcurves/LIGHTCURVE_JACKKNIFE.csv'); repeats=pd.read_csv(ROOT/'lightcurves/REPEAT_SN_SURVEY_SPLITS.csv'); tomo=pd.read_csv(O/'TOMOGRAPHY_ALL_SUBSETS.csv'); fit=read('FITTED_AMPLITUDE_AND_RANK.json'); perms=pd.DataFrame(read('PERMUTATION_SUMMARY.json')); checks=json.loads((ROOT/'VALIDATION.json').read_text()); cells=pd.read_csv(ROOT/'sightlines/LOCAL_MAP_CELLS.csv')
 single=pd.read_csv(ROOT/'lightcurves/SINGLE_EPOCH_JACKKNIFE_SUMMARY.csv'); lc=lc.merge(single[['uid','single_epoch_jackknife_max']],on='uid',how='left',validate='one_to_one')
 coverage=lc.groupby(['dataset','status']).size().unstack(fill_value=0).reset_index(); med=lc.groupby('dataset')[['band_jackknife_max','epoch_jackknife_max','single_epoch_jackknife_max','systematic_common_rms']].median().reset_index()
 sysstatus=sys.groupby(['status','calibration_mapping_complete']).size().rename('rows').reset_index(); sysstatus.to_csv(ROOT/'lightcurves/SYSTEMATICS_COVERAGE_SUMMARY.csv',index=False)
 stability=lc.merge(sight[['uid','primary_selection']],on='uid').groupby(['dataset','IDSURVEY','status'])[['common_delta_mag','band_jackknife_max','epoch_jackknife_max','systematic_common_rms']].agg(['count','median',lambda x:x.quantile(.9)]); stability.to_csv(ROOT/'lightcurves/SURVEY_STABILITY_SUMMARY.csv')
 tags=['P_original_heldout_release_only','P_original_training_to_heldout','P_training_to_heldout_fixed_external_PV_sensitivity','DES_original_heldout_release_only','DES_Pantheon_transfer']
 names=['Pantheon: release only','Pantheon: source/environment','Pantheon: fixed external PV sensitivity','DES: release only','DES: Pantheon nuisance transfer']
 main=ix.loc[tags,['rows','chi2_null','chi2_frozen','logL_null','logL_frozen','delta_chi2','likelihood_ratio','mock_p_frozen_direction']].copy(); main.insert(0,'Test',names); main.columns=['Test','Rows','χ² null','χ² frozen','log L null','log L frozen','Δχ²','L frozen / L null','Gaussian mock p']; main.reset_index(drop=True,inplace=True)
 fig,ax=plt.subplots(1,2,figsize=(12,4.8),layout='constrained'); vals=ix.loc[tags,'delta_chi2'].to_numpy(); ax[0].barh(names,vals,color=['#536878','#197c70','#b18535','#536878','#197c70']); ax[0].axvline(0,color='black',lw=.8); ax[0].invert_yaxis(); ax[0].set_xlabel('Δχ² = frozen − null (negative favors frozen)'); ax[0].set_title('Fixed prediction: held-out comparisons')
 for label,color in [('P_crossfit','#197c70'),('DES_transfer','#b18535')]:
  tt=tomo[tomo.test==label]; ax[1].plot((tt.z_min+tt.z_max)/2,tt.p_max_over_shells,'o-',label=label,color=color)
 ax[1].axhline(.05,color='black',ls='--',lw=.8); ax[1].set(xlabel='Foreground shell midpoint (redshift)',ylabel='Mock p, maximum over shells',ylim=(0,1),title='No aggregate local-cell localization'); ax[1].legend(); fig.savefig(ROOT/'FROZEN_COMPARISONS_AND_TOMOGRAPHY.png',dpi=180); fig.savefig(ROOT/'FROZEN_COMPARISONS_AND_TOMOGRAPHY.pdf'); plt.close(fig)
 fig,ax=plt.subplots(1,2,figsize=(10.5,4),layout='constrained')
 for name,color in [('Pantheon','#197c70'),('DES','#b18535')]:
  q=sight[(sight.dataset==name)&sight.primary_selection]; ax[0].scatter(q.zHD,q.path_map_fraction,s=5,alpha=.3,label=name,color=color)
 ax[0].set(xscale='log',xlabel='Released zHD',ylabel='Fraction of path measured by external map',ylim=(0,1.03),title='Local environment coverage'); ax[0].legend()
 for tag,name in zip(['P_original_training_to_heldout','DES_Pantheon_transfer'],['Pantheon','DES']):
  v=np.load(O/(tag+'_PERMUTATIONS.npz'))['null_calibration_p']; ax[1].hist(v,bins=np.linspace(0,1,11),histtype='step',lw=1.8,label=name)
 ax[1].axhline(20,color='gray',ls=':'); ax[1].set(xlabel='Conditional permutation p under Gaussian null',ylabel='Count / 200 null mocks',title='Permutation calibration is imperfect'); ax[1].legend(); fig.savefig(ROOT/'COVERAGE_AND_PERMUTATION_CALIBRATION.png',dpi=180); plt.close(fig)
 pairs=[]
 for label in ['P_crossfit','DES_transfer']:
  pair=pd.read_csv(O/(label+'_MATCHED_PAIRS.csv'))
  for cat in ['same_z_different_direction','same_direction_different_z','shared_foreground_different_z']:
   key=label+'_'+cat; row=ix.loc[key] if key in ix.index else None
   pairs.append(dict(dataset=label,match=cat,pairs=int((pair.category==cat).sum()),delta_chi2=float(row.delta_chi2) if row is not None else np.nan,mock_p=float(row.mock_p_frozen_direction) if row is not None else np.nan))
 pairtable=pd.DataFrame(pairs); pairtable.to_csv(O/'MATCHED_TEST_COVERAGE.csv',index=False)
 sky=l[l.test.str.startswith('P_leave_sky')][['test','rows','delta_chi2','mock_p_frozen_direction']]; surveys=l[l.test.str.startswith('P_leave_survey')][['test','rows','delta_chi2','mock_p_frozen_direction']]
 density=ix.loc[['P_crossfit_density_low','P_crossfit_density_high'],['rows','delta_chi2','mock_p_frozen_direction']].reset_index(); primary=sight[sight.primary_selection].groupby('dataset').agg(rows=('uid','size'),unique_SN=('sn_block','nunique'),external_PV_rows=('independent_pv_available','sum'),median_map_path_fraction=('path_map_fraction','median')).reset_index()
 rankp=read('P_crossfit_W_RANK.json'); rankd=read('DES_transfer_W_RANK.json'); known=read('FIXED_PV_CORRECTION_PROJECTIONS.json'); extra=tomo[tomo.test.str.startswith('P_leave')].sort_values('p_max_over_shells').iloc[0]
 frozen=json.loads((ROOT/'FROZEN_VALIDATION_RESULT.json').read_text()); state=dict(parameters=frozen['parameters'],amplitude=frozen['training']['A_mag'],amplitude_sigma_from_original_training=frozen['training']['sigma_A_mag'],new_topographic_parameters_fitted=0,secondary_diagnostic='one scalar amplitude only; never used for predictive likelihoods',primary_selection=dict(Pantheon='0.01 <= zHD < 1; original training below .03; held-out .03..1',DES='IDSURVEY 10; .10 <= zHD < 1.20; frozen same-SN overlap exclusions'))
 (ROOT/'FROZEN_SCIENTIFIC_STATE.json').write_text(json.dumps(state,indent=2))
 report=f'''# BHSM supernova sightline tomography V1

## Scientific answer

**A robust residual following the frozen topographic prediction is not established.** With source/environment nuisance relations learned only from the original Pantheon training set, the held-out fixed-template improvement is small: Δχ² = {ix.loc['P_original_training_to_heldout','delta_chi2']:.3f}, likelihood ratio {ix.loc['P_original_training_to_heldout','likelihood_ratio']:.3f}, conditional one-sided Gaussian-mock p = {ix.loc['P_original_training_to_heldout','mock_p_frozen_direction']:.4f}. The overlap-clean DES transfer prefers the null: Δχ² = {ix.loc['DES_Pantheon_transfer','delta_chi2']:.3f}. A fixed independent PV mean produces a stronger Pantheon sensitivity result, but its new field-error covariance is unavailable. It is not a fully recalibrated discovery likelihood.

The aggregate local-cell correlation tests do not localize residual correlations to a radial shell. Pantheon has some directional and high-density-subset preference, but it is uneven across surveys and sky sectors. These outcomes do not determine whether the physical cause is optical history, velocity correction, calibration, or another environment-dependent effect. The complete source-reduced, externally constrained, jointly marginalized scientific question remains **unresolved**, because the necessary joint light-curve/calibration/BBC and external-field error propagation is not supplied by this local reconstruction.

All frozen scientific parameters are preserved. No axis, transfer law, cutoff, cosmology, predictive amplitude, or sample threshold was selected to improve agreement. The exploratory fitted amplitude below is never substituted into a predictive test. Earlier result packages pass their original hash manifests unchanged.

![Fixed comparisons and tomography](FROZEN_COMPARISONS_AND_TOMOGRAPHY.png)

## Frozen state and data ordering

R1: Ωm=0.31, Ωr=0.00009, Ωk=−0.018, λ=1, q=3, V0=2.3986073. Axis (RA, Dec)=(211.48°, −12.81°); A={state['amplitude']:.17g} mag, zref=.02, zc=.03. The local transfer is H_K(χ) min(χ,χc)/[H_K(χref)χref]. The previous global table and 1:5 conditional partition remain bundled unchanged; this analysis tests the local topographic prediction. A is conditional on its previous low-redshift training, not a new independent prior draw.

Pantheon: 1,701 original rows and original STAT+SYS covariance; decimal asymmetry is symmetrized at the documented 3×10⁻⁸ level. DES: 1,820 original rows; the complete released precision is inverted before any row selection. Duplicate Pantheon observations retain their covariance and row identities. DES uses actual SN header positions, and the previous overlap exclusion is retained. No covariance difference is interpreted as an independent PV covariance.

{table(primary)}

The original Pantheon training sample has 357 rows / 274 CIDs; its 1,208 held-out rows contain no training SN. DES contains 1,467 retained rows. All 3,521 raw HD rows have a sightline record even when excluded by the frozen science cuts. `UNIQUE_SN_INDEX.csv` indexes repeated rows; cross-release identities remain explicit through the existing overlap audit.

## Level 1: source reconstruction and stability

Individual epochs and bands were loaded from the pinned public SNANA releases. Nominal SALT2-B21 and SALT3-Dovekie models, released primary-standard/passband tables, and the documented SNANA MAG_OFFSET=0.27 were used. Fits hold redshift fixed and separate log flux normalization, x1, color, and peak epoch. Remaining band contrasts are projected orthogonally to those source columns; phase contrasts are then projected orthogonally to both. Per-epoch decomposed residuals and each four-parameter covariance/estimator are retained.

{table(coverage)}

There are {int((lc.status=='ok').sum())} converged reconstructions and {int((lc.status=='local_fit_not_converged').sum())} non-converged local fits; these numerical failures remain in the package and do not remove any HD row. Exact SNID matching resolves the nominal 2005df row; the two public photometry alternatives are also retained separately. The fixed diagnostic window is −15..45 rest days, with positive finite uncertainties, supported passbands, at least five epochs and two bands. This is a local reconstruction, not an exact SNANA/BBC rerun.

Band jackknifes, three phase-block jackknifes, and individual-MJD epoch deletions use the local linear response at the reconstructed solution. Individual epoch deletion groups simultaneous bands together: 132,494 epoch groups are tested, with rank failures retained. Their uncertainty is the covariance of the **difference of estimators on shared epochs**, not a sum of independent-fit variances. Direct retained-epoch refits verify the local deletion formula to 7.6×10⁻¹⁵ in parameter units. They are not full nonlinear refits. All nine Pantheon and ten DES released spectral/calibration realizations were evaluated as local responses. There are {len(sys):,} nominal-row systematic responses, of which {int((sys.calibration_mapping_complete==False).sum()):,} have incomplete calibration mappings; these are flagged, not silently treated as complete released realizations.

Median absolute maximum jackknife excursions and RMS across the supplied systematic realizations, in magnitudes:

{table(med)}

These excursions are descriptive stability diagnostics, not calibrated additional distance errors. `SURVEY_STABILITY_SUMMARY.csv` preserves survey-specific counts, medians and 90th percentiles. `REPEAT_SN_SURVEY_SPLITS.csv` contains {len(repeats)} cross-survey repeat-SN comparisons; released common-residual differences use Cii+Cjj−2Cij. Raw-refit comparisons have only a clearly labeled diagonal photometric error proxy, because cross-survey raw calibration covariance is not reconstructed.

Critical limits: epoch fitting uses reported photometric variances, not a complete propagated SALT model/intrinsic-scatter/calibration covariance. The nominal DES Dovekie kcor files referenced by its fitter are absent from the downloaded release; corresponding Pantheon passbands are explicitly labeled proxies. Swift UV placeholder passbands are excluded from diagnostics. Source-fit χ² and jackknife errors therefore cannot be interpreted as exact survey-pipeline goodness-of-fit probabilities. Large excursions and non-convergence are reported without cutting the cosmology sample.

**The distance-like common residual for the released full-C likelihood is published standardized MU minus frozen R1 MU.** The epoch common-amplitude departure is a separate column. Replacing released MU by a new raw-refit distance while retaining its old covariance would be invalid, so no such replacement is made. Orthogonal color/shape/phase departures may enter the fixed-feature nuisance test, but their error cross-covariance with HD distances is not supplied; resulting likelihoods are conditional diagnostics.

## Level 2: independent sightlines and environments

The independent external field is [Lilow, Ganeschaiah Veena & Nusser's 2MRS-NeuralNet release](https://github.com/rlilow/2MRS-NeuralNet), trained on galaxy-count mocks rather than these SN distance residuals. Its 128³ Galactic Cartesian grid, 3.125 h⁻¹ Mpc spacing, 3 h⁻¹ Mpc smoothing, CMB-frame velocities and 200 h⁻¹ Mpc valid sphere are kept as released. The simulation/fiducial-cosmology prior is an external assumption, not a retuning of R1. See the [methods paper](https://arxiv.org/abs/2404.02278).

All 3,521 sightlines retain coordinates, redshifts, MU common residual, raw covariance index/hash/semantics, released PV and dust/host information, local density and Green-kernel lensing proxy, external endpoint velocity/components/marginal errors when available, and source-stability columns. Endpoint PV exists for 865 raw rows. Missing distant foregrounds remain missing. Lensing columns with zero placeholders are not promoted to independent lensing measurements. The local lensing proxy is not full convergence, shear, or a measurement of every optical channel.

The [2MRS observed galaxy catalogue](https://vizier.cfa.harvard.edu/viz-bin/VizieR-3?-source=J%2FApJS%2F199%2F26%2Ftable3) supplies 44,599 galaxies, with 43,507 positive barycentric-redshift entries usable for descriptive cell assignment. Its galaxy counts and shell densities are retained separately from the NN's galaxy-conditioned **matter** density. Raw counts are magnitude-limited, mask-uncorrected and in barycentric redshift space; zero count does not establish a void. These descriptive catalogue columns were added without changing the prespecified nuisance design or any likelihood.

![Coverage and permutation calibration](COVERAGE_AND_PERMUTATION_CALIBRATION.png)

## Covariance-correct conditional held-out tests

The nuisance vector uses released x1/color, Milky Way E(B−V), host mass, independent shell-density averages, and identifiable epoch color/shape/phase departures with missingness flags. Means/scales and nuisance coefficients are learned on training rows only; numerical SVD rank tolerance is 10⁻¹⁰. Original redshift-bin monopoles/slopes are the only held-out fitted nuisance terms. The raw epoch common-amplitude departure is never a nuisance predictor. **No external PV normalization, bulk velocity, axis or direction is fitted.** Released-PV means define the primary comparison; the external mean replacement has coefficient exactly one in the separately labeled sensitivity.

For a training estimator H and held-out feature matrix Xv, T=Xv H. For each model independently, the test residual and template are yv−Tyt and fv−Tft. The prediction-error covariance is Cvv+T Ctt Tᵀ−Cvt Tᵀ−T Ctv. This includes the released training/test cross-covariance. Null and frozen tests use the same covariance and design; tabulated log L is the Gaussian prediction-error density evaluated at profiled bin nuisances, not Bayesian evidence. Absolute log L values should only be compared within a row. Cross-survey systematic covariance is unreleased; DES transfer explicitly assumes its cross block is zero, and no combined discovery significance is claimed.

Δχ²=χ²frozen−χ²null; negative favors the frozen prediction. Gaussian p is the one-sided fixed-template score tail under the conditional null, not the probability that a model is true.

{table(main)}

The primary Pantheon result is a weak preference, not a significant held-out detection. The fixed external-PV sensitivity improves agreement, but the missing spatial and source/PV cross-covariance prevents interpreting it as a new fully marginalized likelihood. The stronger all-Pantheon cross-fit score (Δχ²={ix.loc['P_sky_crossfit_joint','delta_chi2']:.3f}, p={ix.loc['P_sky_crossfit_joint','mock_p_frozen_direction']:.4f}) reuses the original amplitude-training population and is **descriptive**, not independent confirmation. Singular prediction-error directions in the aggregate are removed by the documented 10⁻¹⁰ covariance eigenvalue tolerance; positive modes are retained without covariance jitter.

## Matched sightlines, survey and sky checks

Pairs are selected geometrically, preferring the same survey, with no reused SN in a category. Conditions: |Δz|≤.005 and ≥60° separation; ≤5° separation and |Δz|≥.03; or normalized local-W overlap≥.5 and |Δz|≥.03. Pair operators carry the complete covariance and bin-nuisance response.

{table(pairtable)}

DES supplies **no** same-z/≥60° pairs under these fixed conditions. Its other two matched selections are identical; their apparently favorable score is one correlated diagnostic, not two replications. Matched tests using all Pantheon rows also reuse the original amplitude-training population.

Low/high density uses the fixed Pantheon median independent path-density predictor, not a residual-derived threshold:

{table(density)}

Preference is concentrated in the high-density subset; the low-density subset favors the null. This is inconsistent with treating every subset as equally confirming the frozen pattern. These subset p-values are unadjusted for the family of exploratory comparisons.

Four fixed RA sectors [0,90), [90,180), [180,270), [270,360) are held out in turn. Every training observation of a held-out SN is excluded:

{table(sky)}

Independent-survey subset checks, with the same duplicate-SN exclusion:

{table(surveys)}

These folds overlap in their training data; neither likelihood improvements nor p-values may be multiplied as independent evidence. The complete row memberships, estimators, nuisance coefficients, singular spectra and prediction-error covariances are saved.

## Does the standard PV correction remove the same optical history?

Applying the **unchanged saved original low-z estimator** to the released correction gives {known[0]['saved_estimator_projection_mag']:+.9f} mag, reproducing the prior result. Applying it to the independent galaxy-only correction gives {known[1]['saved_estimator_projection_mag']:+.9f} mag. Both quantities describe projections of externally specified correction vectors, not fitted velocity directions. They differ substantially. The earlier released correction moved A from −0.077169950 to −0.043216234 mag; the independent fixed mean would instead give approximately {state['amplitude']-known[0]['saved_estimator_projection_mag']+known[1]['saved_estimator_projection_mag']:.9f} mag on that saved estimator. The predictive A is nevertheless unchanged.

Thus standard PV processing overlaps the frozen optical estimator, and an independent PV prescription does not remove the same projected amount. This is sensitivity to the correction convention, **not proof that PV processing is subtracting actual BHSM optical history**. No unrestricted residual-fitted bulk flow is used, and no matrix subtraction fabricates independent PV covariance.

## Covariance/rank and fitted-only amplitude diagnostic

With the complete fixed nuisance feature space, the one-scalar fitted amplitude is **{fit['amplitude']:.6f} ± {fit['sigma']:.6f} mag**. Direction and redshift transfer remain fixed. This is an in-sample fitted diagnostic, not a prediction or a replacement of A.

The nuisance-projected template retains {100*fit['remaining_template_fraction']:.2f}% of its original bin-profiled Fisher information; its uncertainty inflates by {fit['sigma_inflation']:.3f}. Nuisance rank is {fit['nuisance_rank']}. The severe unrestricted-bulk-flow rank collapse from the old interpretation test is not imposed by this independently constrained construction. That numerical improvement does not by itself establish the signal's physical origin.

## Experimental radial-shell tomography

Wik stores path length through fixed 10 h⁻¹ Mpc Galactic Cartesian cells, integrated by 2 h⁻¹ Mpc midpoint steps, split by the unchanged shell edges [0,.01,.03,.05,.1,.15,.3,.5,1]. `W_FULL_GEOMETRY_PATH_LENGTH.npz` covers this original radial domain through z=1; any remaining geometry is explicitly flagged, not counted as observed. `W_LOCAL_MAP_PATH_LENGTH.npz` contains only sampled path within the external field. It has {len(cells):,} columns, including cells intersecting the map boundary. All foreground cell IDs, shell labels, path lengths, crossings, density proxies and catalogue-count assignments are saved.

Local W ranks at relative singular tolerance 10⁻⁷ are {rankp['rank_at_relative_singular_tolerance_1e_7']} / {rankp['rows']} Pantheon rows and {rankd['rank_at_relative_singular_tolerance_1e_7']} / {rankd['rows']} DES rows. Numerical rank near the Gram-eigenvalue floor is not treated as physical resolution. Many cells and sightlines remain strongly correlated; individual cell amplitudes are not fit to held-out data.

For each shell the row-normalized cell-sharing kernel is W Wᵀ. Its physical perturbation is transformed by the **same training/test residual operator** as the distances before whitening/profiling. Scores test an additional positive covariance component; they are not direct likelihood comparisons between arbitrary mean fields. The null mean and frozen mean are both reported. 2,000 Gaussian mocks calibrate each score and the maximum across four locally covered shells.

{table(tomo[tomo.test.isin(['P_crossfit','DES_transfer'])][['test','shell','score_null_mean','score_after_frozen_mean','p_max_over_shells']])}

Neither aggregate local-shell test is significant. The smallest leave-subset shell-adjusted tail is {extra.p_max_over_shells:.4f} in {extra.test} (shell {int(extra.shell)}); it is one of 22 overlapping leave-subset scans, is not adjusted across that family, and does not establish a localized structure. The tomography therefore **does not localize a common foreground cause**, while the cross-survey results also do not establish a coherent frozen pattern. Both non-identifications are retained.

## Null calibration and limits

Each likelihood comparison has 2,000 Gaussian draws with seed 20260920 and the correctly propagated full covariance. Fixed positions, redshifts, source features and external fields are retained. Additional 2,000 permutations exchange angular-template labels only within survey × original-z-bin × RA-sector strata, with duplicate-SN blocks moved together; each destination retains its exact frozen redshift law. These conditional permutations do not reconstruct new physical skies or preserve an unprovided external-field error covariance.

{table(perms[['test','movable_rows','total_rows','conditional_permutation_p','fraction_null_p_below_05']])}

The last column should be near .05 for a perfectly calibrated nominal .05 test; it is estimated using 200 full-covariance Gaussian mocks. Pantheon gives .09 and DES .035. **Permutation calibration is imperfect**, so raw permutation tails are auxiliary only. Gaussian tails are conditional on fixed measured nuisance features and released covariance; they still exclude feature-estimation, missing cross-release calibration, external spatial-field errors and the original discovery/selection history. No small reported subset tail is advertised as a calibrated global detection.

## Reproduction, validation and deliverables

Run `python run_all.py` from this directory after installing `requirements.txt`. Bundled inputs make numerical reproduction offline. `prepare_inputs.py` and `fetch_foreground.py` record acquisition; `package_outputs.py` refreshes the manifest/archive. Python used: {platform.python_version()}. Input Pantheon commit: c447f0fea703fcd0fff57de5000947b5ca81286b; DES commit: c9a4fcafc4cbd19bd750dee47fc76194a45c181f. Public source links: [Pantheon+ release](https://github.com/PantheonPlusSH0ES/DataRelease/tree/c447f0fea703fcd0fff57de5000947b5ca81286b), [DES release](https://github.com/des-science/DES-SN5YR/tree/c9a4fcafc4cbd19bd750dee47fc76194a45c181f).

Validation passes: {checks['input_hashes_verified']} input hashes; original 20-file and 40-file result manifests unchanged; original fixed-template likelihoods reproduced; {checks['folds_without_same_SN_training_test_overlap']} train/test designs without same-SN overlap; four-parameter epoch injection recovery; orthogonal epoch power accounting; distance/covariance row-permutation invariance; likelihood-ratio identity; projected covariance idempotence; mock-score moments; path-length conservation; and no free PV/common-amplitude nuisance column.

Key files:

- `sightlines/SN_SIGHTLINES_COMPLETE.csv.gz`: every HD row, source diagnostics, environments and covariance provenance.
- `lightcurves/`: all retained epochs, fitted parameters/local covariance, orthogonal components, jackknifes, spectral/calibration responses and survey comparisons.
- `sightlines/W_*`, `FOREGROUND_CELLS_WITH_GALAXIES.csv.gz`, `GALAXY_TO_CELL_ASSIGNMENTS.csv.gz`: experimental foreground geometry and observed-count provenance.
- `interpretation/LIKELIHOODS.csv`, `FOLD_MEMBERSHIP.csv.gz`, `*_OPERATORS.npz`, `*_DESIGN.json`: all predictive tests and reproducible GLS operators.
- `interpretation/GAUSSIAN_NULL_MOCKS.npz`, `*_PERMUTATIONS.npz`, `*_SHELL*_MOCK.npz`: intermediate null distributions.
- `FROZEN_SCIENTIFIC_STATE.json`, `INPUT_PROVENANCE.json`, `VALIDATION.json`, `RESULT_MANIFEST.json`: scientific state and integrity.
- `audit/`: original manifests, methods corrections and an explicitly superseded development run. Its free-PV-feature likelihoods are **not final results**.

Successes are real-data reconstruction, independent galaxy/PV linkage, full released covariance handling, reproducible matched/held-out tests and experimental W geometry. Failures/limits are incomplete exact light-curve calibration/error propagation, local foreground depth, imperfect permutation calibration, uneven subset behavior and lack of DES replication. The frozen model has not been rescued or retuned.
'''
 (ROOT/'FINAL_REPORT.md').write_text(report,encoding='utf-8')
 (ROOT/'README.md').write_text('''# BHSM SN sightline tomography V1

Read **FINAL_REPORT.md** first. The analysis preserves the previous scientific prediction. It does not establish a robust coherent topographic residual.

Install the pinned requirements in Python 3.14, then run:

```powershell
python -m pip install -r requirements.txt
python run_all.py
python package_outputs.py
```

Numerical stages use only bundled data. Prior result packages are checked when present; their original manifests are bundled for portable reproduction. Do not use files in `audit/SUPERSEDED_free_PV_feature.zip` as final results.

Scientific outputs are deterministic at numerical tolerance; gzip/ZIP timestamps and some acquisition metadata need not be byte-identical after a rerun. `RESULT_MANIFEST.json` hashes the delivered snapshot. The archive excludes itself and bytecode caches.
''',encoding='utf-8')
 print('Wrote FINAL_REPORT.md and figures',flush=True)
if __name__=='__main__': main()
