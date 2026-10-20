"""Consumed-data macro revisions; never reads the sealed LCO test rows."""
from __future__ import annotations

import argparse
from dataclasses import replace
from hashlib import sha256
import json
from pathlib import Path
import sys
import time

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path[:0] = [str(ROOT/'code'), str(ROOT/'scripts')]
import strong_lensing_macro_revision as macro
from strong_lensing_lightcurve_revision import RevisionParameters, log_likelihood
from fit_he0435_revision import load_development

OUT = ROOT/'artifacts/strong_lensing_revision_20261020'
DEVELOPMENT = OUT/'revision_development_fit.json'
INITIAL = ROOT/'artifacts/strong_lensing_observed_20261020/conventional_training.json'
PROTOCOL = ROOT/'preregistration/he0435_macro_revision_development_20261020.json'


def digest(path):
    return sha256(Path(path).read_bytes()).hexdigest()


def save(path, value):
    path.parent.mkdir(parents=True, exist_ok=True)
    def convert(x):
        if isinstance(x, np.ndarray):
            return x.tolist()
        if isinstance(x, np.generic):
            return x.item()
        raise TypeError(type(x).__name__)
    path.write_text(json.dumps(value, indent=2, default=convert, allow_nan=False)+'\n', encoding='utf-8', newline='\n')


def setup():
    receipt = json.loads(DEVELOPMENT.read_text())
    if receipt['target_photometry_consumed']:
        raise RuntimeError('Development input must not have consumed the unused targets')
    data = load_development()
    parameters = RevisionParameters(**receipt['map_fit']['parameters'])
    particles = np.array([p['delays_days'][1:] for p in receipt['particles']])
    weights = np.asarray(receipt['weights'])
    mean = weights@particles
    covariance = (particles-mean).T@((particles-mean)*weights[:, None])
    initial = np.asarray(json.loads(INITIAL.read_text())['parameter_vector'])
    def likelihood(delay):
        return log_likelihood(data, replace(parameters, delays_days=tuple(np.r_[0., delay])))
    return parameters, covariance, initial, likelihood


def protocol():
    if PROTOCOL.exists():
        return
    save(PROTOCOL, dict(
        status='PRESPECIFIED_CONSUMED_DATA_MACRO_COMPARISON',
        local_source='rho(r)=sum_j 3 M_j a_j^2/[4 pi (r^2+a_j^2)^(5/2)]; Menc(r)=sum_j M_j r^3/(r^2+a_j^2)^(3/2)',
        scalar_equation='C x^2 p -4 x p^2/lambda = -b Menc_natural(x), regular decaying branch; TOTAL enclosed mass consumed once, no addition of independent nonlinear scalar fields',
        metric_equation='2 x^2 Phi_prime + b x^2 p = 2 Menc_natural(x)',
        compact_mass_candidates_Msun=list(macro.COMPACT_MASS_CANDIDATES_MSUN),
        compact_scale_kpc=macro.COMPACT_SCALE_KPC,
        compact_source_scope='finite positive regular prescribed nuclear matter surrogate, not an inferred SMBH or a solved material EOS',
        continuous_parameters=list(macro.PARAMETER_NAMES),
        lower_bounds=macro.LOWER.tolist(), upper_bounds=macro.UPPER.tolist(),
        conventional_selection='Among candidates passing conservative HST summed-positive-PSF admissibility, minimize the consumed-data astrometry plus 3mas centroid-prior chi2 minus twice exact conditional MAP-nuisance development light log likelihood. All candidates remain visible.',
        likelihood_scope='Two-stage empirical Bayes: revised brightness model MAP hyperparameters fixed during macro optimization; each evaluated objective is its original exact state-space development likelihood, not a quadratic delay fit. Weighted posterior delay covariance only initializes the optimizer.',
        predictive_scope='After selection, finite nuisance-mixture particles transported to each macro delay vector with development-only likelihood-ratio reweighting. No complete macro/posterior integral or unique physical-state inference is claimed.',
        paired_comparison='Select common compact mass using conventional development plus HST; R1 refits nuisance parameters at that same compact mass, fixed state/screen-plane geometry. Separately report scalar activation on the fixed conventional body and maps.',
        image_census='Conventional degree-nine radial polynomial plus logarithmic-seed audit and original 2D residual verification; R1 continuation of all conventional roots plus logarithmic seeds, each verified in original total-source projection. Sum every additional-image PSF before HST comparison.',
        target_photometry_consumed=False,
        development_sha256=digest(DEVELOPMENT),
        source_sha256={p.relative_to(ROOT).as_posix(): digest(p) for p in (Path(__file__), ROOT/'code/strong_lensing_macro_revision.py')},
        physical_admission=False))


def candidate(index, *, exact_maxiter=70):
    protocol()
    mass = macro.COMPACT_MASS_CANDIDATES_MSUN[index]
    parameters, covariance, initial, likelihood = setup()
    started = time.perf_counter()
    result = macro.fit_macro(parameters.delays_days[1:], covariance, compact_mass_Msun=mass,
                             initial=initial, delay_log_likelihood=likelihood,
                             exact_maxiter=exact_maxiter, include_census=True)
    result.update(candidate_index=index, wall_seconds=time.perf_counter()-started,
                  target_photometry_consumed=False,
                  development_sha256=digest(DEVELOPMENT), protocol_sha256=digest(PROTOCOL))
    save(OUT/f'macro_candidate_{index}.json', result)
    hst = {'models': {f'compact_candidate_{index}': {'images': [
        {'position_light_centered_arcsec': row['position_arcsec'],
         'absolute_macro_flux_relative_A': row['flux_ratio_to_image_A_achromatic_no_differential_extinction_or_microlensing']}
        for row in result['additional_images']]}}}
    save(OUT/f'macro_candidate_{index}_hst_input.json', hst)
    print(json.dumps(dict(candidate=index, compact_mass_Msun=mass,
                          objective=result['exact_training_objective'], optimizer_success=result['optimizer_success'],
                          delays=result['predicted_delays_days'], image_count=result['global_image_census']['image_count'],
                          additional_images=result['additional_images'], wall_seconds=result['wall_seconds'])), flush=True)
    return result


def paired(index, *, exact_maxiter=70, warm_start=False):
    conventional = json.loads((OUT/f'macro_candidate_{index}.json').read_text())
    parameters, covariance, _, likelihood = setup()
    mass = conventional['compact_mass_Msun']
    initial = np.array(conventional['parameter_vector'])
    if warm_start and (OUT/'macro_selected_r1.json').exists() and not (OUT/'macro_selected_r1_drift_audit.json').exists():
        save(OUT/'macro_selected_r1_drift_audit.json', json.loads((OUT/'macro_selected_r1.json').read_text()))
    result = macro.fit_macro(parameters.delays_days[1:], covariance, compact_mass_Msun=mass,
                             model='r1', initial=initial, delay_log_likelihood=likelihood,
                             exact_maxiter=exact_maxiter, include_census=True,
                             skip_quadratic_initializer=warm_start,
                             retain_best_original_endpoint=warm_start)
    fixed = macro.macro_application(initial, compact_mass_Msun=mass, model='r1')
    same_body_conventional = macro.macro_application(initial, compact_mass_Msun=mass, model='conventional')
    result['fixed_conventional_body_R1_activation'] = dict(
        image_offsets_microarcsec=((fixed['images_light_centered_arcsec']-same_body_conventional['images_light_centered_arcsec'])*1e6).tolist(),
        delay_offsets_seconds=((fixed['delays_days']-same_body_conventional['delays_days'])*86400).tolist(),
        scope='same total matter body, source, rear screen, maps, clock and physical state; scalar field switched on once')
    result['separately_refitted_R1_minus_conventional_delay_seconds'] = ((np.array(result['predicted_delays_days'])-np.array(conventional['predicted_delays_days']))*86400).tolist()
    result.update(selected_candidate_index=index, target_photometry_consumed=False,
                  development_sha256=digest(DEVELOPMENT), protocol_sha256=digest(PROTOCOL),
                  source_sha256={p.relative_to(ROOT).as_posix(): digest(p) for p in (Path(__file__), ROOT/'code/strong_lensing_macro_revision.py')})
    save(OUT/'macro_selected_r1.json', result)
    save(OUT/'revised_r1_training.json', result)
    save(OUT/'revised_conventional_training.json', conventional)
    save(OUT/'macro_selected_r1_hst_input.json', {'models': {'selected_R1': {'images': [
        {'position_light_centered_arcsec': row['position_arcsec'],
         'absolute_macro_flux_relative_A': row['flux_ratio_to_image_A_achromatic_no_differential_extinction_or_microlensing']}
        for row in result['additional_images']]}}})
    print(json.dumps(dict(selected_candidate=index, model='r1', objective=result['exact_training_objective'],
                          optimizer_success=result['optimizer_success'], delays=result['predicted_delays_days'],
                          additional_images=result['additional_images'])), flush=True)


def select():
    hst = json.loads((OUT/'hst_macro_candidates.json').read_text())
    candidates = []
    for index, mass in enumerate(macro.COMPACT_MASS_CANDIDATES_MSUN):
        fitted = json.loads((OUT/f'macro_candidate_{index}.json').read_text())
        bands = hst['models'][f'compact_candidate_{index}']['bands']
        allowed = all(min(b['worst_psf_centroid_upper_scale_95'], b['worst_adversarial_upper_scale_95']) >= 1. for b in bands.values())
        candidates.append(dict(candidate_index=index, compact_mass_Msun=mass,
                               exact_training_objective=fitted['exact_training_objective'],
                               hst_admissible=allowed,
                               hst_upper_scale95_F814W=bands['F814W']['worst_psf_centroid_upper_scale_95'],
                               hst_stress_upper_scale95_F814W=bands['F814W']['worst_adversarial_upper_scale_95'],
                               hst_upper_scale95_F160W=bands['F160W']['worst_psf_centroid_upper_scale_95'],
                               additional_images=fitted['additional_images'],
                               numerical_optimizer_success=fitted['optimizer_success']))
    allowed = [c for c in candidates if c['hst_admissible']]
    if not allowed:
        raise RuntimeError('No candidate passes the declared conditional HST constraint')
    chosen = min(allowed, key=lambda c: c['exact_training_objective'])
    result = dict(status='CONSUMED_DEVELOPMENT_AND_HST_SELECTION_BEFORE_UNUSED_TARGETS',
                  selected_candidate_index=chosen['candidate_index'], selected_compact_mass_Msun=chosen['compact_mass_Msun'],
                  compact_scale_kpc=macro.COMPACT_SCALE_KPC, candidates=candidates,
                  selection_rule='HST admissibility under every retained band/PSF/centroid stress envelope, then minimum exact conditional-MAP development objective',
                  interpretation='A representative regular compact-source model; tiny development objective differences among allowed grid points do not measure nuclear mass or distinguish their physical composition',
                  hst_photometric_condition=hst['photometric_condition'],
                  target_photometry_consumed=False, physical_admission=False,
                  hst_candidates_sha256=digest(OUT/'hst_macro_candidates.json'),
                  development_sha256=digest(DEVELOPMENT), protocol_sha256=digest(PROTOCOL))
    save(OUT/'macro_selection.json', result)
    save(OUT/'revised_conventional_training.json', json.loads((OUT/f"macro_candidate_{chosen['candidate_index']}.json").read_text()))
    print(json.dumps(result), flush=True)


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('mode', choices=('candidate', 'paired', 'select'))
    parser.add_argument('--index', type=int)
    parser.add_argument('--exact-maxiter', type=int, default=70)
    parser.add_argument('--warm-start', action='store_true')
    args = parser.parse_args()
    if args.mode != 'select' and (args.index is None or not 0 <= args.index < len(macro.COMPACT_MASS_CANDIDATES_MSUN)):
        raise ValueError('candidate index outside prespecified compact-mass grid')
    if args.mode == 'candidate':
        candidate(args.index, exact_maxiter=args.exact_maxiter)
    elif args.mode == 'paired':
        paired(args.index, exact_maxiter=args.exact_maxiter, warm_start=args.warm_start)
    else:
        select()
