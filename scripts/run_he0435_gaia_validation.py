"""Prepare, commit, then separately expose and score Gaia DR3 once."""
from __future__ import annotations
import argparse,csv,hashlib,io,json,subprocess,sys
from datetime import datetime,timezone
from pathlib import Path
from urllib.parse import urlencode
from urllib.request import urlopen
import numpy as np

ROOT=Path(__file__).resolve().parents[1];sys.path.insert(0,str(ROOT/'code'))
from strong_lensing_gaia_validation import (tangent_covariance,relative_astrometry_covariance_mas2,
    match_four_sources,gaussian_comparison_by_floor)

OUT=ROOT/'artifacts/strong_lensing_revision_20261020/gaia'
DATA=ROOT/'data/strong_lensing_revision_20261020/gaia'
CODE=['code/strong_lensing_gaia_validation.py','scripts/run_he0435_gaia_validation.py',
      'tests/test_strong_lensing_gaia_validation.py','tests/test_strong_lensing_gaia_runner.py']
SOURCE_COMMIT='b4cf9a04a7962c03441356ee8f7e18e3da157165'
COLUMNS=['source_id','ra','dec','ra_error','dec_error','ra_dec_corr','ref_epoch',
 'astrometric_params_solved','astrometric_excess_noise','astrometric_excess_noise_sig',
 'ruwe','duplicated_source','ipd_frac_multi_peak','ipd_frac_odd_win',
 'ipd_gof_harmonic_amplitude','visibility_periods_used','astrometric_n_good_obs_al',
 'phot_g_mean_mag','pmra','pmdec','pmra_error','pmdec_error']


def sha(path):return hashlib.sha256(Path(path).read_bytes()).hexdigest()
def dump(path,data):
    path.parent.mkdir(parents=True,exist_ok=True)
    path.write_text(json.dumps(data,indent=2)+'\n',encoding='utf-8')
def committed_bytes(commit,path):
    return subprocess.check_output(['git','show',f'{commit}:{path}'],cwd=ROOT)


def start_query_lock(path,metadata):
    path.parent.mkdir(parents=True,exist_ok=True)
    with path.open('x',encoding='utf-8')as stream:
        json.dump(metadata,stream,indent=2);stream.write('\n')


def prepare(amend=False):
    protocol_path=OUT/'gaia_protocol.json';prediction_path=OUT/'frozen_astrometry_predictions.json'
    if (OUT/'query_started.json').exists():
        raise FileExistsError('Gaia exposure marker exists; preparation cannot change')
    if amend:
        old=protocol_path.read_bytes()
        version=0
        while (OUT/f'gaia_protocol_v{version}_pre_exposure.json').exists():version+=1
        (OUT/f'gaia_protocol_v{version}_pre_exposure.json').write_bytes(old)
    elif any(p.exists()for p in [protocol_path,prediction_path]):
        raise FileExistsError('Existing Gaia preparation/exposure must be preserved')
    models={};source={}
    for name,rel in [('conventional','artifacts/strong_lensing_revision_20261020/revised_conventional_training.json'),
                     ('conditional_R1','artifacts/strong_lensing_revision_20261020/revised_r1_training.json')]:
        raw=committed_bytes(SOURCE_COMMIT,rel);j=json.loads(raw)
        if raw!=(ROOT/rel).read_bytes():raise RuntimeError('Frozen macro source bytes changed')
        source[rel]=hashlib.sha256(raw).hexdigest()
        models[name]={'relative_B_C_D_to_A_east_north_arcsec':j['predicted_relative_astrometry_arcsec'][:6],
                      'images_light_centered_arcsec':j['predicted_images_light_centered_arcsec'],
                      'macro_parameter_vector':j['parameter_vector']}
    prediction={'status':'FROZEN_BEFORE_GAIA_NUMERICAL_POSITION_QUERY','macro_source_commit':SOURCE_COMMIT,
        'source_sha256':source,'models':models,'target_rows_consumed':False,'physical_admission':False,
        'unchanged_physical_guard':2e-10}
    dump(prediction_path,prediction)
    query='SELECT '+','.join(COLUMNS)+" FROM gaiadr3.gaia_source WHERE 1=CONTAINS(POINT('ICRS',ra,dec),CIRCLE('ICRS',69.562083333,-12.287333333,0.001)) ORDER BY source_id"
    protocol={'status':'FIXED_PRIOR_TO_GAIA_NUMERICAL_TARGET_EXPOSURE',
        'chronology':'Added after LCO primary outcome was known; no Gaia numerical position/error/quality rows inspected before this freeze',
        'metadata_count_only_sources':4,'tap_endpoint':'https://gea.esac.esa.int/tap-server/tap/sync',
        'adql':query,'format':'csv','labels':['A','B','C','D'],
        'association':'All24 permutations; minimum unweighted sum squared six conventional relative-coordinate residuals; same association used for R1. Report runner-up gap. No frame/scale/rotation fit.',
        'position_observable':'Fixed ICRS east/north B/C/D offsets at catalogue A, in mas, using the frozen helper tangent-plane convention; quasar macro geometry assumed stationary at Gaia epochJ2016.0; no target-derived proper-motion correction',
        'primary_covariance':'Published source ra_error (already cosdelta),dec_error,ra_dec_corr; common-A full6x6blocks; source-independent small-angle approximation, no unavailable between-source covariance invented',
        'coordinate_floor_controls_mas':[0.,2.,5.],
        'floor_definition':'After propagating published per-source errors and the common-A crossblocks, add floor^2 I6 independently to the six relative-coordinate variances; apply exactly once.',
        'statistics':['all six raw residuals in mas','coordinate RMS','maximum image vector residual','joint Gaussian conditional logdensity in mas^-6','conditional chi2'],
        'quality_policy':'Retain and report all four sources regardless of quality; missing mandatory astrometric errors make formal density/chi2 undefined but raw relative-position residuals are still reported. No reopening or filtering.',
        'limitations':['No full macro posterior is supplied; point-prediction scores are conditional.',
            'Gaia measurements are different2014-2017photons from2003/2012HST; source NIC2pointastrometry2011predatesGaia.',
            '2026HSTcentralpixelmosaicsuseGaiaEDR3fieldframealignment; independence of every calibration is not certified.',
            'Gaia close-image blending/chromatic photocentres may differ from HSTNIC2 astrometry; report diagnostics without post-target filtering.',
            'No new mass/screen/state/clock fit, added optical amplitude, complete-ring reconstruction or unique BHSM claim.'],
        'source_sha256':{p:sha(ROOT/p)for p in CODE},'prediction_sha256':sha(prediction_path),
        'single_exposure_lock':'query_started.json created exclusively before request; never removed or silently retried',
        'pre_exposure_implementation_amendment':('v0 preserved; standardized independent helper API, relative-coordinate floors, stdlibHTTP and raw-residual reporting for null errors before any target query'if amend else None)}
    dump(protocol_path,protocol)
    print(json.dumps({'prepared':str(protocol_path),'protocol_sha256':sha(protocol_path),
                      'prediction_sha256':sha(prediction_path),'targets_retrieved':False},indent=2))


def score_catalogue_rows(rows,models,floors):
    """Pure post-query integration, separately checked with synthetic rows."""
    conventional=models['conventional']['relative_B_C_D_to_A_east_north_arcsec']
    diag=match_four_sources([float(r['ra'])for r in rows],[float(r['dec'])for r in rows],conventional)
    association=diag['permutation_abcd'];observed=np.asarray(diag['relative_arcsec'],float).reshape(3,2)*1000.
    selected=[rows[i]for i in association]
    covariance_error=None
    try:
        source_cov=[tangent_covariance(float(r['ra_error']),float(r['dec_error']),float(r['ra_dec_corr']))for r in selected]
        basecov=relative_astrometry_covariance_mas2(source_cov,coordinate_floor_mas=0.)
    except (ValueError,TypeError)as exc:
        covariance_error=str(exc);basecov=None
    scored={}
    for model,spec in models.items():
        if basecov is not None:
            scored[model]=gaussian_comparison_by_floor(observed.reshape(6)/1000.,spec['relative_B_C_D_to_A_east_north_arcsec'],basecov,floors)
        else:
            residual=observed-np.asarray(spec['relative_B_C_D_to_A_east_north_arcsec']).reshape(3,2)*1000.
            scored[model]={'formal_score_undefined':covariance_error,'residual_mas':residual.tolist(),
                'rms_coordinate_residual_mas':float(np.sqrt(np.mean(residual**2)))}
    return diag,observed,selected,scored,basecov


def evaluate():
    p=OUT/'gaia_protocol.json';f=OUT/'frozen_astrometry_predictions.json'
    protocol=json.loads(p.read_text());pred=json.loads(f.read_text())
    for path in [p,f]+[ROOT/k for k in CODE]:
        rel=path.relative_to(ROOT).as_posix()
        if committed_bytes('HEAD',rel)!=path.read_bytes():raise RuntimeError('Gaia code/protocol/predictions must be committed before target query')
    for rel,expected in protocol['source_sha256'].items():
        if sha(ROOT/rel)!=expected:raise RuntimeError('Frozen Gaia source hash mismatch')
    if sha(f)!=protocol['prediction_sha256']:raise RuntimeError('Frozen Gaia prediction hash mismatch')
    for rel,expected in pred['source_sha256'].items():
        if sha(ROOT/rel)!=expected:raise RuntimeError('Already-frozen macro source hash mismatch')
    stamp={'status':'GAIA_NUMERICAL_TARGET_QUERY_STARTED_PERSISTENT_EXPOSURE',
        'time_utc':datetime.now(timezone.utc).isoformat(),'protocol_sha256':sha(p),
        'prediction_sha256':sha(f),'head_before_query':subprocess.check_output(['git','rev-parse','HEAD'],cwd=ROOT,text=True).strip()}
    start_query_lock(OUT/'query_started.json',stamp)
    queryurl=protocol['tap_endpoint']+'?'+urlencode({'REQUEST':'doQuery','LANG':'ADQL','FORMAT':'csv','QUERY':protocol['adql']})
    with urlopen(queryurl,timeout=120)as response:raw=response.read()
    DATA.mkdir(parents=True,exist_ok=True);target_path=DATA/'GaiaDR3_target_rows.csv';target_path.write_bytes(raw)
    rows=list(csv.DictReader(io.StringIO(raw.decode('utf-8'))))
    if len(rows)!=4:raise RuntimeError('Exactly four sources required; exposure remains marked')
    mandatory=['ra','dec']
    for row in rows:
        for col in mandatory:row[col]=float(row[col])
    diag,observed,selected,scored,basecov=score_catalogue_rows(rows,pred['models'],protocol['coordinate_floor_controls_mas'])
    result={'status':'EXECUTED_SINGLE_PREVIOUSLY_UNUSED_GAIA_ASTROMETRY_APPLICATION',
        'protocol_sha256':sha(p),'prediction_sha256':sha(f),'raw_target_sha256':sha(target_path),
        'source_commit_before_query':stamp['head_before_query'],'association':diag,
        'source_ids_A_B_C_D':[r['source_id']for r in selected],
        'observed_relative_B_C_D_to_A_east_north_mas':observed.tolist(),
        'catalogue_rows_A_B_C_D':selected,'models':scored,
        'conditional_R1_minus_conventional_joint_log_density':({key:scored['conditional_R1'][key]['joint_log_density_per_mas6']-scored['conventional'][key]['joint_log_density_per_mas6']for key in scored['conventional']}if basecov is not None else None),
        'limitations':protocol['limitations'],'numerical_target_rows_consumed':True,'physical_admission':False,
        'unchanged_physical_guard':2e-10}
    dump(OUT/'gaia_astrometry_score.json',result)
    stamp.update({'status':'GAIA_TARGETS_EXPOSED_AND_SCORED_ONCE','target_sha256':sha(target_path),'score_sha256':sha(OUT/'gaia_astrometry_score.json')})
    dump(OUT/'query_started.json',stamp)
    print(json.dumps({'score_file':str(OUT/'gaia_astrometry_score.json'),'models':scored,
        'association':diag},indent=2))


if __name__=='__main__':
    a=argparse.ArgumentParser();a.add_argument('mode',choices=['prepare','amend-before-exposure','evaluate']);args=a.parse_args()
    if args.mode=='evaluate':evaluate()
    else:prepare(amend=args.mode=='amend-before-exposure')
