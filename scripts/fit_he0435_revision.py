"""Development-only revised HE0435 fit and rolling-origin calibration.

Never reads the sealed 2019-onward LCO targets.  The entire already-consumed
Millon2020 table is now development data, including the former test.  LCO R
pre-2019 rows independently establish per-image cross-instrument offsets.
"""
from __future__ import annotations

import argparse
from dataclasses import asdict, replace
from hashlib import sha256
import json
from pathlib import Path
import sys
import time

import numpy as np

ROOT=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT/'code'))
from strong_lensing_lightcurve_revision import (
    RevisionCurves,RevisionParameters,concatenate,fit_revision,
    posterior_particles,mixture_predict,log_likelihood,
)

OLD=ROOT/'data/strong_lensing_observed_20261020/candidates/HE0435-1223/HE0435_all_Millon2020.dat'
LCO=ROOT/'data/strong_lensing_revision_20261020/calibration/he0435r_pre2019.dat'
OUT=ROOT/'artifacts/strong_lensing_revision_20261020'
PROTOCOL=ROOT/'preregistration/he0435_revision_development_20261020.json'
INSTRUMENTS=('Euler','EulerCAM','Maidanak','Mercator','SMARTS','LCO_R')


def digest(path):return sha256(Path(path).read_bytes()).hexdigest()


def encode_json(value):
    if isinstance(value,np.ndarray):return value.tolist()
    if isinstance(value,np.generic):return value.item()
    raise TypeError(type(value).__name__)


def save(path,value):
    path.parent.mkdir(parents=True,exist_ok=True)
    path.write_text(json.dumps(value,indent=2,default=encode_json,allow_nan=False)+'\n',
                    encoding='utf-8',newline='\n')


def load_development():
    old=[];instruments=[]
    for line in OLD.read_text().splitlines():
        tokens=line.split()
        if not tokens:continue
        old.append([float(x) for x in tokens[:9]])
        instruments.append(INSTRUMENTS.index(tokens[9]))
    old=np.asarray(old)
    lco=np.loadtxt(LCO,usecols=range(9))
    if np.max(lco[:,0])>=58484:
        raise RuntimeError('Calibration loader rejects post-cut target times')
    rows=np.vstack([old,lco]);ins=np.r_[instruments,np.full(len(lco),INSTRUMENTS.index('LCO_R'))]
    return RevisionCurves(np.repeat(rows[:,0],4),np.tile(np.arange(4),len(rows)),
        rows[:,1::2].ravel(),rows[:,2::2].ravel(),np.repeat(ins,4),4,INSTRUMENTS)


def write_protocol():
    if PROTOCOL.exists():raise RuntimeError('Preserve existing development protocol')
    save(PROTOCOL,{
        'status':'DECLARED_BEFORE_REVISED_NUMERICAL_SELECTION',
        'data_domain':'All 1002 formerly consumed COSMOGRAIL epochs plus 25 LCO R pre-2019 calibration epochs; NONE is a confirmatory target',
        'source_sha256':{OLD.relative_to(ROOT).as_posix():digest(OLD),LCO.relative_to(ROOT).as_posix():digest(LCO)},
        'mean_model':'m_i,k(t)=m0+b_i,k+q(t-d_i)+e_i(t-d_i)+a_i*(t-d_i-t0)/365.25+w_i(t-d_i)+epsilon_i',
        'state_model':'Shared stationary source OU; independent per-image stationary extrinsic OU; per-image Gaussian random slopes; independent Brownian extrinsic increments starting at t0; independent Gaussian proper per-image/instrument offsets; per-image independent jitter plus quoted marginal photometric errors.',
        'fixed_prior_scales':{'offset_sigma_mag':30.,'slope_sigma_mag_per_year':.03},
        'priors':'Uniform delays [-30,30] days relative to A; independent uniform LOG scales: source sigma [.02,.8]mag, tau[30,30000]day; each extrinsic sigma[.002,.8]mag, shared tau[200,30000]day; each jitter[.0003,.08]mag; each Brownian diffusion amplitude[.001,.2]mag/sqrt(year).',
        'bounds_are_prior_support_not_physical_constraints':True,
        'source_time_reference_rule':'floor(earliest native training epoch-31days); references are transformed consistently within each fit',
        'clock_origin_nuisance':{'old_COSMOGRAIL_offset_days':[-.5,0.,.5],
            'prior_weights':[1/3]*3,'LCO_offset_days':0.,
            'rule':'Apply old-release common offsets to source times. Update discrete weights using development likelihood only. Offset candidates are metadata uncertainty, not a fitted future-target shift.',
            'unresolved_subday_correction':'LCO MJD timescale/heliocentric correction unspecified; heliocentric-vs-topocentric difference up to approximately .006day is retained as a timing limitation.'},
        'initialization':'Deterministic A-referenced delays [0,6,-1,13] days copied from consumed first-experiment neighborhood only as numerical initializer; every fold likelihood re-estimates all delays.',
        'rolling_origin_folds':[{'train_inclusive':c,'test_exclusive_upper':c+500.} for c in (56000.,56800.,57600.)],
        'fold_eligibility':'Only instruments with at least 12 earlier epochs (48 image measurements) may supply validation targets; prevents treating wholly uncalibrated new-instrument zero points as microlensing extrapolation errors.',
        'primary_fold_mode':'All later C measurements in the fold are allowed shared-source driver; A/B/D withheld. This is batch cross-image prediction, not real-time forecasting.',
        'secondary_fold_mode':'No later driver; all four images withheld.',
        'mixture':'Final finite equal-weight random-walk MCMC posterior approximation. Two seeds, adaptation only during burn-in. Time-origin candidates update by development likelihood. Parameter draws, between-component spread and driver posterior weight changes are propagated. Short-chain limitations reported.',
        'calibration_rule':'For each mode and each target image, pool standardized absolute rolling-origin residuals; scale=max(1, empirical higher 97.5th percentile/1.95996398454). This conservative design targets 97.5% development marginal absolute-residual coverage to allow moderate domain shift; no exchangeability guarantee. Primary C is driver, hence identity scale. Per-image scales are frozen before unused targets.',
        'calibration_scoring_rule':'Apply diagonal affine scale to the entire mixture about its mean; normalize joint density with log Jacobian. Report actual unused-target 68.27% and 95% coverage and proper log scores, never assume calibration establishes validity.',
        'model_family_selection':'No post-target model selection. Retain OU+Brownian+random-slope family regardless of development result; remaining undercoverage or excessive width is reported.',
        'not_a_BHSM_state_kernel':True,'physical_admission':False,
    })
    print('PROTOCOL',digest(PROTOCOL),flush=True)


def run_fold(cut,window=500.):
    data=load_development();train=data.subset(data.time_days<=cut)
    eligible=[i for i in range(len(INSTRUMENTS)) if np.sum((train.instrument==i)&(train.image==0))>=12]
    mask=(data.time_days>cut)&(data.time_days<cut+window)&np.isin(data.instrument,eligible)
    later=data.subset(mask)
    initial=RevisionParameters((0.,6.,-1.,13.))
    fitted,receipt=fit_revision(train,initial=initial,maxiter=160)
    particles,diagnostics=posterior_particles(train,fitted,seed=20261020+int(cut),burn=1200,draws=2400,thin=300)
    result={'cut':cut,'window_days':window,'train_observations':len(train.time_days),'later_observations':len(later.time_days),
            'eligible_instruments':[INSTRUMENTS[i] for i in eligible],
            'fit':receipt,'posterior_diagnostics':diagnostics,'particles':particles,'modes':{}}
    for mode in ('primary_C_driver','secondary_forecast'):
        if mode=='primary_C_driver':
            driver=later.subset(later.image==2);targets=later.subset(later.image!=2)
        else:driver=None;targets=later
        prediction=mixture_predict(train,targets,particles,driver=driver)
        standardized=(targets.magnitude-prediction['mean'])/np.sqrt(prediction['variance'])
        result['modes'][mode]={'image':targets.image,'standardized_residual':standardized,
            'time_days':targets.time_days,'instrument':targets.instrument,
            'prediction_mean':prediction['mean'],'prediction_sigma':np.sqrt(prediction['variance']),
            'magnitude':targets.magnitude,'unscaled_95_gaussian_moment_coverage':float(np.mean(abs(standardized)<=1.95996398454))}
    prefix='development_long_fold' if window>500 else 'development_fold'
    save(OUT/f'{prefix}_{int(cut)}.json',result)
    print(json.dumps({'fold':cut,'fit_log_likelihood':receipt['log_likelihood'],
                     'coverage':{k:v['unscaled_95_gaussian_moment_coverage'] for k,v in result['modes'].items()}}),flush=True)


def calibrate():
    paths=[OUT/f'development_fold_{int(c)}.json' for c in (56000.,56800.,57600.)]
    paths += [OUT/f'development_long_fold_{int(c)}.json' for c in (54800.,55600.,56400.)
              if (OUT/f'development_long_fold_{int(c)}.json').exists()]
    folds=[json.loads(path.read_text()) for path in paths]
    result={'development_only':True,'independent_validation':False,'modes':{},
            'protocol_sha256':digest(PROTOCOL),'fold_paths':[p.relative_to(ROOT).as_posix() for p in paths],
            'overlapping_temporal_folds_are_not_independent':True}
    for mode in ('primary_C_driver','secondary_forecast'):
        scales=[];images={}
        for image in range(4):
            z=np.concatenate([np.asarray(f['modes'][mode]['standardized_residual'])[
                np.asarray(f['modes'][mode]['image'])==image] for f in folds])
            if len(z):
                quantile=float(np.quantile(abs(z),.975,method='higher'))
                scale=max(1.,quantile/1.95996398454)
                images[str(image)]={'count':len(z),'quantile_abs_z_97_5':quantile,
                    'unscaled_95_moment_coverage':float(np.mean(abs(z)<=1.95996398454)),
                    'scaled_95_moment_coverage':float(np.mean(abs(z)<=1.95996398454*scale)),
                    'maximum_abs_z':float(max(abs(z))),'scale':scale}
            else:scale=1.;images[str(image)]={'count':0,'scale':1.,'reason':'Image C is the allowed primary driver'}
            scales.append(scale)
        result['modes'][mode]={'per_image_scale':scales,'by_image':images}
    save(OUT/'development_calibration.json',result)
    print(json.dumps(result),flush=True)


def final_fit():
    data=load_development();started=time.perf_counter()
    base=RevisionParameters((0.,6.,-1.,13.))
    fitted,receipt=fit_revision(data,initial=base,maxiter=240)
    save(OUT/'revision_map_fit.json',receipt)
    particles=[];chains=[]
    for seed in (20261020,20261021):
        components,diagnostics=posterior_particles(data,fitted,seed=seed,burn=2400,draws=4800,thin=400)
        particles.extend(components);chains.append(diagnostics)
        print(json.dumps({'chain_seed':seed,'diagnostics':diagnostics}),flush=True)
    # The old-release clock ambiguity is marginalized on the finite retained
    # posterior approximation.  A candidate's likelihood ratio reweights each
    # existing sample, not a future-target fit or unnormalized arbitrary score.
    base_log=np.array([log_likelihood(data,RevisionParameters(**p)) for p in particles])
    expanded=[];logweights=[]
    for offset in (-.5,0.,.5):
        shifts=tuple(0. if name=='LCO_R' else offset for name in INSTRUMENTS)
        for p,oldll in zip(particles,base_log):
            shifted=replace(RevisionParameters(**p),instrument_time_offset_days=shifts)
            expanded.append(asdict(shifted));logweights.append(log_likelihood(data,shifted)-oldll)
    from scipy.special import logsumexp
    weights=np.exp(np.array(logweights)-logsumexp(logweights))
    count=len(particles);origin_weights=[float(np.sum(weights[j*count:(j+1)*count])) for j in range(3)]
    result={'status':'DEVELOPMENT_FIT_READY_FOR_MACRO_AND_PREDICTIVE_FREEZE',
        'data_observations':len(data.time_days),'old_epochs':1002,'LCO_calibration_epochs':25,
        'maximum_native_training_time':float(np.max(data.time_days)),
        'target_photometry_consumed':False,'map_fit':receipt,'particles':expanded,
        'weights':weights,'posterior_diagnostics':chains,
        'time_origin_candidates_days':[-.5,0.,.5],'time_origin_posterior_weights':origin_weights,
        'time_origin_importance_effective_sample_size':float(1/np.sum(weights**2)),
        'calibration':json.loads((OUT/'development_calibration.json').read_text()),
        'source_sha256':{p.relative_to(ROOT).as_posix():digest(p) for p in
            (OLD,LCO,PROTOCOL,Path(__file__),ROOT/'code/strong_lensing_lightcurve_revision.py',
             ROOT/'code/strong_lensing_revision_filter.c')},
        'wall_seconds':time.perf_counter()-started,
        'qualification':'The finite two-chain mixture approximates the declared bounded-prior posterior; no claim of complete global mixing. Calibration and model development reuse consumed data and are not independent confirmation.',
        'physical_admission':False}
    save(OUT/'revision_development_fit.json',result)
    print(json.dumps({'fit_log_likelihood':receipt['log_likelihood'],'time_origin_weights':origin_weights,
                      'particles':len(expanded),'seconds':result['wall_seconds']}),flush=True)


if __name__=='__main__':
    parser=argparse.ArgumentParser();parser.add_argument('mode',choices=('protocol','fold','calibrate','final'))
    parser.add_argument('--cut',type=float);parser.add_argument('--window',type=float,default=500.)
    args=parser.parse_args()
    if args.mode=='protocol':write_protocol()
    else:
        if not PROTOCOL.exists():raise RuntimeError('Declare development protocol first')
        if args.mode=='fold':run_fold(args.cut,args.window)
        elif args.mode=='calibrate':calibrate()
        else:final_fit()
