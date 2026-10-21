"""Classical empirical HE0435 variability model with predictive state uncertainty.

This does not define a BHSM quantum statistical kernel.  The original
observed-lightcurve module and its scored experiment remain immutable.

m_ik(t) = m0 + b_ik + q(t-d_i) + e_i(t) + a_i*(t-d_i-t0)/year
          + w_i(t-d_i) + independent measurement noise.

q and e_i are independent stationary OUs.  a_i are Gaussian random slopes;
w_i are independent Brownian motions initialized to zero at t0.  b_ik are
proper Gaussian per-image/instrument offsets.  Quoted errors and per-image
jitter are independent conditional on the latent state.  No cross-image
measurement covariance is fabricated.  A finite retained MCMC mixture
integrates the declared delay/hyperparameter posterior approximately.

Development-only per-image calibration is an invertible diagonal affine
transform about the *mixture* predictive mean.  Therefore even the joint
calibrated predictive density is normalized (including its Jacobian).
Its scale estimates have no finite-sample exchangeability guarantee for
temporally dependent observations.  Unused target coverage must be measured.
"""
from __future__ import annotations

from dataclasses import asdict, dataclass, replace
import ctypes
from hashlib import sha256
from pathlib import Path
import shutil
import subprocess
import tempfile

import numpy as np
from scipy.linalg import solve
from scipy.optimize import brentq, minimize
from scipy.special import logsumexp, ndtr


YEAR = 365.25
LOG2PI = float(np.log(2*np.pi))


@dataclass(frozen=True)
class RevisionCurves:
    time_days: np.ndarray
    image: np.ndarray
    magnitude: np.ndarray
    error_mag: np.ndarray
    instrument: np.ndarray
    nimage: int = 4
    instrument_names: tuple[str, ...] = ('reference',)

    def __post_init__(self):
        object.__setattr__(self, 'instrument_names', tuple(self.instrument_names))
        for name, dtype in [('time_days', float), ('image', int), ('magnitude', float),
                            ('error_mag', float), ('instrument', int)]:
            value = np.array(getattr(self, name), dtype=dtype, copy=True)
            if value.ndim != 1 or not np.all(np.isfinite(value)):
                raise ValueError('Finite one-dimensional observation arrays required')
            value.setflags(write=False)
            object.__setattr__(self, name, value)
        if len({len(getattr(self, k)) for k in
                ('time_days','image','magnitude','error_mag','instrument')}) != 1:
            raise ValueError('Observation arrays must have the same length')
        if len(self.time_days) == 0 or np.any(self.error_mag <= 0):
            raise ValueError('Nonempty data and positive uncertainties required')
        if (np.any(self.image < 0) or np.any(self.image >= self.nimage)
                or np.any(self.instrument < 0)
                or np.any(self.instrument >= len(self.instrument_names))):
            raise ValueError('Image/instrument identifier outside declared domain')

    def subset(self, mask):
        return RevisionCurves(*(getattr(self,k)[mask] for k in
            ('time_days','image','magnitude','error_mag','instrument')),
            self.nimage, self.instrument_names)


@dataclass(frozen=True)
class RevisionParameters:
    delays_days: tuple[float, ...]
    source_sigma_mag: float = .2
    source_tau_days: float = 2000.
    extrinsic_sigma_mag: tuple[float, ...] = (.1,.1,.1,.1)
    extrinsic_tau_days: float = 5000.
    jitter_mag: tuple[float, ...] = (.01,.01,.01,.01)
    brownian_mag_per_sqrt_year: tuple[float, ...] = (.02,.02,.02,.02)
    slope_sigma_mag_per_year: float = .03
    offset_sigma_mag: float = 30.
    reference_time_days: float = 52000.
    magnitude_zero: float = -10.
    instrument_time_offset_days: tuple[float, ...] = ()

    def __post_init__(self):
        ni=len(self.delays_days)
        for name in ('delays_days','extrinsic_sigma_mag','jitter_mag','brownian_mag_per_sqrt_year'):
            value=tuple(float(x) for x in getattr(self,name))
            if len(value)!=ni or not np.all(np.isfinite(value)):
                raise ValueError('One finite value per image is required')
            if name!='delays_days' and min(value)<=0:
                raise ValueError('Variability scales must be positive')
            object.__setattr__(self,name,value)
        if abs(self.delays_days[0])>1e-12:
            raise ValueError('Reference-image delay must be zero')
        object.__setattr__(self,'instrument_time_offset_days',
                           tuple(float(x) for x in self.instrument_time_offset_days))
        if not np.all(np.isfinite(self.instrument_time_offset_days)):
            raise ValueError('Finite instrument time offsets required')
        if any(not np.isfinite(x) or x<=0 for x in (self.source_sigma_mag,
                self.source_tau_days,self.extrinsic_tau_days,self.offset_sigma_mag,
                self.slope_sigma_mag_per_year)):
            raise ValueError('Proper positive finite model scales required')


def concatenate(first, second):
    if (first.nimage,first.instrument_names)!=(second.nimage,second.instrument_names):
        raise ValueError('Image and instrument domains must match')
    return RevisionCurves(*(np.r_[getattr(first,k),getattr(second,k)] for k in
        ('time_days','image','magnitude','error_mag','instrument')),
        first.nimage,first.instrument_names)


def _setup(data,p):
    ni=data.nimage
    if len(p.delays_days)!=ni:
        raise ValueError('One delay per image is required')
    offsets=p.instrument_time_offset_days or (0.,)*len(data.instrument_names)
    if len(offsets)!=len(data.instrument_names):
        raise ValueError('One time offset per instrument is required')
    source_time=data.time_days+np.asarray(offsets)[data.instrument]-np.asarray(p.delays_days)[data.image]
    if np.min(source_time)<p.reference_time_days:
        raise ValueError('Brownian reference must precede every shifted observation')
    variance=np.r_[p.source_sigma_mag**2,np.asarray(p.extrinsic_sigma_mag)**2,
        np.full(ni,p.slope_sigma_mag_per_year**2),np.zeros(ni),
        np.full(ni*len(data.instrument_names),p.offset_sigma_mag**2)]
    rate=np.r_[1/p.source_tau_days,np.full(ni,1/p.extrinsic_tau_days),
               np.zeros(ni*(2+len(data.instrument_names)))]
    diffusion=np.zeros(len(variance))
    diffusion[1+2*ni:1+3*ni]=np.asarray(p.brownian_mag_per_sqrt_year)**2/YEAR
    # Initialize the nonstationary process at the actual first source time.
    variance[1+2*ni:1+3*ni]=diffusion[1+2*ni:1+3*ni]*(np.min(source_time)-p.reference_time_days)
    return source_time,variance,rate,diffusion


def _row(data,p,original,source_time):
    ni=data.nimage; image=data.image[original]
    selected=np.array([0,1+image,1+ni+image,1+2*ni+image,
                       1+3*ni+data.instrument[original]*ni+image])
    h=np.array([1.,1.,(source_time[original]-p.reference_time_days)/YEAR,1.,1.])
    return selected,h


_ACCELERATOR=None
_ACCELERATOR_TRIED=False


def _accelerator():
    """Compile the optional 80-line accelerator; no generated binary is tracked."""
    global _ACCELERATOR,_ACCELERATOR_TRIED
    if _ACCELERATOR_TRIED:
        return _ACCELERATOR
    _ACCELERATOR_TRIED=True
    compiler=shutil.which('gcc')
    source=Path(__file__).with_name('strong_lensing_revision_filter.c')
    if not compiler or not source.exists():
        return None
    key=sha256(source.read_bytes()).hexdigest()[:16]
    directory=Path(tempfile.gettempdir())/'bhsm_revision_filter'/key
    directory.mkdir(parents=True,exist_ok=True)
    binary=directory/('filter.dll' if __import__('os').name=='nt' else 'filter.so')
    if not binary.exists():
        result=subprocess.run([compiler,'-shared','-O3',str(source),'-o',str(binary),'-lm'],
                              capture_output=True)
        if result.returncode:
            return None
    try:
        function=ctypes.CDLL(str(binary)).revision_filter
        fp=ctypes.POINTER(ctypes.c_double); ip=ctypes.POINTER(ctypes.c_int)
        function.argtypes=[ctypes.c_int,ctypes.c_int,ctypes.c_int,fp,ip,ip,fp,fp,
                           fp,fp,fp,fp,ctypes.c_double,ctypes.c_double,ctypes.c_double]
        function.restype=ctypes.c_double
        _ACCELERATOR=function
    except OSError:
        return None
    return function


def _filter(data,p,assimilate=None,save=False):
    source_time,variance,rate,diffusion=_setup(data,p)
    order=np.argsort(source_time,kind='stable'); ni=data.nimage
    if assimilate is None:
        assimilate=np.ones(len(order),dtype=bool)
    if np.shape(assimilate)!=np.shape(order):
        raise ValueError('Assimilation mask must match the observation domain')
    k=len(variance); count=len(order)
    m=np.zeros(k); covariance=np.diag(variance); ll=0.
    if save:
        trace={key:np.empty((count,k,k)) for key in ('predicted_covariance','filtered_covariance')}
        trace.update({key:np.empty((count,k)) for key in
                      ('predicted_mean','filtered_mean','transitions')})
    for pos,original in enumerate(order):
        dt=0. if pos==0 else source_time[original]-source_time[order[pos-1]]
        a=np.exp(-rate*dt); m*=a
        covariance*=a[:,None]*a[None,:]
        # Brownian initial variance is nonstationary: only rate-zero diffusion
        # contributes after initialization, never an OU stationary reset.
        covariance.flat[::k+1]+=variance*(-np.expm1(-2*rate*dt))+diffusion*dt
        if save:
            trace['predicted_mean'][pos]=m
            trace['predicted_covariance'][pos]=covariance
            trace['transitions'][pos]=a
        if assimilate[original]:
            indices,h=_row(data,p,original,source_time)
            v=covariance[:,indices]@h
            s=h@v[indices]+data.error_mag[original]**2+p.jitter_mag[data.image[original]]**2
            residual=data.magnitude[original]-p.magnitude_zero-h@m[indices]
            ll-=.5*(LOG2PI+np.log(s)+residual**2/s)
            m+=v*(residual/s)
            covariance-=np.outer(v,v)/s
            covariance=.5*(covariance+covariance.T)
        if save:
            trace['filtered_mean'][pos]=m
            trace['filtered_covariance'][pos]=covariance
    result={'log_likelihood':float(ll),'order':order}
    if save:
        result.update(trace)
    return result


def log_likelihood(data,p,*,backend='auto'):
    accelerator=_accelerator() if backend!='python' else None
    if accelerator is None:
        return _filter(data,p)['log_likelihood']
    time,variance,rate,diffusion=_setup(data,p)
    order=np.argsort(time,kind='stable')
    arrays=[np.ascontiguousarray(a,dtype=dtype) for a,dtype in
        [(time[order],np.float64),(data.image[order],np.int32),
         (data.instrument[order],np.int32),(data.magnitude[order],np.float64),
         (data.error_mag[order],np.float64),(variance,np.float64),(rate,np.float64),
         (diffusion,np.float64),(p.jitter_mag,np.float64)]]
    pointers=[a.ctypes.data_as(ctypes.POINTER(ctypes.c_int if a.dtype==np.int32
                                              else ctypes.c_double)) for a in arrays]
    return float(accelerator(len(time),len(variance),data.nimage,*pointers,
                             p.magnitude_zero,p.reference_time_days,YEAR))


def dense_covariance(data,p):
    u,_,_,_=_setup(data,p)
    same_image=data.image[:,None]==data.image[None,:]
    same_instrument=data.instrument[:,None]==data.instrument[None,:]
    sep=np.abs(u[:,None]-u[None,:])
    covariance=p.source_sigma_mag**2*np.exp(-sep/p.source_tau_days)
    covariance+=same_image*np.asarray(p.extrinsic_sigma_mag)[data.image,None]**2*np.exp(-sep/p.extrinsic_tau_days)
    covariance+=same_image*p.slope_sigma_mag_per_year**2*np.outer((u-p.reference_time_days)/YEAR,(u-p.reference_time_days)/YEAR)
    covariance+=same_image*np.asarray(p.brownian_mag_per_sqrt_year)[data.image,None]**2*np.minimum(u[:,None]-p.reference_time_days,u[None,:]-p.reference_time_days)/YEAR
    covariance+=same_image*same_instrument*p.offset_sigma_mag**2
    covariance.flat[::len(u)+1]+=data.error_mag**2+np.asarray(p.jitter_mag)[data.image]**2
    return covariance


def predict(train,test,p):
    """Exact Gaussian state conditioning; test magnitudes are never consumed."""
    data=concatenate(train,test); source_time,_,_,_=_setup(data,p)
    mask=np.arange(len(data.time_days))<len(train.time_days)
    trace=_filter(data,p,mask,save=True)
    mean=trace['filtered_mean'].copy(); covariance=trace['filtered_covariance'].copy()
    for pos in range(len(trace['order'])-2,-1,-1):
        numerator=trace['filtered_covariance'][pos]*trace['transitions'][pos+1][None,:]
        # Brownian states may be exactly zero at t0, hence semidefinite.
        try:
            gain=solve(trace['predicted_covariance'][pos+1],numerator.T,
                       assume_a='pos',check_finite=False).T
        except np.linalg.LinAlgError:
            gain=numerator@np.linalg.pinv(trace['predicted_covariance'][pos+1],rcond=1e-12)
        mean[pos]+=gain@(mean[pos+1]-trace['predicted_mean'][pos+1])
        covariance[pos]+=gain@(covariance[pos+1]-trace['predicted_covariance'][pos+1])@gain.T
        covariance[pos]=.5*(covariance[pos]+covariance[pos].T)
    inverse=np.argsort(trace['order']); positions=inverse[len(train.time_days):]
    out_mean=[]; out_variance=[]
    for j,pos in enumerate(positions):
        original=len(train.time_days)+j
        selected,h=_row(data,p,original,source_time)
        out_mean.append(p.magnitude_zero+h@mean[pos,selected])
        out_variance.append(h@covariance[pos][np.ix_(selected,selected)]@h+
                            test.error_mag[j]**2+p.jitter_mag[test.image[j]]**2)
    if np.min(out_variance)<=0:
        raise ArithmeticError('Nonpositive conditional predictive variance')
    return {'mean':np.array(out_mean),'variance':np.array(out_variance),
            'conditioning_log_likelihood':trace['log_likelihood']}


def mixture_predict(train,test,particles,weights=None,*,driver=None,calibration=None):
    """Condition mixture weights on an allowed driver, then predict targets."""
    particles=[RevisionParameters(**p) if isinstance(p,dict) else p for p in particles]
    weights=np.full(len(particles),1/len(particles)) if weights is None else np.asarray(weights,dtype=float)
    if np.any(weights<0) or not np.isclose(weights.sum(),1.):
        raise ValueError('Mixture weights must be normalized and nonnegative')
    conditioning=train if driver is None else concatenate(train,driver)
    log_weights=np.log(weights)
    if driver is not None:
        log_weights+=np.array([log_likelihood(conditioning,p)-log_likelihood(train,p) for p in particles])
        log_weights-=logsumexp(log_weights)
    effective_weights=np.exp(log_weights)
    components=[predict(conditioning,test,p) for p in particles]
    means=np.array([p['mean'] for p in components]); variances=np.array([p['variance'] for p in components])
    center=effective_weights@means
    scales=np.ones(test.nimage) if calibration is None else np.asarray(calibration,dtype=float)
    if scales.shape!=(test.nimage,) or np.any(scales<=0):
        raise ValueError('One positive predictive calibration scale per image required')
    transform=scales[test.image]
    means=center[None,:]+transform[None,:]*(means-center[None,:])
    variances*=transform[None,:]**2
    variance=effective_weights@(variances+(means-center[None,:])**2)
    quantiles={}
    for probability in (.025,.158655253931457,.841344746068543,.975):
        values=[]
        for j in range(len(center)):
            sig=np.sqrt(variances[:,j]); mu=means[:,j]
            cdf=lambda x: float(effective_weights@ndtr((x-mu)/sig)-probability)
            values.append(brentq(cdf,float(np.min(mu-12*sig)),float(np.max(mu+12*sig))))
        quantiles[str(probability)]=np.array(values)
    return {'mean':center,'variance':variance,'component_mean':means,
            'component_variance':variances,'weights':effective_weights,
            'quantiles':quantiles,'conditioning_log_likelihood':
                np.array([c['conditioning_log_likelihood'] for c in components]),
            'calibration':scales}


def score_mixture(train,test,particles,weights=None,*,driver=None,calibration=None):
    prediction=mixture_predict(train,test,particles,weights,driver=driver,calibration=calibration)
    center=prediction['mean']; residual=test.magnitude-center
    log_weights=np.log(prediction['weights']); means=prediction['component_mean']
    variances=prediction['component_variance']
    marginal=logsumexp(log_weights[:,None]-.5*(LOG2PI+np.log(variances)+(test.magnitude[None,:]-means)**2/variances),axis=0)
    transform=prediction['calibration'][test.image]
    transformed=replace(test,magnitude=center+(test.magnitude-center)/transform)
    conditioning=train if driver is None else concatenate(train,driver)
    parameters=[RevisionParameters(**p) if isinstance(p,dict) else p for p in particles]
    joint=logsumexp(log_weights+np.array([log_likelihood(concatenate(conditioning,transformed),p)
            for p in parameters])-prediction['conditioning_log_likelihood'])-np.sum(np.log(transform))
    def summary(mask):
        lower=prediction['quantiles']['0.025'][mask]
        upper=prediction['quantiles']['0.975'][mask]
        actual=test.magnitude[mask]
        interval_score=(upper-lower)+40*np.maximum(lower-actual,0)+40*np.maximum(actual-upper,0)
        return {'count':int(np.sum(mask)),'rmse_mag':float(np.sqrt(np.mean(residual[mask]**2))),
            'mean_residual_mag':float(np.mean(residual[mask])),
            'coverage_95':float(np.mean((test.magnitude[mask]>=prediction['quantiles']['0.025'][mask])&
                                       (test.magnitude[mask]<=prediction['quantiles']['0.975'][mask]))),
            'coverage_68_27':float(np.mean((test.magnitude[mask]>=prediction['quantiles']['0.158655253931457'][mask])&
                                       (test.magnitude[mask]<=prediction['quantiles']['0.841344746068543'][mask]))),
            'mean_predictive_sigma_mag':float(np.mean(np.sqrt(prediction['variance'][mask]))),
            'mean_95_interval_width_mag':float(np.mean(upper-lower)),
            'mean_95_interval_score_mag':float(np.mean(interval_score)),
            'mean_marginal_negative_log_density':float(-np.mean(marginal[mask]))}
    return {**prediction,'joint_conditional_log_density':float(joint),
        'overall':summary(np.ones(len(center),dtype=bool)),
        'by_image':{str(i):summary(test.image==i) for i in range(test.nimage) if np.any(test.image==i)},
        'residual':residual}


def encode(p):
    return np.r_[p.delays_days[1:],np.log([p.source_sigma_mag,p.source_tau_days,
        *p.extrinsic_sigma_mag,p.extrinsic_tau_days,*p.jitter_mag,*p.brownian_mag_per_sqrt_year])]


def decode(x,template):
    ni=len(template.delays_days); hyper=np.exp(np.asarray(x)[ni-1:]); position=2
    extrinsic=hyper[position:position+ni];position+=ni
    tau=hyper[position];position+=1
    jitter=hyper[position:position+ni];position+=ni
    brownian=hyper[position:position+ni]
    return replace(template,delays_days=tuple(np.r_[0.,x[:ni-1]]),
        source_sigma_mag=float(hyper[0]),source_tau_days=float(hyper[1]),
        extrinsic_sigma_mag=tuple(extrinsic),extrinsic_tau_days=float(tau),
        jitter_mag=tuple(jitter),brownian_mag_per_sqrt_year=tuple(brownian))


def coordinate_bounds(ni=4):
    return [(-30.,30.)]*(ni-1)+[tuple(np.log(pair)) for pair in
        [(.02,.8),(30.,30000.)]+[(.002,.8)]*ni+[(200.,30000.)]+[(.0003,.08)]*ni+[(.001,.2)]*ni]


def fit_revision(data,*,initial=None,maxiter=120):
    template=RevisionParameters(tuple(np.zeros(data.nimage)),
        magnitude_zero=float(np.median(data.magnitude)),
        reference_time_days=float(np.floor(np.min(data.time_days)-31.)))
    if initial is not None:
        template=replace(initial,magnitude_zero=template.magnitude_zero,
                         reference_time_days=template.reference_time_days)
    bounds=coordinate_bounds(data.nimage)
    def objective(x):
        value=-log_likelihood(data,decode(x,template))
        return value if np.isfinite(value) else 1e100
    optimum=minimize(objective,encode(template),method='L-BFGS-B',bounds=bounds,
        options={'maxiter':maxiter,'ftol':1e-9,'gtol':1e-4,'eps':1e-5,'maxls':30})
    parameters=decode(optimum.x,template)
    return parameters,{'parameters':asdict(parameters),'coordinates':optimum.x.tolist(),
        'log_likelihood':float(-optimum.fun),'success':bool(optimum.success),
        'message':str(optimum.message),'iterations':int(optimum.nit),'evaluations':int(optimum.nfev),
        'coordinate_bounds':bounds,'active_bounds':[bool(min(abs(x-lo),abs(x-hi))<1e-4*(hi-lo))
            for x,(lo,hi) in zip(optimum.x,bounds)],
        'likelihood_backend':'C accelerator' if _accelerator() is not None else 'Python reference'}


def posterior_particles(data,parameters,*,seed=20261020,burn=1200,draws=2400,thin=40):
    """Finite adaptive random-walk posterior approximation in log-scale priors.

    Priors are independent uniform on delay and declared log hyperparameter
    bounds.  Adaptation stops at burn-in.  Every retained equal-weight sample
    represents one posterior component.  Diagnostics are reported rather than
    claiming exact posterior integration or convergence from a short chain.
    """
    rng=np.random.default_rng(seed); x=encode(parameters); bounds=np.array(coordinate_bounds(data.nimage))
    scale=np.r_[np.full(data.nimage-1,.25),np.full(len(x)-(data.nimage-1),.12)]
    value=log_likelihood(data,decode(x,parameters)); chain=[]; likelihoods=[]
    accepted=0; blockaccepted=0; scalefactor=1.
    for iteration in range(burn+draws):
        proposal=x+rng.normal(size=len(x))*scale*scalefactor
        if np.all(proposal>=bounds[:,0]) and np.all(proposal<=bounds[:,1]):
            proposed=log_likelihood(data,decode(proposal,parameters))
            if np.isfinite(proposed) and np.log(rng.uniform())<proposed-value:
                x=proposal;value=proposed;accepted+=1;blockaccepted+=1
        if iteration<burn and (iteration+1)%100==0:
            scalefactor*=np.exp(np.clip(blockaccepted/100-.234,-.2,.2))
            blockaccepted=0
        if iteration>=burn:
            chain.append(x.copy());likelihoods.append(value)
    chain=np.array(chain); lag1=[]
    for j in range(chain.shape[1]):
        variance=np.var(chain[:,j])
        lag1.append(float(np.corrcoef(chain[:-1,j],chain[1:,j])[0,1]) if variance>0 else 1.)
    particles=[asdict(decode(x,parameters)) for x in chain[thin//2::thin]]
    return particles,{'seed':seed,'burn':burn,'post_burn_draws':draws,'thin':thin,
        'retained_particles':len(particles),'acceptance_fraction':accepted/(burn+draws),
        'proposal_scale_factor_after_burn':scalefactor,'lag1_coordinate_autocorrelation':lag1,
        'chain_coordinate_mean':chain.mean(axis=0).tolist(),'chain_coordinate_sd':chain.std(axis=0).tolist(),
        'post_burn_log_likelihood_range':[float(min(likelihoods)),float(max(likelihoods))],
        'qualification':'Finite single-chain posterior approximation; lag-one autocorrelation alone does not establish global convergence. Development predictive calibration is separately required.'}
