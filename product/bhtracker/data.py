"""Versioned measurements and full covariance, always in stated physical units."""
from dataclasses import dataclass
from pathlib import Path
import copy
import hashlib
import json
import numpy as np
from scipy.linalg import cholesky, solve_triangular

FIELDS=['x_plus_mas','y_plus_mas','x_minus_mas','y_minus_mas','ln_flux_ratio','delay_seconds']

@dataclass
class Dataset:
    raw: dict
    epochs: np.ndarray
    observations: np.ndarray
    covariance: np.ndarray
    chol: np.ndarray
    sha256: str

    @property
    def n(self): return len(self.epochs)
    @property
    def geometry(self): return self.raw['geometry']
    @property
    def analysis(self): return self.raw.get('analysis',{})
    @property
    def t0(self): return self.raw['reference_epoch_yr']
    def whiten(self,values): return solve_triangular(self.chol,np.asarray(values),lower=True)
    def summary(self):
        return {'rows':self.n,'fields':FIELDS,'reference_epoch_yr':self.t0,
                'span_observer_years':float(np.ptp(self.epochs)),
                'provenance':self.raw['provenance'],'sha256':self.sha256,
                'covariance':'full supplied positive-definite covariance; row order preserved'}

def load_dataset(source):
    if isinstance(source,Dataset): return source
    raw=json.loads(Path(source).read_text(encoding='utf-8-sig')) if isinstance(source,(str,Path)) else copy.deepcopy(source)
    if not isinstance(raw,dict) or raw.get('schema_version')!=1: raise ValueError('schema_version must be 1')
    geometry=raw.get('geometry',{})
    for key in ['D_l_pc','D_s_pc','D_ls_pc']:
        v=geometry.get(key)
        if not isinstance(v,(int,float)) or not np.isfinite(v) or v<=0:raise ValueError(f'geometry.{key} must be positive and finite')
    for key in ['z_l','z_s']:
        v=geometry.get(key)
        if not isinstance(v,(int,float)) or not np.isfinite(v) or v<0:raise ValueError(f'geometry.{key} must be nonnegative and finite')
    if geometry['z_s']<=geometry['z_l'] and not geometry['z_l']==geometry['z_s']==0:
        raise ValueError('source redshift must exceed lens redshift, or both may be zero for supplied local distances')
    epochs=np.asarray(raw.get('epochs_yr',[]),dtype=float)
    obs=np.asarray(raw.get('observations',[]),dtype=float)
    if epochs.ndim!=1 or not 3<=len(epochs)<=80:raise ValueError('provide 3 to 80 epoch/source rows')
    if obs.shape!=(len(epochs),6):raise ValueError('observations must have six fields per row, in the documented order')
    if not np.all(np.isfinite(epochs)) or not np.all(np.isfinite(obs)):raise ValueError('measurements and epochs must be finite')
    if np.any(np.linalg.norm(obs[:,:2]-obs[:,2:4],axis=1)==0):raise ValueError('primary images must be separated')
    t0=raw.get('reference_epoch_yr')
    if not isinstance(t0,(float,int)) or not np.isfinite(t0):raise ValueError('reference_epoch_yr must be finite')
    covariance=np.asarray(raw.get('covariance',[]),dtype=float)
    n=6*len(epochs)
    if covariance.shape!=(n,n) or not np.all(np.isfinite(covariance)):raise ValueError('covariance must be a finite 6N by 6N matrix')
    # Compare symmetry in correlation units, avoiding mixed mas/second scales.
    if np.any(np.diag(covariance)<=0):raise ValueError('covariance variances must be positive')
    scale=np.sqrt(np.diag(covariance));corr=covariance/np.outer(scale,scale)
    if np.max(abs(corr-corr.T))>1e-9:raise ValueError('covariance must be symmetric')
    covariance=(covariance+covariance.T)/2
    try:chol=cholesky(covariance,lower=True)
    except np.linalg.LinAlgError as error:raise ValueError('covariance must be positive definite') from error
    provenance=raw.get('provenance',{})
    if provenance.get('kind') not in ['synthetic','observed']:raise ValueError('provenance.kind must explicitly be synthetic or observed')
    if raw.get('fields',FIELDS)!=FIELDS:raise ValueError('field order differs from the six-field contract')
    analysis=raw.setdefault('analysis',{})
    if 'lens_parallax_factors' in analysis:
        pf=np.asarray(analysis['lens_parallax_factors'],dtype=float)
        if pf.shape!=(len(epochs),2) or not np.all(np.isfinite(pf)):
            raise ValueError('lens_parallax_factors must contain N finite two-vectors')
    if analysis.get('source_policy','free_per_epoch')!='free_per_epoch':raise ValueError('this version supports source_policy=free_per_epoch')
    for key in ['r0_pc','reference_mass_msun']:
        if key in analysis and (not np.isfinite(analysis[key]) or analysis[key]<=0):raise ValueError(f'{key} must be positive and finite')
    if 'mass_bounds_msun' in analysis:
        mb=np.asarray(analysis['mass_bounds_msun'],dtype=float)
        if mb.shape!=(2,) or not np.all(np.isfinite(mb)) or not 0<mb[0]<mb[1]:
            raise ValueError('mass_bounds_msun must be two increasing positive finite masses')
    digest=hashlib.sha256(json.dumps(raw,sort_keys=True,separators=(',',':'),allow_nan=False).encode()).hexdigest()
    return Dataset(raw,epochs,obs,covariance,chol,digest)
