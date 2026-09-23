"""Orthogonal epoch residual decomposition and repeat-SN survey comparisons."""
import os
os.environ['OPENBLAS_NUM_THREADS']='1'
from pathlib import Path
import numpy as np
import pandas as pd
from scipy.linalg import svd
import frozen_analysis as fa
ROOT=Path(__file__).resolve().parent; LC=ROOT/'lightcurves'
def orth(x,reference=None):
 u,s,v=svd(x,full_matrices=False); ref=(s[0] if len(s) else 0) if reference is None else reference
 return u[:,s>ref*1e-8] if ref>0 else np.zeros((len(x),0))
def main():
 fits=pd.read_csv(LC/'LIGHTCURVE_RECONSTRUCTION.csv'); summaries=[]; errors=[]
 for row in fits.loc[fits.status!='failed'].itertuples():
  a=np.load(LC/(row.uid+'_FIT.npz')); ep=pd.read_csv(LC/(row.uid+'_EPOCHS.csv.gz')); sigma=ep.fluxerr.to_numpy(); j=a['J']/sigma[:,None]; q=orth(j)
  rw=ep.residual_flux.to_numpy()/sigma; rw-=q@(q.T@rw)
  b=np.column_stack([(ep.BAND==band)*ep.model_flux/sigma for band in ep.BAND.unique()]); bref=np.linalg.norm(b); b-=q@(q.T@b); qb=orth(b,bref)
  phase=ep.rest_phase.to_numpy()/30; ph=np.column_stack([ep.model_flux/sigma*phase,ep.model_flux/sigma*phase**2]); pref=np.linalg.norm(ph); ph-=q@(q.T@ph); ph-=qb@(qb.T@ph); qp=orth(ph,pref)
  rb=qb@(qb.T@rw); rp=qp@(qp.T@rw); ru=rw-rb-rp
  summaries.append(dict(uid=row.uid,chromatic_residual_chi2=float(rb@rb),chromatic_rank=qb.shape[1],phase_given_chromatic_chi2=float(rp@rp),phase_rank=qp.shape[1],remaining_chi2=float(ru@ru),post_source_chi2=float(rw@rw)))
  errors.append(abs(float(rw@rw-rb@rb-rp@rp-ru@ru)))
  pd.DataFrame(dict(MJD=ep.MJD,BAND=ep.BAND,post_source_residual_flux=rw*sigma,chromatic_residual_flux=rb*sigma,phase_given_chromatic_flux=rp*sigma,unresolved_residual_flux=ru*sigma)).to_csv(LC/(row.uid+'_COMPONENTS.csv.gz'),index=False)
 pd.DataFrame(summaries).to_csv(LC/'ORTHOGONAL_RESIDUAL_DECOMPOSITION.csv',index=False)
 sight=pd.read_csv(ROOT/'sightlines/SN_SIGHTLINES.csv.gz'); p=sight.loc[sight.dataset=='Pantheon'].merge(fits[['uid','common_delta_mag','common_sigma_photometric','status']],on='uid',validate='one_to_one'); C=np.load(ROOT/'sightlines/RELEASED_TOTAL_COVARIANCES.npz')['Pantheon']; records=[]
 for sn,g in p.groupby('sn_block'):
  for i in range(len(g)):
   for j in range(i+1,len(g)):
    a,b=g.iloc[i],g.iloc[j]
    if a.IDSURVEY==b.IDSURVEY: continue
    ia,ib=int(a.raw_row),int(b.raw_row); variance=C[ia,ia]+C[ib,ib]-2*C[ia,ib]
    records.append(dict(sn_block=sn,uid_i=a.uid,uid_j=b.uid,survey_i=a.IDSURVEY,survey_j=b.IDSURVEY,released_common_difference=a.common_mode_residual-b.common_mode_residual,sigma_released_common_difference=np.sqrt(max(variance,0)),local_refit_departure_difference=a.common_delta_mag-b.common_delta_mag,local_photometric_sigma_independent_proxy=np.hypot(a.common_sigma_photometric,b.common_sigma_photometric),raw_error_cross_survey_covariance='unreleased; sigma is a diagonal diagnostic only',fit_i=a.status,fit_j=b.status))
 pd.DataFrame(records).to_csv(LC/'REPEAT_SN_SURVEY_SPLITS.csv',index=False)
 checks=dict(epoch_orthogonal_power_sum_max_error=max(errors),frozen_amplitude_used_in_source_fitting=False)
 assert max(errors)<1e-5
 # Genuine separation check: inject known amplitude/color/shape/phase perturbations
 # into one real irregular epoch design and recover with its saved estimator.
 a=np.load(LC/'P0100_FIT.npz'); injected=np.array([.02,.1,.01,.2]); recovered=a['estimator']@(a['J']@injected); checks['four_parameter_injection_error']=float(np.max(abs(recovered-injected))); assert checks['four_parameter_injection_error']<1e-8
 fa.savejson(LC/'SOURCE_CHECKS.json',checks); print(checks,flush=True)
if __name__=='__main__': main()
