"""Delete each observed MJD epoch, grouping simultaneous bands, in local GLS."""
import os
os.environ['OPENBLAS_NUM_THREADS']='1'
from pathlib import Path
import numpy as np
import pandas as pd
import frozen_analysis as fa
ROOT=Path(__file__).resolve().parent; LC=ROOT/'lightcurves'
def main():
 fits=pd.read_csv(LC/'LIGHTCURVE_RECONSTRUCTION.csv'); records=[]; summary=[]; errors=[]; factor=2.5/np.log(10)
 for row in fits.loc[fits.status!='failed'].itertuples():
  a=np.load(LC/(row.uid+'_FIT.npz')); ep=pd.read_csv(LC/(row.uid+'_EPOCHS.csv.gz')); J=a['J']; H=a['estimator']; residual=ep.residual_flux.to_numpy(); sig=ep.fluxerr.to_numpy(); dp=H@residual; rlin=residual-J@dp; vals=[]
  for mjd,ids in ep.groupby('MJD',sort=True).indices.items():
   ids=np.asarray(ids); D=np.eye(len(ids))-J[ids]@H[:,ids]; rec=dict(uid=row.uid,MJD=mjd,deleted_observations=len(ids),status='rank_failed')
   if len(ep)-len(ids)>=4 and np.linalg.cond(D)<1e10:
    M=np.linalg.solve(D.T,H[:,ids].T).T; shift=-M@rlin[ids]; dh=M[0]@(J[ids]@H); dh[ids]-=M[0]; dm=-factor*shift[0]; se=factor*np.sqrt(np.sum((dh*sig)**2))
    rec.update(status='ok',delta_common_mag=dm,sigma_difference=se,delta_shape=shift[1],delta_color=shift[2],delta_phase=shift[3]); vals.append(abs(dm))
    if len(errors)<25 and row.raw_row%100==0:
     keep=np.ones(len(ep),bool); keep[ids]=False; jw=J[keep]/sig[keep,None]; yr=residual[keep]/sig[keep]; direct=np.linalg.lstsq(jw,yr,rcond=1e-12)[0]-dp; errors.append(float(np.max(abs(direct-shift))))
   records.append(rec)
  summary.append(dict(uid=row.uid,epochs_tested=ep.MJD.nunique(),identifiable_epoch_jackknifes=len(vals),single_epoch_jackknife_max=max(vals) if vals else np.nan))
 pd.DataFrame(records).to_csv(LC/'SINGLE_EPOCH_JACKKNIFE.csv.gz',index=False); pd.DataFrame(summary).to_csv(LC/'SINGLE_EPOCH_JACKKNIFE_SUMMARY.csv',index=False)
 assert max(errors)<1e-7
 fa.savejson(LC/'SINGLE_EPOCH_CHECKS.json',dict(direct_retained_epoch_refit_comparisons=len(errors),max_parameter_difference=max(errors),groups=len(records),successful=int(sum(r['status']=='ok' for r in records)),method='exact local-linear delete-MJD GLS; simultaneous bands grouped; not nonlinear refits'))
 print('Single-epoch groups',len(records),'max direct-refit error',max(errors),flush=True)
if __name__=='__main__': main()
