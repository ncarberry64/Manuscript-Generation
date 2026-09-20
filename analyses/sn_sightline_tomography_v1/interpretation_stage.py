"""Frozen-template cross-validation, exact covariance propagation, mocks and W rank."""
import os
os.environ['OPENBLAS_NUM_THREADS']='1'; os.environ['OMP_NUM_THREADS']='1'
from pathlib import Path
import json
import numpy as np
import pandas as pd
from scipy.linalg import cholesky,solve_triangular,svd,eigh
from scipy.sparse import load_npz
from scipy.stats import chi2,norm
import frozen_analysis as fa
ROOT=Path(__file__).resolve().parent; OUT=ROOT/'interpretation'; OUT.mkdir(exist_ok=True)
A=-0.04321623435997208; SEED=20260920; NM=2000
CHECKS={}; RESULTS=[]; MOCKS={}; OPERATORS={}; FOLDS=[]
def check(name,value,tol=1e-8):
 CHECKS[name]=dict(error=float(value),tolerance=tol,passed=bool(value<tol)); assert value<tol,(name,value)
def white(L,x): return solve_triangular(L,x,lower=True,check_finite=False)
def basis(x):
 if x.shape[1]==0: return np.empty((len(x),0)),np.array([])
 u,s,v=svd(x,full_matrices=False,check_finite=False); keep=s>s[0]*1e-10
 return u[:,keep],s
def projector(C,N):
 C=(C+C.T)/2
 try:
  L=cholesky(C,lower=True); B=white(L,np.eye(len(C)))
  if np.min(np.diag(L))<1e-6*np.max(np.diag(L)): raise np.linalg.LinAlgError('near-singular')
 except np.linalg.LinAlgError:
  e,u=eigh(C); assert e[0]>-1e-8*e[-1]
  keep=e>1e-10*e[-1]; inv=np.zeros(len(e)); inv[keep]=1/np.sqrt(e[keep]); B=inv[:,None]*u.T; L=u*np.sqrt(np.maximum(e,0))
 q,s=basis(B@N); B-=q@(q.T@B)
 return B,L,len(q.T)
def evaluate(label,y,C,f,N,mocks=True):
 B,L,rank=projector(C,N); r=B@y; t=B@f; info=float(t@t)
 q0=float(r@r); q1=float((r-A*t)@(r-A*t)); delta=q1-q0; dof=len(y)-rank
 ce=eigh((C+C.T)/2,eigvals_only=True); keep=ce>1e-10*ce[-1]; covariance_rank=int(keep.sum()); dof=covariance_rank-rank
 lognorm=np.log(ce[keep]).sum()+covariance_rank*np.log(2*np.pi)
 rec=dict(test=label,rows=len(y),profiled_nuisance_rank=rank,chi2_null=q0,chi2_frozen=q1,dof=dof,delta_chi2=delta,logL_null=float(-.5*(q0+lognorm)),logL_frozen=float(-.5*(q1+lognorm)),likelihood_ratio=float(np.exp(np.clip(-delta/2,-700,700))),p_null_gof=float(chi2.sf(q0,dof)),p_frozen_gof=float(chi2.sf(q1,dof)),fitted_amplitude=float(t@r/info) if info>1e-20 else None,fitted_amplitude_sigma=float(info**-.5) if info>1e-20 else None,template_information=info)
 if mocks:
  rng=np.random.default_rng(SEED); z=rng.standard_normal((len(y),NM)); mockr=B@L@z
  md=A*A*info-2*A*(t@mockr); rec['mock_p_frozen_direction']=float((1+np.sum(md<=delta))/(NM+1)); rec['analytic_p_frozen_direction']=float(norm.sf(np.sign(A)*(t@r)/np.sqrt(info))) if info>0 else None
  MOCKS[label]=md
 RESULTS.append(rec); print(label,'n',len(y),'delta',round(delta,5),flush=True)
 rec['covariance_rank']=covariance_rank
 return rec,B,L
def raw_features(df):
 cols=['x1','c','MWEBV','HOST_LOGMASS']+[f'shell{s}_density_proxy' for s in range(8)]+['lc_delta_c','lc_delta_x1','lc_delta_t0']
 # External velocity normalization is FIXED. Never regress its correction
 # difference on the SN residuals; use coefficient one in mean sensitivities.
 x=df.reindex(columns=cols).astype(float).replace([np.inf,-np.inf],np.nan)
 for c in ['HOST_LOGMASS','MWEBV']: x.loc[x[c]<0,c]=np.nan
 # No amplitude departure is used to remove a putative coherent signal.
 return x
def scale_features(train,test):
 means=train.mean().fillna(0); scales=train.std(ddof=0).replace(0,1).fillna(1)
 def trans(x):
  return np.c_[((x-means)/scales).fillna(0).to_numpy(),x.isna().to_numpy().astype(float)]
 return trans(train),trans(test),dict(names=list(train.columns)+[c+'_missing' for c in train.columns],mean=means.to_list(),scale=scales.to_list())
def nuisance_estimator(C,N,X):
 L=cholesky(C,lower=True); Q,_=basis(white(L,N)); B=white(L,np.eye(len(C))); B-=Q@(Q.T@B); ew=B@X
 u,s,vt=svd(ew,full_matrices=False); keep=s>s[0]*1e-10
 H=(vt[keep].T/s[keep])@u[:,keep].T@B
 check('nuisance annihilation '+str(len(C))+'_'+str(len(FOLDS)),np.max(abs(H@N)),1e-7)
 return H,dict(rank=int(keep.sum()),singular=s.tolist(),condition=float(s[0]/s[keep][-1]) if keep.any() else None)
def train_predict(label,tr,te,Ctr,Cte,Ccross,Ftr,Fte,edges):
 Xt,Xv,scale=scale_features(Ftr,Fte); Nt=fa.nuisance(tr.zHD,fa.TB); H,rank=nuisance_estimator(Ctr,Nt,Xt)
 T=Xv@H; yy=te.common_mode_residual.to_numpy()-T@tr.common_mode_residual.to_numpy(); ff=te.frozen_template.to_numpy()-T@tr.frozen_template.to_numpy()
 V=Cte+T@Ctr@T.T-Ccross@T.T-T@Ccross.T; V=(V+V.T)/2
 rec,B,L=evaluate(label,yy,V,ff,fa.nuisance(te.zHD,edges)); rec.update(training_rows=len(tr),environment_rank=rank['rank'],training_SN_blocks=tr.sn_block.nunique())
 rec['heldout_error_semantics']='joint Gaussian prediction-error covariance propagated from released total C; fixed external/source features; no feature-error covariance supplied'
 np.savez_compressed(OUT/(label+'_OPERATORS.npz'),test_uid=te.uid.to_numpy(dtype=str),train_uid=tr.uid.to_numpy(dtype=str),training_estimator=H,train_X=Xt,test_X=Xv,T=T,prediction_error_covariance=V,transformed_template=ff,prediction_error=yy,B=B)
 fa.savejson(OUT/(label+'_DESIGN.json'),dict(scale=scale,rank=rank,beta_null=H@tr.common_mode_residual.to_numpy(),beta_frozen=H@(tr.common_mode_residual.to_numpy()-A*tr.frozen_template.to_numpy())))
 FOLDS.extend(dict(test=label,role='test',uid=u) for u in te.uid); FOLDS.extend(dict(test=label,role='train',uid=u) for u in tr.uid)
 return dict(y=yy,f=ff,C=V,B=B,L=L,T=T,frame=te,record=rec)
def within_fold(label,df,C,F,testmask):
 teidx=np.flatnonzero(testmask); heldblocks=set(df.iloc[teidx].sn_block); tridx=np.flatnonzero(~df.sn_block.isin(heldblocks).to_numpy())
 if len(teidx)<4 or len(tridx)<5: return None
 return train_predict(label,df.iloc[tridx],df.iloc[teidx],C[np.ix_(tridx,tridx)],C[np.ix_(teidx,teidx)],C[np.ix_(teidx,tridx)],F.iloc[tridx],F.iloc[teidx],fa.TB)
def permutations(label,obj):
 df=obj['frame'].reset_index(drop=True); B=obj['B']; y=obj['y']; f=df.frozen_template.to_numpy(); rng=np.random.default_rng(SEED+1)
 transfer=fa.distances()[3](df.zHD,'local'); angular=f/transfer
 # Block permutation retains each duplicate group's complete survey/z/sector signature.
 signatures={}; groups=list(df.groupby('sn_block',sort=True).indices.values())
 zbin=np.searchsorted(fa.TB,df.zHD,side='right'); sector=(df.RA//90).astype(int).to_numpy(); survey=df.IDSURVEY.to_numpy()
 for ids in groups:
  signature=tuple(sorted((int(survey[i]),int(zbin[i]),int(sector[i])) for i in ids))
  signatures.setdefault(signature,[]).append(np.asarray(sorted(ids,key=lambda i:(survey[i],zbin[i],sector[i],i))))
 fp=np.repeat(angular[:,None],NM,axis=1); moved=0
 for blocks in signatures.values():
  if len(blocks)<2: continue
  moved+=sum(len(b) for b in blocks)
  for m in range(NM):
   dest=rng.permutation(len(blocks))
   for j,k in enumerate(dest): fp[blocks[j],m]=angular[blocks[k]]
 fp*=transfer[:,None]  # exact observed redshifts and frozen transfer law retained
 # Training response stays frozen; randomization is held-out conditional-design only.
 fp-= (f-obj['f'])[:,None]; tp=B@fp; r=B@y
 dp=A*A*np.sum(tp*tp,axis=0)-2*A*(r@tp)
 observed=obj['record']['delta_chi2']; tail=float((1+np.sum(dp<=observed))/(NM+1))
 z=rng.standard_normal((len(df),200)); rm=B@obj['L']@z
 null_perm=A*A*np.sum(tp*tp,axis=0)[:,None]-2*A*(tp.T@rm)
 tt=B@obj['f']; null_obs=A*A*(tt@tt)-2*A*(tt@rm)
 nullp=(1+np.sum(null_perm<=null_obs[None,:],axis=0))/(NM+1)
 np.savez_compressed(OUT/(label+'_PERMUTATIONS.npz'),delta_chi2=dp,null_calibration_p=nullp,conditional_template=fp)
 return dict(test=label,draws=NM,movable_rows=moved,total_rows=len(df),conditional_permutation_p=tail,mock_calibration_draws=len(nullp),fraction_null_p_below_05=float(np.mean(nullp<.05)),exchangeability='not guaranteed; use Gaussian mock p for inference')
def matched(label,obj,W):
 df=obj['frame'].reset_index(drop=True); n=len(df); v=fa.vectors(df.RA,df.DEC); angle=np.rad2deg(np.arccos(np.clip(v@v.T,-1,1))); dz=abs(df.zHD.to_numpy()[:,None]-df.zHD.to_numpy()[None,:]); same=df.IDSURVEY.to_numpy()[:,None]==df.IDSURVEY.to_numpy()[None,:]
 gram=(W@W.T).toarray(); normw=np.sqrt(np.diag(gram)); overlap=gram/np.maximum(normw[:,None]*normw[None,:],1e-100)
 types={'same_z_different_direction':(dz<=.005)&(angle>=60),'same_direction_different_z':(angle<=5)&(dz>=.03),'shared_foreground_different_z':(overlap>=.5)&(dz>=.03)}
 records=[]
 for kind,mask in types.items():
  aa,bb=np.where(np.triu(mask,1)); order=sorted(range(len(aa)),key=lambda k:(not same[aa[k],bb[k]],float(dz[aa[k],bb[k]]) if kind.startswith('same_z') else float(angle[aa[k],bb[k]]),aa[k],bb[k])); used=set(); pairs=[]
  for k in order:
   i,j=int(aa[k]),int(bb[k]); si,sj=df.sn_block.iloc[i],df.sn_block.iloc[j]
   if si==sj or si in used or sj in used: continue
   used|={si,sj}; pairs.append((i,j)); records.append(dict(category=kind,uid_i=df.uid.iloc[i],uid_j=df.uid.iloc[j],delta_z=dz[i,j],angle_deg=angle[i,j],W_cosine=overlap[i,j],same_survey=bool(same[i,j])))
  if not pairs: continue
  P=np.zeros((len(pairs),n))
  for k,(i,j) in enumerate(pairs): P[k,i]=1; P[k,j]=-1
  # Profile pair response to ORIGINAL test-bin nuisance, preserving off-diagonal C.
  evaluate(label+'_'+kind,P@obj['y'],P@obj['C']@P.T,P@obj['f'],P@obj.get('N',fa.nuisance(df.zHD,fa.TB)))
  np.savez_compressed(OUT/(label+'_'+kind+'_PAIR_OPERATOR.npz'),P=P,C=P@obj['C']@P.T)
 pd.DataFrame(records).to_csv(OUT/(label+'_MATCHED_PAIRS.csv'),index=False)
def tomography(label,obj,W,cells):
 B=obj['B']; r=B@obj['y']; rf=r-A*(B@obj['f']); n=len(r); rows=[]; scores=[]; rng=np.random.default_rng(SEED+2); zm=B@obj['L']@rng.standard_normal((n,NM))
 rawgram=(W@W.T).toarray(); eig=eigh(rawgram,eigvals_only=True); s=np.sqrt(np.maximum(eig,0)); tol=max(W.shape)*np.finfo(float).eps*s[-1]
 # Raw W rank uses an explicit singular-value tolerance; near duplicates retained.
 rank=int(np.sum(s>max(tol,1e-7*s[-1]))); fa.savejson(OUT/(label+'_W_RANK.json'),dict(rows=W.shape[0],columns=W.shape[1],rank_at_relative_singular_tolerance_1e_7=rank,singular=np.sort(s)[::-1],note='Gram eigenvalue floor limits rank near 1e-8; duplicates and common near-observer paths create degeneracies'))
 for shell in sorted(cells.shell.unique()):
  ws=W[:,cells.shell.to_numpy()==shell]; g=(ws@ws.T).toarray(); scale=np.sqrt(np.diag(g)); k=g/np.maximum(scale[:,None]*scale[None,:],1e-100)
  # A physical cell perturbation acts on original distances BEFORE cross-fitting.
  # Apply exactly the same prediction-error transformation as for y and f.
  if 'R' in obj: k=obj['R']@k@obj['R'].T
  elif 'Wtrain' in obj:
   wt=obj['Wtrain'][:,cells.shell.to_numpy()==shell]; st=np.sqrt(np.asarray(wt.multiply(wt).sum(axis=1)).ravel())
   response=ws.multiply(1/np.maximum(scale,1e-100)[:,None]).toarray()-obj['T']@wt.multiply(1/np.maximum(st,1e-100)[:,None]).toarray(); k=response@response.T
  # Score for a positive covariance component, not fitted cell amplitudes.
  K=B@k@B.T; K=(K+K.T)/2; trace=float(np.trace(K)); sd=float(np.sqrt(2*np.sum(K*K)))
  if sd<=1e-20: continue
  observed=float((r@K@r-trace)/sd); frozen=float((rf@K@rf-trace)/sd); ms=(np.sum(zm*(K@zm),axis=0)-trace)/sd; scores.append(ms)
  rows.append(dict(test=label,shell=int(shell),z_min=fa.EDGES[shell],z_max=fa.EDGES[shell+1],cells=ws.shape[1],score_null_mean=observed,score_after_frozen_mean=frozen,p_mock=float((1+np.sum(ms>=observed))/(NM+1)),kernel_effective_rank=trace*trace/float(np.sum(K*K)),distinct_crossing_SN_rows=int((scale>0).sum())))
  np.savez_compressed(OUT/(label+f'_SHELL{shell}_MOCK.npz'),score=ms)
 maxmock=np.max(scores,axis=0)
 for row in rows: row['p_max_over_shells']=float((1+np.sum(maxmock>=row['score_null_mean']))/(NM+1))
 pd.DataFrame(rows).to_csv(OUT/(label+'_TOMOGRAPHY_SCORES.csv'),index=False)
 return rows
def main():
 sight=pd.read_csv(ROOT/'sightlines/SN_SIGHTLINES.csv.gz',dtype={'CID':str}); lc=pd.read_csv(ROOT/'lightcurves/LIGHTCURVE_RECONSTRUCTION.csv')
 for c in ['delta_c','delta_x1','delta_t0']:
  valid=lc[c].where(lc.status=='ok'); sight['lc_'+c]=sight.uid.map(pd.Series(valid.to_numpy(),index=lc.uid))
 sight['lc_status']=sight.uid.map(lc.set_index('uid').status)
 cov=np.load(ROOT/'sightlines/RELEASED_TOTAL_COVARIANCES.npz'); df=sight.loc[sight.primary_selection].copy(); _,_,_,transfer,_=fa.distances()
 df['frozen_template']=(fa.vectors(df.RA,df.DEC)@fa.vectors(*fa.AXIS)[0])*transfer(df.zHD,'local')
 df.to_csv(OUT/'ANALYSIS_ROWS_WITH_SOURCE_FEATURES.csv.gz',index=False)
 p=df.loc[df.dataset=='Pantheon'].reset_index(drop=True); d=df.loc[df.dataset=='DES'].reset_index(drop=True)
 cp=cov['Pantheon'][np.ix_(p.raw_row,p.raw_row)]; cd=cov['DES'][np.ix_(d.raw_row,d.raw_row)]; fp,fd=raw_features(p),raw_features(d)
 fp.assign(uid=p.uid).to_csv(OUT/'PANTHEON_NUISANCE_FEATURES.csv',index=False); fd.assign(uid=d.uid).to_csv(OUT/'DES_NUISANCE_FEATURES.csv',index=False)
 pi=np.flatnonzero(p.zHD>=.03); tr=np.flatnonzero(p.zHD<.03)
 evaluate('P_original_heldout_release_only',p.common_mode_residual.to_numpy()[pi],cp[np.ix_(pi,pi)],p.frozen_template.to_numpy()[pi],fa.nuisance(p.zHD.to_numpy()[pi],fa.PB))
 evaluate('DES_original_heldout_release_only',d.common_mode_residual.to_numpy(),cd,d.frozen_template.to_numpy(),fa.nuisance(d.zHD,fa.DB))
 baseline=train_predict('P_original_training_to_heldout',p.iloc[tr],p.iloc[pi],cp[np.ix_(tr,tr)],cp[np.ix_(pi,pi)],cp[np.ix_(pi,tr)],fp.iloc[tr],fp.iloc[pi],fa.PB)
 des=train_predict('DES_Pantheon_transfer',p,d,cp,cd,np.zeros((len(d),len(p))),fp,fd,fa.DB)
 # Zero cross-survey block is an explicit unavailable-cross-covariance assumption.
 des['record']['cross_release_covariance']='unreleased; zero assumed for this transfer, never combine surveys into joint significance'
 pe=p.copy(); pe['common_mode_residual']=p.common_mode_residual+p.pv_replacement_delta_mag.fillna(0)
 extbaseline=train_predict('P_training_to_heldout_fixed_external_PV_sensitivity',pe.iloc[tr],pe.iloc[pi],cp[np.ix_(tr,tr)],cp[np.ix_(pi,pi)],cp[np.ix_(pi,tr)],fp.iloc[tr],fp.iloc[pi],fa.PB)
 sky=[]; survey=[]
 for sector in range(4):
  obj=within_fold('P_leave_sky_'+str(sector),p,cp,fp,(p.RA//90).to_numpy()==sector)
  if obj: sky.append(obj)
 for sid in sorted(p.IDSURVEY.unique()):
  obj=within_fold('P_leave_survey_'+str(sid),p,cp,fp,(p.IDSURVEY==sid).to_numpy())
  if obj: survey.append(obj)
 # Assemble out-of-fold errors with their full covariance (folds are correlated).
 n=len(p); R=np.eye(n); covered=np.zeros(n,bool)
 for obj in sky:
  test=p.uid.get_indexer(obj['frame'].uid) if hasattr(p.uid,'get_indexer') else pd.Index(p.uid).get_indexer(obj['frame'].uid)
  held=set(obj['frame'].sn_block); train=np.flatnonzero(~p.sn_block.isin(held).to_numpy()); R[np.ix_(test,train)]-=obj['T']; covered[test]=True
 assert covered.all()
 V=R@cp@R.T; yy=R@p.common_mode_residual.to_numpy(); ff=R@p.frozen_template.to_numpy(); N=R@fa.nuisance(p.zHD,fa.TB)
 rec,B,L=evaluate('P_sky_crossfit_joint',yy,V,ff,N); joint=dict(y=yy,f=ff,C=V,B=B,L=L,frame=p,record=rec,N=N,R=R)
 evaluate('P_sky_crossfit_fixed_external_PV_sensitivity',R@pe.common_mode_residual.to_numpy(),V,ff,N)
 np.savez_compressed(OUT/'P_SKY_CROSSFIT_JOINT_OPERATORS.npz',R=R,C=V,B=B,y=yy,f=ff,uid=p.uid.to_numpy(dtype=str))
 # Conditional amplitude/rank diagnostics after in-sample nuisance projection.
 X,_,sc=scale_features(fp,fp); N0=fa.nuisance(p.zHD,fa.TB); Bt,Lt,rr=projector(cp,np.c_[N0,X]); f=p.frozen_template.to_numpy(); r=p.common_mode_residual.to_numpy(); tt=Bt@f; nn=projector(cp,N0)[0]@f
 info=float(tt@tt); fdiag=dict(label='FITTED one scalar amplitude, not predictive',amplitude=float(tt@(Bt@r)/info),sigma=float(info**-.5),frozen_amplitude=A,information_before=float(nn@nn),information_after=info,sigma_inflation=float(np.sqrt((nn@nn)/info)),nuisance_rank=rr,remaining_template_fraction=float(info/(nn@nn)))
 fa.savejson(OUT/'FITTED_AMPLITUDE_AND_RANK.json',fdiag)
 # Project KNOWN velocity corrections onto the saved optical estimator. No velocity is fitted.
 old=np.load(ROOT/'FROZEN_GLS_OPERATORS.npz'); pfull=sight.loc[sight.dataset=='Pantheon'].set_index('raw_row'); otr=old['training_raw_rows']; oh=old['training_amplitude_estimator'].ravel()
 check('frozen training row identity',np.max(abs(p.raw_row.to_numpy()[tr]-otr)))
 pvproj=[]
 for name,col in [('released','release_pv_delta_mag'),('independent_NN','independent_pv_delta_mag')]:
  value=pfull.loc[otr,col].to_numpy(); pvproj.append(dict(correction=name,rows=len(value),available=int(np.isfinite(value).sum()),saved_estimator_projection_mag=float(oh@value) if np.isfinite(value).all() else None))
 fa.savejson(OUT/'FIXED_PV_CORRECTION_PROJECTIONS.json',pvproj)
 # Fixed external mean replacement: no same-residual velocity inference.
 ids=np.flatnonzero(p.independent_pv_available.to_numpy()); evaluate('P_external_PV_fixed_mean_sensitivity',p.residual_fixed_external_pv.to_numpy()[ids],cp[np.ix_(ids,ids)],p.frozen_template.to_numpy()[ids],fa.nuisance(p.zHD.to_numpy()[ids],fa.TB))
 evaluate('P_external_PV_covered_release_control',p.common_mode_residual.to_numpy()[ids],cp[np.ix_(ids,ids)],p.frozen_template.to_numpy()[ids],fa.nuisance(p.zHD.to_numpy()[ids],fa.TB))
 median=float(p.path_density_proxy.median())
 for label,mask in [('low',p.path_density_proxy<=median),('high',p.path_density_proxy>median)]:
  ix=np.flatnonzero(mask); evaluate('P_crossfit_density_'+label,yy[ix],V[np.ix_(ix,ix)],ff[ix],N[ix])
 fa.savejson(OUT/'DENSITY_SPLIT.json',dict(median_galaxy_conditioned_matter_density=median))
 W=load_npz(ROOT/'sightlines/W_LOCAL_MAP_PATH_LENGTH.npz'); cells=pd.read_csv(ROOT/'sightlines/LOCAL_MAP_CELLS.csv'); wp=W[p.matrix_row.to_numpy()]; wd=W[d.matrix_row.to_numpy()]
 des['Wtrain']=wp
 matched('P_crossfit',joint,wp); matched('DES_transfer',des,wd)
 perms=[permutations('P_original_training_to_heldout',baseline),permutations('DES_Pantheon_transfer',des)]
 fa.savejson(OUT/'PERMUTATION_SUMMARY.json',perms)
 all_tomo=tomography('P_crossfit',joint,wp,cells)+tomography('DES_transfer',des,wd,cells)
 for obj in sky+survey:
  te=pd.Index(p.uid).get_indexer(obj['frame'].uid); held=set(obj['frame'].sn_block); trf=np.flatnonzero(~p.sn_block.isin(held).to_numpy()); obj['Wtrain']=wp[trf]
  all_tomo.extend(tomography(obj['record']['test'],obj,wp[te],cells))
 pd.DataFrame(all_tomo).to_csv(OUT/'TOMOGRAPHY_ALL_SUBSETS.csv',index=False)
 # Sanity checks: amplitude injection and covariance row-order invariance.
 check('one-amplitude injected signal recovery',abs(float(tt@(Bt@(A*f))/info)-A))
 order=np.random.default_rng(SEED).permutation(len(pi)); yc=p.common_mode_residual.to_numpy()[pi]; cc=cp[np.ix_(pi,pi)]; nc=fa.nuisance(p.zHD.to_numpy()[pi],fa.PB)
 BB=projector(cc,nc)[0]; BR=projector(cc[np.ix_(order,order)],nc[order])[0]; check('joint row/covariance permutation invariance',abs(float(np.sum((BB@yc)**2)-np.sum((BR@yc[order])**2))),1e-7)
 fa.savejson(OUT/'CHECKS.json',CHECKS); fa.savejson(OUT/'LIKELIHOODS.json',RESULTS); pd.DataFrame(RESULTS).to_csv(OUT/'LIKELIHOODS.csv',index=False); pd.DataFrame(FOLDS).to_csv(OUT/'FOLD_MEMBERSHIP.csv.gz',index=False); np.savez_compressed(OUT/'GAUSSIAN_NULL_MOCKS.npz',**MOCKS)
 print('Complete',len(RESULTS),'likelihood comparisons',flush=True)
if __name__=='__main__': main()
