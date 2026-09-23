"""Independent galaxy-field sightlines. No SN residual enters environmental construction."""
import os
os.environ['OPENBLAS_NUM_THREADS']='1'
from pathlib import Path
import json,hashlib
import numpy as np
import pandas as pd
from scipy.ndimage import map_coordinates
from scipy.sparse import coo_matrix,save_npz
from astropy.coordinates import SkyCoord
import astropy.units as u
import frozen_analysis as fa
ROOT=Path(__file__).resolve().parent; OUT=ROOT/'sightlines'; OUT.mkdir(exist_ok=True)
SHELLS=fa.EDGES; C100=2997.92458
def sample(field,xyz):
 return map_coordinates(field,(xyz/3.125+63.5).T,order=1,mode='constant',cval=np.nan,prefilter=False)
def main():
 p,cp,d,cd=fa.read_data(); chi,hk,den,transfer,mu=fa.distances()
 meta=pd.read_csv(ROOT/'inputs/des_DES-Dovekie_Metadata.csv',sep=r'\s+',comment='#',dtype={'CID':str})
 for c in ['zCMB','x1','c','MWEBV','HOST_LOGMASS','HOST_LOGSFR','HOST_LOGsSFR','HOST_COLOR','LENSDMU','LENSDMUERR','VPEC','VPECERR']:
  if c not in d and c in meta: d[c]=d.CID.map(meta.set_index('CID')[c])
 frames=[]
 for name,df,C,mucol in [('Pantheon',p,cp,'MU_SH0ES'),('DES',d,cd,'MU')]:
  df=df.copy(); df['dataset']=name; df['uid']=[('P' if name=='Pantheon' else 'D')+f'{i:04d}' for i in df.raw_row]
  df['released_mu']=df[mucol]; df['common_mode_residual']=df[mucol]-mu(df.zHD,df.zHEL)
  df['release_pv_delta_mag']=mu(df.zCMB,df.zHEL)-mu(df.zHD,df.zHEL)
  df['covariance_file']='inputs/pantheon_Pantheon+SH0ES_STAT+SYS.cov' if name=='Pantheon' else 'inputs/des_STAT+SYS.npz'
  df['covariance_sha256']=hashlib.sha256((ROOT/df.covariance_file.iloc[0]).read_bytes()).hexdigest()
  df['covariance_semantics']='released total covariance, symmetrized decimal rounding' if name=='Pantheon' else 'inverse of FULL released precision before any selection'
  df['released_covariance_diagonal']=np.diag(C)
  df['primary_selection']=((df.zHD>=.01)&(df.zHD<1)) if name=='Pantheon' else ((df.IDSURVEY==10)&(df.zHD>=.1)&(df.zHD<1.2)&~df.pantheon_overlap)
  df['original_training']=(name=='Pantheon')&(df.zHD>=.01)&(df.zHD<.03)
  frames.append(df)
 df=pd.concat(frames,ignore_index=True); df['matrix_row']=np.arange(len(df))
 gal=SkyCoord(ra=df.RA.to_numpy()*u.deg,dec=df.DEC.to_numpy()*u.deg,frame='icrs').galactic
 uv=gal.cartesian.xyz.value.T; df['galactic_l']=gal.l.degree; df['galactic_b']=gal.b.degree
 df['galaxy_zone_of_avoidance']=abs(df.galactic_b)<np.where((df.galactic_l<30)|(df.galactic_l>330),8,5)
 fields={k:np.load(ROOT/'inputs/foreground'/f'{k}.npy') for k in ['density','density_error','xVelocity','yVelocity','zVelocity','xVelocity_error','yVelocity_error','zVelocity_error']}
 df['distance_hinv_Mpc']=C100*chi(df.zHD); df['geometry_background_extrapolated']=df.zHD>2.1
 xyz=uv*df.distance_hinv_Mpc.to_numpy()[:,None]
 vv=np.column_stack([sample(fields[k+'Velocity'],xyz) for k in 'xyz']); ve=np.column_stack([sample(fields[k+'Velocity_error'],xyz) for k in 'xyz'])
 vlos=np.sum(uv*vv,axis=1); df['independent_pv_kms']=vlos
 df['pv_sigma_diagonal_component_proxy']=np.sqrt(np.sum((uv*ve)**2,axis=1))
 df['pv_error_cross_components_and_sightlines']='not released; diagonal proxy is not a full PV covariance'
 df['independent_pv_available']=np.isfinite(vlos)&(df.distance_hinv_Mpc<=200)
 df['z_independent_pv']=(1+df.zCMB)/(1+vlos/299792.458)-1
 df.loc[~df.independent_pv_available,'z_independent_pv']=np.nan
 df['independent_pv_delta_mag']=mu(df.zCMB,df.zHEL)-mu(df.z_independent_pv,df.zHEL)
 df['pv_replacement_delta_mag']=df.independent_pv_delta_mag-df.release_pv_delta_mag
 df['residual_fixed_external_pv']=df.common_mode_residual+df.pv_replacement_delta_mag
 df['endpoint_density_proxy']=sample(fields['density'],xyz)
 # Integer Cartesian cell + radial shell makes the foreground volume explicit.
 cells={}; cellrows=[]; ii=[]; jj=[]; vals=[]; localvals=[]; env=[]; redges=C100*chi(SHELLS)
 for i,row in df.iterrows():
  r=float(row.distance_hinv_Mpc); end=min(r,redges[-1]); edges=np.r_[np.arange(0,end,2.),end]; dl=np.diff(edges); mid=(edges[:-1]+edges[1:])/2
  pos=mid[:,None]*uv[i]; shell=np.searchsorted(redges,mid,side='right')-1; index=np.floor(pos/10).astype(int)
  dens=sample(fields['density'],pos); valid=(mid<=200)&np.isfinite(dens); delta=dens-1
  rec={'matrix_row':i,'mapped_path_hinv_Mpc':float(dl[valid].sum()),'path_map_fraction':float(dl[valid].sum()/r),'path_density_proxy':float(np.average(dens[valid],weights=dl[valid])) if valid.any() else np.nan}
  # Closed-space Green weight, dimensionless chi units, no lensing normalization fit.
  x=mid/C100; xs=r/C100; sq=np.sqrt(fa.K)
  green=np.sin(sq*x)*np.sin(sq*(xs-x))/(sq*np.sin(sq*xs))
  rec['local_lensing_green_proxy']=float(np.sum(green[valid]*delta[valid]*dl[valid]/C100))
  for s in range(len(SHELLS)-1):
   a=(shell==s); good=a&valid; path=dl[a].sum()
   rec[f'shell{s}_density_proxy']=float(np.average(dens[good],weights=dl[good])) if good.any() else np.nan
   rec[f'shell{s}_map_fraction']=float(dl[good].sum()/path) if path else np.nan
  env.append(rec); accum={}; localaccum={}
  for ind,ss,w,covered in zip(index,shell,dl,valid):
   key=(*ind,int(ss))
   if key not in cells:
    j=len(cells); cells[key]=j; center=(ind+.5)*10; center_density=sample(fields['density'],center[None])[0]
    cellrows.append(dict(cell_id=j,ix=ind[0],iy=ind[1],iz=ind[2],shell=int(ss),x=center[0],y=center[1],z=center[2],density_proxy=center_density,within_local_map=np.linalg.norm(center)<=200 and np.isfinite(center_density)))
   j=cells[key]; accum[j]=accum.get(j,0.)+float(w)
   if covered: localaccum[j]=localaccum.get(j,0.)+float(w)
  for j,w in accum.items(): ii.append(i); jj.append(j); vals.append(w); localvals.append(localaccum.get(j,0.))
 W=coo_matrix((vals,(ii,jj)),shape=(len(df),len(cells))).tocsr(); ct=pd.DataFrame(cellrows)
 ct['sightline_row_crossings']=np.asarray((W>0).sum(axis=0)).ravel(); ct['total_path_hinv_Mpc']=np.asarray(W.sum(axis=0)).ravel()
 WL=coo_matrix((localvals,(ii,jj)),shape=W.shape).tocsr(); WL.eliminate_zeros(); covered=np.asarray(WL.sum(axis=0)).ravel()>0
 ct['has_sampled_map_path']=covered
 save_npz(OUT/'W_FULL_GEOMETRY_PATH_LENGTH.npz',W); save_npz(OUT/'W_LOCAL_MAP_PATH_LENGTH.npz',WL[:,covered])
 ct.to_csv(OUT/'FOREGROUND_CELLS.csv.gz',index=False); ct.loc[covered].to_csv(OUT/'LOCAL_MAP_CELLS.csv',index=False)
 df=df.merge(pd.DataFrame(env),on='matrix_row',validate='one_to_one')
 df['sn_block']=df.dataset+':'+df.CID.str.lower().str.replace(r'^(sn|at)(?=\d)','',regex=True)
 df['dust_provenance']='released MWEBV; lightcurve MW convention tracked separately'
 df['foreground_provenance']='2MRS-NeuralNet galaxy-conditioned matter density (1+delta), not direct galaxy-count density'
 df['lensing_status']='local Green density proxy only; no independent full-path lensing convergence'
 df.to_csv(OUT/'SN_SIGHTLINES.csv.gz',index=False)
 df.groupby('sn_block').agg(rows=('uid','size'),uids=('uid',lambda x:';'.join(x)),RA=('RA','first'),DEC=('DEC','first'),zHD=('zHD','first'),surveys=('IDSURVEY',lambda x:';'.join(map(str,sorted(set(x)))))).to_csv(OUT/'UNIQUE_SN_INDEX.csv')
 np.savez_compressed(OUT/'RELEASED_TOTAL_COVARIANCES.npz',Pantheon=cp,DES=cd)
 stats=dict(rows=len(df),unique_SN_blocks=df.sn_block.nunique(),full_cells=W.shape[1],local_cells=int(ct.within_local_map.sum()),path_sum_error=float(np.max(abs(np.asarray(W.sum(axis=1)).ravel()-np.minimum(df.distance_hinv_Mpc,redges[-1])))),independent_endpoint_PV_rows=int(df.independent_pv_available.sum()),primary=df.groupby('dataset').primary_selection.sum().to_dict(),checks=fa.CHECKS)
 assert stats['path_sum_error']<1e-8
 fa.savejson(OUT/'BUILD_SUMMARY.json',stats); print(stats,flush=True)
if __name__=='__main__': main()
