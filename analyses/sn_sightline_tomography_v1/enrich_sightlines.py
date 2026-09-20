"""Provenance-complete source/field columns and observed galaxy-count supplement.
The catalogue supplement is NOT a new nuisance or a retuning of the tests.
"""
import os
os.environ['OPENBLAS_NUM_THREADS']='1'
from pathlib import Path
import json,hashlib,urllib.request,io
import numpy as np
import pandas as pd
from scipy.sparse import load_npz
import frozen_analysis as fa
from sightline_stage import sample
ROOT=Path(__file__).resolve().parent; IN=ROOT/'inputs/foreground'; OUT=ROOT/'sightlines'
URL='https://vizier.cfa.harvard.edu/viz-bin/asu-tsv?-source=J/ApJS/199/26/table3&-out.all&-out.max=unlimited'
def main():
 p=IN/'2MRS_GALAXY_CATALOG.tsv'
 if not p.exists(): p.write_bytes(urllib.request.urlopen(URL,timeout=60).read())
 lines=[l for l in p.read_text().splitlines() if l and not l.startswith('#')]; g=pd.read_csv(io.StringIO('\n'.join([lines[0]]+lines[3:])),sep='\t',dtype=str)
 for c in ['GLON','GLAT','cz','RAJ2000','DEJ2000','Kcmag']: g[c]=pd.to_numeric(g[c],errors='coerce')
 assert len(g)==44599
 chi,*_=fa.distances(); g['positive_redshift_available']=g.cz.notna()&(g.cz>0); g['barycentric_z']=g.cz/299792.458; g['distance_hinv_Mpc']=np.where(g.positive_redshift_available,2997.92458*chi(g.barycentric_z),np.nan)
 g['coordinate_status']='barycentric redshift space; no fitted flow correction; not NN real-space positions'
 valid=g.loc[g.positive_redshift_available].copy(); xyz=fa.vectors(valid.GLON,valid.GLAT)*valid.distance_hinv_Mpc.to_numpy()[:,None]; ijk=np.floor(xyz/10).astype(int); shell=np.searchsorted(fa.EDGES,valid.barycentric_z,side='right')-1
 for j,c in enumerate(['ix','iy','iz']): valid[c]=ijk[:,j]
 valid['shell']=shell
 counts=valid.groupby(['ix','iy','iz','shell']).size().rename('observed_2MRS_galaxies').reset_index(); cells=pd.read_csv(OUT/'FOREGROUND_CELLS.csv.gz').merge(counts,on=['ix','iy','iz','shell'],how='left',validate='one_to_one'); cells.observed_2MRS_galaxies=cells.observed_2MRS_galaxies.fillna(0).astype(int)
 # Fixed 5^3 midpoint integration estimates cube/shell intersection volume.
 offsets=np.array(np.meshgrid(*([np.arange(1,10,2)]*3),indexing='ij')).reshape(3,-1).T
 redges=2997.92458*chi(fa.EDGES); vol=[]
 for row in cells.itertuples():
  pts=np.array([row.ix,row.iy,row.iz])*10+offsets; r=np.linalg.norm(pts,axis=1); fraction=np.mean((r>=redges[row.shell])&(r<redges[row.shell+1])); vol.append(1000*fraction)
 cells['shell_cell_volume_midpoint_hinv_Mpc3']=vol; cells['observed_galaxy_density_h3_Mpc3']=np.where(np.array(vol)>0,cells.observed_2MRS_galaxies/np.array(vol),np.nan)
 cells['galaxy_density_status']='observed magnitude-limited counts; no selection-function or sky-mask correction; zero is not evidence of a void'
 cells.to_csv(OUT/'FOREGROUND_CELLS_WITH_GALAXIES.csv.gz',index=False); valid.to_csv(OUT/'GALAXY_TO_CELL_ASSIGNMENTS.csv.gz',index=False); g.to_csv(OUT/'GALAXY_CATALOG_PROCESSED.csv.gz',index=False)
 df=pd.read_csv(OUT/'SN_SIGHTLINES.csv.gz',dtype={'CID':str}); W=load_npz(OUT/'W_FULL_GEOMETRY_PATH_LENGTH.npz'); xy=fa.vectors(df.galactic_l,df.galactic_b)*df.distance_hinv_Mpc.to_numpy()[:,None]
 for c in ['density_error','xVelocity','yVelocity','zVelocity','xVelocity_error','yVelocity_error','zVelocity_error']: df['external_endpoint_'+c]=sample(np.load(IN/(c+'.npy')),xy)
 for s in range(8):
  take=cells.shell.to_numpy()==s; ww=W[:,take]; den=np.asarray(ww.sum(axis=1)).ravel(); density=cells.observed_galaxy_density_h3_Mpc3.to_numpy()[take]; covered=np.isfinite(density); num=ww[:,covered]@density[covered]; actual=np.asarray(ww[:,covered].sum(axis=1)).ravel()
  df[f'shell{s}_observed_galaxy_density_h3_Mpc3']=np.divide(num,actual,out=np.full(len(df),np.nan),where=actual>0)
 df['galaxy_count_density_provenance']='VizieR J/ApJS/199/26/table3; path-weighted observed cell counts; redshift-space and magnitude/footprint selection uncorrected'
 df['tomography_geometry_path_fraction']=np.minimum(df.distance_hinv_Mpc,redges[-1])/df.distance_hinv_Mpc
 df['tomography_geometry_domain']='original shells through z=1; map density/PV limited to 200 h^-1 Mpc'
 lc=pd.read_csv(ROOT/'lightcurves/LIGHTCURVE_RECONSTRUCTION.csv').set_index('uid')
 for c in ['status','common_delta_mag','common_sigma_photometric','delta_c','delta_x1','delta_t0','band_jackknife_max','epoch_jackknife_max','systematic_common_rms','nominal_calibration_status']:
  df['epoch_'+c]=df.uid.map(lc[c])
 if (ROOT/'lightcurves/SINGLE_EPOCH_JACKKNIFE_SUMMARY.csv').exists():
  single=pd.read_csv(ROOT/'lightcurves/SINGLE_EPOCH_JACKKNIFE_SUMMARY.csv').set_index('uid'); df['single_epoch_jackknife_max']=df.uid.map(single.single_epoch_jackknife_max)
 df['common_residual_covariance_status']='released standardized MU only; epoch common delta is separate and cannot inherit released HD covariance'
 df.to_csv(OUT/'SN_SIGHTLINES_COMPLETE.csv.gz',index=False)
 provenance=dict(url=URL,sha256=hashlib.sha256(p.read_bytes()).hexdigest(),rows=len(g),positive_redshift_rows=int(g.positive_redshift_available.sum()),description='descriptive galaxy counts added without changing any nuisance design, fit, threshold or likelihood')
 fa.savejson(IN/'GALAXY_CATALOG_PROVENANCE.json',provenance); print(provenance,flush=True)
if __name__=='__main__': main()
