"""Epoch/band local SALT reconstructions and jackknife/systematic diagnostics.

Distances in the subsequent released-covariance test are NEVER replaced here.
"""
import os
os.environ['OPENBLAS_NUM_THREADS']='1'
os.environ['OMP_NUM_THREADS']='1'
from pathlib import Path
import json,gzip,re,hashlib,time,traceback,warnings
from collections import defaultdict
import numpy as np
import pandas as pd
from astropy.io import fits
from astropy.io.fits.verify import VerifyWarning
warnings.filterwarnings('ignore',category=VerifyWarning)
import sncosmo
from scipy.linalg import cholesky,solve_triangular,svd
from concurrent.futures import ProcessPoolExecutor
ROOT=Path(__file__).resolve().parent; IN=ROOT/'inputs'; OUT=ROOT/'lightcurves'; OUT.mkdir(exist_ok=True)
PR=IN/'pantheon_release'/'Pantheon+_Data'; DR=IN/'des_release'
KC=PR/'2_CALIBRATION'/'SNANA_kcor'; HC=1.9864458571489286e-8
MAP={1:('JLA2014_SDSS_DS17','kcor_SDSS_v6_1.fits','SDSS'),4:('JLA2014_SNLS_DS17','kcor_SNLS_v6_1.fits','SNLS'),5:('CSPDR3_anthony','kcor_CSPDR3_v6_1.fits','CSP'),10:('DES3YR_DES_COMBINED_TEXT','kcor_DES_3yr_v6_1.fits','DES'),15:('Pantheon_PS1MD','kcor_PS1_v6_1.fits','PS1MD'),18:('ASASSN','kcor_ASASSN_v6_1.fits','ASASSN'),50:('PS1s_LOWZ_JRK07_NOCFA_DS15','kcor_PS1_CFA_COMBINED_v6_1.fits','LOWZ'),51:('KAIT_DS15','kcor_KAIT_Mo_v6_1.fits','KAITM'),56:('SWIFT','kcor_SWIFT_v6_1.fits','SWIFT'),57:('LOSS','kcor_KAIT_Stahl_v6_1.fits','KAIT'),61:('PS1_LOWZ_COMBINED_TEXT_DS17','kcor_PS1_CFA_COMBINED_v6_1.fits','PS1_LOWZ_COMBINED'),62:('PS1_LOWZ_COMBINED_TEXT_DS17','kcor_PS1_CFA_COMBINED_v6_1.fits','PS1_LOWZ_COMBINED'),63:('CfA3_DJ20','kcor_PS1_CFA_COMBINED_v6_1.fits','PS1_LOWZ_COMBINED'),64:('CfA3_DJ20','kcor_PS1_CFA_COMBINED_v6_1.fits','PS1_LOWZ_COMBINED'),65:('PS1_LOWZ_COMBINED_TEXT_DS17','kcor_PS1_CFA_COMBINED_v6_1.fits','PS1_LOWZ_COMBINED'),66:('PS1_LOWZ_COMBINED_TEXT_DS17','kcor_PS1_CFA_COMBINED_v6_1.fits','PS1_LOWZ_COMBINED'),100:('PS1_HST_COMBINED_TEXT_DS17','kcor_PS1_HST_COMBINED_ADAMRECAL.fits','PS1MD'),101:('PS1_HST_COMBINED_TEXT_DS17','kcor_PS1_HST_COMBINED_ADAMRECAL.fits','PS1MD'),106:('PS1_HST_COMBINED_TEXT_DS17','kcor_PS1_HST_COMBINED_ADAMRECAL.fits','PS1MD'),150:('Foundation_DJ17','kcor_Foundation_v6_1.fits','FOUNDATION')}
CACHE={}
MAP[51]=('LOSS','kcor_KAIT_Stahl_v6_1.fits','KAIT')
MAP[57]=('KAIT_DS15','kcor_KAIT_Mo_v6_1.fits','KAITM')
def norm(s):
 s=str(s).strip().lower(); s=re.sub(r'^(sn|at)(?=\d)','',s)
 return str(int(s)) if s.isdigit() else s
def text_read(p):
 return gzip.open(p,'rt',errors='replace').read() if p.suffix=='.gz' else p.read_text(errors='replace')
def ascii_lc(p):
 meta={}; names=None; rows=[]
 for line in text_read(p).splitlines():
  s=line.split('#')[0].strip()
  if s.startswith('VARLIST '): s=s.replace('VARLIST ','VARLIST: ',1)
  if not s or ':' not in s: continue
  k,v=s.split(':',1); vals=v.split()
  if k in ['VARLIST','VARNAMES']: names=vals
  elif k=='OBS' and names and len(vals)>=len(names): rows.append(vals[:len(names)])
  elif vals: meta[k]=vals[0]
 if not rows: return None
 d=pd.DataFrame(rows,columns=names); band=next(c for c in ['FLT','BAND','FILTER'] if c in d)
 out=pd.DataFrame(dict(MJD=pd.to_numeric(d.MJD),BAND=d[band].astype(str),FLUXCAL=pd.to_numeric(d.FLUXCAL),FLUXCALERR=pd.to_numeric(d.FLUXCALERR)))
 return meta,out
def index_photometry():
 catalog=[]; cachefiles={}
 for p in sorted((PR/'1_DATA'/'photometry').glob('*/*')):
  if not p.is_file() or p.suffix.lower() not in ['.dat','.txt'] or 'README' in p.name: continue
  q=ascii_lc(p)
  if q is None: continue
  meta,data=q
  catalog.append(dict(source=str(p),group=p.parent.name,SNID=str(meta.get('SNID','')),IAUC=str(meta.get('IAUC','')),kind='ascii',index=-1))
 for head in sorted((PR/'1_DATA'/'photometry').glob('*/*HEAD.FITS.gz'))+sorted((DR/'0_DATA').glob('*/*HEAD.FITS.gz')):
  h=fits.getdata(head,1); phot=head.with_name(head.name.replace('_HEAD.','_PHOT.')); assert phot.exists()
  for i,row in enumerate(h):
   catalog.append(dict(source=str(head),group=head.parent.name,SNID=str(row['SNID']).strip(),IAUC=str(row['IAUC']).strip() if 'IAUC' in h.dtype.names else '',kind='fits',index=i))
 df=pd.DataFrame(catalog); df.to_csv(OUT/'RAW_LIGHTCURVE_INDEX.csv',index=False)
 return df
def load_obs(item):
 p=Path(item['source'])
 if item['kind']=='ascii': return ascii_lc(p)[1]
 if str(p) not in CACHE:
  CACHE[str(p)]=(fits.getdata(p,1),fits.getdata(p.with_name(p.name.replace('_HEAD.','_PHOT.')),1))
 h,phot=CACHE[str(p)]; row=h[int(item['index'])]
 a=int(row['PTROBS_MIN'])-1; b=int(row['PTROBS_MAX']); x=phot[a:b]
 assert (x['MJD']!=-777).all() and len(x)==b-a
 names=x.dtype.names; bn='BAND' if 'BAND' in names else 'FLT'
 return pd.DataFrame({name:np.asarray(x[field]).astype(str if name=='BAND' else float) for name,field in [('MJD','MJD'),('BAND',bn),('FLUXCAL','FLUXCAL'),('FLUXCALERR','FLUXCALERR')]})
def source(dataset,variant=0):
 key=(dataset,variant)
 if key in CACHE: return CACHE[key]
 if dataset=='Pantheon':
  p=PR/'3_SALT2/SALT2_B21trained_withsys'/f'SALT2.MODEL{variant:03d}'
  kwargs={n:stem+'.dat.gz' for n,stem in [('m0file','salt2_template_0'),('m1file','salt2_template_1'),('clfile','salt2_color_correction'),('cdfile','salt2_color_dispersion'),('errscalefile','salt2_lc_dispersion_scaling'),('lcrv00file','salt2_lc_relative_variance_0'),('lcrv11file','salt2_lc_relative_variance_1'),('lcrv01file','salt2_lc_relative_covariance_01')]}
  cache=ROOT/'model_cache'/p.name; cache.mkdir(parents=True,exist_ok=True)
  for k,fn in list(kwargs.items()):
   target=cache/fn.removesuffix('.gz')
   if not target.exists(): target.write_bytes(gzip.decompress((p/fn).read_bytes()))
   kwargs[k]=target.name
  src=sncosmo.SALT2Source(modeldir=str(cache),**kwargs)
 else:
  p=DR/'2_LCFIT_MODEL'/('SALT3.DOVEKIE' if variant==0 else f'SALT3.DOVEKIE-SYS/SALT3.DOV{variant-1}UNC')
  kwargs={n:stem+'.dat.gz' for n,stem in [('m0file','salt3_template_0'),('m1file','salt3_template_1'),('clfile','salt3_color_correction'),('cdfile','salt3_color_dispersion'),('lcrv00file','salt3_lc_variance_0'),('lcrv11file','salt3_lc_variance_1'),('lcrv01file','salt3_lc_covariance_01')]}
  cache=ROOT/'model_cache'/p.name; cache.mkdir(parents=True,exist_ok=True)
  for k,fn in list(kwargs.items()):
   target=cache/fn.removesuffix('.gz')
   if not target.exists(): target.write_bytes(gzip.decompress((p/fn).read_bytes()))
   kwargs[k]=target.name
  src=sncosmo.SALT3Source(modeldir=str(cache),**kwargs)
 model=sncosmo.Model(source=src,effects=[sncosmo.F99Dust()],effect_names=['mw'],effect_frames=['obs'])
 info=(p/('SALT2.INFO' if dataset=='Pantheon' else 'SALT3.INFO')).read_text()
 CACHE[(dataset,variant,'mag_offset')]=float(re.search(r'MAG_OFFSET:\s+([\d.+-]+)',info).group(1))
 CACHE[key]=model; return model
def passbands(kcor):
 if kcor in CACHE: return CACHE[kcor]
 if kcor=='DES_LOWZ_REMAPPED':
  cfa=passbands('kcor_PS1_CFA_COMBINED_v6_1.fits'); csp=passbands('kcor_CSPDR3_v6_1.fits'); out={}
  for new,old in zip('abcdefghijklmno','ugriBVomnYyJjHh'): out[new]=csp[old]
  for new,old in zip('pqrstuvwxy'+'zABCD'+'EFGH','fghijabcde'+'klmno'+'pqrs'): out[new]=cfa[old]
  CACHE[kcor]=out; return out
 f=fits.open(KC/kcor); trans=f['FilterTrans'].data; standard=f['PrimarySED'].data; zp=f['ZPoff'].data
 wave=np.asarray(trans.field(0),float); swave=np.asarray(standard.field(0),float)
 bands={}
 for row in zp:
  name=str(row['Filter Name']).strip(); primary=str(row['Primary Name']).strip(); mag=float(row['Primary Mag'])
  if name.startswith('*'): continue  # SNANA auxiliary wavelength-shift derivative bands
  response=np.asarray(trans[name],float); take=response>1e-8
  ids=np.flatnonzero(take)
  if len(ids)<2: continue
  a=max(0,ids[0]-1); b=min(len(wave),ids[-1]+2); band=sncosmo.Bandpass(wave[a:b],response[a:b],name=kcor+':'+name)
  ref=np.asarray(standard[primary],float)
  refphot=np.trapezoid(np.interp(wave,swave,ref)*wave/HC*response,wave)
  factor=10**(.4*(27.5-mag))/refphot
  short=name.split('/')[-1] if '/' in name else name[-1]
  bands[short]=(band,factor,name,swave,ref,mag)
 f.close(); CACHE[kcor]=bands; return bands
def sys_shifts(dataset,variant,section,bandname,short):
 if variant==0: return 0.,0.,True
 if dataset=='Pantheon':
  key='psys'
  if key not in CACHE:
   current=''; dd={}
   for line in (PR/'3_SALT2/CALIB_fitopts/ALL.fitopts').read_text().splitlines():
    if re.match(r'^\w+:',line): current=line.split(':')[0]
    match=re.search(r'cal_(\d+):',line)
    if not match: continue
    vv=int(match.group(1)); d={}
    for field in ['MAGOBS_SHIFT_ZP','FILTER_LAMSHIFT']:
     mm=re.search(field+r"\s+'([^']+)'",line)
     ls=mm.group(1).split() if mm else []; d[field]={ls[i]:float(ls[i+1]) for i in range(0,len(ls),2)}
    dd[(current,vv)]=d
   CACHE[key]=dd
  d=CACHE[key].get((section,variant))
  if d is None: return 0.,0.,False
  return d['MAGOBS_SHIFT_ZP'].get(short,0.),d['FILTER_LAMSHIFT'].get(short,0.),True
 p=DR/f'2_LCFIT_MODEL/SALT3.DOVEKIE-SYS/SALT3.DOV{variant-1}UNC/SALT3.INFO'
 key=str(p)
 if key not in CACHE:
  dd={}
  for line in p.read_text().splitlines():
   ss=line.split()
   if ss and ss[0] in ['MAGSHIFT:','WAVESHIFT:'] and len(ss)>=4: dd[(ss[0],ss[1],ss[2])]=float(ss[3])
  CACHE[key]=dd
 dd=CACHE[key]; snana_section={'FOUNDATION':'FOUNDATION','PS1_LOWZ_COMBINED':'CFA3','DES':'DES','CSP':'CSP'}.get(section,section)
 bn=bandname if section!='FOUNDATION' else 'PS1-'+short
 if section=='DES_LOWZ_REMAPPED':
  if short in 'abcdefghijklmno':
   snana_section='CSP'; bn={'CSP-g':'CSP-g/A','CSP-r':'CSP-r/L','CSP-i':'CSP-i/C','CSP-B':'CSP-B/u','CSP-o':'CSP-o/v','CSP-m':'CSP-m/w','CSP-n':'CSP-n/x'}.get(bandname,bandname)
  elif short in 'pqrst': snana_section='CFA3'; bn=dict(zip('pqrst',['CFA3K-U','CFA3K-B/h','CFA3K-V/j','CFA3K-r/k','CFA3K-i/l']))[short]
  elif short in 'uvwxy': snana_section='CFA3'
  elif short in 'zABCD': snana_section='CFA4p1'; bn=dict(zip('zABCD',['CFA41-U','CFA41-B/D','CFA41-V/E','CFA41-r/F','CFA41-i/G']))[short]
  else: snana_section='CFA4p2'; bn=dict(zip('EFGH',['CFA42-B/P','CFA42-V/Q','CFA42-r/W','CFA42-i/T']))[short]
 if section!='DES_LOWZ_REMAPPED' and bandname.startswith('CFA41'): snana_section='CFA4p1'
 if section!='DES_LOWZ_REMAPPED' and bandname.startswith('CFA42'): snana_section='CFA4p2'
 return dd.get(('MAGSHIFT:',snana_section,bn),0.),dd.get(('WAVESHIFT:',snana_section,bn),0.),('MAGSHIFT:',snana_section,bn) in dd
def evaluate(model,pars,row,obs,bands,variant=0):
 model.set(z=float(row['zHEL']),t0=pars[3],x0=np.exp(pars[0]),x1=pars[1],c=pars[2],mwebv=max(float(row['MWEBV']),0)*row.get('dust_scale',.86))
 flux=np.empty(len(obs['MJD'])); allbands=[]; factors=np.empty(len(flux)); datafactor=np.ones(len(flux)); complete=True
 for b in np.unique(obs['BAND']):
  band,factor,name,sw,ref,mag=bands[b]; take=obs['BAND']==b
  dm,dw,found=sys_shifts(row['dataset'],variant,row['sys_section'],name,b); complete &= found
  if dw:
   band=sncosmo.Bandpass(band.wave+dw,band.trans,name=name+f'_shift{dw:g}')
   rf=np.trapezoid(np.interp(band.wave,sw,ref)*band.wave/HC*band.trans,band.wave)
   factor=10**(.4*(27.5-mag))/rf
  flux[take]=model.bandflux(band,obs['MJD'][take])*factor*10**(-.4*CACHE[(row['dataset'],variant,'mag_offset')])
  factors[take]=factor; datafactor[take]=10**(-.4*dm)
  allbands.extend([])
 return flux,datafactor,complete
def jacobian(model,pars,row,obs,bands):
 f=evaluate(model,pars,row,obs,bands)[0]; J=np.empty((len(f),4)); J[:,0]=f
 for j,h in [(1,.01),(2,.001),(3,.02)]:
  pp=pars.copy(); pp[j]+=h; pm=pars.copy(); pm[j]-=h
  J[:,j]=(evaluate(model,pp,row,obs,bands)[0]-evaluate(model,pm,row,obs,bands)[0])/(2*h)
 return f,J
def solve(J,y,L):
 if len(y)<4: raise ValueError('Too few retained epochs for four source parameters')
 jw=solve_triangular(L,J,lower=True); yw=solve_triangular(L,y,lower=True)
 norms=np.linalg.norm(jw,axis=0); u,s,vt=svd(jw/norms,full_matrices=False)
 if s[-1]<s[0]*1e-8: raise ValueError('SALT source parameters rank deficient')
 R=(vt.T/s)@u.T/norms[:,None]
 H=solve_triangular(L.T,R.T,lower=False).T
 return H@y,R@R.T,H,float(s[0]/s[-1])
def reconstruct(task):
 row,item=task; uid=row['uid']; base=dict(uid=uid,dataset=row['dataset'],raw_row=row['raw_row'],CID=row['CID'],IDSURVEY=row['IDSURVEY'],source=item.get('source',''),status='failed')
 try:
  if not item: raise ValueError('No unique matching public light curve for release CID/survey')
  obs=load_obs(item); obs['BAND']=obs.BAND.str.strip().str.split('/').str[-1].str.split('-').str[-1]
  if row['dataset']=='Pantheon' and row['IDSURVEY']==50 and set(obs.BAND)<=set('UBVRI'): row['kcor']='kcor_Land2.fits'
  bands=passbands(row['kcor']); model=source(row['dataset'])
  init=np.array([np.log(float(row['x0'])),float(row['x1']),float(row['c']),float(row['PKMJD'])])
  phase=(obs.MJD.to_numpy()-init[3])/(1+float(row['zHEL']))
  usable=(obs.FLUXCALERR>0)&np.isfinite(obs.FLUXCALERR)&np.isfinite(obs.FLUXCAL)&(phase>=-15)&(phase<=45)&obs.BAND.isin(bands)
  if row['IDSURVEY']==56: usable &= ~obs.BAND.isin(['X','N','W'])  # released kcor explicitly supplies placeholder U curves here
  for b,vals in bands.items():
   band=vals[0]
   if band.minwave()/(1+row['zHEL'])<model.source.minwave() or band.maxwave()/(1+row['zHEL'])>model.source.maxwave(): usable &= obs.BAND!=b
  data=obs.loc[usable].reset_index(drop=True); phase=phase[usable]
  if len(data)<5 or data.BAND.nunique()<2: raise ValueError('Insufficient supported epochs/bands in fixed diagnostic window')
  obs={c:data[c].to_numpy() for c in data}; pars=init.copy(); fitflux,J=jacobian(model,pars,row,obs,bands)
  # Measurement covariance only: model and shared calibration covariance are
  # not silently invented. Released full HD covariance is used in stage 2.
  C=np.diag(obs['FLUXCALERR']**2); L=cholesky(C,lower=True)
  converged=False
  for iteration in range(8):
   fitflux,J=jacobian(model,pars,row,obs,bands)
   dp,cov,H,cond=solve(J,obs['FLUXCAL']-fitflux,L)
   step=np.clip(dp,[-.3,-1,-.1,-2],[.3,1,.1,2]); pars+=step
   if np.max(abs(step)/np.array([.001,.01,.001,.02]))<1: converged=True; break
  fitflux,J=jacobian(model,pars,row,obs,bands); residual=obs['FLUXCAL']-fitflux
  dp,cov,H,cond=solve(J,residual,L)
  common=-2.5/np.log(10)*(pars[0]+dp[0]-init[0])
  standardized_chi=float(np.sum((residual/obs['FLUXCALERR'])**2))
  # Residual band/phase components are orthogonal to SALT fitted columns.
  jw=solve_triangular(L,J,lower=True); q=svd(jw,full_matrices=False)[0]
  rw=solve_triangular(L,residual,lower=True); rw=rw-q@(q.T@rw)
  bcols=np.column_stack([(obs['BAND']==b)*fitflux for b in np.unique(obs['BAND'])]); bw=solve_triangular(L,bcols,lower=True); bw-=q@(q.T@bw)
  phcols=np.column_stack([fitflux*phase/30,fitflux*(phase/30)**2]); pw=solve_triangular(L,phcols,lower=True); pw-=q@(q.T@pw)
  def contrast_power(x):
   u,s,v=svd(x,full_matrices=False); ok=s>s[0]*1e-8 if s[0]>0 else s>0
   return float(np.sum((u[:,ok].T@rw)**2)),int(ok.sum())
  bchi,brank=contrast_power(bw); pchi,prank=contrast_power(pw)
  records=[]
  masks=[('band_'+b,obs['BAND']!=b) for b in np.unique(obs['BAND'])]+[(f'phase_{lo}_{hi}',~((phase>=lo)&(phase<hi))) for lo,hi in [(-15,0),(0,15),(15,45.00001)]]
  shifts=[]
  for label,keep in masks:
   rec=dict(uid=uid,kind=label,kept=int(keep.sum()),status='rank_failed')
   try:
    shift,jcov,jH,jcond=solve(J[keep],residual[keep],cholesky(C[np.ix_(keep,keep)],lower=True))
    hh=np.zeros(len(data)); hh[keep]=jH[0]; delta_h=hh-H[0]
    dm=-2.5/np.log(10)*(shift[0]-dp[0]); se=2.5/np.log(10)*np.sqrt(max(float(delta_h@C@delta_h),0))
    rec.update(status='ok',delta_common_mag=dm,sigma_difference=se,color_shift=shift[2]-dp[2],shape_shift=shift[1]-dp[1],phase_shift=shift[3]-dp[3],condition=jcond)
    shifts.append((label,dm))
   except (ValueError,np.linalg.LinAlgError): pass
   records.append(rec)
  # All alternative released realizations: local response, not a new bias fit.
  sysrecords=[]
  for v in range(1,10 if row['dataset']=='Pantheon' else 11):
   sr=dict(uid=uid,variant=v,status='failed')
   try:
    ff,scale,complete=evaluate(source(row['dataset'],v),pars,row,obs,bands,v)
    ds=H@(obs['FLUXCAL']*(scale-1)-(ff-fitflux))
    sr.update(status='ok',delta_common_mag=-2.5/np.log(10)*ds[0],delta_color=ds[2],delta_shape=ds[1],delta_phase=ds[3],calibration_mapping_complete=complete)
   except Exception as e: sr['reason']=str(e)[:200]
   sysrecords.append(sr)
  bsh=[abs(x) for k,x in shifts if k.startswith('band')]; psh=[abs(x) for k,x in shifts if k.startswith('phase')]
  sy=[s['delta_common_mag'] for s in sysrecords if s['status']=='ok']
  base.update(status='ok' if converged else 'local_fit_not_converged',epochs=len(data),bands=data.BAND.nunique(),common_delta_mag=common,common_sigma_photometric=2.5/np.log(10)*np.sqrt(cov[0,0]),delta_x1=pars[1]-init[1],delta_c=pars[2]-init[2],delta_t0=pars[3]-init[3],chi2_photometric=standardized_chi,dof_photometric=len(data)-4,band_residual_chi2=bchi,band_residual_rank=brank,phase_residual_chi2=pchi,phase_residual_rank=prank,band_jackknife_max=max(bsh) if bsh else None,epoch_jackknife_max=max(psh) if psh else None,systematic_common_rms=float(np.std(sy)) if sy else None,systematic_variants=len(sy),condition=cond,nominal_calibration_status=row['calibration_status'])
  np.savez_compressed(OUT/(uid+'_FIT.npz'),initial=init,parameters=pars,parameter_covariance=cov,J=J,estimator=H,epoch_covariance=C)
  pd.DataFrame(dict(MJD=obs['MJD'],BAND=obs['BAND'],rest_phase=phase,flux=obs['FLUXCAL'],fluxerr=obs['FLUXCALERR'],model_flux=fitflux,residual_flux=residual)).to_csv(OUT/(uid+'_EPOCHS.csv.gz'),index=False)
  return base,records,sysrecords
 except Exception as e:
  base['reason']=str(e)[:250]; return base,[],[]
def build_tasks():
 idx=index_photometry(); byname=defaultdict(list)
 for item in idx.to_dict('records'):
  for k in set([norm(item['SNID']),norm(item['IAUC'])]):
   if k and k not in ['unknown','null']: byname[k].append(item)
 p=pd.read_csv(IN/'pantheon_Pantheon+SH0ES.dat',sep=r'\s+',dtype={'CID':str}); p['dataset']='Pantheon'; p['raw_row']=np.arange(len(p))
 d=pd.read_csv(IN/'des_DES-Dovekie_HD.csv',sep=r'\s+',comment='#',dtype={'CID':str}); m=pd.read_csv(IN/'des_DES-Dovekie_Metadata.csv',sep=r'\s+',comment='#',dtype={'CID':str}); d['raw_row']=np.arange(len(d)); d=d.merge(m[['CID','x0','x1','c','PKMJD','MWEBV']],on='CID',how='left',validate='one_to_one',sort=False); d['dataset']='DES'
 tasks=[]
 for row in pd.concat([p,d],ignore_index=True).to_dict('records'):
  row['uid']=('P' if row['dataset']=='Pantheon' else 'D')+f"{row['raw_row']:04d}"
  rule=MAP.get(int(row['IDSURVEY']))
  if rule is None: tasks.append((row,{})); continue
  group,kcor,sec=rule
  if row['dataset']=='DES':
   group='DES-SN5YR_DES' if row['IDSURVEY']==10 else ('DES-SN5YR_Foundation' if row['IDSURVEY']==150 else 'DES-SN5YR_LOWZ')
   if row['IDSURVEY']==10: kcor='kcor_DES_5yr_v6_1.fits'
   elif group=='DES-SN5YR_LOWZ': kcor='DES_LOWZ_REMAPPED'; sec='DES_LOWZ_REMAPPED'
  row.update(kcor=kcor,sys_section=sec,calibration_status='Pantheon released kcor' if row['dataset']=='Pantheon' else 'Pantheon kcor proxy; nominal Dovekie kcor not released here')
  candidates=[a for a in byname[norm(row['CID'])] if a['group']==group]
  if not candidates and row['IDSURVEY']==5: candidates=[a for a in byname[norm(row['CID'])] if a['group']=='CSP_data2']
  exact=[a for a in candidates if norm(a['SNID'])==norm(row['CID'])]
  if len(exact)==1: candidates=exact
  if len(candidates)>1:
   preferred=[]
   for candidate in candidates:
    pp=Path(candidate['source']); listfile=pp.parent/(pp.parent.name+'.LIST')
    if listfile.exists() and pp.name in listfile.read_text().split(): preferred.append(candidate)
   if len(preferred)==1: candidates=preferred
  item=candidates[0] if len(candidates)==1 else {}
  row['dust_scale']=.86
  if item and item['kind']=='ascii':
   mwline=next((l for l in text_read(Path(item['source'])).splitlines() if l.startswith('MWEBV:')),'')
   if 'SF11' in mwline: row['dust_scale']=1.
  tasks.append((row,item))
 pd.DataFrame([{**{k:r.get(k) for k in ['uid','dataset','raw_row','CID','IDSURVEY','kcor','sys_section','calibration_status']},**item} for r,item in tasks]).to_csv(OUT/'ROW_TO_LIGHTCURVE.csv',index=False)
 return tasks
def main():
 tasks=build_tasks()
 limit=int(os.environ.get('LC_LIMIT','0'))
 if limit: tasks=tasks[:limit]
 results=[]; jack=[]; systems=[]; start=time.time()
 # Populate derived model cache before workers can read partially written files.
 for dataset,count in [('Pantheon',10),('DES',11)]:
  for variant in range(count): source(dataset,variant)
 with ProcessPoolExecutor(max_workers=4) as ex:
  for i,(b,j,s) in enumerate(ex.map(reconstruct,tasks,chunksize=5)):
   results.append(b); jack.extend(j); systems.extend(s)
   if (i+1)%100==0: print(i+1,'/',len(tasks),'ok',sum(r['status']=='ok' for r in results),'elapsed',round(time.time()-start),flush=True)
 pd.DataFrame(results).to_csv(OUT/'LIGHTCURVE_RECONSTRUCTION.csv',index=False)
 pd.DataFrame(jack).to_csv(OUT/'LIGHTCURVE_JACKKNIFE.csv',index=False)
 pd.DataFrame(systems).to_csv(OUT/'LIGHTCURVE_SYSTEMATIC_REALIZATIONS.csv',index=False)
 print(pd.DataFrame(results).groupby(['dataset','status']).size().to_string(),flush=True)
 print('Failures',pd.DataFrame([r for r in results if r['status']=='failed']).get('reason',pd.Series(dtype=str)).value_counts().head(12).to_string(),flush=True)
if __name__=='__main__': main()
