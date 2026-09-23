from pathlib import Path
import json,hashlib,urllib.request,urllib.parse,shutil,time,zipfile
from concurrent.futures import ThreadPoolExecutor,as_completed
ROOT=Path(__file__).resolve().parent; IN=ROOT/'inputs'
P='c447f0fea703fcd0fff57de5000947b5ca81286b'; D='c9a4fcafc4cbd19bd750dee47fc76194a45c181f'
DESLOCAL=ROOT.parent/'DES-SN5YR-main'/'DES-SN5YR-main'
provenance=[]
def sha(p): return hashlib.sha256(p.read_bytes()).hexdigest()
def blob(p):
 b=p.read_bytes(); return hashlib.sha1(f'blob {len(b)}\0'.encode()+b).hexdigest()
def fetch(url,p,gitsha=None):
 p.parent.mkdir(parents=True,exist_ok=True)
 if not p.exists():
  for attempt in range(3):
   try:
    with urllib.request.urlopen(urllib.request.Request(url,headers={'User-Agent':'BHSM-scientific-reproducibility'}),timeout=180) as f: p.write_bytes(f.read())
    break
   except Exception:
    if attempt==2: raise
    time.sleep(2)
 if gitsha and blob(p)!=gitsha: raise ValueError('Git hash mismatch '+str(p))
 return dict(path=str(p.relative_to(ROOT)),url=url,git_blob=gitsha,sha256=sha(p),bytes=p.stat().st_size)
jobs=[]
tree=json.loads((IN/'pantheon_tree.json').read_text(encoding='utf-8-sig'))['tree']
for entry in tree:
 path=entry['path']
 if entry['type']!='blob': continue
 if path.startswith(('Pantheon+_Data/1_DATA/photometry/','Pantheon+_Data/2_CALIBRATION/SNANA_kcor/','Pantheon+_Data/3_SALT2/')) and '/.#' not in path:
  jobs.append((f'https://raw.githubusercontent.com/PantheonPlusSH0ES/DataRelease/{P}/'+urllib.parse.quote(path),IN/'pantheon_release'/path,entry['sha']))
with ThreadPoolExecutor(max_workers=16) as pool:
 fs=[pool.submit(fetch,*job) for job in jobs]
 for i,f in enumerate(as_completed(fs)):
  provenance.append(f.result())
  if (i+1)%500==0: print('Pantheon files',i+1,'/',len(jobs),flush=True)
tree=json.loads((IN/'des_tree.json').read_text(encoding='utf-8-sig'))['tree']
for entry in tree:
 path=entry['path']
 if entry['type']!='blob': continue
 if path.startswith(('0_DATA/','2_LCFIT_MODEL/','7_PIPPIN_FILES/base_files/lcfit/')) and '/plots/' not in path:
  p=IN/'des_release'/path; source=DESLOCAL/path
  if source.exists() and blob(source)==entry['sha']:
   p.parent.mkdir(parents=True,exist_ok=True); shutil.copy2(source,p)
  provenance.append(fetch(f'https://raw.githubusercontent.com/des-science/DES-SN5YR/{D}/'+urllib.parse.quote(path),p,entry['sha']))
base=ROOT.parent/'BHSM_PV_INTERPRETATION_TESTS'
for name in ['pantheon_Pantheon+SH0ES.dat','pantheon_Pantheon+SH0ES_STAT+SYS.cov','des_DES-Dovekie_HD.csv','des_DES-Dovekie_Metadata.csv','des_STAT+SYS.npz','DES-SN5YR_DES_HEAD.FITS.gz','DES-SN5YR_Foundation_HEAD.FITS.gz','DES-SN5YR_LOWZ_HEAD.FITS.gz','R1_SI_preregistered_redshift_table.csv','background.py','bhsm_optical_background_reconstruction.py']:
 src=base/'inputs'/name; p=IN/name; shutil.copy2(src,p); provenance.append(dict(path=str(p.relative_to(ROOT)),source=str(src),sha256=sha(p),bytes=p.stat().st_size))
for name in ['frozen_analysis.py','FROZEN_VALIDATION_RESULT.json','FROZEN_GLS_OPERATORS.npz']:
 src=base/name; p=ROOT/name; shutil.copy2(src,p); provenance.append(dict(path=name,source=str(src),sha256=sha(p),bytes=p.stat().st_size))
(ROOT/'INPUT_PROVENANCE.json').write_text(json.dumps(provenance,indent=2))
print('Completed',len(provenance),'input files',flush=True)
