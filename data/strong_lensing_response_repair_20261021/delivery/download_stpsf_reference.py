from pathlib import Path
import hashlib,urllib.request,tarfile,json
ROOT=Path(__file__).resolve().parent/'stpsf_reference'
ROOT.mkdir(parents=True,exist_ok=True)
url='https://stsci.box.com/shared/static/kqfolg2bfzqc4mjkgmujo06d3iaymahv.gz'
path=ROOT/'stpsf-data-LATEST.tar.gz'
sha=hashlib.sha256();n=0
with urllib.request.urlopen(url,timeout=90) as response,path.open('wb') as f:
 while True:
  b=response.read(1024*1024)
  if not b:break
  n+=len(b)
  if n>250*1024**2:raise RuntimeError('Unexpected reference archive size')
  f.write(b);sha.update(b)
with tarfile.open(path) as tf:tf.extractall(ROOT,filter='data')
receipt={'url':url,'sha256':sha.hexdigest(),'bytes':n,'scope':'Official STScI reference pupil/filter/instrument response; no quasar target observations'}
(ROOT/'download_receipt.json').write_text(json.dumps(receipt,indent=2)+'\n')
print(json.dumps(receipt),flush=True)
print([str(p) for p in ROOT.iterdir() if p.is_dir()],flush=True)

