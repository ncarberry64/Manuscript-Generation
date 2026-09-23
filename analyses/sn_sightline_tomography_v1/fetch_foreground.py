from pathlib import Path
import urllib.request,zipfile,hashlib,json
ROOT=Path(__file__).resolve().parent; IN=ROOT/'inputs'/'foreground'; IN.mkdir(exist_ok=True)
url='https://www.dropbox.com/scl/fo/wb8iyg113hyin4ni7srkg/h?rlkey=bfry3x0s612qtnmgb6n82rnny&dl=1'
p=IN/'2MRS_NeuralNet_release.zip'
if not p.exists():
 with urllib.request.urlopen(url,timeout=180) as r, p.open('wb') as f:
  while chunk:=r.read(1024*1024): f.write(chunk)
print('Downloaded',p.stat().st_size,'bytes',flush=True)
records=[]
with zipfile.ZipFile(p) as z:
 for info in z.infolist():
  if not info.filename.endswith('.npy'): continue
  target=IN/Path(info.filename).name
  target.write_bytes(z.read(info)); records.append(dict(file=target.name,bytes=target.stat().st_size,sha256=hashlib.sha256(target.read_bytes()).hexdigest()))
(IN/'PROVENANCE.json').write_text(json.dumps(dict(url=url,archive_sha256=hashlib.sha256(p.read_bytes()).hexdigest(),reference='https://github.com/rlilow/2MRS-NeuralNet',files=records),indent=2))
print(records,flush=True)
