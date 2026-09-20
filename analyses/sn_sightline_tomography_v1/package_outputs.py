"""Hash the delivered snapshot and build a self-contained reproducible archive."""
from pathlib import Path
import hashlib,json,zipfile
ROOT=Path(__file__).resolve().parent
ARCHIVE='BHSM_SN_SIGHTLINE_TOMOGRAPHY_V1.zip'
def digest(p):
 h=hashlib.sha256()
 with p.open('rb') as f:
  for b in iter(lambda:f.read(1024*1024),b''): h.update(b)
 return h.hexdigest()
def main():
 excluded={ARCHIVE,'RESULT_MANIFEST.json','ARCHIVE_SHA256.txt','PACKAGING_LOG.txt'}
 files=sorted(p for p in ROOT.rglob('*') if p.is_file() and '__pycache__' not in p.parts and p.name not in excluded)
 manifest={p.relative_to(ROOT).as_posix():digest(p) for p in files}
 m=ROOT/'RESULT_MANIFEST.json'; m.write_text(json.dumps(manifest,indent=2),encoding='utf-8'); files.append(m)
 with zipfile.ZipFile(ROOT/ARCHIVE,'w',compression=zipfile.ZIP_DEFLATED,compresslevel=3,allowZip64=True) as z:
  for p in files: z.write(p,Path(ROOT.name)/p.relative_to(ROOT))
 with zipfile.ZipFile(ROOT/ARCHIVE) as z:
  bad=z.testzip(); assert bad is None,bad
  assert len(z.infolist())==len(files)
 h=digest(ROOT/ARCHIVE); (ROOT/'ARCHIVE_SHA256.txt').write_text(h+'  '+ARCHIVE+'\n')
 print(json.dumps(dict(files=len(files),archive_bytes=(ROOT/ARCHIVE).stat().st_size,sha256=h,crc_check='passed'),indent=2),flush=True)
if __name__=='__main__': main()
