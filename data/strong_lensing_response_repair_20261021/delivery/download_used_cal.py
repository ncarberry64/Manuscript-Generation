from pathlib import Path
import json,hashlib,urllib.request,concurrent.futures
ROOT=Path(__file__).resolve().parent/'repo'
receipt=json.loads((ROOT/'artifacts/strong_lensing_geometry_20261021/jwst/native_acquisition_complete.json').read_text())
OUT=ROOT/'data/strong_lensing_geometry_20261021/native';OUT.mkdir(parents=True,exist_ok=True)
def get(item):
    p=OUT/item['filename']
    if p.exists() and hashlib.sha256(p.read_bytes()).hexdigest()==item['sha256']:
        return {'filename':p.name,'verified':True,'cached':True}
    uri=item['dataURI'];url='https://mast.stsci.edu/api/v0.1/Download/file?uri='+uri
    digest=hashlib.sha256();n=0;tmp=p.with_suffix('.download')
    with urllib.request.urlopen(url,timeout=90) as response,tmp.open('wb') as f:
        while True:
            block=response.read(1024*1024)
            if not block:break
            n+=len(block)
            if n>item['bytes']+1024:raise RuntimeError('Unexpected archive product size')
            digest.update(block);f.write(block)
    if digest.hexdigest()!=item['sha256']:raise RuntimeError('Archive checksum changed '+p.name)
    tmp.replace(p)
    return {'filename':p.name,'bytes':n,'verified':True}
with concurrent.futures.ThreadPoolExecutor(max_workers=4) as pool:
    for result in pool.map(get,receipt['products']):print(json.dumps(result),flush=True)

