"""Acquire only an already-consumed F115W exposure and its recorded references."""
from pathlib import Path
import hashlib, json, urllib.request, concurrent.futures

ROOT = Path(__file__).resolve().parent
OUT = ROOT / 'ramp_diagnostic'
OUT.mkdir(exist_ok=True)
names = ['jwst_nircam_linearity_0054.fits', 'jwst_nircam_superbias_0180.fits',
         'jwst_nircam_saturation_0100.fits']
items = [(n, 'https://jwst-crds.stsci.edu/unchecked_get/references/jwst/' + n) for n in names]
raw = 'jw01198007001_02101_00001_nrcb3_uncal.fits'
items.append((raw, 'https://mast.stsci.edu/api/v0.1/Download/file?uri=mast:JWST/product/' + raw))

def get(item):
    name, url = item
    dest = OUT / name
    if not dest.exists():
        part = dest.with_suffix('.download')
        size = 0
        with urllib.request.urlopen(url, timeout=90) as response, part.open('wb') as f:
            while block := response.read(1024 * 1024):
                size += len(block)
                if size > 300 * 1024**2:
                    raise RuntimeError('Unexpected product size')
                f.write(block)
        part.replace(dest)
    digest = hashlib.sha256()
    with dest.open('rb') as f:
        while block := f.read(1024 * 1024):
            digest.update(block)
    row = {'filename': name, 'url': url, 'bytes': dest.stat().st_size,
           'sha256': digest.hexdigest()}
    print(json.dumps(row), flush=True)
    return row

with concurrent.futures.ThreadPoolExecutor(max_workers=3) as pool:
    rows = list(pool.map(get, items))
(OUT / 'download_receipt.json').write_text(json.dumps({
    'purpose': 'Development diagnostic of existing F115W photons; no unused band acquired.',
    'products': rows}, indent=2) + '\n')
