"""Acquire only the preregistered native files after a committed physical freeze.

This writes bytes without opening any numerical FITS plane. Field calibration
and target measurement have separate guards and exposure receipts.
"""
from __future__ import annotations

import argparse
from concurrent.futures import ThreadPoolExecutor
from datetime import datetime, timezone
import hashlib
import json
from pathlib import Path
import subprocess

import requests

ROOT = Path(__file__).resolve().parents[1]


def digest(path):
    h = hashlib.sha256()
    with Path(path).open('rb') as stream:
        for block in iter(lambda: stream.read(1024*1024), b''):
            h.update(block)
    return h.hexdigest()


def verify_freeze(freeze_path, commit):
    freeze_path = Path(freeze_path).resolve()
    freeze = json.loads(freeze_path.read_text(encoding='utf-8'))
    inputs = dict(freeze['source_sha256'])
    inputs[freeze_path.relative_to(ROOT).as_posix()] = digest(freeze_path)
    for relative, expected in inputs.items():
        path = ROOT/relative
        if path.resolve().is_relative_to(ROOT) is False:
            raise ValueError('Frozen input leaves workspace')
        committed = subprocess.check_output(['git', 'show', f'{commit}:{relative}'], cwd=ROOT)
        if hashlib.sha256(committed).hexdigest() != expected or digest(path) != expected:
            raise RuntimeError(f'Frozen input changed: {relative}')
    return freeze


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--freeze', required=True)
    parser.add_argument('--freeze-commit', required=True)
    parser.add_argument('--download', action='store_true')
    args = parser.parse_args()
    freeze = verify_freeze(args.freeze, args.freeze_commit)
    protocol = json.loads((ROOT/freeze['protocol']).read_text(encoding='utf-8'))
    filenames = protocol['primary_observation']['products']
    metadata_path = ROOT/'data/strong_lensing_geometry_20261021/metadata/jwst_nircam_f115w_products_metadata.json'
    rows = json.loads(metadata_path.read_text(encoding='utf-8'))['response']['data']
    products = {r['productFilename']: r for r in rows
                if r.get('productSubGroupDescription') == 'CAL' and r['productFilename'] in filenames}
    if len(products) != 4 or set(products) != set(filenames):
        raise RuntimeError('Exactly the four declared native CAL exposures are required')
    if not args.download:
        print(json.dumps({'status': 'FREEZE_VERIFIED_NO_DOWNLOAD', 'products': filenames}))
        return
    receipt_dir = ROOT/'artifacts/strong_lensing_geometry_20261021/jwst'
    receipt_dir.mkdir(parents=True, exist_ok=True)
    marker = receipt_dir/'native_acquisition_started.json'
    with marker.open('x', encoding='utf-8') as stream:
        json.dump({'created_utc': datetime.now(timezone.utc).isoformat(),
                   'freeze_commit': args.freeze_commit, 'freeze_sha256': digest(args.freeze),
                   'products': filenames, 'science_arrays_opened': False}, stream, indent=2)
    destination = ROOT/'data/strong_lensing_geometry_20261021/native'
    destination.mkdir(parents=True, exist_ok=True)

    def fetch(filename):
        row = products[filename]
        path = destination/filename
        url = 'https://mast.stsci.edu/api/v0.1/Download/file'
        with requests.get(url, params={'uri': row['dataURI']}, stream=True, timeout=(30, 120)) as response:
            response.raise_for_status()
            with path.open('xb') as stream:
                for block in response.iter_content(1024*1024):
                    stream.write(block)
        if path.stat().st_size != row['size']:
            raise RuntimeError(f'Native file size mismatch: {filename}')
        return {'filename': filename, 'dataURI': row['dataURI'], 'bytes': path.stat().st_size,
                'sha256': digest(path), 'science_arrays_opened': False}

    with ThreadPoolExecutor(max_workers=2) as pool:
        receipts = list(pool.map(fetch, filenames))
    (receipt_dir/'native_acquisition_complete.json').write_text(json.dumps({
        'completed_utc': datetime.now(timezone.utc).isoformat(),
        'freeze_commit': args.freeze_commit, 'products': receipts,
        'science_arrays_opened': False}, indent=2)+'\n', encoding='utf-8')
    print(json.dumps({'status': 'FOUR_FROZEN_NATIVE_FILES_ACQUIRED_WITHOUT_OPENING_ARRAYS',
                      'bytes': sum(row['bytes'] for row in receipts)}))


if __name__ == '__main__':
    main()
