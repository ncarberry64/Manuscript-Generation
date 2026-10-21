"""Finish the exact official dark reference with verified HTTP byte ranges."""
from pathlib import Path
import concurrent.futures
import hashlib
import json
import urllib.request
import time

OUT = Path(__file__).resolve().parent
REFS = OUT / "refs"
NAME = "jwst_nircam_dark_0337.fits"
URL = "https://jwst-crds.stsci.edu/unchecked_get/references/jwst/" + NAME
TOTAL = 6291483840
PART = REFS / (NAME + ".download")
CHUNKS = OUT / "dark_ranges"
CHUNKS.mkdir(exist_ok=True)
START = PART.stat().st_size
BLOCK = 128 * 1024**2
EXPECTED_ETAG = '"177006cc0-5ce5252c09169"'


def fetch(bounds):
    first, last = bounds
    path = CHUNKS / f"{first}_{last}.bin"
    expected = last - first + 1
    if not path.exists() or path.stat().st_size != expected:
        with urllib.request.urlopen(urllib.request.Request(URL, headers={
            "Range": f"bytes={first}-{last}",
            "If-Range": EXPECTED_ETAG,
        }), timeout=120) as response:
            assert response.status == 206
            assert response.headers["Content-Range"] == f"bytes {first}-{last}/{TOTAL}"
            assert response.headers["ETag"] == EXPECTED_ETAG
            with path.open("wb") as stream:
                while block := response.read(4 * 1024**2):
                    stream.write(block)
        assert path.stat().st_size == expected
    digest = hashlib.sha256(path.read_bytes()).hexdigest()
    print(json.dumps({"first": first, "last": last, "bytes": expected, "sha256": digest}), flush=True)
    return dict(first=first, last=last, bytes=expected, sha256=digest)


if __name__ == "__main__":
    started = time.time()
    ranges = [(first, min(first + BLOCK, TOTAL) - 1) for first in range(START, TOTAL, BLOCK)]
    with concurrent.futures.ThreadPoolExecutor(max_workers=8) as pool:
        receipt = list(pool.map(fetch, ranges))
    assert PART.stat().st_size == START
    with PART.open("ab") as stream:
        for row in receipt:
            path = CHUNKS / f"{row['first']}_{row['last']}.bin"
            with path.open("rb") as chunk:
                while block := chunk.read(4 * 1024**2):
                    stream.write(block)
    assert PART.stat().st_size == TOTAL
    digest = hashlib.sha256()
    with PART.open("rb") as stream:
        while block := stream.read(8 * 1024**2):
            digest.update(block)
    destination = REFS / NAME
    PART.replace(destination)
    (OUT / "dark_range_receipt.json").write_text(json.dumps({
        "url": URL, "etag": EXPECTED_ETAG, "serial_prefix_bytes": START,
        "range_receipts": receipt, "bytes": TOTAL, "sha256": digest.hexdigest(),
        "elapsed_seconds": time.time() - started,
    }, indent=2) + "\n")
    print("EXACT_DARK_COMPLETE", TOTAL, digest.hexdigest(), time.time() - started, flush=True)
