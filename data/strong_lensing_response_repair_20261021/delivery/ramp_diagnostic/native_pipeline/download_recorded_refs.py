"""Acquire only references named by the already consumed first F115W CAL."""
from pathlib import Path
import concurrent.futures
import hashlib
import json
import urllib.request

OUT = Path(__file__).resolve().parent
REFS = OUT / "refs"
REFS.mkdir(exist_ok=True)
NAMES = [
    "jwst_nircam_sirskernel_0002.asdf",
    "jwst_nircam_mask_0075.fits",
    "jwst_nircam_gain_0094.fits",
    "jwst_nircam_readnoise_0282.fits",
    "jwst_nircam_flat_0667.fits",
    "jwst_nircam_dark_0337.fits",
]


def acquire(name):
    url = "https://jwst-crds.stsci.edu/unchecked_get/references/jwst/" + name
    destination = REFS / name
    if not destination.exists():
        temporary = destination.with_suffix(destination.suffix + ".download")
        size = 0
        digest = hashlib.sha256()
        with urllib.request.urlopen(url, timeout=120) as response:
            expected = int(response.headers["Content-Length"])
            if expected > 7 * 1024**3:
                raise RuntimeError(f"Unexpected reference size: {name}: {expected}")
            with temporary.open("wb") as stream:
                while block := response.read(4 * 1024**2):
                    size += len(block)
                    digest.update(block)
                    stream.write(block)
        if size != expected:
            raise RuntimeError(f"Incomplete reference {name}: {size}/{expected}")
        temporary.replace(destination)
    else:
        size = destination.stat().st_size
        digest = hashlib.sha256()
        with destination.open("rb") as stream:
            while block := stream.read(4 * 1024**2):
                digest.update(block)
    row = dict(filename=name, url=url, bytes=size, sha256=digest.hexdigest())
    print(json.dumps(row), flush=True)
    return row


if __name__ == "__main__":
    with concurrent.futures.ThreadPoolExecutor(max_workers=3) as pool:
        rows = list(pool.map(acquire, NAMES))
    (OUT / "reference_receipt.json").write_text(json.dumps({
        "scope": "Exact named references for existing F115W photons; no new science observation",
        "products": rows,
    }, indent=2) + "\n")
