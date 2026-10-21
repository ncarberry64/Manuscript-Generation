"""Independent raw-dark covariance source trail, with exact matched read support.

Only public detector calibration darks are acquired.  No HE0435 residual or
reserved science observation is used to select data or estimate this noise.
"""
from __future__ import annotations
from pathlib import Path
from concurrent.futures import ThreadPoolExecutor, as_completed
import argparse
import hashlib
import importlib.util
import json
import time
import numpy as np

ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT/"artifacts/strong_lensing_response_completion_20261021/dark_calibration"
DATA = ROOT/"data/strong_lensing_response_completion_20261021/dark_calibration"
RAW = ROOT/"data/strong_lensing_response_completion_20261021/native_cache/dark_calibration"
READS = np.array([[5*j, 5*j+1] for j in range(9)])
SOURCES = (
    "jw01453022001_02301_00001_nrcb3_uncal.fits",
    "jw01453020001_02201_00001_nrcb3_uncal.fits",
)


def sha(path):
    h = hashlib.sha256()
    with Path(path).open("rb") as f:
        for block in iter(lambda: f.read(8*1024**2), b""):
            h.update(block)
    return h.hexdigest()


def dump(path, value):
    Path(path).write_text(json.dumps(value, indent=2, allow_nan=False)+"\n",
                         encoding="utf-8", newline="\n")


def acquire():
    OUT.mkdir(parents=True, exist_ok=True)
    DATA.mkdir(parents=True, exist_ok=True)
    RAW.mkdir(parents=True, exist_ok=True)
    spec = importlib.util.spec_from_file_location("_native_range_acquisition",
                     ROOT/"scripts/acquire_he0435_target_histories.py")
    acquisition = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(acquisition)
    plan = {"selection": "metadata-only earliest and latest public program1453 NRCB3 FULL RAPID within actual0282 pedigree",
        "sources": list(SOURCES), "individual_RAPID_frame_indices_zero_based": READS.tolist(),
        "integrations_zero_based": [0, 1], "pixels": "all2048x2048 science and reference pixels, raw ADC before IPC/refpix/linearity",
        "program4442": "in-pedigree NRCB3 RAPID archive entries are subarrays; first FULL entry is outside pedigree",
        "independent_of_HE0435_residuals": True, "no_reserved_science_observations": True,
        "source_script_sha256_before_acquisition": sha(Path(__file__)),
        "intended_estimator": "raw spatial/temporal summaries and full science-reference cross terms; finite4-realization rank explicit; any stationary pooling must be separately stated"}
    plan_path = OUT/"dark_selection_before_pixels.json"
    if not plan_path.exists():
        dump(plan_path, plan)
    elif json.loads(plan_path.read_text())["sources"] != list(SOURCES):
        raise RuntimeError("Existing source selection differs; never silently replace it")
    jobs, records = [], []
    for source in SOURCES:
        url = "https://mast.stsci.edu/api/v0.1/Download/file?uri=mast:JWST/product/"+source
        layout, meta = acquisition.remote_layout(url)
        primary = layout[0][0]
        header, offset, size = next(row for row in layout if row[0].get("EXTNAME") == "SCI")
        expected = {"NAXIS1": 2048, "NAXIS2": 2048, "NAXIS3": 150, "NAXIS4": 2,
                    "BITPIX": 16, "BSCALE": 1, "BZERO": 32768}
        if any(header.get(k) != v for k, v in expected.items()):
            raise RuntimeError("Unexpected raw dark storage/read support")
        if any(primary.get(k) != v for k, v in {"DETECTOR": "NRCB3", "READPATT": "RAPID",
                    "NFRAMES": 1, "SUBARRAY": "FULL", "NINTS": 2, "NGROUPS": 150}.items()):
            raise RuntimeError("Source dark configuration does not match the selected metadata")
        record = {"source": source, "remote": meta,
            "primary_header": {k: primary.get(k) for k in ("DETECTOR", "READPATT", "NFRAMES", "GROUPGAP", "NGROUPS", "NINTS", "SUBARRAY", "TFRAME", "EXPSTART", "EXPEND", "CAL_VER", "CRDS_CTX")},
            "primary_header_cards": primary.tostring(sep="\n", padding=False),
            "SCI_header_cards": header.tostring(sep="\n", padding=False),
            "source_SCI_offset": offset, "source_SCI_bytes": size,
            "range_chunks": [], "subsets": []}
        for integration in (0, 1):
            target = RAW/(source.replace("_uncal.fits", f"_integration{integration}_matched_reads.npy"))
            temporary = RAW/(target.name+".partial.npy")
            subset = {"integration_zero_based": integration, "file": str(target.relative_to(ROOT)).replace('\\', '/')}
            record["subsets"].append(subset)
            if target.exists():
                old = np.load(target, mmap_mode="r")
                if old.shape != (18, 2048, 2048) or old.dtype != np.uint16:
                    raise RuntimeError("Existing exact raw subset has wrong shape/dtype")
                subset.update(bytes=target.stat().st_size, sha256=sha(target), reused=True)
                # The first transfer completed before a FITS-commentary JSON
                # serialization error. Recover exact source-range hashes from
                # reversible uint16 storage, without repeating any data transfer.
                for group in range(9):
                    raw = (np.asarray(old[2*group:2*group+2], np.int32)-32768).astype(">i2").tobytes()
                    record["range_chunks"].append({"integration_zero_based": integration,
                        "group_zero_based": group,
                        "first_source_byte": offset+(integration*150+int(READS[group, 0]))*2048*2048*2,
                        "bytes": len(raw), "raw_range_sha256": hashlib.sha256(raw).hexdigest(),
                        "ETag": meta.get("etag"), "range_digest_recovered_from_reversible_exact_uint16_subset": True})
                continue
            np.lib.format.open_memmap(temporary, mode="w+", dtype=np.uint16,
                                     shape=(18, 2048, 2048)).flush()
            for group in range(9):
                frame = int(READS[group, 0])
                first = offset+(integration*150+frame)*2048*2048*2
                nbytes = 2*2048*2048*2
                jobs.append((url, first, nbytes, temporary, group, subset, record))
        records.append(record)
    def fetch(job):
        url, first, nbytes, target, group, subset, record = job
        raw, headers = acquisition.get_range(url, first, first+nbytes-1)
        signed = np.frombuffer(raw, dtype=">i2").reshape(2, 2048, 2048)
        values = (signed.astype(np.int32)+32768).astype(np.uint16)
        mem = np.load(target, mmap_mode="r+")
        mem[2*group:2*group+2] = values
        mem.flush()
        return record, {"integration_zero_based": subset["integration_zero_based"],
            "group_zero_based": group, "first_source_byte": first, "bytes": nbytes,
            "raw_range_sha256": hashlib.sha256(raw).hexdigest(),
            "ETag": next((v for k, v in headers.items() if k.lower() == "etag"), None)}
    started = time.perf_counter()
    with ThreadPoolExecutor(max_workers=4) as pool:
        for index, future in enumerate(as_completed([pool.submit(fetch, job) for job in jobs]), 1):
            record, chunk = future.result()
            record["range_chunks"].append(chunk)
            print("EXACT_DARK_RANGES", index, "/", len(jobs), "bytes", index*16777216, flush=True)
    for record in records:
        tags = {x["ETag"] for x in record["range_chunks"] if x["ETag"] is not None}
        if len(tags) > 1:
            raise RuntimeError("Source ETag changed during raw range acquisition")
        for subset in record["subsets"]:
            target = ROOT/subset["file"]
            temporary = RAW/(target.name+".partial.npy")
            if not target.exists():
                temporary.replace(target)
                subset.update(bytes=target.stat().st_size, sha256=sha(target), reused=False)
        record["range_chunks"].sort(key=lambda x: (x["integration_zero_based"], x["group_zero_based"]))
    report = {"status": "exact public raw dark matched-read subsets acquired",
        "records": records, "transferred_pixel_bytes": len(jobs)*16777216,
        "seconds": time.perf_counter()-started, "source_selection_sha256": sha(plan_path),
        "full_5GB_raw_files_downloaded": False, "no_science_reserved_data": True}
    dump(OUT/"raw_dark_acquisition.json", report)
    return report


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--acquire-independent-dark-subsets", action="store_true")
    args = parser.parse_args()
    if not args.acquire_independent_dark_subsets:
        parser.error("Choose the explicit calibration-data acquisition stage")
    acquire()


if __name__ == "__main__":
    main()
