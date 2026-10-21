"""Query the recorded context's IPC best reference without running a pipeline.

Only metadata from the already-consumed four F115W raw files is consulted.
CRDS cache changes are confined to the ignored native_cache directory.
"""
from pathlib import Path
import hashlib
import json
import os

ROOT = Path(__file__).resolve().parents[1]
BASE = ROOT / "data/strong_lensing_response_completion_20261021"
OUT = ROOT / "artifacts/strong_lensing_response_completion_20261021/ipc"
CACHE = BASE / "native_cache/crds_ipc"
os.environ["CRDS_PATH"] = str(CACHE)
os.environ["CRDS_SERVER_URL"] = "https://jwst-crds.stsci.edu"
os.environ["CRDS_MODE"] = "remote"
os.environ["CRDS_CONTEXT"] = "jwst_1535.pmap"


def sha(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def main():
    import crds
    from jwst.stpipe.core import JwstStep
    from astropy.io import fits
    OUT.mkdir(parents=True, exist_ok=True)
    CACHE.mkdir(parents=True, exist_ok=True)
    rows = []
    for exposure in range(1, 5):
        path = BASE / f"native_cache/jw01198007001_02101_0000{exposure}_nrcb3_uncal.fits"
        parameters, observatory = JwstStep._get_crds_parameters(path)
        selected = {key: parameters.get(key) for key in (
            "meta.instrument.name", "meta.instrument.detector", "meta.observation.date",
            "meta.observation.time", "meta.exposure.type", "meta.instrument.filter")}
        print("IPC_QUERY", path.name, selected, flush=True)
        recommendation = crds.getrecommendations(parameters, reftypes=["ipc"],
            context="jwst_1535.pmap", observatory=observatory)
        print("IPC_RECOMMENDATION", recommendation, flush=True)
        row = {"science_file": str(path.relative_to(ROOT)).replace("\\", "/"),
            "science_sha256": sha(path), "context": "jwst_1535.pmap",
            "parameters": parameters, "recommendation": recommendation,
            "metadata_only_science_access": True}
        # Direct public reference retrieval avoids CRDS's Windows260-character
        # instrument-cache-directory limitation; selection is still CRDS's.
        import urllib.request
        refs = {}
        for kind, filename in recommendation.items():
            if filename == "N/A":
                refs[kind] = filename
                continue
            destination = BASE / "native_cache" / Path(filename).name
            if not destination.exists():
                url = "https://jwst-crds.stsci.edu/unchecked_get/references/jwst/" + Path(filename).name
                with urllib.request.urlopen(url, timeout=180) as response:
                    destination.write_bytes(response.read())
            refs[kind] = str(destination)
        row["cached_reference"] = refs
        for kind, filename in refs.items():
            if filename == "N/A":
                continue
            ref = Path(filename)
            row["reference_sha256"] = sha(ref)
            row["reference_bytes"] = ref.stat().st_size
            with fits.open(ref, memmap=False) as hdul:
                row["reference_primary_header_cards"] = [str(c) for c in hdul[0].header.cards]
                row["reference_extensions"] = [
                    {"name": h.name, "shape": list(h.data.shape) if h.data is not None else None}
                    for h in hdul]
        rows.append(row)
        (OUT / "archived_bestref_query.json").write_text(json.dumps({
            "crds_version": crds.__version__, "query_context": "jwst_1535.pmap",
            "records": rows, "no_pipeline_execution": True,
            "no_reserved_observations": True, "script_sha256": sha(__file__)},
            indent=2, default=str, allow_nan=False) + "\n", encoding="utf-8")
        print("IPC_RESULT", recommendation, refs, flush=True)


if __name__ == "__main__":
    main()
