"""Exact temporal subsets of the official 187-frame RAPID dark reference.

The nine SHALLOW2 groups consume frames[0:2],[5:7],...,[40:42]. The
two reconstructed RAPID reads consume frames[0:2]. No pixel or uncertainty
value is approximated; omitted frames are never consumed by either replay.
"""
from pathlib import Path
import hashlib
import json
import numpy as np
from astropy.io import fits

OUT = Path(__file__).resolve().parent
SOURCE = OUT / "refs/jwst_nircam_dark_0337.fits"


def sha256(path):
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        while block := stream.read(8 * 1024**2):
            digest.update(block)
    return digest.hexdigest()


def main():
    source_hash = sha256(SOURCE)
    receipt = {"source_file": str(SOURCE), "source_sha256": source_hash,
               "purpose": "Exact used-frame subsets to bound memory, not an alternative dark model",
               "source_group_indices_used_by_SHALLOW2": [[5*j, 5*j+1] for j in range(9)],
               "products": []}
    with fits.open(SOURCE, memmap=False, lazy_load_hdus=True) as hd:
        assert hd[0].header["NGROUPS"] == 187
        assert hd[0].header["NFRAMES"] == 1
        assert hd[0].header["GROUPGAP"] == 0
        for count in [42, 2]:
            primary = hd[0].header.copy()
            primary["NGROUPS"] = count
            filename = f"derived_jwst_nircam_dark_0337_first{count}.fits"
            primary["FILENAME"] = filename
            primary.add_history(f"Exact first{count}frames of jwst_nircam_dark_0337.fits")
            primary.add_history("Original reference SHA256=" + source_hash)
            items = [fits.PrimaryHDU(header=primary)]
            for source_hdu in hd[1:]:
                if source_hdu.name == "ASDF":
                    # Reconstructed from identical FITS metadata; original ASDF
                    # dimensions describe187 frames and must not be carried over.
                    continue
                if source_hdu.name in ["SCI", "ERR"]:
                    array = np.array(source_hdu.section[:count, :, :])
                    items.append(fits.ImageHDU(data=array, header=source_hdu.header.copy(), name=source_hdu.name))
                else:
                    items.append(source_hdu.copy())
            path = OUT / "refs" / filename
            fits.HDUList(items).writeto(path, overwrite=True)
            # Re-read a small row at each retained frame; the arrays were sliced,
            # never averaged, interpolated or assigned replacement values.
            with fits.open(path, memmap=False) as subset:
                for name in ["SCI", "ERR"]:
                    assert np.array_equal(subset[name].section[:, 868, 892:909],
                                          hd[name].section[:count, 868, 892:909], equal_nan=True)
                assert np.array_equal(subset["DQ"].data, hd["DQ"].data)
            receipt["products"].append({"filename": filename, "frames": count,
                                        "bytes": path.stat().st_size,
                                        "sha256": sha256(path),
                                        "data_subset_check": True})
            del items
            print("EXACT_DARK_SUBSET", filename, flush=True)
    (OUT / "dark_subset_receipt.json").write_text(json.dumps(receipt, indent=2) + "\n")


if __name__ == "__main__":
    main()
