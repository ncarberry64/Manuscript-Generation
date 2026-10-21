"""Archive a small, exact field-only stamp for continued read-history research."""
from pathlib import Path
import json
import hashlib
import numpy as np
from astropy.io import fits

OUT = Path(__file__).resolve().parent
ROOT = OUT.parents[1]
SL = np.s_[843:909, 867:933]


def collect(path, definitions, arrays):
    with fits.open(path, memmap=False) as hd:
        for name, key, leading in definitions:
            index = (*leading, SL[0], SL[1])
            arrays[key] = np.array(hd[name].section[index])


def main():
    arrays = {}
    paths = [OUT / "first_exposure_native_corrected_ramp.fits",
             OUT / "first_two_reads_corrected_ramp.fits",
             ROOT / "repo/data/strong_lensing_geometry_20261021/native/jw01198007001_02101_00001_nrcb3_cal.fits"]
    collect(paths[0], [("SCI", "native_corrected_groups", (0, slice(None))),
                      ("GROUPDQ", "native_groupdq", (0, slice(None))),
                      ("PIXELDQ", "native_pixeldq", ()),
                      ("ZEROFRAME", "native_corrected_zeroframe", (0,))], arrays)
    collect(paths[1], [("SCI", "reconstructed_corrected_reads", (0, slice(None))),
                      ("GROUPDQ", "reconstructed_groupdq", (0, slice(None))),
                      ("PIXELDQ", "reconstructed_pixeldq", ())], arrays)
    collect(paths[2], [("SCI", "archived_cal_MJy_sr", ()), ("ERR", "archived_cal_ERR", ()),
                      ("DQ", "archived_cal_DQ", ())], arrays)
    for filename, prefix in [("first_exposure_native_rate.fits", "rate_stcal1191"),
                              ("first_exposure_stcal118_rate.fits", "rate_stcal118")]:
        collect(OUT / filename, [("SCI", prefix, ()), ("ERR", prefix + "_ERR", ()),
                                 ("DQ", prefix + "_DQ", ())], arrays)
    for filename, definitions in [
        ("jwst_nircam_gain_0094.fits", [("SCI", "gain_e_per_DN", ())]),
        ("jwst_nircam_readnoise_0282.fits", [("SCI", "readnoise_CDS_DN", ())]),
        ("jwst_nircam_flat_0667.fits", [("SCI", "flat", ()), ("ERR", "flat_ERR", ()), ("DQ", "flat_DQ", ())]),
        ("derived_jwst_nircam_dark_0337_first2.fits", [("SCI", "dark_first2", (slice(None),)), ("ERR", "dark_first2_ERR", (slice(None),))]),
    ]:
        collect(OUT / "refs" / filename, definitions, arrays)
    arrays["x_pixel"], arrays["y_pixel"] = np.meshgrid(np.arange(867, 933), np.arange(843, 909))
    arrays["frame_time_s"] = 10.73677
    arrays["native_group_time_header_s"] = 53.684
    arrays["native_read_times_s"] = np.array([[5*j+1, 5*j+2] for j in range(9)]) * 10.73677
    output = OUT / "native_corrected_field_ramp.npz"
    np.savez_compressed(output, **arrays)
    receipt = {"status": "EXACT_CORRECTED_FIELD_STAMP",
               "detector_bounds_xyxy": [867, 843, 933, 909],
               "sources": [str(p) for p in paths],
               "arrays": {key: list(np.shape(value)) for key, value in arrays.items()},
               "file": output.name, "bytes": output.stat().st_size,
               "sha256": hashlib.sha256(output.read_bytes()).hexdigest(),
               "science_scope": "Already consumed isolated field star; no target fit or unused observation",
               "zero_frame_note": "Native zero frame includes saturation/superbias/SIRS/linearity; native dark step leaves it unchanged"}
    (OUT / "native_corrected_field_stamp_receipt.json").write_text(json.dumps(receipt, indent=2) + "\n")
    print(json.dumps(receipt, indent=2), flush=True)


if __name__ == "__main__":
    main()
