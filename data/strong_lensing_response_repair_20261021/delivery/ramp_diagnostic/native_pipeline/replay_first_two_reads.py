"""Calibrate reconstructed first and second reads using official detector steps.

This is a consumed-data development experiment, not an original native product.
Raw read2 = 2 * first averaged group - ZEROFRAME, before nonlinear correction.
On-board integer averaging prevents claiming bit-exact recovery of raw read2.
"""
from pathlib import Path
import json
import time
import numpy as np
from astropy.io import fits
from stdatamodels.jwst import datamodels
from replay_detector1 import OUT, RAW, ROOT, configured_pipeline


def main():
    started = time.time()
    pipeline, steps = configured_pipeline()
    pipeline.dark_current.override_dark = str(OUT / "refs/derived_jwst_nircam_dark_0337_first2.fits")
    steps["dark_current"]["override_dark"] = pipeline.dark_current.override_dark
    model = datamodels.RampModel(RAW)
    assert model.meta.exposure.nframes == 2
    assert model.meta.exposure.ngroups == 9
    assert model.meta.exposure.nints == 1
    assert model.meta.exposure.zero_frame
    raw_zero = np.array(model.zeroframe, dtype=np.float32)
    raw_second = 2.0 * np.array(model.data[:, 0], dtype=np.float32) - raw_zero
    model.data = np.stack([raw_zero, raw_second], axis=1)
    model.groupdq = np.zeros(model.data.shape, dtype=np.uint8)
    # ZEROFRAME zeros are invalid sentinels, not count measurements.
    bad_reconstruction = (raw_zero == 0) | (raw_second < 0) | (raw_second > 65535)
    model.groupdq[:, 0][bad_reconstruction] |= 1
    model.groupdq[:, 1][bad_reconstruction] |= 1
    interval = float(model.meta.exposure.frame_time)
    model.meta.exposure.nframes = 1
    model.meta.exposure.frame_divisor = 1
    model.meta.exposure.ngroups = 2
    model.meta.exposure.groupgap = 0
    model.meta.exposure.group_time = interval
    model.meta.exposure.readpatt = "RAPID"
    model.meta.exposure.zero_frame = False
    model.meta.exposure.integration_time = 2 * interval
    model.meta.exposure.exposure_time = 2 * interval
    operations = ["group_scale", "dq_init", "saturation", "ipc", "superbias", "refpix", "linearity", "persistence", "dark_current"]
    field = np.s_[868:885, 892:909]
    checkpoints = {}
    for name in operations:
        model = getattr(pipeline, name).run(model)
        if name in ["refpix", "linearity", "dark_current"]:
            checkpoints[name + "_field_reads"] = np.array(model.data[0, :, field[0], field[1]])
        print("COMPLETED_TWO_READ_STEP", name, time.time() - started, flush=True)
    # Preserve the original science-ramp cosmic-ray/static exclusions at group0.
    native = OUT / "first_exposure_native_corrected_ramp.fits"
    with fits.open(native, memmap=False) as hd:
        original_first_dq = np.asarray(hd["GROUPDQ"].data[:, 0], dtype=np.uint32)
        first_bad = (original_first_dq & (1 | 4)) != 0
    model.groupdq[:, 0][first_bad] |= 1
    model.groupdq[:, 1][first_bad] |= 1
    destination = OUT / "first_two_reads_corrected_ramp.fits"
    model.save(destination)
    np.savez_compressed(OUT / "first_two_reads_field_checkpoints.npz", **checkpoints)
    record = {
        "status": "RECONSTRUCTED_COMMON_TIME_READ_PAIR_DEVELOPMENT",
        "input": str(RAW), "output": str(destination),
        "construction": "Read1=raw ZEROFRAME; read2=2*raw group0-read1, before all calibrations",
        "limitations": ["Read2 contains on-board group-averaging rounding, not bit-exact original ADC samples", "No new independent observation or lens fit", "Pixel-noise assumptions and spatial covariance require separate audit"],
        "interval_s": interval,
        "synthetic_read_pattern": {"nframes": 1, "frame_divisor": 1, "ngroups": 2, "groupgap": 0, "zero_frame": False},
        "steps": operations, "reference_overrides": steps,
        "original_first_group_DNU_or_JUMP_excluded": int(first_bad.sum()),
        "raw_invalid_reconstruction_pixel_count": int(bad_reconstruction.sum()),
        "elapsed_seconds": time.time() - started,
    }
    (OUT / "first_two_reads_receipt.json").write_text(json.dumps(record, indent=2) + "\n")
    model.close()
    print("TWO_READ_REPLAY_COMPLETE", record["elapsed_seconds"], flush=True)


if __name__ == "__main__":
    main()
