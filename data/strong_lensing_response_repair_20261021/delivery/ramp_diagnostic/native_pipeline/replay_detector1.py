"""Replay recorded detector1 settings on the same consumed F115W exposure.

Uses public, released jwst2.0.1 steps and the exact named reference files.
No lens fit or previously unused observation enters this calculation.
"""
from pathlib import Path
import os
import json
import time
import hashlib
import importlib.metadata

OUT = Path(__file__).resolve().parent
ROOT = OUT.parents[1]
RAW = OUT.parent / "jw01198007001_02101_00001_nrcb3_uncal.fits"
CAL = ROOT / "repo/data/strong_lensing_geometry_20261021/native" / RAW.name.replace("_uncal", "_cal")
os.environ["CRDS_PATH"] = str(OUT / "crds_cache")
os.environ["CRDS_CONTEXT"] = "jwst_1535.pmap"
os.environ["CRDS_SERVER_URL"] = "https://jwst-crds.stsci.edu"
os.environ["OMP_NUM_THREADS"] = "1"
os.environ["OPENBLAS_NUM_THREADS"] = "1"
os.environ["MPLCONFIGDIR"] = str(OUT / "matplotlib_cache")

import numpy as np
import yaml
from astropy.io import fits
from stdatamodels.jwst import datamodels
from jwst.pipeline import Detector1Pipeline


class TaggedLoader(yaml.SafeLoader):
    pass


def tagged(loader, suffix, node):
    if isinstance(node, yaml.MappingNode):
        return loader.construct_mapping(node)
    if isinstance(node, yaml.SequenceNode):
        return loader.construct_sequence(node)
    return loader.construct_scalar(node)


TaggedLoader.add_multi_constructor("", tagged)


def configuration():
    with fits.open(CAL) as hd:
        text = hd["ASDF"].data["ASDF_METADATA"][0].tobytes().decode("latin1")
    tree = yaml.load(text[:text.index("\n...\n") + 5], Loader=TaggedLoader)
    logs = tree["cal_logs"]["calwebb_detector1"]
    setting = next(line for line in logs if "Step Detector1Pipeline parameters are:" in line)
    config = yaml.safe_load(setting.split("parameters are:\n", 1)[1])
    (OUT / "native_recorded_detector1_config.json").write_text(json.dumps(config, indent=2) + "\n")
    return config


def configured_pipeline():
    config = configuration()
    framework = {"pre_hooks", "post_hooks", "output_file", "output_dir", "output_ext",
                 "output_use_model", "output_use_index", "save_results", "suffix",
                 "search_output_file", "input_dir"}
    common = OUT.parent
    refs = OUT / "refs"
    override = {
        "dq_init": {"mask": refs / "jwst_nircam_mask_0075.fits"},
        "saturation": {"saturation": common / "jwst_nircam_saturation_0100.fits",
                       "superbias": common / "jwst_nircam_superbias_0180.fits"},
        "superbias": {"superbias": common / "jwst_nircam_superbias_0180.fits"},
        "refpix": {"refpix": "N/A", "sirskernel": refs / "jwst_nircam_sirskernel_0002.asdf"},
        "linearity": {"linearity": common / "jwst_nircam_linearity_0054.fits"},
        "dark_current": {"dark": refs / "derived_jwst_nircam_dark_0337_first42.fits"},
        "jump": {"gain": refs / "jwst_nircam_gain_0094.fits",
                 "readnoise": refs / "jwst_nircam_readnoise_0282.fits"},
        "ramp_fit": {"gain": refs / "jwst_nircam_gain_0094.fits",
                     "readnoise": refs / "jwst_nircam_readnoise_0282.fits"},
        "gain_scale": {"gain": refs / "jwst_nircam_gain_0094.fits"},
    }
    steps = {}
    for name, values in config["steps"].items():
        steps[name] = {key: (None if value == "None" else value)
                       for key, value in values.items() if key not in framework}
        steps[name].update({"override_" + kind: str(path)
                            for kind, path in override.get(name, {}).items()})
    return Detector1Pipeline(steps=steps), steps


def main():
    started = time.time()
    pipeline, steps = configured_pipeline()
    model = datamodels.RampModel(RAW)
    operations = ["group_scale", "dq_init", "saturation", "ipc", "superbias", "refpix",
                  "linearity", "persistence", "dark_current", "charge_migration", "jump",
                  "picture_frame", "clean_flicker_noise"]
    for name in operations:
        model = getattr(pipeline, name).run(model)
        print("COMPLETED_NATIVE_STEP", name, time.time() - started, flush=True)
    pre = OUT / "first_exposure_native_corrected_ramp.fits"
    model.save(pre)
    print("SAVED_PRE_RAMP", str(pre), time.time() - started, flush=True)
    rate, rateints = pipeline.ramp_fit.run(model)
    rate = pipeline.gain_scale.run(rate)
    rate.save(OUT / "first_exposure_native_rate.fits")
    if rateints is not None:
        rateints.close()
    model.close()
    rate.close()
    record = {"status": "SAME_EXPOSURE_NATIVE_DETECTOR_REPLAY",
              "jwst_version": importlib.metadata.version("jwst"),
              "stcal_version": importlib.metadata.version("stcal"),
              "native_stcal_version_recorded": False,
              "reference_context_recorded": "jwst_1535.pmap",
              "science_parameters": steps,
              "input": str(RAW), "output_ramp": str(pre),
              "elapsed_seconds": time.time() - started,
              "new_observations_accessed": False,
              "physical_lens_parameters_fitted": False}
    (OUT / "replay_receipt.json").write_text(json.dumps(record, indent=2) + "\n")
    print("NATIVE_REPLAY_COMPLETE", record["elapsed_seconds"], flush=True)


if __name__ == "__main__":
    main()
