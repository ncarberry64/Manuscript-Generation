"""Resolve unrecorded native stcal version using unchanged pre-ramp data.

stcal1.19 changed segment-slope weights. The original1.19.1 replay and its
failed numerical gate are retained. This only reruns the official ramp fit
with the earlier compatible1.18.0 release; no PSF, noise or lens fit changes.
"""
import importlib.metadata
import json
from replay_detector1 import OUT, configured_pipeline
from stdatamodels.jwst import datamodels


def main():
    version = importlib.metadata.version("stcal")
    assert version == "1.18.0", version
    pipeline, steps = configured_pipeline()
    path = OUT / "first_exposure_native_corrected_ramp.fits"
    model = datamodels.RampModel(path)
    rate, integrations = pipeline.ramp_fit.run(model)
    rate = pipeline.gain_scale.run(rate)
    destination = OUT / "first_exposure_stcal118_rate.fits"
    rate.save(destination)
    record = {"scope": "Numerical provenance reconciliation of missing native stcal version",
              "input_corrected_ramp": str(path), "output_rate": str(destination),
              "jwst_version": importlib.metadata.version("jwst"), "stcal_version": version,
              "ramp_fit_settings": steps["ramp_fit"],
              "native_version_inferred_not_recorded": True,
              "original_stcal1191_replay_retained": True,
              "physical_models_changed": False, "errors_rescaled": False}
    (OUT / "reramp_stcal118_receipt.json").write_text(json.dumps(record, indent=2) + "\n")
    model.close()
    rate.close()
    if integrations is not None:
        integrations.close()
    print(json.dumps(record, indent=2), flush=True)


if __name__ == "__main__":
    main()
