"""Apply the conditional direct-group diagnostic to the delivered consumed field."""
from __future__ import annotations
import argparse
import hashlib
import json
from pathlib import Path
import platform
import sys
import time
import numpy as np
import scipy
from astropy.io import fits

ROOT=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT/"code"))
from strong_lensing_read_history import fit_history


def main():
    parser=argparse.ArgumentParser()
    parser.add_argument("--input",type=Path,default=ROOT/"data/strong_lensing_response_repair_20261021/delivery/ramp_diagnostic/native_pipeline/native_corrected_field_ramp.npz")
    parser.add_argument("--output",type=Path,default=ROOT/"artifacts/strong_lensing_response_repair_20261021/read_history")
    args=parser.parse_args()
    args.output.mkdir(parents=True,exist_ok=True)
    start=time.perf_counter()
    source=np.load(args.input,allow_pickle=False)
    kwargs=dict(groups=source["native_corrected_groups"],groupdq=source["native_groupdq"],
                pixeldq=source["native_pixeldq"],gain=source["gain_e_per_DN"],
                readnoise_cds=source["readnoise_CDS_DN"],read_times=source["native_read_times_s"])
    result={"input":str(args.input.relative_to(ROOT)),"input_sha256":hashlib.sha256(args.input.read_bytes()).hexdigest(),
            "runtime":{"python":platform.python_version(),"numpy":np.__version__,"scipy":scipy.__version__},
            "new_native_pipeline_processing":False,"reserved_observation_access":False,
            "new_independent_observations":False,
            "array_schema":{"train_records":["flat index in 66x66 context","zero-based group index","zero-based segment index","measured DN","fitted group mean DN","raw residual DN","sqrt frozen covariance diagonal DN"],
                "held_records":["flat index in 66x66 context","measured ninth-group DN","conditional forecast DN","conditional raw residual DN","unconditioned mean forecast DN","shared-read conditional correction DN"],
                "parameter_covariance":"Unscaled inverse profile information in parameter_coordinates; all segment intercepts profiled, with their conditional uncertainty restored in held predictive covariance",
                "held_predictive_covariance":"Conditional read innovation plus profiled-intercept and full joint optical/rate/kappa covariance; no Q/dof multiplier"},
            "source_sha256":{str(p.relative_to(ROOT)):hashlib.sha256(p.read_bytes()).hexdigest() for p in [ROOT/"code/strong_lensing_read_history.py",Path(__file__)]},
            "equations":{"read_support":"t_jr=(5*j+r)*10.73677 s, r=1,2; exact delivered times used",
                         "mean":"Q_pj=a_p,segment + m_j f_p + kappa v_j/2 L(f^2)_p; m=mean(t), v=mean(t^2)",
                         "leading_mean_hypothesis":"dot Q_p=f_p+kappa*t*sum_neighbours(f_q^2-f_p^2); antisymmetric edge flux, not a unique detector derivation",
                         "covariance":"C_jk=max(f_train,0)/gain * mean_rs min(t_jr,t_ks) + delta_jk READNOISE_CDS^2/4",
                         "holdout":"mu_h+C_ht C_tt^-1(y_t-mu_t); V=C_hh-C_ht C_tt^-1 C_th + nuisance and joint fitted-parameter propagation",
                         "slope_relation":"For normalized supplied weights a_j: E[r_p]=f_p+kappa/2 sum_j a_j v_j L(f^2)_p. Exact archived OLS_C weights are not reconstructed."},
            "meaning":"Leading positive, charge-conserving stationary-source mean hypothesis; conditional field diagnostic only."}
    delivered=ROOT/"data/strong_lensing_response_repair_20261021/delivery"
    psf_path=delivered/"stpsf_equal_sampling/jw01198007001_02101_00001_nrcb3_FIELD_stpsf_epsf.fits"
    static_receipt=json.loads((delivered/"stpsf_equal_sampling/field_trial.json").read_text(encoding="utf-8"))
    init=static_receipt["records"][0]
    cal_path=ROOT/"data/strong_lensing_geometry_20261021/native/jw01198007001_02101_00001_nrcb3_cal.fits"
    # Header-only access: all science values in this application come exclusively
    # from the consumed delivered field stamp, never from quasar/reserved arrays.
    photmjsr=float(fits.getheader(cal_path,"SCI")["PHOTMJSR"])
    psf=np.asarray(fits.getdata(psf_path),float)
    optical=dict(psf_grid=psf,oversampling=5,flat=source["flat"],photmjsr=photmjsr,
                 initial_center_xy=np.asarray(init["centroid_xy_pixel"])-[867,843],
                 initial_flux_SCI=init["flux_in_SCI_pixel_units"],
                 initial_background_SCI=init["background_plane_SCI"])
    result["fixed_optical_input"]={"PSF":str(psf_path.relative_to(ROOT)),"PSF_sha256":hashlib.sha256(psf_path.read_bytes()).hexdigest(),
        "oversampling":5,"origin_grid_xy":[(psf.shape[1]-1)/2,(psf.shape[0]-1)/2],
        "PHOTMJSR":photmjsr,"PHOTMJSR_source":"SCI header only of the already archived native CAL",
        "spectrum":"Prescribed 4000 K blackbody; no later shared-temperature/diffusion or phase optimization reused",
        "wavefront":"Delivered dated OPD R2022101203-NRCA3_FP1-1.fits; fixed STPSF2.2.0",
        "initial_parameters":"Earlier consumed static fit supplies optimizer initial values only",
        "pixel_integration":"Delivered oversampled ePSF is already integrated; sampled once with cubic ImagePSF-convention spline"}
    for name,enabled,optic in [("constant_response",False,None),("history_dependent",True,None),
                               ("optical_constant_response",False,optical),("optical_history_dependent",True,optical)]:
        print(f"Fitting {name}",flush=True)
        summary,arrays=fit_history(**kwargs,redistribution=enabled,optical=optic)
        result[name]=summary
        np.savez_compressed(args.output/f"{name}_arrays.npz",**arrays)
        (args.output/f"{name}.json").write_text(json.dumps(summary,indent=2)+"\n",encoding="utf-8")
        print(json.dumps({k:summary[k] for k in ["optimizer_success","nfev","training_Q","training_Q_per_dof","kappa","kappa_local_curvature_standard_error","held_Q","held_Q_per_dof","held_raw_RMS_DN"]}),flush=True)
    old,new=result["constant_response"],result["history_dependent"]
    result["comparison"]={"training_Q_reduction":old["training_Q"]-new["training_Q"],
        "held_Q_reduction":old["held_Q"]-new["held_Q"],
        "held_raw_RMS_fractional_reduction":1-new["held_raw_RMS_DN"]/old["held_raw_RMS_DN"],
        "physical_calibration_pass":False,"lensing_validation_pass":False,
        "notes":"No error rescaling. Training and held groups share Poisson reads; conditional innovations used. Working covariance remains incomplete."}
    old,new=result["optical_constant_response"],result["optical_history_dependent"]
    result["optical_comparison"]={"training_Q_reduction":old["training_Q"]-new["training_Q"],
        "held_Q_reduction":old["held_Q"]-new["held_Q"],
        "held_raw_RMS_fractional_reduction":1-new["held_raw_RMS_DN"]/old["held_raw_RMS_DN"],
        "physical_calibration_pass":False,"lensing_validation_pass":False,
        "meaning":"Same fixed optical PSF coupled to the temporal charge mean; absolute group residuals remain the adequacy test."}
    result["elapsed_seconds"]=time.perf_counter()-start
    (args.output/"results.json").write_text(json.dumps(result,indent=2)+"\n",encoding="utf-8")
    print(json.dumps(result["comparison"],indent=2),flush=True)


if __name__=="__main__":main()
