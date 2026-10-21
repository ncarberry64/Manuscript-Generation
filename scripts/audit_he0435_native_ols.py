"""Compare piecewise source equations to native C on saved consumed histories.

This executes compact OLS_C value comparisons only, without native detector
processing, a numerical fit, or opening any unused observation.
"""
from pathlib import Path
import hashlib,json,sys
import numpy as np
from astropy.io import fits
from stcal.ramp_fitting import ramp_fit
from stcal.ramp_fitting.ramp_fit_class import RampData
from stdatamodels.jwst.datamodels import dqflags
ROOT=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT/"code"))
from strong_lensing_native_ols import piecewise_ols_arrays

BASE=ROOT/"data/strong_lensing_response_completion_20261021"
OUT=ROOT/"artifacts/strong_lensing_response_completion_20261021/native_ols"


def sha(p): return hashlib.sha256(Path(p).read_bytes()).hexdigest()


def c_value(groups,zero,gdq,pdq,dark,gain,rn,header):
    r=RampData()
    r.set_arrays(np.asarray(groups[None,:,None],np.float32).copy(),
        np.asarray(gdq[None,:,None],np.uint8).copy(),np.asarray(pdq[None],np.uint32).copy(),
        np.asarray(dark[None],np.float32).copy(),zeroframe=np.asarray(zero[None,None],np.float32).copy())
    r.set_meta("NIRCAM",header["TFRAME"],header["TGROUP"],header["GROUPGAP"],header["NFRAMES"],0)
    r.set_dqflags(dqflags.pixel);r.start_row=0;r.num_rows=1
    result,_,_=ramp_fit.ramp_fit_data(r,False,np.asarray(rn[None],np.float32).copy(),
        np.asarray(gain[None],np.float32).copy(),"OLS_C","optimal","none")
    return result["slope"].ravel()


def main():
    OUT.mkdir(exist_ok=True,parents=True);rows=[];records={}
    for e in range(1,5):
        stem=f"jw01198007001_02101_0000{e}_nrcb3"
        hp=BASE/f"histories/{stem}_histories.npz"
        header=fits.getheader(BASE/f"native_cache/{stem}_uncal.fits")
        with np.load(hp) as h:
            for p,label in enumerate(h["labels"]):
                roi=np.s_[24:41,24:41];gain=h["gain_e_per_DN"][p][roi].ravel();rn=h["readnoise_CDS_DN"][p][roi].ravel()
                dark=h["average_dark_current_DN_s"][p][roi].ravel();pdq=h["pixeldq"][p][roi].ravel()
                gdq=h["groupdq"][p,:,24:41,24:41].reshape(9,-1)
                mp=ROOT/f"artifacts/strong_lensing_response_completion_20261021/bounded_native/training_plus_held/{stem}_{label}_storage_ceiling_cap1_k0.npz"
                with np.load(mp) as m:
                    cases={"observed":(h["corrected_group_DN"][p,:,24:41,24:41].reshape(9,-1),h["native_corrected_zeroframe_DN"][p][roi].ravel()),
                        "canonical_conditional_model":(m["native_groups_DN"][:,24:41,24:41].reshape(9,-1),m["native_zeroframe_DN"][roi].ravel())}
                for kind,(y,z) in cases.items():
                    value,jac,branches=piecewise_ols_arrays(y,gdq,gain=gain,readnoise_cds_dn=rn,
                        dark_dn_s=dark,pixeldq=pdq,zeroframe=z,group_time=header["TGROUP"],
                        frame_time=header["TFRAME"],nframes=header["NFRAMES"],cast_input=True)
                    native=c_value(y,z,gdq,pdq,dark,gain,rn,header)
                    finite=np.isfinite(value)&np.isfinite(native)
                    delta=abs(value[finite]-native[finite]);rounded=np.float32(value[finite])
                    key=f"{stem}_{label}_{kind}"
                    records[key+"_jacobian_group_and_zeroframe_per_s"]=jac
                    records[key+"_continuous_source_slope_DN_s"]=value
                    records[key+"_native_C_slope_DN_s"]=native
                    rows.append({"exposure":stem,"label":str(label),"case":kind,"finite_pixels":int(finite.sum()),
                        "finite_pattern_equal":bool(np.array_equal(np.isfinite(value),np.isfinite(native))),
                        "max_abs_double_source_vs_published_C_DN_s":float(delta.max(initial=0.)),
                        "max_abs_float32_source_vs_C_DN_s":float(abs(rounded-native[finite]).max(initial=0.)),
                        "bitwise_equal_float32_values":int(np.sum(rounded==native[finite])),
                        "SNR_boundary_pixels":sum(bool(b["snr_boundary"]) for b in branches),
                        "median_tie_pixels":sum(bool(b["median_tie"]) for b in branches),
                        "nonunique_local_derivative_pixels":sum(b.get("locally_unique_derivative") is False for b in branches),
                        "nonzero_median_variance_weight_derivative_pixels":sum(b.get("median_weight_derivative_coefficient",0.)!=0. for b in branches)})
                print("OLS_SOURCE_COMPARED",stem,label,flush=True)
    array=OUT/"piecewise_source_native_comparison.npz";np.savez_compressed(array,**records)
    result={"status":"actual saved observed and conditional-model source equations compared to original native C",
        "rows":rows,"array_file":str(array.relative_to(ROOT)).replace("\\","/"),"array_sha256":sha(array),
        "native_value_not_replaced":True,"native_pipeline_rerun":False,"model_fit":False,
        "continuous_extension_not_quantization_derivative":True,
        "derivative_relation":"dR=sum(alpha_s W_s)+dmedian*sum(alpha_s a_s(R-s_s)/v_s), v_s=a_s(max(median,0)+dark)+varread_s; fixed native segments and SNR-power branches",
        "source_sha256":{str(p.relative_to(ROOT)).replace("\\","/"):sha(p) for p in (
            Path(__file__),ROOT/"code/strong_lensing_native_ols.py",ROOT/"tests/test_strong_lensing_native_ols.py")},
        "physical_validation":False,"physical_guard":2e-10}
    (OUT/"source_equation_comparison.json").write_text(json.dumps(result,indent=2,allow_nan=False)+"\n",encoding="utf-8")
    print("OLS_COMPARISON_DONE",max(r["max_abs_double_source_vs_published_C_DN_s"] for r in rows),flush=True)


if __name__=="__main__":main()
