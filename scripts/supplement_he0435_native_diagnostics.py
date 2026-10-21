"""Archive existing official segment/intercept and reference-control diagnostics.

This reads already computed files only. It never reruns detector processing or
changes the immutable compact history files.
"""
from pathlib import Path
import hashlib,json
import numpy as np
from astropy.io import fits
ROOT=Path(__file__).resolve().parents[1]
BASE=ROOT/"data/strong_lensing_response_completion_20261021"
CACHE=BASE/"native_cache";OUT=ROOT/"artifacts/strong_lensing_response_completion_20261021"


def sha(p):
    h=hashlib.sha256()
    with p.open("rb") as f:
        for b in iter(lambda:f.read(8*1024**2),b""):h.update(b)
    return h.hexdigest()


def ref_values(name,ext,x,y):
    with fits.open(CACHE/name,memmap=False) as h:return np.array(h[ext].data[y,x])


def main():
    rows=[]
    for exp in range(1,5):
        stem=f"jw01198007001_02101_0000{exp}_nrcb3"
        original=BASE/"histories"/f"{stem}_histories.npz"
        output=BASE/"histories"/f"{stem}_native_diagnostics.npz"
        receipt=OUT/f"{stem}_native_diagnostics.json"
        if output.exists() and receipt.exists():rows.append(json.loads(receipt.read_text()));continue
        base=np.load(original);bounds=base["bounds_xyxy"];x,y=base["reference_control_native_xy"].T
        arrays={"labels":base["labels"],"bounds_xyxy":bounds,"reference_control_native_xy":base["reference_control_native_xy"]}
        # stpipe appends its fitopt suffix and uses output_dir for location.
        filename=f"{stem}_fitopt_fitopt.fits";old=ROOT/filename;opt=CACHE/filename
        if old.exists():
            if opt.exists():raise RuntimeError("refuse to overwrite existing optional archive")
            old.rename(opt)
        if not opt.exists():raise FileNotFoundError(opt)
        with fits.open(opt,memmap=False) as h:
            for hd in h[1:]:
                if hd.name=="ASDF" or hd.header.get("NAXIS",0)<2:continue
                if hd.header.get("NAXIS1")==2048 and hd.header.get("NAXIS2")==2048:
                    arrays["native_optional_"+hd.name]=np.stack([np.array(hd.section[...,y0:y1,x0:x1]) for x0,y0,x1,y1 in bounds])
        for key,name,ext in [("reference_control_mask_DQ","jwst_nircam_mask_0075.fits","DQ"),
            ("reference_control_saturation_DQ","jwst_nircam_saturation_0100.fits","DQ"),
            ("reference_control_superbias_DQ","jwst_nircam_superbias_0180.fits","DQ"),
            ("reference_control_readnoise_CDS_DN","jwst_nircam_readnoise_0282.fits","SCI"),
            ("reference_control_gain_e_per_DN","jwst_nircam_gain_0094.fits","SCI"),
            ("reference_control_superbias_ERR_DN","jwst_nircam_superbias_0180.fits","ERR")]:
            arrays[key]=ref_values(name,ext,x,y)
        with fits.open(CACHE/f"{stem}_uncal.fits",memmap=False) as h:
            arrays["reference_control_raw_pixeldq"]=np.array(h["PIXELDQ"].data[y,x]) if "PIXELDQ" in h else np.zeros(len(x),np.uint32)
        pre=arrays["reference_control_raw_pixeldq"].astype(np.uint32)
        for key in ["reference_control_mask_DQ","reference_control_saturation_DQ","reference_control_superbias_DQ"]:pre|=arrays[key].astype(np.uint32)
        arrays["reference_control_pre_sirs_pixeldq"]=pre
        with fits.open(CACHE/f"{stem}_corrected_ramp.fits",memmap=False) as h:
            arrays["reference_control_final_pixeldq"]=np.array(h["PIXELDQ"].data[y,x])
        j1=base["native_rate_jacobian_group_and_zeroframe_h1_per_s"][:,:,24:41,24:41]
        j2=base["native_rate_jacobian_group_and_zeroframe_h2_per_s"][:,:,24:41,24:41]
        np.savez_compressed(output,**arrays)
        row={"status":"EXISTING_NATIVE_SEGMENT_AND_REFERENCE_CONTROL_SUPPLEMENT","file":str(output.relative_to(ROOT)),"sha256":sha(output),
            "main_history_file":str(original.relative_to(ROOT)),"main_history_sha256":sha(original),
            "official_optional_file":str(opt.relative_to(ROOT)),"official_optional_sha256":sha(opt),
            "arrays":{k:list(v.shape) for k,v in arrays.items()},
            "pre_sirs_control_DNU_disagrees_final_count":int(np.sum((pre&1)!=(arrays["reference_control_final_pixeldq"]&1))),
            "readnoise_reference_control_DN_min_max":[float(np.min(arrays["reference_control_readnoise_CDS_DN"])),float(np.max(arrays["reference_control_readnoise_CDS_DN"]))],
            "h1_h2_likelihood17_max_abs_J_difference_per_s":float(np.max(abs(j1-j2))),
            "optional_intercept_scope":"Official post-correction segment YINT/SIGYINT retained with native segment axes. Not an exact noiseless pre-linearity reset state or its full covariance.",
            "native_processing_rerun":False,"main_history_overwritten":False,"physical_admission":False}
        receipt.write_text(json.dumps(row,indent=2)+"\n");rows.append(row)
        print(stem,'SUPPLEMENT_READY',output.stat().st_size,row['readnoise_reference_control_DN_min_max'],flush=True)
    (OUT/"native_diagnostic_supplements.json").write_text(json.dumps({"records":rows},indent=2)+"\n")


if __name__=="__main__":main()
