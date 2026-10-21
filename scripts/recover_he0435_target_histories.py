"""Execute recorded native detector1 operations on four consumed F115W raws.

Full-detector reference correction precedes compact A/B/C/D/FIELD extraction.
The earlier completed field-only pipeline evidence is retained and compared,
not regenerated as a new scientific result. New processing is required for
the missing target histories and remaining three calibration histories.
"""
from pathlib import Path
import os,sys,json,time,copy,hashlib,importlib.metadata,argparse

ROOT=Path(__file__).resolve().parents[1]
BASE=ROOT/"data/strong_lensing_response_completion_20261021"
CACHE=BASE/"native_cache"
COMPACT=BASE/"histories"
OUT=ROOT/"artifacts/strong_lensing_response_completion_20261021"
os.environ["CRDS_PATH"]=str(CACHE/"crds")
os.environ["CRDS_CONTEXT"]="jwst_1535.pmap"
os.environ["CRDS_SERVER_URL"]="https://jwst-crds.stsci.edu"
os.environ["OMP_NUM_THREADS"]="1"
os.environ["OPENBLAS_NUM_THREADS"]="1"
os.environ["MPLCONFIGDIR"]=str(CACHE/"mpl")
import numpy as np
import yaml
from astropy.io import fits
from stdatamodels.jwst import datamodels
from jwst.pipeline import Detector1Pipeline
from stcal.ramp_fitting import ramp_fit as official_ramp_fit
from stdatamodels.jwst.datamodels import dqflags


def dump(path,obj):
    def convert(x):
        if isinstance(x,np.ndarray):return x.tolist()
        if isinstance(x,np.generic):return x.item()
        raise TypeError(type(x).__name__)
    path.write_text(json.dumps(obj,indent=2,default=convert)+"\n",encoding="utf-8")


def sha(path):
    h=hashlib.sha256()
    with path.open("rb") as f:
        for b in iter(lambda:f.read(8*1024**2),b""):h.update(b)
    return h.hexdigest()


class TaggedLoader(yaml.SafeLoader):pass


def tagged(loader,suffix,node):
    if isinstance(node,yaml.MappingNode):return loader.construct_mapping(node)
    if isinstance(node,yaml.SequenceNode):return loader.construct_sequence(node)
    return loader.construct_scalar(node)


TaggedLoader.add_multi_constructor("",tagged)


def configure(cal,stem):
    with fits.open(cal,memmap=False) as hd:
        text=hd["ASDF"].data["ASDF_METADATA"][0].tobytes().decode("latin1")
    tree=yaml.load(text[:text.index("\n...\n")+5],Loader=TaggedLoader)
    logs=tree["cal_logs"]["calwebb_detector1"]
    config=yaml.safe_load(next(x for x in logs if "Step Detector1Pipeline parameters are:" in x).split("parameters are:\n",1)[1])
    dump(OUT/f"{stem}_recorded_config.json",config)
    framework={"pre_hooks","post_hooks","output_file","output_dir","output_ext","output_use_model","output_use_index","save_results","suffix","search_output_file","input_dir"}
    overrides={"dq_init":{"mask":"jwst_nircam_mask_0075.fits"},
        "saturation":{"saturation":"jwst_nircam_saturation_0100.fits","superbias":"jwst_nircam_superbias_0180.fits"},
        "superbias":{"superbias":"jwst_nircam_superbias_0180.fits"},
        "refpix":{"refpix":"N/A","sirskernel":"jwst_nircam_sirskernel_0002.asdf"},
        "linearity":{"linearity":"jwst_nircam_linearity_0054.fits"},
        "dark_current":{"dark":"derived_jwst_nircam_dark_0337_first42.fits"},
        "jump":{"gain":"jwst_nircam_gain_0094.fits","readnoise":"jwst_nircam_readnoise_0282.fits"},
        "ramp_fit":{"gain":"jwst_nircam_gain_0094.fits","readnoise":"jwst_nircam_readnoise_0282.fits"},
        "gain_scale":{"gain":"jwst_nircam_gain_0094.fits"}}
    steps={}
    for name,values in config["steps"].items():
        steps[name]={key:(None if value=="None" else value) for key,value in values.items() if key not in framework}
        steps[name].update({"override_"+kind:(name if name=="N/A" else str(CACHE/name)) for kind,name in overrides.get(name,{}).items()})
    # Diagnostic storage only; numerical replay checks retain the original gates.
    steps["ramp_fit"].update(save_opt=True,opt_name=str(CACHE/f"{stem}_fitopt.fits"))
    return Detector1Pipeline(steps=steps),steps


def patch_info(cal_name):
    old=ROOT/"artifacts/strong_lensing_geometry_20261021"
    item=json.loads((old/"jwst_target_diagnostic/exposures"/cal_name.replace(".fits",".json")).read_text())
    fieldrows=json.loads((old/"standard_psf_convention_diagnostic.json").read_text())["records"]
    field=next(x for x in fieldrows if x["exposure"]==cal_name)
    labels=[p["label"] for p in item["points"]]+["FIELD"]
    guesses=np.asarray([p["centroid_xy_pixel"] for p in item["points"]]+[field["reproduced_global_centroid"]],float)
    centers=np.rint(guesses).astype(int)
    bounds=np.column_stack((centers[:,0]-32,centers[:,1]-32,centers[:,0]+33,centers[:,1]+33))
    return labels,guesses,bounds


def patches(array,bounds):
    return np.stack([np.array(array[...,y0:y1,x0:x1]) for x0,y0,x1,y1 in bounds])


def reference_patches(filename,ext,bounds):
    with fits.open(CACHE/filename,memmap=False) as h:
        return np.stack([np.array(h[ext].section[...,y0:y1,x0:x1]) for x0,y0,x1,y1 in bounds])


def collect_references(arrays,bounds):
    sources={"gain_e_per_DN":("jwst_nircam_gain_0094.fits","SCI"),
       "readnoise_CDS_DN":("jwst_nircam_readnoise_0282.fits","SCI"),
       "flat":("jwst_nircam_flat_0667.fits","SCI"),"flat_ERR":("jwst_nircam_flat_0667.fits","ERR"),"flat_DQ":("jwst_nircam_flat_0667.fits","DQ"),
       "superbias_DN":("jwst_nircam_superbias_0180.fits","SCI"),"superbias_ERR_DN":("jwst_nircam_superbias_0180.fits","ERR"),
       "linearity_coefficients":("jwst_nircam_linearity_0054.fits","COEFFS"),"linearity_DQ":("jwst_nircam_linearity_0054.fits","DQ"),
       "saturation_threshold_DN":("jwst_nircam_saturation_0100.fits","SCI"),"saturation_DQ":("jwst_nircam_saturation_0100.fits","DQ"),
       "dark_RAPID42_DN":("derived_jwst_nircam_dark_0337_first42.fits","SCI"),"dark_RAPID42_ERR_DN":("derived_jwst_nircam_dark_0337_first42.fits","ERR"),
       "dark_DQ":("derived_jwst_nircam_dark_0337_first42.fits","DQ")}
    for key,(name,ext) in sources.items():arrays[key]=reference_patches(name,ext,bounds)
    arrays["readnoise_CDS_electron"]=arrays["gain_e_per_DN"]*arrays["readnoise_CDS_DN"]


def collect_optional(arrays,path,bounds):
    record={"file":str(path.relative_to(ROOT)),"exists":path.exists(),"extensions":{}}
    if path.exists():
        with fits.open(path,memmap=False) as h:
            for ext in h[1:]:
                if ext.name=="ASDF" or ext.data is None:continue
                if getattr(ext.data,"ndim",0)>=2 and ext.data.shape[-2:]==(2048,2048):
                    key="native_optional_"+ext.name
                    arrays[key]=patches(ext.data,bounds);record["extensions"][ext.name]=list(arrays[key].shape)
    record["interpretation"]="Official optional segment diagnostics; segment weights are NOT asserted to be read/group estimator coefficients."
    return record


def native_compact_mapping(model,arrays):
    """Actual official OLS_C application, plus a declared local finite Jacobian.

    The native data-dependent SNR/segment weights are reevaluated, not replaced
    with an ordinary/uniform OLS formula. DQ is held at its observed value.
    """
    template=official_ramp_fit.create_ramp_fit_class(model,"OLS_C",dqflags.pixel,False)
    npatch,ng,ny,nx=arrays["corrected_group_DN"].shape
    slope=[];j1=[];j2=[]
    for p in range(npatch):
        def apply(group_values,zero_values):
            local=copy.copy(template)
            local.set_arrays(np.array(group_values[None],dtype=np.float32,copy=True),
                np.array(arrays["groupdq"][p,None],copy=True),np.array(arrays["pixeldq"][p],copy=True),
                np.array(arrays["average_dark_current_DN_s"][p],copy=True),
                zeroframe=np.array(zero_values[None],dtype=np.float32,copy=True))
            local.start_row=0;local.num_rows=ny
            local.zframe_mat=None;local.zframe_locs=None;local.zframe_cnt=0
            local.one_groups_locs=None;local.one_groups_time=None
            image,_,_=official_ramp_fit.ramp_fit_data(local,False,
                np.array(arrays["readnoise_CDS_DN"][p],copy=True),np.array(arrays["gain_e_per_DN"][p],copy=True),
                "OLS_C","optimal","none")
            return np.array(image["slope"],copy=True)
        values=arrays["corrected_group_DN"][p];zero=arrays["native_corrected_zeroframe_DN"][p]
        base=apply(values,zero);slope.append(base)
        steps=[]
        for h in [1.,2.]:
            jac=[]
            for j in range(ng+1):
                plus=values.copy();minus=values.copy();zp=zero.copy();zm=zero.copy()
                if j<ng:plus[j]+=h;minus[j]-=h
                else:zp+=h;zm-=h
                jac.append((apply(plus,zp)-apply(minus,zm))/(2*h))
            steps.append(np.asarray(jac))
        j1.append(steps[0]);j2.append(steps[1])
        print("NATIVE_GROUP_RATE_MAPPING",p,arrays["labels"][p],flush=True)
    slope=np.asarray(slope);j1=np.asarray(j1);j2=np.asarray(j2)
    arrays["native_compact_OLS_rate_DN_s"]=slope
    arrays["native_rate_jacobian_group_and_zeroframe_h1_per_s"]=j1
    arrays["native_rate_jacobian_group_and_zeroframe_h2_per_s"]=j2
    original=arrays["native_rate_pre_gain_scale_DN_s"]
    keep=np.isfinite(slope)&np.isfinite(original)&((arrays["native_rate_DQ"]&1)==0)
    difference=np.abs(slope-original)
    change=np.abs(j1-j2);use=np.maximum(np.abs(j1),np.abs(j2))>1e-5
    return {"application":"Official stcal1.18 ramp_fit_data OLS_C on actual compact histories; native metadata/selection/zeroframe/adaptive weights retained",
       "base_rate_max_abs_difference_from_full_detector_DN_s":float(np.max(difference[keep])),
       "base_rate_max_relative_difference_from_full_detector":float(np.max(difference[keep]/np.maximum(np.abs(original[keep]),1e-10))),
       "jacobian_columns":"groups0..8 then native corrected ZEROFRAME; output is pre gain_scale slope; units s^-1",
       "finite_difference_h_DN":[1.,2.],"jacobian_step_comparison_max_abs_per_s":float(np.nanmax(change)),
       "jacobian_step_comparison_max_relative_where_abs_gt_1e_5":float(np.nanmax(change[use]/np.maximum(np.abs(j1[use]),np.abs(j2[use])))),
       "limits":"Observed GROUPDQ fixed; adaptive native weights reevaluated. This finite Jacobian is a local conditional numerical linearization, not an exact global read-weight operator or calibrated covariance."}


def process(cal):
    stem=cal.stem.replace("_cal","");labels,guesses,bounds=patch_info(cal.name)
    output=COMPACT/f"{stem}_histories.npz";receipt=OUT/f"{stem}_history_receipt.json"
    if output.exists() and receipt.exists():
        print("REUSE_COMPLETED_HISTORY",output.name,flush=True);return json.loads(receipt.read_text())
    started=time.time();pipeline,settings=configure(cal,stem)
    raw=CACHE/f"{stem}_uncal.fits";arrays={"labels":np.asarray(labels),"bounds_xyxy":bounds,"prior_guess_native_xy":guesses}
    with fits.open(raw,memmap=False) as h:
        arrays["raw_group_DN"]=patches(h["SCI"].data[0],bounds)
        arrays["raw_zeroframe_DN"]=patches(h["ZEROFRAME"].data[0],bounds)
        header=h[0].header;tf=float(header["TFRAME"]);nf=int(header["NFRAMES"]);gap=int(header["GROUPGAP"]);ng=int(header["NGROUPS"])
        arrays["native_read_times_s"]=(np.arange(ng)[:,None]*(nf+gap)+np.arange(1,nf+1)[None,:])*tf
        arrays["native_dark_RAPID_frame_indices"]=np.arange(ng)[:,None]*(nf+gap)+np.arange(nf)[None,:]
        arrays["frame_time_s"]=tf;arrays["group_time_header_s"]=float(header["TGROUP"])
    arrays["x_pixel"],arrays["y_pixel"]=np.stack([np.meshgrid(np.arange(x0,x1),np.arange(y0,y1)) for x0,y0,x1,y1 in bounds],axis=1)
    collect_references(arrays,bounds)
    model=datamodels.RampModel(raw)
    operations=["group_scale","dq_init","saturation","ipc","superbias","refpix","linearity","persistence","dark_current","charge_migration","jump","picture_frame","clean_flicker_noise"]
    for name in operations:
        if name=="refpix":
            arrays["pre_refpix_group_DN"]=patches(model.data[0],bounds)
            ry,rx=np.where((model.pixeldq.astype(np.uint32)&np.uint32(2147483648))!=0)
            arrays["reference_control_native_xy"]=np.column_stack((rx,ry))
            arrays["pre_refpix_reference_control_DN"]=np.array(model.data[0][:,ry,rx])
            arrays["pre_refpix_zeroframe_reference_control_DN"]=np.array(model.zeroframe[0,ry,rx])
        if name=="linearity":
            arrays["pre_linearity_group_DN"]=patches(model.data[0],bounds)
            arrays["pre_linearity_zeroframe_DN"]=patches(model.zeroframe[0],bounds)
            x=arrays["pre_linearity_group_DN"].astype(float);coeff=arrays["linearity_coefficients"].astype(float)
            first=np.zeros_like(x);second=np.zeros_like(x)
            for power in range(1,coeff.shape[1]):
                first+=power*coeff[:,power,None]*x**(power-1)
                if power>=2:second+=power*(power-1)*coeff[:,power,None]*x**(power-2)
            arrays["linearity_polynomial_first_derivative"]=first
            arrays["linearity_polynomial_second_derivative"]=second
        if name=="dark_current":arrays["pre_dark_group_DN"]=patches(model.data[0],bounds)
        model=getattr(pipeline,name).run(model)
        if name=="refpix":arrays["post_refpix_group_DN"]=patches(model.data[0],bounds)
        if name=="linearity":arrays["post_linearity_group_DN"]=patches(model.data[0],bounds)
        if name=="dark_current":arrays["native_subtracted_dark_group_DN"]=arrays["pre_dark_group_DN"]-patches(model.data[0],bounds)
        print("NATIVE_COMPLETED",stem,name,round(time.time()-started,2),flush=True)
    arrays["corrected_group_DN"]=patches(model.data[0],bounds)
    arrays["groupdq"]=patches(model.groupdq[0],bounds)
    arrays["pixeldq"]=patches(model.pixeldq,bounds)
    arrays["native_corrected_zeroframe_DN"]=patches(model.zeroframe[0],bounds)
    # Current native RampModel has no ERR array. It is not an absent covariance
    # relation and must never be fabricated as zeros or copied from final CAL.
    if model.hasattr("err"):arrays["native_group_ERR_DN"]=patches(model.err[0],bounds)
    arrays["average_dark_current_DN_s"]=patches(model.average_dark_current,bounds) if model.hasattr("average_dark_current") else np.zeros_like(arrays["gain_e_per_DN"])
    pre=CACHE/f"{stem}_corrected_ramp.fits";model.save(pre)
    np.savez_compressed(CACHE/f"{stem}_pre_ramp_arrays.npz",**arrays)
    rate,rateints=pipeline.ramp_fit.run(model)
    arrays["native_rate_pre_gain_scale_DN_s"]=patches(rate.data,bounds)
    rate=pipeline.gain_scale.run(rate)
    ratepath=CACHE/f"{stem}_native_rate.fits";rate.save(ratepath)
    arrays["native_rate_DN_s"]=patches(rate.data,bounds)
    arrays["native_rate_ERR_DN_s"]=patches(rate.err,bounds)
    arrays["native_rate_DQ"]=patches(rate.dq,bounds)
    arrays["native_rate_VAR_POISSON"]=patches(rate.var_poisson,bounds)
    arrays["native_rate_VAR_RNOISE"]=patches(rate.var_rnoise,bounds)
    actual_mapping=native_compact_mapping(model,arrays)
    with fits.open(cal,memmap=False) as h:
        arrays["archived_cal_MJy_sr"]=patches(h["SCI"].data,bounds);arrays["archived_cal_ERR"]=patches(h["ERR"].data,bounds);arrays["archived_cal_DQ"]=patches(h["DQ"].data,bounds)
        arrays["photmjsr"]=float(h["SCI"].header["PHOTMJSR"])
    optional=collect_optional(arrays,CACHE/f"{stem}_fitopt.fits",bounds)
    # Native per-pixel group validity and jump/gap segment identifiers are data,
    # not uniform-OLS coefficients. Keep the actual GROUPDQ bit pattern too.
    valid=(arrays["groupdq"]&3)==0
    seg=np.full(valid.shape,-1,np.int16)
    for k in range(len(labels)):
        current=np.full(valid.shape[-2:],-1,np.int16);previous=np.zeros_like(current,bool)
        for j in range(ng):
            begin=valid[k,j]&((~previous)|((arrays["groupdq"][k,j]&4)!=0))
            current[begin]+=1;seg[k,j]=np.where(valid[k,j],current,-1);previous=valid[k,j]
    arrays["valid_group_DNU_SAT_only"]=valid;arrays["jump_gap_segment_index"]=seg
    equivalent=arrays["archived_cal_MJy_sr"]*arrays["flat"]/arrays["photmjsr"]
    good=((arrays["archived_cal_DQ"]&1)==0)&np.isfinite(equivalent)&np.isfinite(arrays["native_rate_DN_s"])
    rel=np.abs(arrays["native_rate_DN_s"]-equivalent)/np.maximum(np.abs(equivalent),1e-10)
    expected_var=arrays["photmjsr"]**2*((arrays["native_rate_VAR_POISSON"]+arrays["native_rate_VAR_RNOISE"])/arrays["flat"]**2+arrays["native_rate_DN_s"]**2*arrays["flat_ERR"]**2/arrays["flat"]**4)
    errrel=np.sqrt(expected_var)/arrays["archived_cal_ERR"]-1
    comparisons=[]
    for k,label in enumerate(labels):
        center=np.zeros((65,65),bool);center[24:41,24:41]=True;keep=good[k]&center
        comparisons.append({"label":label,"pixels":int(keep.sum()),"rate_max_relative_difference":float(np.max(rel[k,keep])),"ERR_max_relative_difference":float(np.max(abs(errrel[k,keep]))),
            "native_replay_gate_pass":bool(np.max(rel[k,keep])<1e-4 and np.max(abs(errrel[k,keep]))<1e-3),"rate_tolerance":1e-4,"ERR_tolerance":1e-3})
    np.savez_compressed(output,**arrays)
    record={"status":"ACTUAL_NATIVE_TARGET_AND_FIELD_HISTORIES_RECOVERED","file":str(output.relative_to(ROOT)),"sha256":sha(output),"bytes":output.stat().st_size,
        "cal_file":str(cal.relative_to(ROOT)),"raw_file":str(raw.relative_to(ROOT)),"labels":labels,"bounds_xyxy":bounds,"arrays":{k:list(np.shape(v)) for k,v in arrays.items()},
        "jwst":importlib.metadata.version("jwst"),"stcal":importlib.metadata.version("stcal"),"native_stcal_version_inferred_not_recorded":True,
        "full_detector_SIRS_before_crop":True,"settings":settings,"native_replay_checks":comparisons,"optional_slope_output":optional,
        "actual_native_slope_mapping":actual_mapping,
        "native_group_weights":"Actual optional segment quantities retained. Exact group estimator Jacobian not asserted by these segment weights; no arbitrary OLS replacement.",
        "read_level_availability":"Raw9onboard2-read averages plus trueZEROFRAME; individual later18readvalues were not downlinked. Exact supports retained; no invented read arrays.",
        "corrected_group_ERR":"Native RampModel supplies no ERR extension; no zero/independent noise array substituted. References and native operation Jacobians retained.",
        "reference_covariance_limits":["Exact gain/readnoise, biasERR, flatERR, dark42ERR and nonlinearity coefficients retained; reference ERR alone does not define their joint spatial/time covariance.",
            "Full SIRS kernel reference and actual input reference controls retained; actual data-dependent refpixel rejection/operator covariance must be derived separately.",
            "Native ZEROFRAME is not dark-subtracted by the dark step; do not call it fully dark-corrected.",
            "Polynomial Lprime/Lsecond are supplied at actual pre-linearity means; SAT/NO_LIN correction branches must be respected by consumers."],
        "new_filters_or_observations":False,"physical_admission":False,"elapsed_s":time.time()-started}
    dump(receipt,record);model.close();rate.close()
    if rateints is not None:rateints.close()
    print("COMPACT_HISTORY_READY",output.name,record["bytes"],comparisons,flush=True)
    return record


def main():
    parser=argparse.ArgumentParser();parser.add_argument("--exposure",type=int,choices=[1,2,3,4]);args=parser.parse_args()
    OUT.mkdir(parents=True,exist_ok=True);COMPACT.mkdir(parents=True,exist_ok=True)
    assert importlib.metadata.version("jwst")=="2.0.1"
    assert importlib.metadata.version("stcal").split("+")[0]=="1.18.0"
    native=ROOT/"data/strong_lensing_geometry_20261021/native"
    files=sorted(native.glob("jw01198007001_02101_0000[1-4]_nrcb3_cal.fits"))
    if args.exposure:files=[files[args.exposure-1]]
    records=[process(p) for p in files]
    dump(OUT/("history_recovery.json" if args.exposure is None else f"history_recovery_{args.exposure}.json"),{"records":records,"reserved_filters_accessed":False,"physical_admission":False})


if __name__=="__main__":main()
