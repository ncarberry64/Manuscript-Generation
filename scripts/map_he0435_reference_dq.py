"""Resolve reference DQ_DEF encodings with the official native data models.

Preserve the earlier diagnostic sidecars, whose pre-SIRS OR used raw reference
bits, and produce explicitly versioned corrected sidecars. No native rerun.
"""
from pathlib import Path
import hashlib,json
import numpy as np
from astropy.io import fits
from stdatamodels.jwst import datamodels
ROOT=Path(__file__).resolve().parents[1];BASE=ROOT/"data/strong_lensing_response_completion_20261021"
CACHE=BASE/"native_cache";OUT=ROOT/"artifacts/strong_lensing_response_completion_20261021"


def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()


def main():
    refs={"mask":"jwst_nircam_mask_0075.fits","saturation":"jwst_nircam_saturation_0100.fits",
          "superbias":"jwst_nircam_superbias_0180.fits","linearity":"jwst_nircam_linearity_0054.fits","flat":"jwst_nircam_flat_0667.fits"}
    mapped={};defs={}
    classes={'mask':datamodels.MaskModel,'saturation':datamodels.SaturationModel,
             'superbias':datamodels.SuperBiasModel,'linearity':datamodels.LinearityModel,'flat':datamodels.FlatModel}
    for key,name in refs.items():
        with classes[key](CACHE/name) as m:mapped[key]=np.array(m.dq)
        with fits.open(CACHE/name,memmap=False) as h:
            defs[key]=[{"bit":int(row['BIT']),"value":int(row['VALUE']),"name":str(row['NAME']),"description":str(row['DESCRIPTION'])} for row in h['DQ_DEF'].data]
    rows=[]
    olddir=BASE/"histories/superseded_reference_dq_v1";olddir.mkdir(exist_ok=True)
    oldreports=OUT/"superseded_reference_dq_v1";oldreports.mkdir(exist_ok=True)
    for exp in range(1,5):
        stem=f"jw01198007001_02101_0000{exp}_nrcb3";base=np.load(BASE/"histories"/f"{stem}_histories.npz")
        original=BASE/"histories"/f"{stem}_native_diagnostics.npz"
        if original.exists():original.rename(olddir/original.name)
        oldreceipt=OUT/f"{stem}_native_diagnostics.json"
        if oldreceipt.exists():oldreceipt.rename(oldreports/(oldreceipt.name+".txt"))
        old=np.load(olddir/f"{stem}_native_diagnostics.npz")
        arrays={k:np.array(old[k]) for k in old.files};bounds=base['bounds_xyxy'];x,y=base['reference_control_native_xy'].T
        for key,dq in mapped.items():
            arrays[key+"_DQ_mapped"]=np.stack([dq[y0:y1,x0:x1] for x0,y0,x1,y1 in bounds])
            arrays['reference_control_'+key+'_DQ_mapped']=dq[y,x]
        pre=arrays['reference_control_raw_pixeldq'].astype(np.uint32).copy()
        for key in ['mask','saturation','superbias']:pre|=arrays['reference_control_'+key+'_DQ_mapped'].astype(np.uint32)
        arrays['reference_control_pre_sirs_pixeldq']=pre
        delta=base['native_rate_jacobian_group_and_zeroframe_h1_per_s'][:,:,24:41,24:41]-base['native_rate_jacobian_group_and_zeroframe_h2_per_s'][:,:,24:41,24:41]
        keep=np.broadcast_to(((base['archived_cal_DQ'][:,24:41,24:41]&1)==0)[:,None],delta.shape)&np.isfinite(delta)
        output=BASE/"histories"/f"{stem}_native_diagnostics_v2.npz"
        if output.exists():raise RuntimeError('refuse existing v2 overwrite')
        np.savez_compressed(output,**arrays)
        row={'status':'OFFICIAL_REFERENCE_DQ_MAPPING_AMENDMENT','file':str(output.relative_to(ROOT)),'sha256':sha(output),
             'source':'stdatamodels5.0.2 native dynamic_mask through datamodels.open; exact same native reference files',
             'superseded_diagnostic_file':str((olddir/f'{stem}_native_diagnostics.npz').relative_to(ROOT)),
             'main_history_unchanged':True,'main_reference_DQ_encoding':'Main reference *_DQ fields retain original FITS DQ_DEF codes; use new *_DQ_mapped for global JWST bits. Science GROUPDQ/PIXELDQ were already mapped by the native pipeline.',
             'pre_sirs_control_DNU_count':int(np.sum((pre&1)!=0)),
             'pre_sirs_control_DNU_disagrees_final_count':int(np.sum((pre&1)!=(arrays['reference_control_final_pixeldq']&1))),
             'h1_h2_likelihood17_usable_max_abs_J_difference_per_s':float(np.max(abs(delta[keep]))),'excluded_nonfinite_or_DNU_J_entries':int(np.sum(~keep)),
             'arrays':{k:list(v.shape) for k,v in arrays.items()},'native_processing_rerun':False,'physical_admission':False}
        (OUT/f'{stem}_native_diagnostics_v2.json').write_text(json.dumps(row,indent=2,allow_nan=False)+'\n');rows.append(row)
        print(stem,row['pre_sirs_control_DNU_count'],row['pre_sirs_control_DNU_disagrees_final_count'],flush=True)
    (OUT/'reference_dq_mapping_amendment.json').write_text(json.dumps({'records':rows,'DQ_DEF_tables':defs,'error_corrected':'Earlier diagnostic pre-SIRS DQ OR incorrectly used raw reference bit codes. Native science processing and replay were correct; all original arrays and receipts are preserved.','native_processing_rerun':False},indent=2,allow_nan=False)+'\n')


if __name__=='__main__':main()
