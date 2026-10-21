"""Acquire exact existing F115W raw exposures and recorded detector references.

The huge RAPID dark is transferred only for its exact first42 SCI/ERR frames
and full DQ. FITS pixels are copied byte-for-byte; gaps are retained, not filled.
"""
from pathlib import Path
from concurrent.futures import ThreadPoolExecutor,as_completed
import hashlib,json,time,urllib.request,threading
import numpy as np
from astropy.io import fits

ROOT=Path(__file__).resolve().parents[1]
BASE=ROOT/"data/strong_lensing_response_completion_20261021"
CACHE=BASE/"native_cache"
OUT=ROOT/"artifacts/strong_lensing_response_completion_20261021"
CRDS="https://jwst-crds.stsci.edu/unchecked_get/references/jwst/"


def digest(path):
    h=hashlib.sha256()
    with path.open("rb") as f:
        for b in iter(lambda:f.read(8*1024**2),b""):h.update(b)
    return h.hexdigest()


def request(url,first=None,last=None):
    headers={"User-Agent":"HE0435-same-photon-history/1"}
    if first is not None:headers["Range"]=f"bytes={first}-{last}"
    return urllib.request.urlopen(urllib.request.Request(url,headers=headers),timeout=180)


def get_range(url,first,last):
    for attempt in range(5):
        try:
            with request(url,first,last) as r:
                if r.status!=206:raise RuntimeError(f"Range not honored {r.status}")
                actual=r.headers.get("Content-Range","")
                if not actual.startswith(f"bytes {first}-{last}/"):raise RuntimeError(actual)
                b=r.read()
                if len(b)!=last-first+1:raise RuntimeError("short range")
                return b,dict(r.headers)
        except Exception:
            if attempt==4:raise
            time.sleep(1+attempt)


def acquire(name,url,expected_hash=None):
    p=CACHE/name
    if not p.exists():
        tmp=p.with_suffix(p.suffix+".download")
        for attempt in range(5):
            try:
                with request(url) as r,tmp.open("wb") as f:
                    total=int(r.headers.get("Content-Length",0));done=0
                    while b:=r.read(4*1024**2):f.write(b);done+=len(b)
                if total and done!=total:raise RuntimeError("incomplete download")
                tmp.replace(p);break
            except Exception:
                if attempt==4:raise
                time.sleep(1+attempt)
    sha=digest(p)
    if expected_hash and sha!=expected_hash:raise RuntimeError(f"Reference identity mismatch {name}")
    row={"file":str(p.relative_to(ROOT)),"url":url,"bytes":p.stat().st_size,"sha256":sha}
    print("ACQUIRED",name,row["bytes"],flush=True)
    return row


def remote_layout(url):
    first,headers=get_range(url,0,2879)
    size=int(headers["Content-Range"].split("/")[-1]);offset=0;layout=[]
    while offset<size:
        chunks=[];ended=False
        for block in range(100):
            b,_=get_range(url,offset+block*2880,offset+(block+1)*2880-1)
            chunks.append(b)
            if any(b[i:i+8]==b"END     " for i in range(0,len(b),80)):
                ended=True;break
        if not ended:raise RuntimeError("FITS header END missing")
        header=fits.Header.fromstring(b"".join(chunks).decode("ascii"),sep="")
        headbytes=len(chunks)*2880;dataoffset=offset+headbytes
        count=int(np.prod([header.get(f"NAXIS{i}",0) for i in range(1,header.get("NAXIS",0)+1)])) if header.get("NAXIS",0) else 0
        datasize=count*abs(header.get("BITPIX",8))//8+header.get("PCOUNT",0)
        layout.append((header,dataoffset,datasize))
        offset=dataoffset+((datasize+2879)//2880)*2880
    return layout,{"source_size":size,"etag":headers.get("ETag",headers.get("Etag",headers.get("etag"))),"url":url}


def dark_subset():
    name="derived_jwst_nircam_dark_0337_first42.fits";dest=CACHE/name
    receipt=OUT/"dark_subset_acquisition.json"
    if dest.exists() and receipt.exists():return json.loads(receipt.read_text())
    url=CRDS+"jwst_nircam_dark_0337.fits"
    layout,meta=remote_layout(url);pieces=[];offset=0
    for i,(h,src,nbytes) in enumerate(layout):
        ext=h.get("EXTNAME","PRIMARY")
        if ext=="ASDF":continue
        # DQ_DEF and any additional metadata extension are copied unchanged.
        h=h.copy()
        for key in ("CHECKSUM","DATASUM"):
            if key in h:del h[key]
        if i==0:
            h["NGROUPS"]=42;h["FILENAME"]=name
            h.add_history("Exact first42 RAPID dark frames for native SHALLOW2 groups; no new detector model")
        if ext in ("SCI","ERR"):
            assert h["NAXIS3"]==187
            nbytes=nbytes//187*42;h["NAXIS3"]=42
        encoded=h.tostring(endcard=True,padding=True).encode("ascii")
        pieces.append({"extension":ext,"header":encoded,"output_header_offset":offset,"output_data_offset":offset+len(encoded),"source_data_offset":src,"bytes":nbytes})
        offset+=len(encoded)+((nbytes+2879)//2880)*2880
    tmp=dest.with_suffix(".fits.download")
    with tmp.open("wb") as f:
        f.truncate(offset)
        for item in pieces:f.seek(item["output_header_offset"]);f.write(item["header"])
    jobs=[]
    for item in pieces:
        for delta in range(0,item["bytes"],32*1024**2):
            n=min(32*1024**2,item["bytes"]-delta)
            jobs.append((item["source_data_offset"]+delta,item["output_data_offset"]+delta,n,item["extension"]))
    lock=threading.Lock();rows=[];done=0
    def fetch(job):
        src,dst,n,ext=job;b,headers=get_range(url,src,src+n-1)
        with lock:
            with tmp.open("r+b") as f:f.seek(dst);f.write(b)
        return {"extension":ext,"source_first_byte":src,"bytes":n,"sha256":hashlib.sha256(b).hexdigest(),"etag":headers.get("ETag",headers.get("Etag",headers.get("etag")))}
    with ThreadPoolExecutor(max_workers=4) as pool:
        for future in as_completed([pool.submit(fetch,job) for job in jobs]):
            row=future.result();rows.append(row);done+=row["bytes"]
            print("DARK_EXACT_BYTES",done,"/",sum(x[2] for x in jobs),flush=True)
    tmp.replace(dest)
    record={**meta,"status":"EXACT_USED_TEMPORAL_DARK_SUBSET","file":str(dest.relative_to(ROOT)),"sha256":digest(dest),"bytes":dest.stat().st_size,
            "retained_RAPID_frame_indices":list(range(42)),"native_group_frame_indices":[[5*j,5*j+1] for j in range(9)],"ranges":sorted(rows,key=lambda x:x["source_first_byte"]),
            "ASDF":"Omitted original187frame dimensions; FITS model metadata rebuilds exact42frame reference",
            "full_6GB_dark_downloaded":False}
    receipt.write_text(json.dumps(record,indent=2)+"\n")
    print("DARK_SUBSET_READY",record["bytes"],flush=True)
    return record


def main():
    CACHE.mkdir(parents=True,exist_ok=True);OUT.mkdir(parents=True,exist_ok=True)
    native=ROOT/"data/strong_lensing_geometry_20261021/native"
    files=sorted(native.glob("jw01198007001_02101_0000[1-4]_nrcb3_cal.fits"))
    if len(files)!=4:raise RuntimeError("exact four consumed CAL exposures required")
    keys=["R_SIRSKL","R_DARK","R_FLAT","R_GAIN","R_LINEAR","R_MASK","R_READNO","R_SATURA","R_SUPERB"]
    records=[];refs=set()
    for p in files:
        h=fits.getheader(p,0);record={"cal_file":str(p.relative_to(ROOT)),"detector_header":{k:h.get(k) for k in ["CAL_VER","READPATT","NINTS","NGROUPS","NFRAMES","GROUPGAP","TFRAME","TGROUP","CRDS_CTX",*keys]}}
        records.append(record)
        refs.update(h[k].replace("crds://","") for k in keys)
    (OUT/"native_source_headers.json").write_text(json.dumps(records,indent=2)+"\n")
    manifest=json.loads((ROOT/"data/strong_lensing_response_repair_20261021/delivery/ramp_diagnostic/native_pipeline/native_input_manifest.json").read_text())
    hashes={Path(x["path_relative_to_REPAIR"]).name:x["sha256"] for x in manifest["inputs"]}
    jobs=[(p.name.replace("_cal","_uncal"),"https://mast.stsci.edu/api/v0.1/Download/file?uri=mast:JWST/product/"+p.name.replace("_cal","_uncal")) for p in files]
    jobs.extend((name,CRDS+name) for name in sorted(refs) if name!="jwst_nircam_dark_0337.fits")
    with ThreadPoolExecutor(max_workers=4) as pool:
        pending=[pool.submit(acquire,name,url,hashes.get(name)) for name,url in jobs]
        rows=[future.result() for future in as_completed(pending)]
    dark=dark_subset()
    (OUT/"acquisition.json").write_text(json.dumps({"same_four_consumed_F115W_only":True,"new_filters_or_observations":False,"products":rows,"dark_subset":dark["file"]},indent=2)+"\n")


if __name__=="__main__":main()
