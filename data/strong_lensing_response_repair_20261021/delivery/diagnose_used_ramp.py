"""Limited group-level diagnostic, not a replacement JWST calibration pipeline.

Uses the same consumed F115W exposure. Applies its recorded superbias and
linearity reference to a star stamp, then checks group-difference profiles.
No lens parameters enter this calculation; native CAL values are unchanged.
"""
from pathlib import Path
import json
import numpy as np
from astropy.io import fits

ROOT=Path(__file__).resolve().parent
OUT=ROOT/'ramp_diagnostic'
RAW=OUT/'jw01198007001_02101_00001_nrcb3_uncal.fits'
CAL=ROOT/'repo/data/strong_lensing_geometry_20261021/native'/RAW.name.replace('_uncal','_cal')
cx,cy=899.7549609755204,875.9549994961739
x0,x1,y0,y1=867,933,843,909
sl=np.s_[y0:y1,x0:x1]
with fits.open(RAW) as h:
    raw=np.array(h['SCI'].data[0,:,y0:y1,x0:x1],dtype=float)
    zero=np.array(h['ZEROFRAME'].data[0,y0:y1,x0:x1],dtype=float)
    head=dict(h[0].header)
    reference=np.concatenate([np.asarray(h['SCI'].data[0,:,:4,512:1024],float),
                              np.asarray(h['SCI'].data[0,:,-4:,512:1024],float)],axis=1)
    zero_reference=np.concatenate([np.asarray(h['ZEROFRAME'].data[0,:4,512:1024],float),
                                   np.asarray(h['ZEROFRAME'].data[0,-4:,512:1024],float)])
with fits.open(OUT/'jwst_nircam_superbias_0180.fits') as h:
    bias=np.array(h['SCI'].data[sl],dtype=float)
    bq=np.array(h['DQ'].data[sl],dtype=np.uint64)
    reference_bias=np.concatenate([np.asarray(h['SCI'].data[:4,512:1024],float),
                                   np.asarray(h['SCI'].data[-4:,512:1024],float)])
with fits.open(OUT/'jwst_nircam_linearity_0054.fits') as h:
    coeff=np.array(h['COEFFS'].data[:,y0:y1,x0:x1],dtype=float)
    lq=np.array(h['DQ'].data[sl],dtype=np.uint64)
with fits.open(OUT/'jwst_nircam_saturation_0100.fits') as h:
    sat=np.array(h['SCI'].data[sl],dtype=float)
with fits.open(CAL) as h:
    cq=np.array(h['DQ'].data[sl],dtype=np.uint64)
    sci=np.array(h['SCI'].data[sl],dtype=float)
    err=np.array(h['ERR'].data[sl],dtype=float)

reference_offsets=np.stack([np.median((reference-reference_bias)[:,:,parity::2],axis=(1,2))
                            for parity in range(2)],axis=1)
zero_offsets=np.array([np.median((zero_reference-reference_bias)[:,parity::2]) for parity in range(2)])
col_parity=np.arange(x0,x1)%2

def linearize(data,offsets):
    v=data-bias-offsets
    out=np.zeros_like(v)
    for c in coeff[::-1]:out=out*v+c
    return out

corrected=linearize(raw,reference_offsets[:,col_parity][:,None,:])
zc=linearize(zero,zero_offsets[col_parity][None,:])
yy,xx=np.mgrid[y0:y1,x0:x1]
r=np.hypot(xx-cx,yy-cy)
# A common mask avoids comparing different spatial samples between groups.
valid=np.all(np.isfinite(corrected),axis=0)&((cq&3)==0)&((bq&1)==0)&((lq&1)==0)
valid&=np.all(raw<sat,axis=0)
bg=valid&(r>20)&(r<30)
background=np.median(corrected[:,bg],axis=1)
corrected=corrected-background[:,None,None]
rates=np.diff(corrected,axis=0)/float(head['TGROUP'])
rows=[]
for lo,hi in [(0,1),(1,2),(2,3),(3,4),(4,6),(6,8),(8,12)]:
    m=valid&(r>=lo)&(r<hi)
    values=rates[:,m].sum(axis=1)
    rows.append({'radius_pixel':[lo,hi],'n_common_pixels':int(m.sum()),
                 'group_difference_sum_dn_s':values.tolist(),
                 'late_minus_early_dn_s':float(values[4:].mean()-values[:4].mean()),
                 'relative_late_minus_early':float((values[4:].mean()-values[:4].mean())/values.mean()) if values.mean()!=0 else None})
early=rates[:4].mean(axis=0)
late=rates[4:].mean(axis=0)
np.savez_compressed(OUT/'field_ramp_stamp.npz',raw=raw,zero=zero,superbias=bias,
    coefficients=coeff,threshold=sat,cal_dq=cq,cal_sci=sci,cal_err=err,
    reference_offsets_dn=reference_offsets,
    linearized_background_subtracted=corrected,group_differences_dn_s=rates,
    common_valid_mask=valid,early=early,late=late,x=xx,y=yy)
report={'status':'LIMITED_RAW_RAMP_DEVELOPMENT_DIAGNOSTIC',
    'exposure':RAW.name,'group_time_s':head['TGROUP'],'n_groups':head['NGROUPS'],
    'operations':['recorded superbias subtraction',
                  'diagnostic median upper/lower reference-pixel subtraction, same amplifier and odd/even column',
                  'recorded polynomial linearity correction',
                  'common CAL non-DNU/non-SAT mask','per-group local annular median subtraction'],
    'limitations':['Not the complete detector1 pipeline: approximate median reference-pixel offsets instead of native SIRS; no dark/cosmic-ray reprocessing.',
                   'No propagated group covariance or formal significance test.',
                   'This check cannot isolate charge migration from wavefront or pointing changes.'],
    'reference_offsets_dn':reference_offsets.tolist(),
    'annulus_background_dn':background.tolist(),'radial_bins':rows,
    'native_CAL_changed':False,'target_parameters_fitted':False,'physical_validation':False}
(OUT/'ramp_diagnostic.json').write_text(json.dumps(report,indent=2)+'\n')
print(json.dumps(report,indent=2))
