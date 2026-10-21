"""Actual common two-group rate products and explicitly conditional noise.

These are new temporal projections of consumed histories, not synthetic data
or an error rescaling of the original full-duration native CAL products.
"""
from __future__ import annotations
import numpy as np


def first_two_group_product(groups,groupdq,times,group_time,flat,photmjsr,
                            gain,cds_noise_dn,positive_pilot_rate_dn_s):
    """Return the exact native length-two slope and a named pilot covariance.

    The covariance assumes independent post-correction per-read CDS marginals
    plus conditional independent Poisson carriers, with no added SIRS factor.
    That assumption does not determine the missing spatial/time read covariance.
    Its illumination is a fixed positive native pilot, not an updated likelihood
    variance or a fitted error multiplier.
    """
    g=np.asarray(groups,float);dq=np.asarray(groupdq);t=np.asarray(times,float)[:2]
    slope=(g[1]-g[0])/group_time
    valid=((dq[0]&3)==0)&((dq[1]&7)==0)&np.isfinite(slope)&(flat>0)&(gain>0)
    rate=np.asarray(positive_pilot_rate_dn_s,float)
    valid &= np.isfinite(rate)&(rate>=0)&np.isfinite(cds_noise_dn)&(cds_noise_dn>=0)
    w=np.array([-.5,-.5,.5,.5])/group_time
    tt=t.ravel();time_factor=float(w@np.minimum(tt[:,None],tt[None,:])@w)
    # The archived fitter uses rounded TGROUP, while physical frame times
    # retain their recorded support. Convert its pilot slope to physical rate.
    physical_rate=rate*group_time/(float(t[1].mean()-t[0].mean()))
    photon_variance=physical_rate/gain*time_factor
    read_variance=np.asarray(cds_noise_dn,float)**2/(2*group_time**2)
    scale=photmjsr/flat
    return {'rate_SCI':slope*scale,'variance_SCI2':(photon_variance+read_variance)*scale**2,
            'photon_variance_SCI2':photon_variance*scale**2,'CDS_marginal_variance_SCI2':read_variance*scale**2,
            'valid':valid,'poisson_time_factor_s_inverse':time_factor,'SCI_per_DN_s':scale}


def adapt_record(record,history,patch):
    """Restrict existing native rows to the same actual early-read domain."""
    bounds=history['bounds_xyxy'][patch];ix=record['x'].astype(int)-bounds[0];iy=record['y'].astype(int)-bounds[1]
    product=first_two_group_product(history['corrected_group_DN'][patch],history['groupdq'][patch],
          history['native_read_times_s'],float(history['group_time_header_s']),history['flat'][patch],float(history['photmjsr']),
          history['gain_e_per_DN'][patch],history['readnoise_CDS_DN'][patch],history['native_rate_DN_s'][patch])
    keep=product['valid'][iy,ix]&(product['variance_SCI2'][iy,ix]>0)
    for key in ['x','y','data','error','detector_background','scene_matrix','fine_ix','fine_iy']:
        if key in record:record[key]=record[key][keep]
    ix,iy=ix[keep],iy[keep]
    record['original_full_native_SCI']=record['data'].copy();record['original_full_native_ERR']=record['error'].copy()
    record['data']=product['rate_SCI'][iy,ix];record['error']=np.sqrt(product['variance_SCI2'][iy,ix])
    record['early_product_receipt']={'original_rows':len(keep),'selected_rows':int(keep.sum()),
         'exact_two_group_weights_per_s':[-1/float(history['group_time_header_s']),1/float(history['group_time_header_s'])],
         'poisson_time_factor_s_inverse':product['poisson_time_factor_s_inverse'],
         'conditional_CDS_stage':'Independent post-correction per-read CDS marginals; no SIRS addition; full cross covariance unresolved.',
         'negative_pilot_not_clipped':True,'source_complete_noise':False}
    return product
