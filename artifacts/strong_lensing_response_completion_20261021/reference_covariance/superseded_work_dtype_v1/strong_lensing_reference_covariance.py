"""Conditional NIRCam SIRS reference-control operator and named noise component.

This is a reconstruction of the native refpix relation, not a detector replay.
Observed rejection sets are frozen. Raw-ADC covariance is NOT identified by a
post-calibration READNOISE reference; covariance amplitudes below are hypotheses.
"""
from __future__ import annotations
import numpy as np
from scipy import sparse,linalg,stats


def _clipped(indices,values,dq):
    idx=indices[(dq[indices]&1)==0]
    if not len(idx):return idx,None
    sample=values[idx]
    if np.std(sample,dtype=np.float64)!=0:
        clipped,low,high=stats.sigmaclip(sample,3.,3.)
        idx=idx[(sample>=low)&(sample<=high)]
        return idx,clipped.mean()
    return idx,sample.mean(dtype=np.float64)


def conditional_sirs(control_values,control_xy_native,control_dq,output_xy_native,
                     kernels,variance,size=2048,side=4,sigreject=4.):
    """Return c=Jref r=(D+UT)r and the conditional white-component Gram.

    NRCB3 FASTAXIS=1,SLOWAXIS=-2 is explicit. T comprises eight independent
    clipped top/bottom odd/even offset rows. D is direct side convolution; U
    subtracts the same top/bottom offsets from the side controls before SIRS.
    Retaining DΣT' cross terms is essential. Side MAD rejection ignores DQ,
    exactly as the native SIRS implementation; top/bottom clipping uses DNU.
    """
    r=np.asarray(control_values,dtype=np.float32);xy=np.asarray(control_xy_native,int)
    out=np.asarray(output_xy_native,int);dq=np.asarray(control_dq,np.uint32)
    variance=np.asarray(variance,float)
    if np.any(~np.isfinite(variance)) or np.any(variance<0):raise ValueError('finite nonnegative conditional variances required')
    detx=xy[:,0];dety=size-1-xy[:,1]
    index=np.full((size,size),-1,np.int32);index[dety,detx]=np.arange(len(r))
    ampwidth=size//4
    offsets=np.zeros((4,2),dtype=np.float32);tr=[];tc=[];tv=[];clip_counts=np.zeros((4,2,2),int)
    for a in range(4):
        for parity in range(2):
            selected=[];means=[]
            for k,ys in enumerate([slice(size-side,size),slice(0,side)]):
                ids=index[ys,a*ampwidth+parity:(a+1)*ampwidth:2].ravel()
                if np.any(ids<0):raise ValueError('missing native top/bottom reference control')
                ids,mean=_clipped(ids,r,dq);selected.append(ids);means.append(mean);clip_counts[a,parity,k]=len(ids)
            available=[k for k,m in enumerate(means) if m is not None]
            if not available:raise ValueError('no native top/bottom controls for an output')
            level=means[available[0]] if len(available)==1 else .5*(means[0]+means[1])
            offsets[a,parity]=level
            for k in available:
                ids=selected[k];tr.extend([2*a+parity]*len(ids));tc.extend(ids);tv.extend([1/(len(available)*len(ids))]*len(ids))
    t=sparse.csr_matrix((tv,(tr,tc)),shape=(8,len(r)))
    # Exactly the native float32 top/bottom subtraction, before float64 SIRS.
    after=np.asarray(r-offsets[detx//ampwidth,detx%2],dtype=np.float32)
    left_idx=index[:,:side];right_idx=index[:,-side:]
    if np.any(left_idx<0) or np.any(right_idx<0):raise ValueError('missing native side reference control')
    lv=after[left_idx].astype(float);rv=after[right_idx].astype(float)
    masks=[]
    for values in [lv,rv]:
        median=np.median(values);sigma=1.48*np.median(abs(values-median))
        masks.append(abs(values-median)<=sigreject*sigma)
    gl,gr=masks;nl=gl.sum(axis=1);nr=gr.sum(axis=1)
    l=(lv*gl).sum(axis=1)/side;rr=(rv*gr).sum(axis=1)/side
    keys=np.column_stack((size-1-out[:,1],out[:,0]//ampwidth,out[:,0]%2))
    features,inverse=np.unique(keys,axis=0,return_inverse=True)
    dr=[];dc=[];dv=[];u=np.zeros((len(features),8));native_correction=np.zeros(len(features))
    norm_left=[];norm_right=[];templates=[]
    for a in range(4):
        kl=np.asarray(kernels[0][a],float);kr=np.asarray(kernels[1][a],float)
        norml=np.convolve(np.ones(size),kl,mode='same')/np.convolve(nl/side,kl,mode='same')
        normr=np.convolve(np.ones(size),kr,mode='same')/np.convolve(nr/side,kr,mode='same')
        if np.any(~np.isfinite(norml)) or np.any(~np.isfinite(normr)):raise ValueError('undefined native SIRS normalization')
        norm_left.append(norml);norm_right.append(normr)
        templates.append(np.convolve(l,kl,mode='same')*norml+np.convolve(rr,kr,mode='same')*normr)
    for p,(row,a,parity) in enumerate(features):
        u[p,2*a+parity]=1
        native_correction[p]=float(offsets[a,parity])+templates[a][row]
        for which,(idx,good,kernel,norm,side_amp) in enumerate([
            (left_idx,gl,kernels[0][a],norm_left[a],0),(right_idx,gr,kernels[1][a],norm_right[a],3)]):
            kernel=np.asarray(kernel);half=len(kernel)//2
            rows=np.arange(max(0,row-half),min(size,row+half+1))
            coeff=(norm[row]*kernel[row-rows+half])[:,None]*good[rows]/side
            ids=idx[rows].ravel();coeff=coeff.ravel();keep=coeff!=0;ids=ids[keep];coeff=coeff[keep]
            dr.extend([p]*len(ids));dc.extend(ids);dv.extend(coeff)
            for par in range(2):u[p,2*side_amp+par]-=coeff[detx[ids]%2==par].sum()
    d=sparse.csr_matrix((dv,(dr,dc)),shape=(len(features),len(r)))
    ds=d.multiply(variance[None]);ts=t.multiply(variance[None])
    cross=(ds@t.T).toarray();small=(ts@t.T).toarray()
    covariance=(ds@d.T).toarray()+u@small@u.T+cross@u.T+u@cross.T
    covariance=(covariance+covariance.T)/2
    eig,vec=linalg.eigh(covariance);tol=np.finfo(float).eps*len(features)*max(float(eig[-1]),1.)
    if eig[0]<-tol:raise ValueError('non-PSD derived reference covariance')
    keep=eig>tol;factor=vec[:,keep]*np.sqrt(eig[keep])
    linear=d@r.astype(float)+u@(t@r.astype(float))
    return dict(T=t,D=d,U=u,features=features,pixel_feature_index=inverse,
                native_correction=native_correction,native_offsets=offsets,
                side_template=np.asarray(templates),side_good_left=gl,side_good_right=gr,
                top_bottom_clip_counts=clip_counts,covariance=covariance,factor=factor,
                eigenvalues=eig,eigenvalue_roundoff_tolerance=tol,
                linear_operator_reconstruction=linear)


def predicted_native_post(pre,xy_native,operator,size=2048):
    xy=np.asarray(xy_native,int);amp=xy[:,0]//(size//4);par=xy[:,0]%2;row=size-1-xy[:,1]
    first=np.asarray(np.asarray(pre,dtype=np.float32)-operator['native_offsets'][amp,par],dtype=np.float32)
    return np.asarray(first-operator['side_template'][amp,row],dtype=np.float32)


def export_sparse(arrays,prefix,matrix):
    m=matrix.tocsr();arrays[prefix+'_data']=m.data;arrays[prefix+'_indices']=m.indices;arrays[prefix+'_indptr']=m.indptr;arrays[prefix+'_shape']=np.array(m.shape)
