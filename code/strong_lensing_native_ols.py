"""Double continuous extension/Jacobian of the recorded stcal1.18 OLS_C.

This implements its source equations for ONE integration and fixed GROUPDQ,
SNR-power and median-order branches. It does not replace the native float32
value application, differentiate quantization, or invent fixed universal OLS
weights. Median ties and SNR boundaries are individually reported.
"""
import numpy as np


def _segments(dq):
    """Native tests dq==0, admitting a preceding PURE JUMP as next start."""
    segments = []
    start = None
    for j, flag in enumerate(dq):
        if flag == 0:
            if start is None:
                start = j-1 if j > 0 and dq[j-1] == 4 else j
        elif start is not None:
            segments.append((start, j)); start = None
    if start is not None:
        segments.append((start, len(dq)))
    if segments:
        longest = max(b-a for a,b in segments)
        segments = [s for s in segments if s[1]-s[0] > 1] if longest > 1 else segments[:1]
    return segments


def piecewise_ols_pixel(groups_dn, groupdq, *, gain, readnoise_cds_dn,
                        group_time, frame_time, nframes=2, dark_dn_s=0.,
                        zeroframe=None, pixeldq=0, suppress_one_group=False,
                        cast_input=False):
    """Return slope, d(slope)/d(groups,zeroframe), and exact branch diagnostics.

    Input units are native ramp DN, CDS-reference DN, e/DN, and seconds.
    Native C converts times to float32 and readnoise to float32 electron
    sigma=CDS_DN*gain/sqrt(2*NFRAMES); those source constants are reproduced.
    cast_input=True mirrors native float32 data before its double arithmetic;
    False evaluates a continuous extension at a model mean. Native published
    values still come from the original C estimator.
    """
    y = np.asarray(groups_dn, dtype=np.float32 if cast_input else float).astype(float).copy()
    dq = np.asarray(groupdq, np.uint8).copy()
    if y.ndim != 1 or dq.shape != y.shape or len(y) < 1:
        raise ValueError("One integration's group vector and matching DQ required")
    n = len(y); jac = np.zeros(n+1)
    tg = float(np.float32(group_time)); tf = float(np.float32(frame_time))
    gain = float(np.float32(gain)); dark = float(np.float32(dark_dn_s))
    rn = float(np.float32(float(np.float32(readnoise_cds_dn))*gain/np.sqrt(2.*nframes)))
    if min(tg,tf) <= 0 or nframes < 1:
        raise ValueError("Positive timing and NFRAMES required")
    branch = {"cast_input": bool(cast_input), "native_group_time_s": tg,
        "native_frame_time_s": tf, "native_rnoise_e": rn,
        "native_value_application_required": True, "single_integration": True,
        "median_tie": False, "snr_boundary": False, "median_zero_boundary": False,
        "zeroframe_branch": False, "one_good_first_group_branch": False}
    if (pixeldq & 3) or not np.isfinite(gain) or gain <= 0:
        return {"slope_double": np.nan, "jacobian_group_and_zeroframe": jac,
            "segments": [], "branch": {**branch,"status":"invalid_pixel"}}
    y[(dq & 2) != 0] = np.nan
    origin = np.eye(n+1)[:n]
    good = int(np.sum(dq == 0))
    timing_one = (nframes+1.)*tf/2.
    if not suppress_one_group:
        if good == 1 and dq[0] == 0:
            branch["one_good_first_group_branch"] = True
        elif good == 0 and zeroframe is not None and zeroframe != 0:
            y[0] = float(np.float32(zeroframe)) if cast_input else float(zeroframe)
            dq[0] = 0; origin[0] = 0.; origin[0,-1] = 1.
            branch["zeroframe_branch"] = True
    segs = [] if dq[0] & 2 else _segments(dq)
    if suppress_one_group and sum(dq == 0) == 1:
        segs = []
    if not segs:
        return {"slope_double": np.nan, "jacobian_group_and_zeroframe": jac,
            "segments": [], "branch": {**branch,"status":"no_segment"}}
    # Native median uses all adjacent differences except DNU/NaN endpoints
    # and those whose SECOND group has JUMP. SAT data were made NaN above.
    if branch["one_good_first_group_branch"]:
        median = y[0]/timing_one; dm = origin[0]/timing_one
    elif branch["zeroframe_branch"]:
        median = y[0]/tf; dm = origin[0]/tf
    elif n == 1:
        median = y[0]/timing_one; dm = origin[0]/timing_one
    else:
        candidates = []
        for j in range(n-1):
            if dq[j] & 1 or dq[j+1] & (1|4) or not np.isfinite(y[j:j+2]).all():
                continue
            candidates.append(((y[j+1]/tg-y[j]/tg), (origin[j+1]-origin[j])/tg, j))
        if not candidates and not (dq[0] & 1) and np.isfinite(y[0]):
            candidates = [(y[0]/tg,origin[0]/tg,-1)]
        if not candidates:
            median = 0.; dm = np.zeros(n+1)
        else:
            candidates.sort(key=lambda v:v[0])
            k = len(candidates)//2
            chosen = [k] if len(candidates)%2 else [k-1,k]
            median = sum(candidates[k][0] for k in chosen)/len(chosen)
            dm = sum((candidates[k][1] for k in chosen),np.zeros(n+1))/len(chosen)
            # Ties across the selected order statistic make the median map
            # nondifferentiable; report instead of pretending unique weights.
            for k in chosen:
                branch["median_tie"] |= any(j != k and candidates[j][0] == candidates[k][0]
                    for j in range(len(candidates)))
    branch["median_rate_dn_s"] = float(median)
    branch["median_zero_boundary"] = bool(median == 0.)
    segment_records=[]; slopes=[]; rows=[]; variance=[]; variance_derivative=[]
    for start,end in segs:
        length=end-start; w=np.zeros(n+1); power=None; snr=None
        if length == 1:
            timing = timing_one if branch["one_good_first_group_branch"] else tf if branch["zeroframe_branch"] else tg
            w=origin[start]/timing; slope=y[start]/timing
            a=1./(timing*gain); vr=2.*rn**2/(timing**2*gain**2)
        elif length == 2:
            w=(origin[start+1]-origin[start])/tg; slope=(y[start+1]-y[start])/tg
            a=1./(tg*gain); vr=2.*rn**2/(tg**2*gain**2)
        else:
            difference=y[end-1]-y[start]; denominator=rn**2+difference*gain
            snr=max(difference*gain/np.sqrt(denominator),0.) if denominator>0 else 0.
            thresholds=np.array([5.,10.,20.,50.,100.]);powers=[0.,.4,1.,3.,6.,10.]
            power=powers[int(np.searchsorted(thresholds,snr,side="right"))]
            branch["snr_boundary"] |= bool(np.any(snr==thresholds))
            x=np.arange(start,end,dtype=float); mid=(length-1)/2.
            # C uses powf(float,float), then double accumulation. A common
            # inverse-rnoise² factor cancels from the weighted slope.
            weights=np.power(np.asarray(abs((np.arange(length)-mid)/mid),np.float32),
                             np.float32(power),dtype=np.float32).astype(float)
            sw=float(weights.sum());sx=float(weights@x);sxx=float(weights@(x*x));den=sw*sxx-sx*sx
            local=weights*(sw*x-sx)/(den*tg)
            w=local@origin[start:end]
            slope=float(local@np.nan_to_num(y[start:end],nan=0.))
            a=1./(tg*gain*(length-1));vr=12.*rn**2/(tg**2*gain**2*(length**3-length))
        v=a*(max(median,0.)+dark)+vr
        if not np.isfinite(v) or v <= 0:
            raise ValueError("Positive finite native segment variance required")
        slopes.append(slope);rows.append(w);variance.append(v)
        variance_derivative.append(a if median>0 else 0.)
        segment_records.append({"start":start,"end":end,"length":length,
            "slope_dn_s":float(slope),"variance_native_working":float(v),
            "snr":float(snr) if snr is not None else None,"power":power,
            "snr_distance_to_threshold":float(np.min(abs(np.array([5.,10.,20.,50.,100.])-snr))) if snr is not None else None})
    inv=1./np.asarray(variance);alpha=inv/inv.sum(); slopes=np.asarray(slopes)
    value=float(alpha@slopes)
    # Derivative of alpha_s=var_s^-1/sum(var^-1), including the common
    # observed/model median-rate dependence of ALL segment variances.
    coupling=float(np.sum(alpha*np.asarray(variance_derivative)*(value-slopes)/variance))
    jac=alpha@np.asarray(rows)+coupling*dm
    branch["median_weight_derivative_coefficient"] = coupling
    branch["locally_unique_derivative"] = bool(not branch["snr_boundary"]
        and (coupling==0. or not(branch["median_tie"] or branch["median_zero_boundary"])))
    branch["status"]="fixed_branch_continuous_extension"
    return {"slope_double":value,"jacobian_group_and_zeroframe":jac,
        "segments":segment_records,"branch":branch}


def piecewise_ols_arrays(groups_dn, groupdq, *, gain, readnoise_cds_dn,
                         group_time, frame_time, nframes=2, dark_dn_s=0.,
                         zeroframe=None, pixeldq=0, cast_input=False):
    """Selected-pixel/image wrapper: groups[ngroup,...], Jacobian[ngroup+1,...]."""
    y=np.asarray(groups_dn);dq=np.asarray(groupdq)
    if y.shape!=dq.shape or y.ndim<2:
        raise ValueError("Group-first arrays with one or more pixel dimensions required")
    shape=y.shape[1:];n=y.shape[0];npixel=int(np.prod(shape))
    flatten=lambda a:np.broadcast_to(np.asarray(a),shape).ravel()
    g=flatten(gain);rn=flatten(readnoise_cds_dn);dark=flatten(dark_dn_s);pdq=flatten(pixeldq)
    z=None if zeroframe is None else flatten(zeroframe)
    yy=y.reshape(n,npixel);dd=dq.reshape(n,npixel)
    values=np.empty(npixel);jac=np.empty((n+1,npixel));diagnostics=[]
    for p in range(npixel):
        r=piecewise_ols_pixel(yy[:,p],dd[:,p],gain=g[p],readnoise_cds_dn=rn[p],
            group_time=group_time,frame_time=frame_time,nframes=nframes,dark_dn_s=dark[p],
            zeroframe=None if z is None else z[p],pixeldq=int(pdq[p]),cast_input=cast_input)
        values[p]=r["slope_double"];jac[:,p]=r["jacobian_group_and_zeroframe"]
        diagnostics.append(r["branch"])
    return values.reshape(shape),jac.reshape((n+1,*shape)),diagnostics
