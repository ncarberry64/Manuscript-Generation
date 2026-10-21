"""Source-fixed stationary READ+DARK covariance, including cross-amplifier modes.

An independent dark integration difference is an observed filter template z.
The additional conditional law is epsilon_c(p)=sum_u z_c(u)*eta(p-u)/sqrt(N),
where eta is one common unit-variance stationary innovation field for all
amplifier/group channels c.  Its covariance is the measured circular spatial
cross-correlation.  Shared innovations retain cross-stamp and cross-amplifier
correlations.  No diagonal white/reference noise is added to this template.

This is an explicit periodic/stationary calibration model, not an identified
unrestricted detector covariance or a BHSM optical-state statistical closure.
The fixed calibration mask is a spectral data window; its convolution is not
silently deconvolved. Finite epoch and mask/state uncertainty remain separate.
"""
from __future__ import annotations
import numpy as np
from scipy import fft


def template_spectra(templates, valid):
    """Orthonormal Fourier filters, with explicit masked-window normalization.

    templates: amplifier,group,row,clock-column; valid: amplifier,row,column.
    Zeroing unobserved samples here is an explicit periodogram data window, not
    a physical zero-flux or zero-noise declaration at those detector pixels.
    Divide by sqrt(mean(mask)) so lag-zero variances equal valid-pixel moments.
    """
    z = np.asarray(templates, float)
    m = np.asarray(valid, bool)
    if z.ndim != 4 or m.shape != (z.shape[0], *z.shape[2:]) or np.any(m.mean(axis=(1, 2)) == 0):
        raise ValueError("Amplifier/group field and nonempty explicit calibration windows required")
    if any(not np.isfinite(z[a, :, m[a]]).all() for a in range(len(m))):
        raise ValueError("Finite actual calibration values inside every window required")
    coverage = m.mean(axis=(1, 2))
    field = np.where(m[:, None], z, 0.)/np.sqrt(coverage)[:, None, None, None]
    return fft.rfft2(field, axes=(-2, -1), norm="ortho"), coverage


class SpectralTemplateModel:
    """Fixed measured filter; only the caller's observation weights vary.

    Gaussian scoring is a conditional electronics/noise approximation. It does
    not factorize or replace the BHSM prepared optical state's correlations.
    """
    def __init__(self, spectra, spatial_shape):
        self.spectra = spectra
        self.shape = tuple(map(int, spatial_shape))
        if (len(self.shape) != 2 or min(self.shape) < 1 or len(spectra.shape) != 4
                or spectra.shape[-2:] != (self.shape[0], self.shape[1]//2+1)):
            raise ValueError("Amplifier/group orthonormal rFFT filters and exact spatial cell required")
        self.namp, self.ngroup = spectra.shape[:2]

    def cross_kernel(self, amp, group, other_amp, other_group):
        return fft.irfft2(np.asarray(self.spectra[amp, group])*
                          np.conj(self.spectra[other_amp, other_group]), s=self.shape)

    def project(self, requests, progress=None):
        """Project multiple output sets while reusing each source cross spectrum.

        Each request has amp[n], clock_xy[n,2], group_weights[n,ngroup].
        Weights carry actual native L', estimator, flat and photometric units.
        Requests are distinct exposures: no noise realization is shared across
        requests. Each returned block retains all within-exposure cross-stamp
        and cross-amplifier terms. Common calibration uncertainty across requests
        is NOT thereby known and is not silently added to them.
        """
        prepared = []
        for request in requests:
            amp = np.asarray(request["amp"], int)
            xy = np.asarray(request["clock_xy"], int)
            weight = np.asarray(request["group_weights"], float)
            if (amp.ndim != 1 or xy.shape != (len(amp), 2)
                    or weight.shape != (len(amp), self.ngroup) or np.any(amp < 0)
                    or np.any(amp >= self.namp) or not np.isfinite(weight).all()
                    or np.any(xy < 0) or np.any(xy[:, 0] >= self.shape[1])
                    or np.any(xy[:, 1] >= self.shape[0])):
                raise ValueError("Explicit in-cell output coordinates/amplifiers and native group weights required")
            groups = []
            for a in range(self.namp):
                idx = np.flatnonzero(amp == a)
                groups.append({"idx": idx, "xy": xy[idx], "weights": weight[idx]})
            prepared.append(groups)
        answer = [np.zeros((len(r["amp"]),)*2) for r in requests]
        # Cross kernels are computed only when actually consumed; requests may
        # restrict to early groups without manufacturing nine-group support.
        channels = [(a, g) for a in range(self.namp) for g in range(self.ngroup)]
        count = 0
        for c, (a, g) in enumerate(channels):
            for d in range(c, len(channels)):
                b, h = channels[d]
                useful = [(k, rows[a], rows[b]) for k, rows in enumerate(prepared)
                          if len(rows[a]["idx"]) and len(rows[b]["idx"])
                          and np.any(rows[a]["weights"][:, g] != 0)
                          and np.any(rows[b]["weights"][:, h] != 0)]
                if not useful:
                    continue
                kernel = self.cross_kernel(a, g, b, h)
                for k, left, right in useful:
                    dy = (left["xy"][:, 1, None]-right["xy"][None, :, 1]) % self.shape[0]
                    dx = (left["xy"][:, 0, None]-right["xy"][None, :, 0]) % self.shape[1]
                    block = (left["weights"][:, g, None]*right["weights"][None, :, h]*kernel[dy, dx])
                    answer[k][np.ix_(left["idx"], right["idx"])] += block
                    if c != d:
                        answer[k][np.ix_(right["idx"], left["idx"])] += block.T
                count += 1
                if progress is not None and count % 30 == 0:
                    progress(count)
        return [(v+v.T)/2 for v in answer]


def logarithmic_frequency_bins(length, *, real_half=False):
    """Fixed signed powers-of-two frequency cells; DC/Nyquist kept separate."""
    if length < 2 or int(length) != length:
        raise ValueError("Integer Fourier cell length required")
    freq = np.rint(np.fft.rfftfreq(length)*length if real_half
                   else np.fft.fftfreq(length)*length).astype(int)
    keys = []
    for value in freq:
        if value == 0:
            keys.append((0, 0))
        elif length % 2 == 0 and abs(value) == length//2:
            keys.append((2, 0))
        else:
            keys.append((1 if value > 0 else -1, int(np.floor(np.log2(abs(value))))))
    unique = sorted(set(keys))
    return np.array([unique.index(k) for k in keys]), unique


def pool_spectral_matrices(spectra):
    """Compact PSD cross spectra from prespecified frequency-bin pooling.

    Every bin is an average of actual Fourier-vector outer products across
    amplifier/group channels. Thus PSD is constructive, without eigenvalue
    clipping; no diagonal variance changes because frequency sums are retained.
    Bins are an explicit stationarity/smoothing model, not instrument truth.
    """
    f = np.asarray(spectra)
    if f.ndim != 4:
        raise ValueError("Amplifier/group orthonormal rFFT template array required")
    namp, ngroup, ny, nxhalf = f.shape
    nx = 2*(nxhalf-1)
    ybin, ykeys = logarithmic_frequency_bins(ny)
    xbin, xkeys = logarithmic_frequency_bins(nx, real_half=True)
    vectors = f.reshape(namp*ngroup, ny, nxhalf)
    pooled = np.empty((len(ykeys), len(xkeys), namp*ngroup, namp*ngroup), complex)
    counts = np.empty((len(ykeys), len(xkeys)), int)
    for iy in range(len(ykeys)):
        rows = np.flatnonzero(ybin == iy)
        for ix in range(len(xkeys)):
            cols = np.flatnonzero(xbin == ix)
            sample = vectors[:, rows[:, None], cols].reshape(namp*ngroup, -1)
            counts[iy, ix] = sample.shape[1]
            cross = sample@sample.conj().T/sample.shape[1]
            pooled[iy, ix] = (cross+cross.conj().T)/2
    return {"spectral_matrices": pooled, "y_frequency_bin": ybin,
            "x_frequency_bin": xbin, "bin_frequency_counts": counts,
            "spatial_shape": np.array([ny, nx]), "namp": namp, "ngroup": ngroup,
            "y_bin_keys": ykeys, "x_bin_keys": xkeys}


class BinnedSpectralModel(SpectralTemplateModel):
    """Same full projection using compact measured36-channel PSD matrices."""
    def __init__(self, spectral_matrices, y_frequency_bin, x_frequency_bin,
                 spatial_shape, namp, ngroup):
        self.shape = tuple(map(int, spatial_shape))
        self.namp, self.ngroup = int(namp), int(ngroup)
        self.matrices = np.asarray(spectral_matrices)
        self.ybin = np.asarray(y_frequency_bin, int)
        self.xbin = np.asarray(x_frequency_bin, int)
        if (self.ybin.shape != (self.shape[0],) or self.xbin.shape != (self.shape[1]//2+1,)
                or self.matrices.shape != (self.ybin.max()+1, self.xbin.max()+1,
                                           self.namp*self.ngroup, self.namp*self.ngroup)):
            raise ValueError("Exact compact spectral calibration dimensions required")

    def cross_kernel(self, amp, group, other_amp, other_group):
        c, d = amp*self.ngroup+group, other_amp*self.ngroup+other_group
        spectrum = self.matrices[self.ybin[:, None], self.xbin[None, :], c, d]
        return fft.irfft2(spectrum, s=self.shape)


def anchored_group0_residual_jacobian(observed_rate_jacobian, predicted_rate_jacobian,
                                      predicted_group_intercept_jacobian):
    """PIXELWISE delta-method map when noisy G0 fixes a local intercept.

    dRpred/dG0=(sum_h Wpred_h*dGpred_h/dbeta)/(dGpred_0/dbeta).
    Residual R_obs-R_pred(theta,G0_obs) has J=Wobs-e0*dRpred/dG0.
    This is NOT a full conditional read likelihood or a noiseless boundary.
    Observation ADC-to-corrected-group L' is applied after this group-DN map.
    With IPC the actual anchor may have H*diag(ratio)*G spatial derivatives;
    that nonlocal map MUST use full_anchored_residual_maps instead. A diagonal
    substitute is not implied by this helper or by fixed STATIC covariance.
    """
    obs, pred, beta = (np.asarray(x, float) for x in (
        observed_rate_jacobian, predicted_rate_jacobian, predicted_group_intercept_jacobian))
    if (obs.ndim != 2 or pred.shape != obs.shape or beta.shape != obs.shape
            or any(not np.isfinite(x).all() for x in (obs, pred, beta)) or np.any(beta[:, 0] == 0)):
        raise ValueError("Actual observed/predicted estimator and group-intercept derivatives required")
    result = obs.copy()
    result[:, 0] -= np.sum(pred*beta, axis=1)/beta[:, 0]
    return result


def full_anchored_residual_maps(observed_noise_maps, prediction_boundary_jacobian,
                                boundary_noise_maps):
    """Preserve the actual NONLOCAL boundary/IPC derivative of a cut estimator.

    J_res,g = J_observed,g - J_prediction,boundary * J_boundary,noise,g.
    All matrices must be supplied from the actual original application, e.g.
    the H diag(<inverseL'>_g/<inverseL'>_0) G relation. No diagonal, Gaussian
    optical-state closure or inverse boundary map is silently substituted.
    """
    prediction = np.asarray(prediction_boundary_jacobian, float)
    observed = [np.asarray(v, float) for v in observed_noise_maps]
    boundary = [np.asarray(v, float) for v in boundary_noise_maps]
    if (prediction.ndim != 2 or len(observed) != len(boundary)
            or not np.isfinite(prediction).all()):
        raise ValueError("Actual prediction/boundary derivative and matched noise maps required")
    result = []
    for obs, anchor in zip(observed, boundary):
        if (obs.ndim != 2 or anchor.ndim != 2 or obs.shape != (prediction.shape[0], anchor.shape[1])
                or anchor.shape[0] != prediction.shape[1]
                or not np.isfinite(obs).all() or not np.isfinite(anchor).all()):
            raise ValueError("Nonlocal anchor and observed noise-map domains must match")
        result.append(obs-prediction@anchor)
    return result


def native_to_clock_coordinates(native_xy):
    """NRCB3 native coordinates to the fixed2040x512 read-clock cell.

    Native FASTAXIS=1 alternates direction across four512-column amplifiers;
    SLOWAXIS=-2 reverses row order. Only the actual outer detector reference
    pixels are excluded; INTERNAL amplifier boundaries contain science pixels.
    This convention must match spectra, including patches crossing a boundary.
    """
    xy = np.asarray(native_xy, int)
    if xy.ndim != 2 or xy.shape[1] != 2:
        raise ValueError("Native output x/y pairs required")
    amp = xy[:, 0]//512
    local_x = xy[:, 0] % 512
    clock_x = np.where(amp % 2 == 0, local_x, 511-local_x)
    clock_y = 2047-xy[:, 1]-4
    if (np.any(xy[:, 0] < 4) or np.any(xy[:, 0] >= 2044) or np.any(amp < 0)
            or np.any(amp > 3) or np.any(clock_x < 0)
            or np.any(clock_x >= 512) or np.any(clock_y < 0) or np.any(clock_y >= 2040)):
        raise ValueError("Output is outside the declared calibrated read-clock interior")
    return amp, np.c_[clock_x, clock_y]


def dark_field_to_clock_templates(field, valid):
    """Same native-to-clock convention for actual full-frame calibration data."""
    z = np.asarray(field)
    m = np.asarray(valid, bool)
    if z.ndim != 3 or z.shape[1:] != (2048, 2048) or m.shape != (2048, 2048):
        raise ValueError("Group/full-native-detector calibration field and mask required")
    templates, masks = [], []
    for amp in range(4):
        a = z[:, 4:2044, amp*512:(amp+1)*512][:, ::-1]
        b = m[4:2044, amp*512:(amp+1)*512][::-1].copy()
        if amp == 0: b[:, :4] = False
        if amp == 3: b[:, -4:] = False
        if amp % 2:
            a = a[:, :, ::-1]; b = b[:, ::-1]
        templates.append(a); masks.append(b)
    return np.asarray(templates), np.asarray(masks)
