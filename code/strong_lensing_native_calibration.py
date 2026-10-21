"""Conditional native history with the source-calibrated IPC/linearity order.

This supplies an executable mean and analytic directional chain, including the
implicitly fitted first-group electronic intercept. It identifies neither a
physical post-full-well law nor photon/background partition. No likelihood,
noise scaling, optical-PSF convolution or physical admission occurs here.
"""
from __future__ import annotations
import numpy as np
from scipy import sparse
from scipy.sparse.linalg import splu
from strong_lensing_stochastic_response import collection_boundary_generator,nearest_neighbor_edges
from strong_lensing_bounded_collection import bounded_leading_collection
from strong_lensing_group_calibration import polynomial_derivatives,inverse_linearity
from strong_lensing_calibration_order import first_group_preipc_intercept,ipc_averaged_read_observation


class FixedHarmonicContext:
    """A declared fixed-mask positive Dirichlet trial and its exact linear action.

    Trusted values remain measured/model values. Unknown values are harmonic
    continuation on the finite graph; they are not zero photons or measured
    full wells. The mask is fixed at initialization. A changed trusted value
    leaving its positive domain is rejected rather than silently reselected.
    """
    def __init__(self,shape,trusted):
        self.shape=tuple(shape);self.known=np.asarray(trusted,bool).ravel()
        if len(self.known)!=int(np.prod(shape)) or not self.known.any():
            raise ValueError('At least one trusted fixed context value required')
        self.k=np.flatnonzero(self.known);self.u=np.flatnonzero(~self.known)
        self.cross=None;self.factor=None
        if len(self.u):
            lap=collection_boundary_generator(np.ones(len(self.known)),.5,nearest_neighbor_edges(shape))
            self.cross=lap[self.u][:,self.k]
            self.factor=splu(lap[self.u][:,self.u].tocsc())

    def action(self,values):
        value=np.asarray(values,float);answer=value.copy()
        if value.shape[0]!=len(self.known) or not np.isfinite(value[self.k]).all():
            raise ValueError('Matching finite trusted context values required')
        if len(self.u):answer[self.u]=self.factor.solve(-(self.cross@value[self.k]))
        return answer

    def positive(self,values):
        value=np.asarray(values,float)
        if np.any(value[self.k]<0):
            raise ValueError('Fixed trusted photon/ceiling context left its positive domain')
        result=self.action(value)
        if not np.isfinite(result).all() or np.any(result<0):
            raise ValueError('Harmonic context continuation left its positive domain')
        return result


class CalibrationOrderedHistory:
    """One actual native patch, fixed observed branches and conditional boundary.

    Effective photon rate is in e/s AFTER the low-charge IPC mean already in
    the empirical ePSF. Its preIPC rate is gain*G*(effective/gain). Source and
    detector-plane rates MUST be passed separately; this class consumes only
    the stated photon source. Physical caps are explicit translated trials.
    Native GROUPDQ/selection and deposited cosmic-ray history are not solved.
    """
    def __init__(self,history,patch,deconvolution,forward_ipc,
                 initial_effective_rate_e_s,*,ceiling_multiplier=1.,
                 prescription='storage_ceiling'):
        h=history;p=patch
        self.shape=h['gain_e_per_DN'][p].shape;self.n=int(np.prod(self.shape))
        self.g=sparse.csr_matrix(deconvolution);self.h=sparse.csr_matrix(forward_ipc)
        if self.g.shape!=(self.n,self.n) or self.h.shape!=self.g.shape:
            raise ValueError('IPC maps and actual padded context must match')
        self.gain=np.asarray(h['gain_e_per_DN'][p],float).ravel()
        self.coeff=np.asarray(h['linearity_coefficients'][p],float).reshape((-1,self.n)).copy()
        flags=np.asarray(h['pixeldq'][p]).ravel()
        self.identity=((flags&1048576)!=0)|~np.isfinite(self.coeff).all(axis=0)
        self.coeff[:,self.identity]=0.;self.coeff[1,self.identity]=1.
        self.times=np.asarray(h['native_read_times_s'],float)
        self.frame=float(h['frame_time_s'])
        self.alltimes=np.r_[self.frame,self.times.ravel()]
        self.order=np.argsort(self.alltimes,kind='stable');self.inverse=np.argsort(self.order)
        self.dark=np.asarray(h['dark_RAPID42_DN'][p][h['native_dark_RAPID_frame_indices']],float).reshape((9,2,self.n))
        self.darkgroup=np.asarray(h['native_subtracted_dark_group_DN'][p],float).reshape((9,self.n))
        native_adc=np.asarray(h['pre_linearity_group_DN'][p],float).reshape((9,self.n))
        self.anchor=(self.g@native_adc.T).T
        self.target0=self.anchor[0]
        self.initial_linear0=polynomial_derivatives(self.coeff,self.target0)[0]
        rawoffset=np.asarray(h['raw_group_DN'][p,0],float).ravel()-native_adc[0]
        threshold=np.asarray(h['saturation_threshold_DN'][p],float).ravel()-rawoffset
        self.threshold_pre=self.g@threshold
        self.threshold_linear=polynomial_derivatives(self.coeff,self.threshold_pre)[0]
        self.nativezero_present=np.asarray(h['native_corrected_zeroframe_DN'][p]).ravel()!=0
        self.edges=nearest_neighbor_edges(self.shape)
        self.cap_multiplier=float(ceiling_multiplier);self.prescription=prescription
        if not np.isfinite(self.cap_multiplier) or self.cap_multiplier<=0:
            raise ValueError('Positive declared ceiling multiplier required')
        initial=np.asarray(initial_effective_rate_e_s,float).ravel()
        initial_pre=self.gain*(self.g@(initial/self.gain))
        self.rate_context=FixedHarmonicContext(self.shape,np.isfinite(initial_pre)&(initial_pre>=0))
        initial_pre=self.rate_context.positive(initial_pre)
        caps=self._raw_caps(initial_pre)
        self.cap_context=FixedHarmonicContext(self.shape,np.isfinite(caps)&(caps>0))
        self.cap_context.positive(caps)

    def _initial_beta(self,rate):
        return self.initial_linear0-rate/self.gain*self.times[0].mean()-self.dark[0].mean(axis=0)

    def _raw_caps(self,rate):
        return self.cap_multiplier*self.gain*(self.threshold_linear-self._initial_beta(rate)-self.dark[0].mean(axis=0))

    def evaluate(self,effective_rate_e_s,kappa_e_inverse=0.,*,derivatives=False):
        effective=np.asarray(effective_rate_e_s,float).ravel()
        raw_pre=self.gain*(self.g@(effective/self.gain))
        rate=self.rate_context.positive(raw_pre)
        caps=self.cap_context.positive(self._raw_caps(rate))
        physical=bounded_leading_collection(rate,self.alltimes[self.order],kappa_e_inverse,caps,self.edges,
             prescription=self.prescription,ceiling_source='Nominal IPC-deconvolved initial-intercept SAT-cutoff trial, fixed harmonic ceiling mask; not identified physical full well.',derivatives=derivatives)
        mu=physical['observed_read_mean_e'][self.inverse]
        signal=mu[1:].reshape((9,2,self.n))/self.gain+self.dark
        fit=first_group_preipc_intercept(signal[0],self.coeff,self.target0,self._initial_beta(rate))
        beta=fit['beta_preipc_DN']
        observed=ipc_averaged_read_observation(signal+beta[None,None],self.coeff,
             np.repeat(self.anchor[:,None],2,axis=1),self.h,dark_group_dn=self.darkgroup)
        qzero=beta+mu[0]/self.gain+self.dark[0,0]
        zero_adc=inverse_linearity(qzero,self.coeff,self.target0)
        _,zero_pre_d,_=polynomial_derivatives(self.coeff,zero_adc)
        zero,zero_native_d,_=polynomial_derivatives(self.coeff,self.h@zero_adc)
        zero=np.where(self.nativezero_present,zero,0.)
        state={'groups_DN':observed['mean_DN'].reshape((9,*self.shape)),
               'zeroframe_DN':zero.reshape(self.shape),'rate_preipc_e_s':rate,
               'uncontinued_rate_preipc_e_s':raw_pre,'ceiling_e':caps,
               'beta_preipc_DN':beta,'physical':physical,'observation':observed,
               'intercept':fit,'inverse_zero_linearity_derivative':1/zero_pre_d,
               'native_zero_linearity_derivative':zero_native_d,
               'derivatives_available':bool(derivatives),
               'diagnostics':{'prescription':self.prescription,'ceiling_multiplier':self.cap_multiplier,
                   'ceiling_from_initial_intercept':True,'rate_context_unknown_count':len(self.rate_context.u),
                   'ceiling_context_unknown_count':len(self.cap_context.u),
                   'identity_linearity_context_pixels':int(self.identity.sum()),
                   'readwise_dark_stage':'Nominal RAPID reference assigned to preIPC linearized read; undocumented for nonzero hot reference values.',
                   'group0_inverse_ADC_match_DN':fit['first_group_inverse_ADC_match_max_abs_DN'],
                   'fixed_observed_DQ_and_jump_selection':True,'cosmic_ray_deposits_reconstructed':False,
                   'physical_full_well_identified':False,'extra_PSF_convolution':False,
                   'photon_partition':'Caller photon source only; additive electronic rate plane excluded.',
                   'physical_validation':False,'error_rescaling':False}}
        return state

    def project_directions(self,state,effective_rate_directions,kappa_directions,
                           selected_pixels,native_group_weights,native_zero_weights):
        """Exact local chain to caller-supplied current-branch native slope map.

        d rate has padded-context-pixel,parameter shape. Current OLS derivatives
        are consumed, never fabricated uniform weights. This differentiates
        preIPC harmonic continuations, translated caps, charge mean, implicit
        beta, inverse L, H, native L and ZEROFRAME in their actual order.
        Native adaptive branch boundary behavior remains the estimator's domain.
        """
        if not state['derivatives_available']:raise ValueError('Evaluate with derivatives first')
        direction=np.asarray(effective_rate_directions,float)
        dk=np.asarray(kappa_directions,float)
        selected=np.asarray(selected_pixels,int)
        weights=np.asarray(native_group_weights,float);wzero=np.asarray(native_zero_weights,float)
        if (direction.ndim!=2 or direction.shape[0]!=self.n or dk.shape!=(direction.shape[1],)
                or weights.shape!=(9,len(selected)) or wzero.shape!=(len(selected),)):
            raise ValueError('Matching parameter directions and actual selected native weights required')
        pre=self.gain[:,None]*(self.g@(direction/self.gain[:,None]))
        dr=self.rate_context.action(pre)
        dc=self.cap_context.action(self.cap_multiplier*self.times[0].mean()*dr)
        physical=state['physical']
        dm=[]
        for i in self.inverse:
            dm.append(physical['observed_rate_jacobian_s'][i]@dr
                      +physical['observed_cap_jacobian'][i]@dc
                      +physical['observed_kappa_jacobian_e2'][i,:,None]*dk[None,:])
        beta_weights=state['intercept']['beta_signal_derivative']
        dbeta=(beta_weights[0,:,None]*dm[1]+beta_weights[1,:,None]*dm[2])/self.gain[:,None]
        hselected=self.h[selected]
        result=np.zeros((len(selected),direction.shape[1]))
        for j in range(9):
            da=0.
            for r in range(2):
                dq=dbeta+dm[1+2*j+r]/self.gain[:,None]
                da=da+state['observation']['inverse_linearity_derivative'][j,r,:,None]*dq/2
            derivative=(hselected@da)*state['observation']['native_linearity_derivative'][j,selected,None]
            result+=weights[j,:,None]*derivative
        dz=dbeta+dm[0]/self.gain[:,None]
        derivative=(hselected@(state['inverse_zero_linearity_derivative'][:,None]*dz))*state['native_zero_linearity_derivative'][selected,None]
        derivative*=self.nativezero_present[selected,None]
        result+=wzero[:,None]*derivative
        return result
