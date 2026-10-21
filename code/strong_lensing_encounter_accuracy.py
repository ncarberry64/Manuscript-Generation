"""Finite-window variational control of the retained frozen R1 encounter.

The background generating functions come separately from closed-FLRW
observer/lens/source boundary data. The local action is the same leading
weak-field, paraxial optical action used by the thin plane, with the
force distributed along the ray. This is a thickness/path control within
that action, not a complete scalar/GR higher-order field calculation.
"""
from __future__ import annotations

from dataclasses import dataclass, replace
from functools import lru_cache
import math

import numpy as np
from scipy.integrate import quad, solve_bvp
from scipy.optimize import brentq

C_KM_S = 299792.458
MPC_KM = 3.085677581491367e19
MPC_OVER_C_SECONDS = MPC_KM/C_KM_S
ARCSEC = np.pi/(180*3600)


@lru_cache(maxsize=8)
def _radial_rule(order):
    return np.polynomial.legendre.leggauss(order)


@dataclass(frozen=True)
class PhysicalBoundaryGeometry:
    a_l: float
    a_s: float
    curvature_per_Mpc2: float
    observer_length_Mpc: float
    source_length_Mpc: float
    source_lens_scaled_Mpc: float
    D_ol_Mpc: float
    D_os_Mpc: float
    D_ls_Mpc: float

    @classmethod
    def from_saved(cls, report, Omega_k0):
        experiment=report['experiment']
        al=1/(1+experiment['z_l']);a_s=1/(1+experiment['z_s'])
        Dol,Dos,Dls=[float(report['physical_maps_Mpc'][k][0][0])
                     for k in ('D_ol','D_os','D_ls')]
        K=-float(Omega_k0)/(C_KM_S/experiment['H0_km_s_Mpc'])**2
        if K<=0:
            raise ValueError('this study owns the declared closed reference')
        root=math.sqrt(K)
        chi_l=math.asin(root*Dol/al)/root
        chi_s=math.asin(root*Dos/a_s)/root
        inferred_Dls=a_s*math.sin(root*(chi_s-chi_l))/root
        if not math.isclose(inferred_Dls,Dls,rel_tol=3e-12):
            raise ValueError('saved maps do not share this closed reference')
        source=experiment['source_transverse_Mpc']
        if abs(source[1])>1e-25:
            raise ValueError('this control owns the three declared axial rays')
        return cls(al,a_s,K/al**2,al*chi_l,al*(chi_s-chi_l),
                   al/a_s*source[0],Dol,Dos,Dls)

    def S(self, distance):
        d=np.asarray(distance)
        return d*np.sinc(np.sqrt(self.curvature_per_Mpc2)*d/np.pi)

    def C(self, distance):
        return np.cos(np.sqrt(self.curvature_per_Mpc2)*np.asarray(distance))

    def exterior(self, half_window_Mpc):
        left=self.observer_length_Mpc-half_window_Mpc
        right=self.source_length_Mpc-half_window_Mpc
        if min(left,right)<=0:
            raise ValueError('encounter window must stay inside both endpoints')
        return (float(self.C(left)/self.S(left)),
                float(self.C(right)/self.S(right)),
                float(self.source_lens_scaled_Mpc/self.S(right)))

    @property
    def plane_coefficients(self):
        ko,ks,j=self.exterior(0.)
        return ko+ks,j

    def boundary_receipt(self):
        ko,ks,j=self.exterior(0.)
        return dict(observer_wavefront_curvature_per_Mpc=ko,
                    source_wavefront_curvature_per_Mpc=ks,
                    source_linear_term=j,
                    sum_from_wavefronts_per_Mpc=ko+ks,
                    sum_from_saved_maps_per_Mpc=self.D_os_Mpc/(self.D_ol_Mpc*self.D_ls_Mpc),
                    linear_term_from_saved_maps=self.source_lens_scaled_Mpc*self.a_s/self.a_l/self.D_ls_Mpc,
                    D_ol_from_observer_propagator_Mpc=float(self.S(self.observer_length_Mpc)),
                    source_lens_scaled_Mpc=self.source_lens_scaled_Mpc,
                    observer_length_Mpc=self.observer_length_Mpc,
                    source_length_Mpc=self.source_length_Mpc,
                    local_curvature_per_Mpc2=self.curvature_per_Mpc2)


@dataclass(frozen=True)
class RetainedEncounterField:
    mass_length_Mpc: float
    radius_Mpc: float
    rV_Mpc: float
    epsilon_infinity: float
    radial_quadrature_order: int = 32

    @classmethod
    def from_saved(cls, report, radial_quadrature_order=32):
        # The original numerical mass length uses the same declared G conversion.
        mass=report['body']['mass_Msun']
        m=4.300917270036279e-9*mass/C_KM_S**2
        return cls(m,report['body']['radius_Mpc'],report['screening']['rV_Mpc'],
                   report['screening']['unscreened_extra_fraction'],radial_quadrature_order)

    def force(self,r,scalar_on):
        r=np.asarray(r,dtype=float)
        x=r/self.radius_Mpc
        h=2.5-1.5*x*x
        with np.errstate(divide='ignore',invalid='ignore'):
            newton=np.where(x<1,self.mass_length_Mpc*r*h/self.radius_Mpc**3,
                            self.mass_length_Mpc/r**2)
            mu_over_r3=np.where(x<1,h/self.radius_Mpc**3,1/r**3)
        if scalar_on:
            extra=2*self.epsilon_infinity/(1+np.sqrt(1+self.rV_Mpc**3*mu_over_r3))
            return newton*(1+extra)
        return newton

    def extra_force(self,r):
        r=np.asarray(r,dtype=float);x=r/self.radius_Mpc
        with np.errstate(divide='ignore',invalid='ignore'):
            newton=np.where(x<1,self.mass_length_Mpc*r*(2.5-1.5*x*x)/self.radius_Mpc**3,
                            self.mass_length_Mpc/r**2)
            ratio=np.where(x<1,(2.5-1.5*x*x)/self.radius_Mpc**3,1/r**3)
        return newton*2*self.epsilon_infinity/(1+np.sqrt(1+self.rV_Mpc**3*ratio))

    def transverse_gradient(self,x,z,scalar_on):
        x=np.asarray(x,dtype=float);r=np.hypot(x,z);y=r/self.radius_Mpc
        with np.errstate(divide='ignore',invalid='ignore'):
            coefficient=np.where(y<1,self.mass_length_Mpc*(2.5-1.5*y*y)/self.radius_Mpc**3,
                                 self.mass_length_Mpc/r**3)
            ratio=np.where(y<1,(2.5-1.5*y*y)/self.radius_Mpc**3,1/r**3)
        if scalar_on:
            coefficient*=1+2*self.epsilon_infinity/(1+np.sqrt(1+self.rV_Mpc**3*ratio))
        return coefficient*x

    def _extra_radial_integral(self,r0,dr):
        nodes,weights=_radial_rule(self.radial_quadrature_order)
        r0,dr=np.broadcast_arrays(r0,dr)
        def segment(a,length):
            radii=a[None,...]+.5*length[None,...]*(nodes[:,None]+1)
            values=self.extra_force(radii)
            return .5*length*np.sum(weights[:,None]*values,axis=0)
        # Split a crossing of the material surface; both pieces retain dr's sign.
        r1=r0+dr
        crossing=(r0-self.radius_Mpc)*(r1-self.radius_Mpc)<0
        first=np.where(crossing,self.radius_Mpc-r0,dr)
        return segment(r0,first)+segment(r0+first,dr-first)

    def potential_difference(self,x1,x0,z,scalar_on):
        """Phi(hypot(x1,z))-Phi(hypot(x0,z)), without distant-tail cancellation."""
        x1,x0,z=np.broadcast_arrays(np.asarray(x1,float),np.asarray(x0,float),np.asarray(z,float))
        shape=x1.shape;x1=x1.reshape(-1);x0=x0.reshape(-1);z=z.reshape(-1)
        r1=np.hypot(x1,z);r0=np.hypot(x0,z)
        squared=(x1-x0)*(x1+x0)
        dr=np.divide(squared,r1+r0,out=np.zeros_like(squared),where=r1+r0>0)
        R=self.radius_Mpc;m=self.mass_length_Mpc
        inside=(r1<=R)&(r0<=R);outside=(r1>=R)&(r0>=R)
        answer=np.empty_like(r1)
        answer[inside]=m/R*(1.25/R**2-.375*(r1[inside]**2+r0[inside]**2)/R**4)*squared[inside]
        answer[outside]=m*dr[outside]/(r1[outside]*r0[outside])
        cross=~(inside|outside)
        def relative_R(r):
            u=(r/R)**2
            return np.where(r<R,m/R*((r-R)*(r+R)/R**2)*(.875-.375*u),
                            m*(r-R)/(R*r))
        answer[cross]=relative_R(r1[cross])-relative_R(r0[cross])
        if scalar_on:
            answer+=self._extra_radial_integral(r0,dr)
        result=answer.reshape(shape)
        return float(result) if result.ndim==0 else result

    def extra_potential_difference(self,x1,x0,z):
        x1,x0,z=np.broadcast_arrays(np.asarray(x1,float),np.asarray(x0,float),np.asarray(z,float))
        shape=x1.shape;x1=x1.reshape(-1);x0=x0.reshape(-1);z=z.reshape(-1)
        r1=np.hypot(x1,z);r0=np.hypot(x0,z)
        dr=np.divide((x1-x0)*(x1+x0),r1+r0,out=np.zeros_like(r0),where=r1+r0>0)
        result=self._extra_radial_integral(r0,dr).reshape(shape)
        return float(result) if result.ndim==0 else result

    def alpha(self,impact,half_window,scalar_on):
        b=abs(float(impact))
        if b==0:return 0.
        upper=math.atan(half_window/b)
        points=[math.acos(b/r) for r in (self.radius_Mpc,self.rV_Mpc)
                if b<r<math.hypot(b,half_window)]
        # z=b tan(t): 4 b int Phi'(r)/r dz =4 int Phi'(r) b/cos(t) dt.
        return 4*quad(lambda t: float(self.force(b/math.cos(t),scalar_on))*b/math.cos(t),
                      0,upper,points=points,epsabs=1e-21,epsrel=3e-13,limit=160)[0]

    def tail_bounds(self,b1,b0,half_window):
        """Rigorous isolated extra-field Born tail bounds, not global-model errors."""
        b1=abs(float(b1));b0=abs(float(b0));L=float(half_window)
        h1=math.hypot(L,b1);h0=math.hypot(L,b0)
        change=(b1-b0)*(b1+b0)/(h1+h0)
        potential=4*self.mass_length_Mpc*self.epsilon_infinity*abs(math.log1p(change/(L+h0)))
        alpha=4*self.mass_length_Mpc*self.epsilon_infinity*b1/(h1*(h1+L))
        return dict(alpha_radians=alpha,potential_difference_Mpc=potential,
                    quadratic_potential_bound_Mpc=self.mass_length_Mpc*self.epsilon_infinity*abs((b1-b0)*(b1+b0))/L**2)


def _positive_breaks(field,half_window,impact_scale):
    small=max(min(impact_scale,field.radius_Mpc),field.radius_Mpc*1e-3)
    if half_window<=small:return []
    return sorted(set(float(v) for v in np.geomspace(small,half_window,23)[:-1]))


def _integrate_symmetric(function,field,half_window,impact_scale,epsabs=2e-23):
    points=_positive_breaks(field,half_window,impact_scale)
    return quad(function,0,half_window,points=points,epsabs=epsabs,epsrel=2e-12,limit=300)


def thin_window_solution(geometry,field,half_window,scalar_on,seeds_rad):
    K,j=geometry.plane_coefficients
    positions=[]
    for theta in seeds_rad:
        centre=geometry.D_ol_Mpc*float(theta)
        def residual(x):return K*x-j-np.sign(x)*field.alpha(abs(x),half_window,scalar_on)
        width=geometry.D_ol_Mpc*.01*ARCSEC
        positions.append(brentq(residual,centre-width,centre+width,xtol=3e-19,rtol=1e-14))
    positions=np.asarray(positions);x0=positions[0]
    relative=[];errors=[];linear_extra=[];extra_errors=[]
    for x in positions[1:]:
        integral,error=_integrate_symmetric(lambda z:field.potential_difference(x,x0,z,scalar_on),
                                          field,half_window,min(abs(x),abs(x0)))
        action=.5*K*(x-x0)*(x+x0)-j*(x-x0)-4*integral
        relative.append(action/geometry.a_l*MPC_OVER_C_SECONDS)
        errors.append(4*error/geometry.a_l*MPC_OVER_C_SECONDS)
        extra,extra_error=_integrate_symmetric(lambda z:field.extra_potential_difference(x,x0,z),field,half_window,min(abs(x),abs(x0)),epsabs=1e-29)
        linear_extra.append(-4*extra/geometry.a_l*MPC_OVER_C_SECONDS)
        extra_errors.append(4*extra_error/geometry.a_l*MPC_OVER_C_SECONDS)
    return dict(linear_extra_response_seconds=np.array(linear_extra),linear_extra_quadrature_error_seconds=np.array(extra_errors),images_rad=positions/geometry.D_ol_Mpc,plane_positions_Mpc=positions,
                relative_delays_seconds=np.array(relative),quadrature_error_seconds=np.array(errors))


def distributed_ray(geometry,field,half_window,scalar_on,theta_seed,*,tol=2e-11):
    R=field.radius_Mpc;eps=max(field.mass_length_Mpc/R,1e-8)
    ko,ks,j=geometry.exterior(half_window)
    small=min(R*.02,half_window/20)
    positive=np.unique(np.r_[np.geomspace(small,half_window,160),R,.2*R,field.rV_Mpc])
    positive=positive[(positive>0)&(positive<=half_window)]
    z=np.r_[-positive[::-1],0.,positive]
    xplane=geometry.D_ol_Mpc*theta_seed
    _,kplane,jplane=geometry.exterior(0.)
    pout=jplane-kplane*xplane
    left=z<=0
    x=np.where(left,geometry.S(geometry.observer_length_Mpc+z)*theta_seed,
               geometry.C(z)*xplane+geometry.S(z)*pout)
    slope=np.where(left,geometry.C(geometry.observer_length_Mpc+z)*theta_seed,
                   -geometry.curvature_per_Mpc2*geometry.S(z)*xplane+geometry.C(z)*pout)
    def rhs(xi,y):
        physical_x=R*y[0];physical_z=R*xi
        return np.vstack([eps*y[1],-geometry.curvature_per_Mpc2*R**2/eps*y[0]
                          -2*R/eps*field.transverse_gradient(physical_x,physical_z,scalar_on)])
    def boundary(yl,yr):
        return np.array([yl[1]-ko*R/eps*yl[0],yr[1]+ks*R/eps*yr[0]-j/eps])
    sol=solve_bvp(rhs,boundary,z/R,np.vstack([x/R,slope/eps]),tol=tol,
                  bc_tol=tol*.1,max_nodes=80000)
    if not sol.success:raise ArithmeticError(sol.message)
    xleft=R*sol.y[0,0]
    theta=xleft/float(geometry.S(geometry.observer_length_Mpc-half_window))
    breaks=_positive_breaks(field,half_window,R*.1)
    def force_balance(z):
        xx=R*sol.sol(z/R)[0]
        return geometry.curvature_per_Mpc2*xx+2*field.transverse_gradient(xx,z,scalar_on)
    left_integral=quad(force_balance,-half_window,0,points=[-v for v in breaks],
                       epsabs=1e-22,epsrel=2e-12,limit=300)[0]
    right_integral=quad(force_balance,0,half_window,points=breaks,
                        epsabs=1e-22,epsrel=2e-12,limit=300)[0]
    momentum_residual=eps*(sol.y[1,-1]-sol.y[1,0])+left_integral+right_integral
    return dict(solution=sol,images_rad=theta,nodes=sol.x.size,
                integrated_momentum_residual_rad=float(momentum_residual),
                max_collocation_relative_residual=float(np.max(sol.rms_residuals)),
                max_scaled_boundary_residual=float(np.max(abs(boundary(sol.y[:,0],sol.y[:,-1])))))


def distributed_window_solution(geometry,field,half_window,scalar_on,seeds_rad,*,tol=2e-11):
    rays=[distributed_ray(geometry,field,half_window,scalar_on,float(t),tol=tol) for t in seeds_rad]
    R=field.radius_Mpc;eps=max(field.mass_length_Mpc/R,1e-8)
    ko,ks,j=geometry.exterior(half_window)
    reference=rays[0]['solution']
    relative=[];errors=[];linear_extra=[];extra_errors=[]
    for ray in rays[1:]:
        solution=ray['solution']
        def interior(z):
            yi=solution.sol(z/R);y0=reference.sol(z/R)
            xi=R*yi[0];x0=R*y0[0];pi=eps*yi[1];p0=eps*y0[1]
            return .5*((pi-p0)*(pi+p0)-geometry.curvature_per_Mpc2*(xi-x0)*(xi+x0))-2*field.potential_difference(xi,x0,z,scalar_on)
        points=_positive_breaks(field,half_window,R*.1)
        value1,error1=quad(interior,-half_window,0,points=[-v for v in points],epsabs=3e-23,epsrel=3e-12,limit=300)
        value2,error2=quad(interior,0,half_window,points=points,epsabs=3e-23,epsrel=3e-12,limit=300)
        xl=R*solution.y[0,0];x0l=R*reference.y[0,0]
        xr=R*solution.y[0,-1];x0r=R*reference.y[0,-1]
        outside=.5*ko*(xl-x0l)*(xl+x0l)+.5*ks*(xr-x0r)*(xr+x0r)-j*(xr-x0r)
        relative.append((outside+value1+value2)/geometry.a_l*MPC_OVER_C_SECONDS)
        errors.append((error1+error2)/geometry.a_l*MPC_OVER_C_SECONDS)
        def extra_integrand(z):
            return field.extra_potential_difference(R*solution.sol(z/R)[0],R*reference.sol(z/R)[0],z)
        extra1,extra_error1=quad(extra_integrand,-half_window,0,points=[-v for v in points],epsabs=1e-29,epsrel=2e-12,limit=300)
        extra2,extra_error2=quad(extra_integrand,0,half_window,points=points,epsabs=1e-29,epsrel=2e-12,limit=300)
        linear_extra.append(-2*(extra1+extra2)/geometry.a_l*MPC_OVER_C_SECONDS)
        extra_errors.append(2*(extra_error1+extra_error2)/geometry.a_l*MPC_OVER_C_SECONDS)
    return dict(linear_extra_response_seconds=np.asarray(linear_extra),linear_extra_quadrature_error_seconds=np.asarray(extra_errors),images_rad=np.array([r['images_rad'] for r in rays]),
                relative_delays_seconds=np.asarray(relative),quadrature_error_seconds=np.asarray(errors),
                nodes=[r['nodes'] for r in rays],
                max_integrated_momentum_residual_rad=max(abs(r['integrated_momentum_residual_rad']) for r in rays),
                max_collocation_relative_residual=max(r['max_collocation_relative_residual'] for r in rays),
                max_scaled_boundary_residual=max(r['max_scaled_boundary_residual'] for r in rays))


def scalar_response_homotopy(geometry,field,half_window,seeds_rad,*,orders=(3,5),tol=2e-11):
    """Integrate exact on-shell optical response between the fixed two fields.

    eta parametrizes Phi_N+eta delta_Phi only as a numerical path integral;
    no intermediate profile is asserted to solve another scalar field theory.
    Stationarity removes the ray variation. This avoids subtraction of the
    large ordinary-gravity arrival terms when computing the finite response.
    """
    cache={};answers=[]
    for order in orders:
        nodes,weights=np.polynomial.legendre.leggauss(order)
        thin_total=np.zeros(2);distributed_total=np.zeros(2)
        records=[]
        for node,weight in zip(nodes,weights):
            eta=float((node+1)/2);weight=float(weight/2)
            if eta not in cache:
                current=replace(field,epsilon_infinity=eta*field.epsilon_infinity)
                thin=thin_window_solution(geometry,current,half_window,True,seeds_rad)
                ray=distributed_window_solution(geometry,current,half_window,True,thin['images_rad'],tol=tol)
                cache[eta]={'eta_optical_integration_coordinate':eta,
                            'thin_on_shell_derivative_seconds':thin['linear_extra_response_seconds']/eta,
                            'distributed_on_shell_derivative_seconds':ray['linear_extra_response_seconds']/eta,
                            'max_integrated_momentum_residual_rad':ray['max_integrated_momentum_residual_rad']}
            item=cache[eta]
            thin_total+=weight*item['thin_on_shell_derivative_seconds']
            distributed_total+=weight*item['distributed_on_shell_derivative_seconds']
            records.append(dict(item,weight=weight))
        answers.append({'Gauss_order':order,'nodes':records,
                        'thin_finite_offsets_seconds':thin_total,
                        'distributed_finite_offsets_seconds':distributed_total,
                        'finite_thickness_change_seconds':distributed_total-thin_total})
    return {'definition':'Integrate stationary derivative for Phi_N+eta delta_Phi from eta0 toeta1, fixed physical endpoints/clock and fixed given delta_Phi; eta is not a fitted theory or state parameter.',
            'orders':answers,
            'change_between_orders_seconds':answers[-1]['finite_thickness_change_seconds']-answers[0]['finite_thickness_change_seconds'],
            'physical_intermediate_states_claimed':False}
