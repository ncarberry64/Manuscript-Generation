"""Reproduce the first matched strong-lensing research checkpoint."""
from __future__ import annotations
import os
os.environ.setdefault('OMP_NUM_THREADS','1')
os.environ.setdefault('OPENBLAS_NUM_THREADS','1')
from pathlib import Path
import argparse
import hashlib
import json
import subprocess
import sys
import numpy as np

ROOT=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT/'code'))
from strong_lensing_transport import ActionNativeAxialTransport,ANCHOR_NAMES
from strong_lensing_comparison import (
    ARCSEC,C_KM_S,LensGeometry,GalaxyPotential,solve_images,
    complete_map_response,observable_vector,response_vector,nuisance_columns,
    nuisance_projected_information,
)
from strong_lensing_envelope import checkpoint_snapshot


def serial(value):
    if isinstance(value,np.ndarray):return value.tolist()
    if isinstance(value,np.generic):return value.item()
    if isinstance(value,dict):return {k:serial(v) for k,v in value.items()}
    if isinstance(value,(tuple,list)):return [serial(v) for v in value]
    return value


def dump(path,data):
    path.parent.mkdir(parents=True,exist_ok=True)
    path.write_text(json.dumps(serial(data),indent=2,sort_keys=True,allow_nan=False)+'\n',encoding='utf8',newline='\n')


def geometry_from_scalars(scalars,hubble_distance,z_l=.3):
    return LensGeometry(*(np.eye(2)*d*hubble_distance for d in scalars),z_l)


def main(output,make_plots=True):
    output.mkdir(parents=True,exist_ok=True)
    provider=ActionNativeAxialTransport()
    transport=provider.response()
    seven_state_control=provider.seven_state_control(transport)
    baseline=provider.maps()
    H0=67.4
    hubble_distance=C_KM_S/H0
    geometry=geometry_from_scalars(baseline.scalar_maps,hubble_distance)
    galaxy=GalaxyPotential(4*np.pi*(220/C_KM_S)**2)
    theta_E=float(geometry.G[0,0]*galaxy.b)
    beta=np.array([.2*theta_E,0.])
    source=geometry.D_os@beta
    solution=solve_images(geometry,galaxy,source)
    coarse_maps=np.array(transport['d_maps_d_anchor'])
    fine_maps=[];responses=[];checks=[]
    sigma=np.r_[np.full(solution['images_rad'].size,.001),
        np.full(len(solution['relative_delays_days']),.1)]
    for i,trial in enumerate(transport['finite_trials']):
        h=trial['step']/2
        plus=provider.maps(np.eye(6)[i]*h)
        minus=provider.maps(-np.eye(6)[i]*h)
        maps=(plus.scalar_maps-minus.scalar_maps)/(2*h)
        fine_maps.append(maps)
        d=np.array([np.eye(2)*v*hubble_distance for v in maps])
        response=complete_map_response(solution['images_rad'],geometry,galaxy,source,d)
        responses.append(response)
        p=solve_images(geometry_from_scalars(plus.scalar_maps,hubble_distance),galaxy,source)
        m=solve_images(geometry_from_scalars(minus.scalar_maps,hubble_distance),galaxy,source)
        finite=(observable_vector(p)-observable_vector(m))/(2*h)
        prediction=response_vector(response)
        coarse_p=solve_images(geometry_from_scalars(trial['plus_maps'],hubble_distance),galaxy,source)
        coarse_m=solve_images(geometry_from_scalars(trial['minus_maps'],hubble_distance),galaxy,source)
        coarse=(observable_vector(coarse_p)-observable_vector(coarse_m))/(4*h)
        checks.append(dict(coordinate=ANCHOR_NAMES[i],half_step=h,
            map_derivative_half_step_absolute_max=float(np.max(abs(maps-coarse_maps[:,i]))),
            independently_resolved_observable_derivative=finite,
            analytic_observable_derivative=prediction,
            analytic_vs_half_step_whitened_relative=float(np.linalg.norm((finite-prediction)/sigma)/max(np.linalg.norm(prediction/sigma),1e-12)),
            coarse_vs_half_step_whitened_relative=float(np.linalg.norm((finite-coarse)/sigma)/max(np.linalg.norm(finite/sigma),1e-12)),
            measured_redshift_errors=[plus.lens['observed_z']-.3,plus.source['observed_z']-1.,minus.lens['observed_z']-.3,minus.source['observed_z']-1.],
            plus_reciprocity=plus.reciprocity,minus_reciprocity=minus.reciprocity))
    fine_maps=np.array(fine_maps).T
    B=np.column_stack([response_vector(r) for r in responses])
    N,nuisance_names=nuisance_columns(geometry,galaxy,source,solution)
    covariance=np.diag(sigma**2)
    free=nuisance_projected_information(B,N,covariance)
    raw=nuisance_projected_information(B,np.empty((len(sigma),0)),covariance)
    # Synthetic prior control: extra independent measurements on nuisance
    # coordinates, not information attributed to these lens observables.
    prior_sigma=np.array([.01,.01,.02,.005,.005,.05,.02,.02])
    B_aug=np.vstack((B,np.zeros((len(nuisance_names),6))))
    N_aug=np.vstack((N,np.diag(1/prior_sigma)))
    C_aug=np.zeros((len(sigma)+len(nuisance_names),)*2)
    C_aug[:len(sigma),:len(sigma)]=covariance
    C_aug[len(sigma):,len(sigma):]=np.eye(len(nuisance_names))
    prior=nuisance_projected_information(B_aug,N_aug,C_aug)
    delay_only=nuisance_projected_information(B[-1:],N[-1:],covariance[-1:,-1:])

    smooth_galaxy=GalaxyPotential(galaxy.b,core_Mpc=5e-5)
    smooth=solve_images(geometry,smooth_galaxy,source)
    smooth_B=np.column_stack([response_vector(complete_map_response(smooth['images_rad'],geometry,
        smooth_galaxy,source,np.array([np.eye(2)*v*hubble_distance for v in fine_maps[:,i]]))) for i in range(6)])
    smooth_N,smooth_names=nuisance_columns(geometry,smooth_galaxy,source,smooth)
    smooth_sigma=np.r_[np.full(smooth['images_rad'].size,.001),np.full(len(smooth['relative_delays_days']),.1)]
    smooth_free=nuisance_projected_information(smooth_B,smooth_N,np.diag(smooth_sigma**2))
    delays=smooth['relative_delays_days']
    ratio=float(delays[1]/delays[0])
    ratio_response=(smooth_B[-1]*delays[0]-delays[1]*smooth_B[-2])/delays[0]**2
    ratio_nuisance=(smooth_N[-1]*delays[0]-delays[1]*smooth_N[-2])/delays[0]**2

    # A declared small input, used only to display a finite paired example.
    scenario_anchor=np.array([1e-3,0.,0.,0.,0.,0.])
    scenario_maps=provider.maps(scenario_anchor)
    scenario_geometry=geometry_from_scalars(scenario_maps.scalar_maps,hubble_distance)
    scenario=solve_images(scenario_geometry,galaxy,source)
    envelope=checkpoint_snapshot(.3,provider.background)
    projection_fields=['signal_rank','nuisance_rank','raw_signal_singular_values',
        'projected_signal_singular_values','relative_rank_tolerance','absolute_rank_threshold',
        'Fisher','null_directions','null_direction_interpretation']
    def compact(info):return {k:info[k] for k in projection_fields}
    state_response=dict(anchor_coordinates=ANCHOR_NAMES,
        d_physical_maps_Mpc_d_anchor=fine_maps*hubble_distance,
        d_observable_d_anchor=B,
        observable_order=['theta_plus_x_arcsec','theta_plus_y_arcsec','theta_minus_x_arcsec',
            'theta_minus_y_arcsec','delay_minus_relative_plus_days'],
        relative_delay_response_days=B[-1],
        relative_delay_fractional_response=B[-1]/solution['relative_delays_days'][0],
        complete_arrival_response_terms=[r['arrival_response_terms_days'] for r in responses],
        arrival_term_order=responses[0]['terms_order'])
    report=dict(schema='strong-lensing-first-matched-comparison-v1',
        status='EXECUTED_ILLUSTRATIVE_BENCHMARK_AND_CONDITIONAL_SIX_STATE_RESPONSE',
        physical_BHSM_amplitude_or_detection=False,
        experiment=dict(z_l=.3,z_s=1.,H0_km_s_Mpc=H0,galaxy_sigma_km_s=220.,
            physical_deflection_strength_rad=galaxy.b,
            source_transverse_Mpc=source,reference_beta_arcsec=beta/ARCSEC,
            reference_Einstein_angle_arcsec=theta_E/ARCSEC,
            observer_wavelength_nm=550.,source_rest_wavelength_nm=275.,
            light='common narrowband transient point source in geometric optics; same emission phase and observer proper-clock protocol',
            incoming_light_distinct_from='cosmological delta_r/v_r and muon/reference covariance',
            fixed='physical source displacement, physical galaxy potential, observed lens/source redshifts, observer event and conditional normalized profile axis',
            varied='six harmonic anchor amplitudes; inferred angular source position and image positions re-solved',
            envelope_role='conditional smooth metric correction plus a separately derived static-phase encounter null; a galaxy-replacement amplitude is not established',
            double_counting='smooth environment enters Jacobi maps once; main galaxy enters encounter potential once; external convergence used only as a nuisance alternative'),
        transport=transport,
        physical_maps_Mpc={k:getattr(geometry,k) for k in ('D_ol','D_os','D_ls')},
        time_delay_matrix_Mpc=geometry.T,
        time_delay_matrix_symmetry_relative=float(np.max(abs(geometry.T-geometry.T.T))/np.max(abs(geometry.T))),
        benchmark=solution,state_response=state_response,
        finite_conditional_example=dict(anchor=scenario_anchor,maps=scenario_maps.as_dict(),
            solution=scenario,observable_difference=observable_vector(scenario)-observable_vector(solution),
            first_order_difference=B@scenario_anchor,
            classification='specified numerical probe, not an action-derived or fitted amplitude'),
        verification=dict(derivative_checks=checks,seven_state_control=seven_state_control,
            analytic_SIS_delay_days=float(2*geometry.T[0,0]*theta_E*beta[0]*3.0856775814913673e19/C_KM_S/86400),
            baseline_reciprocity=baseline.reciprocity,
            lens_to_source_not_observer_difference=float(baseline.D_ls[0,0]-(baseline.D_os[0,0]-baseline.D_ol[0,0]))),
        identifiability=dict(covariance=covariance,
            covariance_status='declared illustrative measurement covariance, no observed-system likelihood',
            nuisance_columns=N,nuisance_names=nuisance_names,
            unrestricted=compact(free),fixed_nuisance_control=compact(raw),
            delay_only_unrestricted=compact(delay_only),
            synthetic_prior_control=dict(prior_standard_deviations=prior_sigma,
                source='explicit synthetic independent nuisance measurements only',information=compact(prior))),
        smooth_check=dict(galaxy_core_kpc=.05,solution=smooth,
            response=smooth_B,nuisance_names=smooth_names,unrestricted=compact(smooth_free),
            delay_ratio=ratio,delay_ratio_response=ratio_response,
            delay_ratio_nuisance_response=ratio_nuisance,
            note='a ratio removes the common mass-sheet scaling but can remain degenerate with source and galaxy structure'),
        envelope_action_projection=envelope,
        dependencies=dict(derived=['current reduced ADM Hessian and complete six-anchor evolution',
            'axial null-geodesic/angular Jacobi equations and measured-redshift event remapping',
            'original lens/arrival equations and complete stationary-image response',
            'static spatial phase mixed-current null under ji=0 or constant texture'],
            specified=['frozen R1 background and H0 conversion','normalized n2 profile and positive axial ray',
                'illustrative galaxy and physical source','light/clock protocol and covariance'],
            conditional=['six unit response columns','small finite zeta probe','synthetic nuisance-prior control'],
            unresolved=['realized six-coordinate state and independently selected spatial profile',
                'same-domain local Newtonian metric/phase Hessian and normalized envelope source',
                'seam state/source normalization and action-derived moving-seam response',
                'observational lens/environment likelihood and independent degeneracy-breaking information']),
        approximation='linear cosmological amplitudes; exact ray and paraxial angular integration in declared first-order-metric completion for derivative checks; thin static galaxy, axial isotropic screen, geometric optics; no nonlinear/UV or arbitrary-sky claim')
    dump(output/'results.json',report)
    np.savez_compressed(output/'response_arrays.npz',B=B,N=N,covariance=covariance,
        d_maps_d_anchor=fine_maps,anchor_to_canonical=provider.anchor_to_canonical,
        smooth_B=smooth_B,smooth_N=smooth_N)
    input_paths=['preregistration/prediction_manifest.json','code/background.py','code/closed_horndeski_n2.py',
        'code/action_native_matter_n2.py','code/bhsm_optical_transfer.py','code/strong_lensing_transport.py',
        'code/strong_lensing_comparison.py','code/strong_lensing_envelope.py','scripts/run_strong_lensing_comparison.py']
    dump(output/'input_manifest.json',dict(base_revision='3427c5297fac71faa501685bd38e5d420619c761',
        historical_recovery_revision='5fd65ea18a6055dbc6028839c2918415e3c00f71',
        historical_recovery_is_reset_target=False,
        input_sha256={p:hashlib.sha256((ROOT/p).read_bytes()).hexdigest() for p in input_paths},
        command='python scripts/run_strong_lensing_comparison.py',
        numerical_steps=[t['step'] for t in transport['finite_trials']],
        current_provider='action-native six-coordinate ADM+Schutz-Sorkin; standard Phi+Psi',
        normalized_profile='sqrt(3)/(pi R^(3/2)) sin(2 chi/R) mu; integral Q² dV_gamma=1; mu=1 conditional axis'))
    write_table(report)
    if make_plots:plot_results(report,output)
    print(json.dumps(serial(dict(images_arcsec=solution['images_arcsec'],delay_days=solution['relative_delays_days'],
        fractional_delay_response=state_response['relative_delay_fractional_response'],
        unrestricted_rank=free['signal_rank'],prior_assisted_rank=prior['signal_rank'],
        maximum_independent_response_relative_error=max(c['analytic_vs_half_step_whitened_relative'] for c in checks))),indent=2))
    return report


def write_table(r):
    g=r['physical_maps_Mpc'];s=r['benchmark'];i=r['identifiability']
    rows=[('Lens/source redshift','$0.3 / 1.0$'),
        ('$D_{ol},D_{os},D_{ls}$ [Mpc]',', '.join(f'{g[k][0,0]:.6f}' for k in ('D_ol','D_os','D_ls'))),
        ('SIS image $x$ positions [arcsec]',', '.join(f'{x:.7f}' for x in s['images_arcsec'][:,0])),
        ('SIS relative delay [days]',f"{s['relative_delays_days'][0]:.8f}"),
        ('Free-nuisance surviving state rank',str(i['unrestricted']['signal_rank'])),
        ('Fixed-nuisance control rank',str(i['fixed_nuisance_control']['signal_rank'])),
        ('Synthetic-prior-assisted rank',str(i['synthetic_prior_control']['information']['signal_rank'])),
        ('Static spatial phase $(u,v)$','$ (0,0)$ under the stated current gate')]
    table=['\\begin{tabular}{p{0.47\\linewidth}p{0.46\\linewidth}}','\\toprule','Quantity & Executed result \\\\','\\midrule']
    table += [label+' & '+value+' \\\\' for label,value in rows]
    table += ['\\bottomrule','\\end{tabular}']
    target=ROOT/'manuscript/tables/strong_lensing_results.tex';target.parent.mkdir(parents=True,exist_ok=True)
    target.write_text('\n'.join(table)+'\n',encoding='utf8',newline='\n')


def plot_results(r,output):
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    plt.rcParams.update({'font.size':10,'axes.spines.top':False,'axes.spines.right':False})
    fig=plt.figure(figsize=(7.3,6.8),layout='constrained')
    grid=fig.add_gridspec(2,2,height_ratios=[1.,1.05])
    axes=[fig.add_subplot(grid[0,0]),fig.add_subplot(grid[0,1]),fig.add_subplot(grid[1,:])]
    s=r['benchmark'];positions=np.array(s['images_arcsec'])
    axes[0].scatter(positions[:,0],positions[:,1],s=55,c=['#19677d','#a45b30'],label='SIS images')
    axes[0].scatter([0],[0],marker='+',s=90,c='black',label='Galaxy centre')
    axes[0].set(xlabel='Screen angle x [arcsec]',ylabel='Screen angle y [arcsec]',title=f"Benchmark: {s['relative_delays_days'][0]:.3f} days",ylim=(-.15,.15))
    axes[0].legend(fontsize=8,loc='upper left')
    f=np.array(r['state_response']['relative_delay_fractional_response'])
    axes[1].barh(np.arange(6),f,color=['#19677d' if a>=0 else '#a45b30' for a in f])
    axes[1].set(yticks=np.arange(6),yticklabels=['zeta','Pi2','delta_m','delta_r','v_m','v_r'],
        xlabel='Fractional delay derivative / anchor unit',title='Conditional six-state response')
    axes[1].axvline(0,color='#777777',lw=.7)
    info=r['identifiability'];raw=np.array(info['fixed_nuisance_control']['raw_signal_singular_values'])
    proj=np.array(info['unrestricted']['projected_signal_singular_values'])
    axes[2].semilogy(np.arange(1,len(raw)+1),np.maximum(raw,1e-17),'o-',color='#19677d',label='Fixed nuisance control')
    axes[2].semilogy(np.arange(1,len(proj)+1),np.maximum(proj,1e-17),'s--',color='#a45b30',label='Free nuisance projection')
    axes[2].axhline(info['unrestricted']['absolute_rank_threshold'],ls=':',c='#666666',label='Declared rank threshold')
    axes[2].set(xlabel='Singular-value index',ylabel='Whitened response singular value',title='No surviving free-nuisance state direction')
    axes[2].legend(fontsize=7)
    figure=ROOT/'manuscript/figures/strong_lensing_comparison.pdf';figure.parent.mkdir(parents=True,exist_ok=True)
    fig.savefig(figure,metadata={'CreationDate':None,'ModDate':None})
    fig.savefig(output/'comparison.png',dpi=180)
    plt.close(fig)


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--output',type=Path,default=ROOT/'artifacts/strong_lensing_first_comparison_20261020')
    p.add_argument('--no-plots',action='store_true')
    args=p.parse_args()
    main(args.output,not args.no_plots)
