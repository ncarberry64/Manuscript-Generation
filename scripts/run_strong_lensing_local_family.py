"""Reproduce the conditional finite local-envelope and matched optical encounter."""
from __future__ import annotations
import os
os.environ.setdefault('OMP_NUM_THREADS','1')
os.environ.setdefault('OPENBLAS_NUM_THREADS','1')
from pathlib import Path
import argparse
import hashlib
import json
import sys
import numpy as np

ROOT=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT/'code'))
from strong_lensing_local_family import (
    TolmanVII,TolmanVIILens,saved_geometry,phase_family_receipt,compute_local_comparison,
)
from strong_lensing_local_hessian import evaluate_local_hessian


def serial(value):
    if isinstance(value,np.ndarray):return value.tolist()
    if isinstance(value,np.generic):return value.item()
    if isinstance(value,dict):return {k:serial(v) for k,v in value.items()}
    if isinstance(value,(tuple,list)):return [serial(v) for v in value]
    return value


def dump(path,data):
    path.parent.mkdir(parents=True,exist_ok=True)
    path.write_text(json.dumps(serial(data),indent=2,sort_keys=True,allow_nan=False)+'\n',
        encoding='utf8',newline='\n')


def write_table(report):
    checks=report['fluid_checks'];encounter=report['optical_encounter'];solution=encounter['solution']
    responses=encounter['intrinsic_responses'];ident=encounter['identifiability']
    rows=[
        ('Finite mass / areal radius',r'$10^{11}M_\odot$ / $3$ kpc'),
        ('$GM/(Rc^2)$',f"{checks['compactness_GM_over_Rc2']:.8g}"),
        ('Independent TOV residual (scaled)',f"{checks['max_independent_TOV_relative_scale_residual']:.3g}"),
        ('Maximum $c_s^2/c^2$',f"{checks['max_sound_speed_squared_over_c_squared']:.8g}"),
        ('Surface pressure / lapse-log residual',f"{checks['surface_pressure_residual']:.1g} / {checks['surface_lapse_log_residual']:.1g}"),
        ('Image $x$ positions [arcsec]',', '.join(f'{v:.7f}' for v in solution['images_arcsec'][:,0])),
        ('Relative delays [days]',', '.join(f'{v:.7f}' for v in solution['relative_delays_days'])),
        (r'$\partial_{\ln M}\Delta t$ [days]',', '.join(f'{v:.7f}' for v in responses['log_mass']['relative_delays_days'])),
        (r'$\partial_{\ln R}\Delta t$ [days]',', '.join(f'{v:.7f}' for v in responses['log_radius']['relative_delays_days'])),
        ('Maximum image-equation residual [rad]',f"{solution['max_lens_residual_rad']:.3g}"),
        ('Cosmological fixed/free-nuisance rank',f"{ident['cosmological_fixed_nuisance']['signal_rank']} / {ident['cosmological_unrestricted']['signal_rank']}"),
        (r'Local $q_0$ response / free-nuisance rank',r'$0$ / $0$'),
    ]
    lines=[r'\begin{tabular}{p{0.49\linewidth}p{0.44\linewidth}}',r'\toprule',
        r'Quantity & Evaluated conditional finite reference \\',r'\midrule']
    lines += [label+' & '+value+r' \\' for label,value in rows]
    lines += [r'\bottomrule',r'\end{tabular}']
    target=ROOT/'manuscript/tables/strong_lensing_local_family_results.tex'
    target.write_text('\n'.join(lines)+'\n',encoding='utf8',newline='\n')


def main(output):
    baseline_path=ROOT/'artifacts/strong_lensing_first_comparison_20261020/results.json'
    baseline,geometry,source=saved_geometry(baseline_path)
    fluid=TolmanVII();lens=TolmanVIILens()
    checks=fluid.checks()
    optical=compute_local_comparison(baseline,geometry,source,lens)
    hessian=evaluate_local_hessian(fluid.compactness)
    grid=np.r_[np.linspace(0,fluid.radius_Mpc,401),np.linspace(1.005,2,200)*fluid.radius_Mpc]
    profile=fluid.application(grid)
    report=dict(schema='strong-lensing-local-regular-reference-v1',
        status='EXECUTED_CONDITIONAL_PATH_B_CONTROL_AND_NONZERO_MATCHED_FINITE_GRAVITATIONAL_ENCOUNTER',
        physical_state_uniquely_selected=False,BHSM_extra_deflection_claim=False,
        physical_inputs=dict(mass_Msun=fluid.mass_Msun,areal_radius_Mpc=fluid.radius_Mpc,
            mass_length_Mpc=fluid.mass_length_Mpc,source_transverse_Mpc=source,
            physical_maps_Mpc={k:getattr(geometry,k) for k in ('D_ol','D_os','D_ls')},
            z_l=geometry.z_l,z_s=baseline['experiment']['z_s'],
            incident_light=baseline['experiment']['light'],
            observer_wavelength_nm=baseline['experiment']['observer_wavelength_nm'],
            matter_profile='rho=15M/(8piR³)(1-r²/R²) for r<=R, zero outside; total finite physical mass M',
            units='Einstein solve uses lengths in Mpc and geometrized rho=G rho_mass/c², p=G p_physical/c⁴; lens psi is Mpc, alpha radians; fluid EOS returned by TOV',
            local_action_assumptions=dict(
                C_partial='c^4/(16 pi G), specified effective four-dimensional Einstein-Hilbert coefficient per physical coordinate time',
                additional_cap_traction='excluded from this conditional local action block',
                extra_local_vacuum_energy='excluded from this conditional local action block',
                eta_weight_and_phase_scale='w>0 and kappa1>0; phase Hessian is reported per k_eta=w*kappa1, not calibrated to a lens',
                sigma_and_cap='dynamical sigma/cap reaction and global closed-cosmology matching are not solved; zero eta adds no eta-dependent force, while the changed metric can still alter other-sector equations'),
            source_owner='same fixed physical source displacement and observer proper-clock protocol as saved first checkpoint'),
        original_equations=dict(
            mass='mprime=4pi r² rho_g',
            radial_Einstein='nuprime=(m+4pi r³ p_g)/(r²(1-2m/r))',
            TOV='p_gprime=-(rho_g+p_g)nuprime',
            angular_Einstein='f(nusecond+nuprime²+nuprime/r)+(fprime/2)(nuprime+1/r)=8pi p_g',
            momentum='Kij=0 and matter momentum=0',
            eta='tangent projection of div[w(kappa1+X³)Deta]=0; eta.eta=1',
            color='Dnu Fnumu/g3²=Jeta^mu; Deta=F=0 gives both sides zero',
            matching='p(R)=rho(R)=0; nu(R)=log(1-2m/R)/2; f(R)=1-2m/R; same extrinsic curvature on opposite sides'),
        phase_reference=phase_family_receipt(),fluid_checks=checks,
        exterior_physical_hessian=hessian,optical_encounter=optical,
        approximations=dict(
            exact='static spherical perfect-fluid Einstein/TOV, unit-target eta/color null rows and vacuum annulus Hessian',
            optical='leading weak-field projected finite density; order compactness relative corrections from pressure and nonlinear metric are not included in quoted optical values',
            cosmological_local_scale='local star and exterior annulus use a conditional effective Einstein-Hilbert sector. Pure background-curvature terms scale as O((H_l R/c)^2); additional rolling-scalar/cap mixed responses and global closed-FLRW matching are not bounded or solved here.',
            orientation='axial incoming ray with isotropic screen; spherical local reference has no preferred interface orientation',
            kinetic_domain='positive phase kinetic coefficient w kappa1 and usual reduced Einstein tensor kinetic form; indefinite static scalar metric Hessian is a constraint form, not a ghost test'),
        dependencies=dict(derived=['finite regular constant-q family from original eta and gauge Euler equations',
            'nonzero finite-mass Einstein/TOV gravitational field and zero-layer matching',
            'same-background physical annulus action Hessian and phase null',
            'finite-profile optical potential, gradient and fixed-impact mass/radius derivatives',
            'new images, delays and nuisance analysis on this matching encounter'],
            specified=['conditional total mass, radius, EOS profile, q0 Dirichlet state',
                'same saved propagation geometry and incoming light/clock preparation'],
            conditional=['M4 Path-B low-energy sector and ordinary perfect-fluid matter',
                'C_partial=c^4/(16 pi G) is the specified effective four-dimensional Einstein-Hilbert truncation; no extra local cap reaction or vacuum-energy term is included',
                'dynamical sigma/cap and global closed-cosmology matching are not solved; zero eta adds no eta-dependent force, but a changed metric can affect the other-sector equations',
                'absence of an added phase interface charge in this source-free problem'],
            unresolved=['selection of physical q0 and nonconstant sourced states',
                'common-action normalization and derivative mapping from exported scalar T/psi to eta phase',
                'nonzero seam source/orientation and moving-seam dynamics',
                'global finite-lens/closed-cosmology metric matching beyond the declared local expansion']))
    dump(output/'results.json',report)
    np.savez_compressed(output/'profiles.npz',**profile)
    input_paths=['code/strong_lensing_local_family.py','code/strong_lensing_local_hessian.py',
        'code/strong_lensing_comparison.py','scripts/run_strong_lensing_local_family.py',
        'artifacts/strong_lensing_local_family_source_receipt_20261020.json',
        'artifacts/strong_lensing_first_comparison_20261020/results.json',
        'artifacts/strong_lensing_first_comparison_20261020/response_arrays.npz']
    dump(output/'input_manifest.json',dict(reviewed_checkpoint='97bf5571335cf60e15e2e486e1197ca226198ed2',
        checkpoint_is_reset_target=False,
        command='python scripts/run_strong_lensing_local_family.py',
        baseline_reused_without_rerun=True,
        input_sha256={p:hashlib.sha256((ROOT/p).read_bytes()).hexdigest() for p in input_paths},
        action_provenance=['BHSM_PATH_B_FOUNDATIONAL_ACTION_V14_31.md',
            'BHSM_PATH_B_PHYSICAL_TOPOLOGY_AND_FR_GATE_V14_32.md',
            'BHSM_STATIC_ETA_METRIC_SPIN4_SOURCE_V14_39.md'],
        action_location='C:/Users/carbe/Downloads/BHSM_muon_parent_maxwell_worktree_20261001/docs',
        source_versions='source hashes and revision pinned in the separate local action recovery receipt'))
    write_table(report)
    print(json.dumps(serial(dict(compactness=fluid.compactness,
        images_arcsec=optical['solution']['images_arcsec'],
        relative_delays_days=optical['solution']['relative_delays_days'],
        intrinsic_delay_derivatives={k:v['relative_delays_days'] for k,v in optical['intrinsic_responses'].items()},
        q_response_rank=optical['identifiability']['q_only']['signal_rank'],
        cosmological_raw_rank=optical['identifiability']['cosmological_fixed_nuisance']['signal_rank'],
        cosmological_free_rank=optical['identifiability']['cosmological_unrestricted']['signal_rank'],
        own_factorization_error=optical['finite_profile_factorization']['relative_whitened_discrepancy'])),indent=2))
    return report


if __name__=='__main__':
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output',type=Path,default=ROOT/'artifacts/strong_lensing_local_family_20261020')
    main(parser.parse_args().output)
