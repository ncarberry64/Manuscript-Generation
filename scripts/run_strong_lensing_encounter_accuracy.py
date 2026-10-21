"""Run the new finite-window transfer control without modifying saved predictions."""
from __future__ import annotations
import argparse
from dataclasses import asdict, replace
import hashlib
import json
from pathlib import Path
import sys
import time

ROOT=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT/'code'))
import numpy as np
from strong_lensing_encounter_accuracy import (
    ARCSEC,MPC_OVER_C_SECONDS,PhysicalBoundaryGeometry,RetainedEncounterField,
    thin_window_solution,distributed_window_solution,scalar_response_homotopy,
)


def encode(value):
    if isinstance(value,np.ndarray):return value.tolist()
    if isinstance(value,np.generic):return value.item()
    raise TypeError(type(value).__name__)


def evaluate(*,half_windows=None,convergence=True):
    source=ROOT/'artifacts/strong_lensing_r1_screened_20261020/results.json'
    manifest=ROOT/'preregistration/prediction_manifest.json'
    report=json.loads(source.read_text(encoding='utf-8'))
    frozen=json.loads(manifest.read_text(encoding='utf-8'))
    geometry=PhysicalBoundaryGeometry.from_saved(report,frozen['background_R1']['Omega_k0'])
    field=RetainedEncounterField.from_saved(report)
    seeds=np.asarray(report['conventional']['images_rad'])[:,0]
    windows=half_windows or [field.rV_Mpc,4*field.rV_Mpc,20*field.rV_Mpc,20.]
    records=[]
    for L in windows:
        start=time.perf_counter()
        row={'half_window_Mpc':float(L),'half_window_over_rV':float(L/field.rV_Mpc)}
        for on,name in [(False,'fixed_body'),(True,'R1_response')]:
            thin=thin_window_solution(geometry,field,L,on,seeds)
            thick=distributed_window_solution(geometry,field,L,on,thin['images_rad'])
            row[name]={'thin_window':thin,'distributed_window':thick,
                       'thickness_delay_correction_seconds':thick['relative_delays_seconds']-thin['relative_delays_seconds'],
                       'thickness_image_correction_microarcsec':(thick['images_rad']-thin['images_rad'])/ARCSEC*1e6}
        row['thin_reference_relative_delay_offsets_seconds']=row['R1_response']['thin_window']['relative_delays_seconds']-row['fixed_body']['thin_window']['relative_delays_seconds']
        row['distributed_reference_relative_delay_offsets_seconds']=row['R1_response']['distributed_window']['relative_delays_seconds']-row['fixed_body']['distributed_window']['relative_delays_seconds']
        row['thickness_change_of_delay_offsets_seconds']=row['distributed_reference_relative_delay_offsets_seconds']-row['thin_reference_relative_delay_offsets_seconds']
        row['thin_reference_relative_image_offsets_microarcsec']=(row['R1_response']['thin_window']['images_rad']-row['fixed_body']['thin_window']['images_rad'])/ARCSEC*1e6
        row['distributed_reference_relative_image_offsets_microarcsec']=(row['R1_response']['distributed_window']['images_rad']-row['fixed_body']['distributed_window']['images_rad'])/ARCSEC*1e6
        row['thickness_change_of_image_offsets_microarcsec']=row['distributed_reference_relative_image_offsets_microarcsec']-row['thin_reference_relative_image_offsets_microarcsec']
        impacts=geometry.D_ol_Mpc*row['fixed_body']['thin_window']['images_rad']
        row['extra_field_fixed_impact_Born_tail_bounds']=[dict(field.tail_bounds(b,impacts[0],L),
            delay_difference_seconds=field.tail_bounds(b,impacts[0],L)['potential_difference_Mpc']/geometry.a_l*MPC_OVER_C_SECONDS) for b in impacts[1:]]
        row['linear_extra_thickness_change_seconds']=row['fixed_body']['distributed_window']['linear_extra_response_seconds']-row['fixed_body']['thin_window']['linear_extra_response_seconds']
        row['finite_minus_linear_thickness_change_seconds']=row['thickness_change_of_delay_offsets_seconds']-row['linear_extra_thickness_change_seconds']
        row['elapsed_seconds']=time.perf_counter()-start
        records.append(row)
        print('window',L,'offset thickness corrections seconds',row['thickness_change_of_delay_offsets_seconds'],flush=True)
    final=records[-1]
    convergence_record=None
    if convergence:
        high=replace(field,radial_quadrature_order=64)
        L=final['half_window_Mpc'];solutions={}
        for on,name in [(False,'fixed_body'),(True,'R1_response')]:
            thin=thin_window_solution(geometry,high,L,on,seeds)
            thick=distributed_window_solution(geometry,high,L,on,thin['images_rad'],tol=8e-12)
            solutions[name]={'thin_window':thin,'distributed_window':thick}
        difference=solutions['R1_response']['distributed_window']['relative_delays_seconds']-solutions['fixed_body']['distributed_window']['relative_delays_seconds']
        thin_difference=solutions['R1_response']['thin_window']['relative_delays_seconds']-solutions['fixed_body']['thin_window']['relative_delays_seconds']
        convergence_record={'variation':'collocation tolerance 2e-11 to8e-12; radial potential quadrature32 to64',
            'higher_accuracy_solutions':solutions,
            'distributed_delay_offsets_seconds':difference,
            'thin_delay_offsets_seconds':thin_difference,
            'change_of_distributed_delay_offsets_seconds':difference-final['distributed_reference_relative_delay_offsets_seconds'],
            'change_of_thickness_delay_correction_seconds':(difference-thin_difference)-final['thickness_change_of_delay_offsets_seconds'],
            'linear_extra_thickness_change_seconds':solutions['fixed_body']['distributed_window']['linear_extra_response_seconds']-solutions['fixed_body']['thin_window']['linear_extra_response_seconds'],
            'change_of_linear_extra_thickness_change_seconds':(solutions['fixed_body']['distributed_window']['linear_extra_response_seconds']-solutions['fixed_body']['thin_window']['linear_extra_response_seconds'])-final['linear_extra_thickness_change_seconds']}
    homotopy=scalar_response_homotopy(geometry,replace(field,radial_quadrature_order=64),final['half_window_Mpc'],seeds) if convergence else None
    return {
      'status':'COMPUTED_CONDITIONAL_THICKNESS_AND_PATH_CONTROL_OF_RETAINED_OPTICAL_ACTION',
      'experiment':'Same fixed physical source, observer clock, measured-redshift distance maps, finite density, frozen R1 coefficients and three image branches. Both calculations use the same field window and separately derived FLRW exterior wavefronts.',
      'action':'c t=(1/a_l)[S_background(left)+integral(0.5(x_prime^2-kappa x^2)-2Phi)dz+S_background(right)]',
      'physical_boundary_geometry':asdict(geometry),
      'wavefront_identity_checks':geometry.boundary_receipt(),
      'field':asdict(field),
      'windows':records,
      'numerical_convergence':convergence_record,
      'finite_response_stationary_homotopy':homotopy,
      'saved_infinite_thin_plane_delay_offsets_seconds':report['changes']['relative_delays_seconds'],
      'unchanged_original_evidence_sha256':{str(source.relative_to(ROOT)):hashlib.sha256(source.read_bytes()).hexdigest(),str(manifest.relative_to(ROOT)):hashlib.sha256(manifest.read_bytes()).hexdigest()},
      'qualification':[
        'This distributes the retained leading weak-field paraxial action; it is not an exact null-ray solve of a completed higher-order scalar metric.',
        'The FLRW exterior wavefront split is fixed by individual Jacobi boundary data, not chosen only to reproduce their sum.',
        'Tail bounds concern the declared isolated extra field at fixed impact within Born projection, not an arbitrary environmental or global matching error.',
        'Scalar-dependent field-equation higher orders, pressure/EOS evolution, time-dependent background and global matching remain outside this numerical control.',
        'No state parameter, amplitude, matter profile, boundary normalization or physical guard is retuned; this is not physical admission.'
      ],
      'physical_admission':False
    }


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output',type=Path,default=ROOT/'artifacts/strong_lensing_acceptance_20261020/encounter_accuracy.json')
    parser.add_argument('--windows',type=float,nargs='+',help='Optional half-window sizes in Mpc.')
    parser.add_argument('--skip-convergence',action='store_true')
    args=parser.parse_args()
    result=evaluate(half_windows=args.windows,convergence=not args.skip_convergence)
    result['source_sha256']={p:hashlib.sha256((ROOT/p).read_bytes()).hexdigest() for p in ['code/strong_lensing_encounter_accuracy.py','scripts/run_strong_lensing_encounter_accuracy.py']}
    args.output.parent.mkdir(parents=True,exist_ok=True)
    with args.output.open('w',encoding='utf-8',newline='\n') as f:json.dump(result,f,indent=2,default=encode);f.write('\n')
    print('output',args.output)

if __name__=='__main__':main()
