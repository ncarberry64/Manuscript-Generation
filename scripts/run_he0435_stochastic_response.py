"""Reproducible synthetic/source application; never consumes native target data.

Native application must supply actual observation/estimator Jacobians and an
admissible unsaturated context.  This executable tests the derived operator and
records the unresolved instrument relations rather than inventing those inputs.
"""
from pathlib import Path
import hashlib
import json
import sys
import numpy as np
from scipy import sparse

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT/"code"))
from strong_lensing_stochastic_response import (
    BoundaryMomentModel, nearest_neighbor_edges, pixel_read_projection,
    project_collection_covariance, independent_adc_read_covariance,
    shared_calibration_covariance, fixed_residual_covariance_application,
)


def sha(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def main():
    output = ROOT/"artifacts/strong_lensing_response_completion_20261021/stochastic"
    output.mkdir(parents=True, exist_ok=True)
    rate = np.array([7., 12., 8., 10., 120., 13., 8., 11., 6.])
    gain = np.array([1.9, 2.1, 2., 2.05, 1.95, 2.2, 2., 1.8, 2.1])
    read_times = np.array([1., 2., 6., 7., 11., 12.])*10.73677
    group_weights = np.zeros((27, 6))
    for group in range(3):
        group_weights[group*9:(group+1)*9, group*2:(group+1)*2] = .5/gain[:, None]
    pixels = np.tile(np.arange(9), 3)
    model = BoundaryMomentModel(rate, 2e-7, nearest_neighbor_edges((3, 3)),
        max_time_s=float(read_times[-1]), charge_limit_e=np.full(9, 80000.),
        domain_source="synthetic declared unsaturated 3x3 no-exterior-edge context")
    result = project_collection_covariance(model, read_times,
        pixel_read_projection(pixels, group_weights, 9),
        provenance="synthetic equal two-read group averages / individual pixel electronic gain")
    zero = BoundaryMomentModel(rate, 0., nearest_neighbor_edges((3, 3)))
    zero_result = project_collection_covariance(zero, read_times,
        pixel_read_projection(pixels, group_weights, 9),
        provenance="identical synthetic map with kappa=0")
    adc_weights = group_weights*gain[pixels, None]
    read = independent_adc_read_covariance(pixel_read_projection(pixels, adc_weights, 9),
        np.full(9, 12.), stage_assumption="synthetic independently declared raw ADC CDS law")
    bias_j = -np.tile(np.eye(9), (3, 1))
    # The diagonal prior below is declared synthetic, not inferred from CRDS ERR.
    reference = shared_calibration_covariance(bias_j, np.eye(9)*.3**2,
                    source="synthetic full covariance specified as diagonal shared bias-DN errors")
    covariance = result["covariance"]+read+reference
    residual = np.linspace(-2., 3., 27)
    original_error = np.sqrt(np.diag(zero_result["covariance"]+read))
    score = fixed_residual_covariance_application(residual, covariance, original_error,
                covariance_scope="synthetic conditional boundary LNA + independent ADC + shared bias")
    count_total = project_collection_covariance(model, [read_times[-1]],
        [sparse.csr_matrix(np.ones((1, 9)))], provenance="total physical electron count at final read")
    np.savez_compressed(output/"synthetic_application.npz", incident_rate_e_s=rate,
        gain_e_per_DN=gain, read_times_s=read_times, group_weights=group_weights,
        mean_DN=result["mean"], original_leading_mean_DN=result["original_leading_mean"],
        poisson_boundary_covariance_DN2=result["covariance"],
        kappa_zero_poisson_covariance_DN2=zero_result["covariance"],
        adc_read_covariance_DN2=read, shared_bias_covariance_DN2=reference,
        complete_conditional_covariance_DN2=covariance, fixed_synthetic_residual_DN=residual)
    report = {
        "status": "executed synthetic/source application; no native data fitted",
        "scope": "conditional detector counting/readout law; no BHSM dressed-state statistical closure",
        "source_mean": "G_j=a_segment+m_j*f_DN+(kappa_e/2)*v_j*diag(1/g)*L((g*f_DN)^2)",
        "same_leading_mean_derivation": {
            "incident_electron_rate": "lambda=g*f_DN, strictly nonnegative physical carrier rate; detector plane offset is not a rate",
            "boundary": "A_pq=kappa_e*(lambda_p+lambda_q) for adjacent p!=q; A_pp=-sum_q A_pq",
            "identity": "A*lambda=kappa_e*L(lambda^2); row/column sums of A are zero",
            "positive_completion": "a_p(N)=Lambda*softmax(log(lambda_p)+(A*N)_p/lambda_p); lambda_p>0; value/Jacobian at N=0 are lambda/A",
            "completion_domain": "positive rates; higher orders are a new explicit conditional assumption, not fixed by the mean",
            "local_moments": "mu'=lambda+A*mu; C'=A*C+C*A^T+diag(lambda+A*mu); mu(0)=C(0)=0",
            "cross_time": "K(t,s)=exp(A*(t-s))*C(s), t>=s",
            "projected": "Cov(sum_t H_t*N(t))=integral_0^T B(u)*diag(exp(A*u)*lambda)*B(u)^T du; B=sum_{t>=u}H_t exp(A*(t-u))",
            "approximation": "resummed local affine/linear-noise moments, not exact nonlinear softmax CME; agrees with leading original mean and leading covariance in kappa",
            "non_identification": "A2=kappa_e*L*diag(lambda) also has A2*lambda=kappa_e*L(lambda^2), but its covariance differs; the source mean does not determine the boundary law",
            "naive_hopping": "propensity kappa*N*(N-1) would first alter the mean at order t^3 and does not reproduce the prescribed t^2 term",
            "gain": "charge conservation is in electron units, not a sum of nonuniform-gain DN",
            "dark": "physical dark births before the map require their own readwise rate/law; native subtracted group-dark and its reference ERR are not that stochastic law"
        },
        "observation_chain": {
            "mean": "H_j=L(mean_r inverseL(q_jr))-dark_j, with q_jr including physical dark and intercept before inverseL",
            "photon_map": "H_t contains native dH/dq_read divided by gain, then actual estimator/flat/photometric derivatives",
            "ADC_noise_map": "read noise enters raw ADC after inverseL: native L'(mean ADC)/N, not dH/dq_read",
            "CDS": "CONDITIONAL independent single ADC read variance=CDS_DN^2/2; actual0282 CDS marginals were measured after reference/IPC correction and do not identify raw ADC independence",
            "shared_references": "sum_named J_ref*C_ref*J_ref^T, preserving same physical reference columns across groups/exposures",
            "parameter_score": "full Gaussian score includes -.5*r^T C^-1 C_theta C^-1*r + .5*tr(C^-1 C_theta), in addition to mean derivative",
            "no_estimator_substitution": "native rate covariance requires its producing adaptive derivative, including actual DQ/segmentation/ZEROFRAME branch"
        },
        "synthetic_parameters": {"shape": [3, 3], "kappa_e_inverse": 2e-7,
            "electronic_gain_range_e_per_DN": [float(gain.min()), float(gain.max())],
            "maximum_time_s": float(read_times[-1]), "charge_ceiling_e": 80000.,
            "readnoise_CDS_DN": 12., "shared_bias_sigma_DN": .3},
        "synthetic_result": {"number_of_group_pixels": 27,
            "poisson_covariance_min_eigenvalue_DN2": float(np.linalg.eigvalsh(result["covariance"]).min()),
            "largest_interpixel_boundary_covariance_DN2": float(result["covariance"][22, 23]),
            "same_pixel_cross_group_covariance_DN2": float(result["covariance"][4, 22]),
            "max_resummed_minus_original_mean_DN": float(np.max(abs(result["resummed_minus_original_mean"]))),
            "total_electron_mean": float(count_total["mean"][0]),
            "total_electron_variance": float(count_total["covariance"][0, 0]),
            "expected_conserved_Poisson_total": float(rate.sum()*read_times[-1]),
            "projection_diagnostics": result["diagnostics"], "fixed_residual_application": score},
        "actual_available_input_relations": [
            "native CAL metadata: SHALLOW2, 9 groups, 2 frames, gap3, TFRAME10.73677s; actual five-patch histories recovered separately",
            "gain0094, CDS-readnoise0282, linearity0054, saturation0100, bias0180, flat0667, dark0337, SIRSkernel0002 native identities",
            "per-read native polynomial derivative and actual adaptive slope-map Jacobian are recovered by separate source-trail applications"
        ],
        "individually_unresolved_relations": [
            "actual NRCB3 wavelength/illumination/position-dependent collection-boundary Jacobian; kappa and isotropy not uniquely fixed by one mean fit",
            "IPC/NLIPC charge-to-voltage spatial kernel and its signal dependence, separated from already IPC-broadened effective optical PSF",
            "post-saturation stored-charge/overflow/collection response of influential cores; no constant-lambda extrapolation physically admitted",
            "initial/reset charge state and persistence/trapping history law",
            "full raw-ADC temporal/interpixel read covariance and SIRS reference/control covariance/Jacobian",
            "actual0282 HISTORY says CDS was measured after reference/bias/IPC corrections: its marginals cannot uniquely identify the raw covariance needed before another SIRS projection",
            "full shared spatial/temporal calibration covariance: CRDS ERR provides marginals only",
            "physical dark-current photon/carrier counting law before native readwise/group dark mapping",
            "signal-dependent masking and adaptive rate-estimator selection law beyond a local fixed-branch Jacobian",
            "BHSM prepared U,V contractions and dressed optical statistical kernel remain separate from this detector counting approximation"
        ],
        "domain": "time/charge-ceiling checks include unobserved context cores; mean below ceiling is necessary, not a random-overflow law",
        "physical_guard": 2e-10, "physical_validation": False,
        "no_reserved_observations_consumed": True, "no_native_fit_by_this_script": True,
        "sources": [
            {"url": "https://arxiv.org/abs/1906.01846", "use": "primary infrared moving-boundary counting, unequal-time moments, IPC/BFE distinction; NOT an NRCB3 calibration"},
            {"url": "https://jwst-pipeline.readthedocs.io/en/2.0.1/jwst/references_general/readnoise_reffile.html", "use": "native READNOISE is CDS-DN noise"},
            {"url": "https://raw.githubusercontent.com/spacetelescope/stcal/1.18.0/src/stcal/ramp_fitting/ramp_fit.py", "use": "producing ramp-fit noise scaling and adaptive-estimator source; actual recovered mapping must be used"},
            {"url": "https://jwst-docs.stsci.edu/jwst-near-infrared-camera/nircam-instrumentation/nircam-detector-overview/nircam-detector-readout", "use": "nondestructive NIRCam readout with accumulated charge in place"}
        ],
        "file_sha256": {str(p.relative_to(ROOT)).replace('\\', '/'): sha(p) for p in (
            ROOT/"code/strong_lensing_stochastic_response.py", Path(__file__),
            ROOT/"tests/test_strong_lensing_stochastic_response.py", output/"synthetic_application.npz")},
    }
    (output/"derivation_and_application.json").write_text(
        json.dumps(report, indent=2, allow_nan=False)+"\n", encoding="utf-8", newline="\n")
    print(json.dumps(report["synthetic_result"], indent=2))


if __name__ == "__main__":
    main()
