"""Render completed frozen-test receipts; never opens raw/sealed photometry.

Run only after the parent researcher announces that the scored receipt is
complete.  This script does no fitting, state filtering or model selection.
It checks receipt ordering against reported metrics and writes the separate
results section, a three-panel prediction figure and a compact metric receipt.
"""
from __future__ import annotations

from hashlib import sha256
import json
from pathlib import Path

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np


ROOT=Path(__file__).resolve().parents[1]
OUT=ROOT/'artifacts/strong_lensing_revision_20261020'
SCORES=OUT/'revision_unused_scores.json'
CONVENTIONAL=OUT/'revised_conventional_training.json'
R1=OUT/'revised_r1_training.json'
SELECTION=OUT/'macro_selection.json'
SECTION=ROOT/'manuscript/sections/he0435_revision_results_20261020.tex'
FIGURE=ROOT/'manuscript/figures/he0435_revision_unused_prediction'
MODEL_ORDER=('empirical_development_delays','revised_conventional_macro',
             'revised_conditional_R1_macro','fixed_body_R1_activation','zero_delay_control')
LABELS={'empirical_development_delays':'Empirical delay ensemble',
        'revised_conventional_macro':'Conventional macro',
        'revised_conditional_R1_macro':'Conditional R1 macro',
        'fixed_body_R1_activation':'Fixed-body R1 activation',
        'zero_delay_control':'Zero-delay control'}


def digest(path):
    return sha256(Path(path).read_bytes()).hexdigest()


def save(path,value):
    path.parent.mkdir(parents=True,exist_ok=True)
    path.write_text(json.dumps(value,indent=2,allow_nan=False)+'\n',
                    encoding='utf-8',newline='\n')


def scientific(value,digits=4):
    mantissa,exponent=f'{float(value):.{digits-1}e}'.split('e')
    return rf'\({mantissa}\times10^{{{int(exponent)}}}\)'


def primary_design(scores):
    design=scores['target_design']
    arrays={k:np.asarray(design[k]) for k in ('time_days','image','magnitude','error_mag')}
    if len({a.size for a in arrays.values()})!=1:
        raise ValueError('Scored design arrays have inconsistent sizes')
    mask=arrays['image']!=2
    result={k:v[mask] for k,v in arrays.items()}
    if len(result['time_days'])!=scores['primary_target_count']:
        raise ValueError('Primary C-driver target count disagrees with scored design')
    for model in MODEL_ORDER:
        receipt=scores['primary'][model]
        means=np.asarray(receipt['mean']); q=receipt['quantiles']
        if means.shape!=result['magnitude'].shape or any(np.asarray(x).shape!=means.shape for x in q.values()):
            raise ValueError('Scored prediction shape/order is inconsistent')
        rmse=float(np.sqrt(np.mean((means-result['magnitude'])**2)))
        if not np.isclose(rmse,receipt['overall']['rmse_mag'],rtol=1e-10,atol=1e-12):
            raise ValueError('Prediction order does not reproduce the stored RMSE')
        inside=(result['magnitude']>=np.asarray(q['0.025']))&(result['magnitude']<=np.asarray(q['0.975']))
        if not np.isclose(np.mean(inside),receipt['overall']['coverage_95'],atol=1e-12):
            raise ValueError('Prediction quantiles do not reproduce stored coverage')
    return result


def plot_predictions(scores,design):
    plt.rcParams.update({'font.size':9,'axes.labelsize':9,'legend.fontsize':8,
                         'pdf.fonttype':42,'ps.fonttype':42})
    fig,axes=plt.subplots(3,1,figsize=(7.1,7.4),sharex=True,constrained_layout=True)
    for ax,image,label in zip(axes,(0,1,3),('A','B','D')):
        selected=np.flatnonzero(design['image']==image)
        selected=selected[np.argsort(design['time_days'][selected],kind='stable')]
        time=design['time_days'][selected]-58000.
        ax.errorbar(time,design['magnitude'][selected],yerr=design['error_mag'][selected],
                    fmt='o',ms=3.,color='black',ecolor='.45',elinewidth=.6,
                    alpha=.8,label='Unused LCO R observations',zorder=4)
        for model,color,style,name in (
            ('revised_conventional_macro','#2468a0','-','Conventional'),
            ('revised_conditional_R1_macro','#ce6c18','--','Conditional R1')):
            prediction=scores['primary'][model]
            mean=np.asarray(prediction['mean'])[selected]
            lower=np.asarray(prediction['quantiles']['0.025'])[selected]
            upper=np.asarray(prediction['quantiles']['0.975'])[selected]
            ax.fill_between(time,lower,upper,color=color,alpha=.13,
                            label=f'{name}: marginal 95% interval')
            ax.plot(time,mean,style,color=color,lw=1.35,label=f'{name}: mean',zorder=3)
        ax.invert_yaxis();ax.set_ylabel(f'Image {label}: R magnitude')
        ax.grid(alpha=.2,lw=.5)
    axes[-1].set_xlabel('Native LCO MJD $-$ 58000 (days)')
    handles,labels=axes[0].get_legend_handles_labels()
    axes[0].legend(handles,labels,ncol=2,loc='best',framealpha=.95)
    FIGURE.parent.mkdir(parents=True,exist_ok=True)
    fig.savefig(FIGURE.with_suffix('.pdf'),bbox_inches='tight')
    fig.savefig(FIGURE.with_suffix('.png'),dpi=180,bbox_inches='tight')
    plt.close(fig)


def build_section(scores,conventional,r1,selection):
    primary=scores['primary'];secondary=scores['secondary']
    conv=primary['revised_conventional_macro'];bhsm=primary['revised_conditional_R1_macro']
    convs=conv['overall'];bs=bhsm['overall']
    sec=secondary['revised_conventional_macro'];secs=sec['overall']
    secb=secondary['revised_conditional_R1_macro']
    delta=bhsm['joint_conditional_log_density']-conv['joint_conditional_log_density']
    fixed_delta=primary['fixed_body_R1_activation']['joint_conditional_log_density']-conv['joint_conditional_log_density']
    secdelta=secb['joint_conditional_log_density']-sec['joint_conditional_log_density']
    mass=float(conventional['compact_mass_Msun'])
    if not np.isclose(mass,float(r1['compact_mass_Msun'])):
        raise ValueError('Reported physical models must share the selected compact source')
    if not np.isclose(mass,float(selection['selected_compact_mass_Msun'])):
        raise ValueError('Macro selection and reported paired models disagree')
    sections=r1['original_section_source_applications']
    row=min(sections,key=lambda item:abs(item['radius_kpc']-r1['compact_scale_kpc']))
    scalar_res=max(x['scalar_relative_residual'] for x in sections)
    gauss_res=max(x['metric_relative_residual'] for x in sections)
    activation=r1['fixed_conventional_body_R1_activation']
    activation_delays=activation['delay_offsets_seconds']
    activation_images=np.asarray(activation['image_offsets_microarcsec'])
    nconv=conventional['global_image_census']['image_count']
    nr1=r1['global_image_census']['image_count']
    covered=int(round(convs['coverage_95']*convs['count']))
    lines=[
        r'\subsection{Revised source and frozen predictions on unused observations}',
        r'\label{sec:he0435-revised-results}',
        rf'The development selection retains candidate {selection["selected_candidate_index"]} with a regular compact matter component of {scientific(mass)}\,$M_\odot$ and scale {r1["compact_scale_kpc"]:.3g}\,kpc.',
        rf'The conventional and conditional R1 applications have {nconv} and {nr1} verified stationary images, respectively; all additional images enter the summed central-light constraint.',
        r'The nearly degenerate development objectives do not measure a nuclear mass. The selected positive component is a conditional matter surrogate, not an inferred black hole or a solved stellar equation of state.',
        r'',
        r'The original scalar and metric equations are applied to the total enclosed matter once:',
        r'\[ Cx^2p-\frac{4}{\widehat\lambda}xp^2=-\widehat b\,m(x),\qquad 2x^2\frac{d\Phi}{dx}+\widehat b x^2p=2m(x). \]',
        r'Here $x=rH_0$, $m(x)=GM(<r)H_0$, $\widehat b=b_T/M_{\rm Pl}$, $p=\pi_T^{\prime}/(M_{\rm Pl}H_0)$ and $c=1$, with the previously fixed $\widehat\lambda=\lambda/(M_{\rm Pl}H_0^2)$.',
        rf'At $r={row["radius_kpc"]:.4g}$\,kpc, the computed section has $\widehat\rho=\rho/(M_{{\rm Pl}}^2H_0^2)$ equal to {scientific(row["density_over_Mpl2H02"])}, reduced source $\widehat J_{{\rm eff}}=J_{{\rm eff}}/(M_{{\rm Pl}}H_0^2)=-\widehat b\widehat\rho/2$ equal to {scientific(row["consumed_reduced_source_J_eff"])}, and $p$ equal to {scientific(row["scalar_gradient_over_MplH0"])}.',
        rf'The scalar integrated left side and source are {scientific(row["scalar_integrated_lhs"])} and {scientific(row["scalar_integrated_source"])}, while the metric Gauss left side and source are {scientific(row["metric_Gauss_lhs"])} and {scientific(row["metric_Gauss_source"])}.',
        rf'Across the evaluated sections the maximum relative scalar and Gauss residuals are {scientific(scalar_res)} and {scientific(gauss_res)}.',
        rf'Activating R1 on the fixed conventional body, source, maps and clock gives B/C/D-minus-A delay offsets $({activation_delays[0]:.5f},{activation_delays[1]:.5f},{activation_delays[2]:.5f})$\,s and a maximum coordinate image displacement {scientific(np.max(abs(activation_images)))}\,$\mu$as. Separately refitted model predictions below need not retain those fixed-reference offsets.',
        r'',
        rf'The final test opens {scores["primary_target_count"]} previously unused A/B/D measurements, conditional on {scores["driver_count"]} allowed C measurements from the same later LCO epochs. The no-driver secondary forecast opens all {scores["secondary_target_count"]} measurements. No brightness, delay, offset, calibration or source parameter is refitted to these target values.',
        r'Table~\ref{tab:he0435-revised-primary} reports the normalized joint conditional density and marginal predictive diagnostics. The joint log density is a density in the declared magnitude units; epochs are correlated, and its value is not a chi-square significance or an independent-event count.',
        r'\begin{table}[t]',
        r'\caption{Frozen revised models on unused LCO R A/B/D observations, with image C as the allowed driver. Coverage is the empirical fraction within marginal 95\% mixture intervals; width and the proper 95\% interval score are in magnitudes. Larger joint log density and smaller interval score are preferred.}\label{tab:he0435-revised-primary}',
        r'\centering\small',
        r'\begin{tabular}{lrrrrr}\toprule',
        r'Model & Joint $\log p$ & RMSE & Coverage (\%) & Width & Interval score\\\midrule',
    ]
    for model in MODEL_ORDER:
        result=primary[model];s=result['overall']
        lines.append(f'{LABELS[model]} & {result["joint_conditional_log_density"]:.3f} & {s["rmse_mag"]:.4f} & {100*s["coverage_95"]:.2f} & {s["mean_95_interval_width_mag"]:.4f} & {s["mean_95_interval_score_mag"]:.4f}'+r'\\')
    lines += [r'\bottomrule\end{tabular}',r'\end{table}',r'',
        rf'The conventional intervals contain {covered} of {convs["count"]} primary targets ({100*convs["coverage_95"]:.2f}\%). The conditional R1-minus-conventional joint log-score difference is {scientific(delta)}; this is a finite-model predictive comparison, not a Bayes factor from a completed physical-state posterior.',
        rf'The fixed-body R1 activation has a primary log-score difference of {scientific(fixed_delta)} relative to conventional, with the same optical nuisance parameters. It isolates the original scalar response; the separately refitted R1 comparison also includes development nuisance and optimizer movement and cannot be attributed entirely to that activation.',
        r'The matched conventional effective-mass representation has the same retained optical transfer as conditional R1 and therefore the same score by construction; it is not an independent confirmation.',
        r'\begin{table}[t]',
        r'\caption{Image-resolved primary diagnostics. Each image contributes 54 unused observations. These are correlated marginal diagnostics, not 54 independent coverage trials.}\label{tab:he0435-revised-images}',
        r'\centering\small\begin{tabular}{lrrr|rrr}\toprule',
        r'&\multicolumn{3}{c}{Conventional}&\multicolumn{3}{c}{Conditional R1}\\',
        r'Image & RMSE & Coverage (\%) & Width & RMSE & Coverage (\%) & Width\\\midrule',
    ]
    for image,label in ((0,'A'),(1,'B'),(3,'D')):
        a=conv['by_image'][str(image)];b=bhsm['by_image'][str(image)]
        lines.append(f'{label} & {a["rmse_mag"]:.4f} & {100*a["coverage_95"]:.2f} & {a["mean_95_interval_width_mag"]:.4f} & {b["rmse_mag"]:.4f} & {100*b["coverage_95"]:.2f} & {b["mean_95_interval_width_mag"]:.4f}'+r'\\')
    lines += [r'\bottomrule\end{tabular}',r'\end{table}',r'',
        rf'Without the C driver, the conventional all-image forecast has joint log density {sec["joint_conditional_log_density"]:.3f}, RMSE {secs["rmse_mag"]:.4f}\,mag, 95\% coverage {100*secs["coverage_95"]:.2f}\%, mean interval width {secs["mean_95_interval_width_mag"]:.4f}\,mag and mean interval score {secs["mean_95_interval_score_mag"]:.4f}\,mag. The secondary R1-minus-conventional log-score difference is {scientific(secdelta)}. Primary and secondary targets overlap and cannot be counted as separate confirmations.',
        r'\begin{figure}[t]',
        r'\centering\includegraphics[width=\textwidth]{figures/he0435_revision_unused_prediction.pdf}',
        r'\caption{Unused LCO R observations and the frozen C-driven conventional and conditional R1 predictions. Shading shows marginal 95\% mixture intervals, including the development-only image calibration. The curves and bands can nearly coincide. Straight segments join predictions only at observed epochs and are not an additional fitted variability curve.}\label{fig:he0435-revised-unused}',
        r'\end{figure}',r'',
        r'The experiment tests a fixed, calibrated finite nuisance ensemble. Its short chains have maximum basic $\widehat R=5.19$, so the mixture is not a demonstrated converged posterior and supplies no complete physical-state marginalization. The independent LCO exposures share a retrospectively constructed difference-imaging reference; common reduction errors and prior qualitative publication awareness limit blinding. The central-image HST result is conditional on the empirical PSFs, adopted marginal noise and macro-to-observed flux relation. These qualifications, and any remaining undercoverage in the actual test, remain part of the result. Neither the statistical ensemble nor a close conventional/R1 score establishes a unique envelope state or a BHSM detection; the physical guard remains $2\times10^{-10}$ and no physical admission is asserted.',
        '',
    ]
    return '\n'.join(lines),row


def main():
    if not SCORES.exists():
        raise RuntimeError('Completed scored receipt required; no target parser is available here')
    scores=json.loads(SCORES.read_text());conventional=json.loads(CONVENTIONAL.read_text())
    r1=json.loads(R1.read_text());selection=json.loads(SELECTION.read_text())
    if scores.get('status')!='ONE_FROZEN_INDEPENDENT_EXPOSURE_RETROSPECTIVE_TEST':
        raise RuntimeError('Unexpected or incomplete scored receipt')
    for mode in ('primary','secondary'):
        if any(model not in scores[mode] for model in MODEL_ORDER):
            raise RuntimeError('Scored receipt lacks a prespecified model')
    design=primary_design(scores)
    section,row=build_section(scores,conventional,r1,selection)
    plot_predictions(scores,design)
    SECTION.write_text(section,encoding='utf-8',newline='\n')
    metrics={'status':'GENERATED_FROM_COMPLETED_FROZEN_SCORE_RECEIPTS_NO_REFIT',
        'freeze_commit':scores['freeze_commit'],'source_section_application':row,
        'primary':{m:{'joint_log_density':scores['primary'][m]['joint_conditional_log_density'],
                      'overall':scores['primary'][m]['overall'],'by_image':scores['primary'][m]['by_image']}
                   for m in MODEL_ORDER},
        'secondary':{m:{'joint_log_density':scores['secondary'][m]['joint_conditional_log_density'],
                        'overall':scores['secondary'][m]['overall']} for m in MODEL_ORDER},
        'selected_compact_mass_Msun':conventional['compact_mass_Msun'],
        'fixed_conventional_body_R1_activation':r1['fixed_conventional_body_R1_activation'],
        'input_sha256':{p.relative_to(ROOT).as_posix():digest(p) for p in
                        (SCORES,CONVENTIONAL,R1,SELECTION,Path(__file__))},
        'output_sha256':{p.relative_to(ROOT).as_posix():digest(p) for p in
                         (SECTION,FIGURE.with_suffix('.pdf'),FIGURE.with_suffix('.png'))},
        'raw_or_sealed_target_file_read':False,'physical_admission':False,'physical_guard':2e-10}
    save(OUT/'revision_report_metrics.json',metrics)
    print(json.dumps({'section':SECTION.relative_to(ROOT).as_posix(),
                     'primary':{m:metrics['primary'][m]['overall'] for m in MODEL_ORDER}},indent=2))


if __name__=='__main__':
    main()
