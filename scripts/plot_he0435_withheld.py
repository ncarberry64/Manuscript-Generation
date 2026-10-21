"""Plot the previously frozen/scored observations; no fits or model selection."""
from pathlib import Path
import json

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np

ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT/'artifacts/strong_lensing_observed_20261020'


def main():
    result = json.loads((OUT/'withheld_scores.json').read_text())
    design = result['target_design']
    t = np.asarray(design['time_days'])-57447.03918
    image = np.asarray(design['image'])
    value, error = np.asarray(design['magnitude']), np.asarray(design['error_mag'])
    fig, axes = plt.subplots(3, 1, figsize=(7.1, 6.7), sharex=True,
                             gridspec_kw={'hspace': .13})
    conv = result['primary']['conventional_macro']
    r1 = result['primary']['conditional_R1_macro']
    for ax, i in zip(axes, [0, 2, 3]):
        keep = image == i
        mean, sigma = np.asarray(conv['predicted_magnitude'])[keep], np.asarray(conv['predictive_sigma_mag'])[keep]
        order = np.argsort(t[keep]); x = t[keep][order]
        ax.fill_between(x, (mean-1.95996398454*sigma)[order],
                        (mean+1.95996398454*sigma)[order], color='#246c96', alpha=.17,
                        label='Conventional 95% conditional interval')
        ax.plot(x, mean[order], color='#246c96', lw=1.1, label='Conventional mean')
        ax.plot(x, np.asarray(r1['predicted_magnitude'])[keep][order], color='#a54e00',
                ls='--', lw=1.0, label='Conditional R1 mean')
        ax.errorbar(t[keep], value[keep], error[keep], fmt='.', ms=3.3,
                    color='#242424', alpha=.85, elinewidth=.55, label='Withheld observations')
        ax.set_ylabel('Image '+ 'ABCD'[i]+'\nR magnitude + const.')
        ax.invert_yaxis(); ax.grid(alpha=.2)
        ax.text(.985, .94, 'No target refit', transform=ax.transAxes, ha='right', va='top', fontsize=8)
    axes[0].legend(fontsize=7, loc='lower right', ncol=2, framealpha=.95)
    axes[-1].set_xlabel('Native observation day minus training cutoff')
    fig.suptitle('HE 0435−1223: later B conditions predictions of withheld A, C, D', fontsize=10.5, y=.965)
    fig.subplots_adjust(left=.13, right=.98, top=.925, bottom=.085)
    figdir = ROOT/'manuscript/figures'
    for suffix in ['pdf', 'png', 'svg']:
        fig.savefig(figdir/('he0435_withheld_prediction.'+suffix), dpi=180,
                    metadata={'Creator': 'Frozen HE0435 retrospective application'} if suffix=='pdf' else None)
    plt.close(fig)
    summary = {'figure_interpretation': 'Point-parameter conditional predictions; batch later-B conditioning, no target fit. Original COSMOGRAIL/CDS data. Conventional/R1 curves may overlap; no separation is visually exaggerated.',
               'plotted_source': 'withheld_scores.json', 'models': {}}
    for name, score in result['primary'].items():
        summary['models'][name] = {'joint_log_density': score['joint_conditional_test_log_density'],
                                   **score['overall']}
    (OUT/'plot_summary.json').write_text(json.dumps(summary,indent=2)+'\n',encoding='utf-8',newline='\n')


if __name__ == '__main__':
    main()
