"""Portable, self-contained research exports with traceable assumptions."""
from pathlib import Path
import csv
import hashlib
import html
import json
import numpy as np

def report_html(result):
    r=result;e=r['estimates'];u=r['uncertainty']
    escape=lambda value:html.escape(str(value))
    def fmt(v):
        if v is None:return 'unresolved'
        if isinstance(v,list):return ', '.join(fmt(x) for x in v)
        if isinstance(v,(float,int)):return f'{v:.7g}'
        return escape(v)
    items=[('Mass (solar masses)',e['mass_msun'],u['mass_msun']),
           ('Center at reference epoch (mas)',e['center_mas'],u['center_mas']),
           ('Proper motion (mas / observer year)',e['proper_motion_mas_yr'],u['proper_motion_mas_yr'])]
    rows=''.join(f'<tr><th>{label}</th><td>{fmt(value)}</td><td>{fmt(err)}</td></tr>' for label,value,err in items)
    diagnostics=''.join(f"<li><strong>{escape(item['code'])}</strong>: {escape(item['message'])}</li>" for item in r['diagnostics'])
    return f'''<!doctype html><html lang="en"><meta charset="utf-8"><meta name="viewport" content="width=device-width">
<title>Blackhole Tracker research report</title><style>body{{font:16px/1.55 system-ui;max-width:1050px;margin:3rem auto;padding:0 1.5rem;color:#172b40}}h1{{font-size:2.2rem}}table{{border-collapse:collapse;width:100%}}th,td{{padding:.7rem;border-bottom:1px solid #dce3e9;text-align:left}}pre{{background:#f1f5f8;padding:1rem;overflow:auto;white-space:pre-wrap}}small{{color:#4f657a}}.status{{padding:1rem;background:#edf6fa;border-left:4px solid #167e9b}}</style>
<h1>Blackhole Tracker</h1><p class="status">{escape(r['status'])} · {escape(r['data_summary']['provenance']['kind'])} measurements</p>
<p>Reference epoch {fmt(r['reference_epoch_yr'])}; family {escape(r['model']['family'])}; nuisance {escape(r['model']['nuisance'])}.</p>
<table><tr><th>Quantity</th><th>Estimate</th><th>Local 1σ</th></tr>{rows}</table>
<p>χ² {fmt(r['chi2'])} for {r['dof']} degrees of freedom; goodness-of-fit p={fmt(r['goodness_of_fit_p'])}. Uncertainties are local known-covariance approximations; active bounds and nulls require the stated qualification.</p>
<h2>Interpretation</h2><ul>{diagnostics}</ul><h2>Conditional envelope</h2><pre>{escape(json.dumps(r['conditional_envelope'],indent=2))}</pre>
<h2>Model and domain</h2><pre>{escape(json.dumps(r['model'],indent=2))}</pre><pre>{escape(json.dumps(r.get('physics_admission'),indent=2))}</pre>
<h2>Reproduction</h2><pre>{escape(json.dumps(r['data_summary'],indent=2))}</pre><pre>{escape(json.dumps(r['runtime'],indent=2))}</pre>
<p><small>Product {escape(r['product_version'])}. A compact exterior fit does not identify a horizon. Synthetic examples and supplied uncertainties are not empirical detections or instrument forecasts. Full measurements, covariance, predictions and parameter covariance accompany the JSON export.</small></p></html>'''

def export_result(result,directory):
    dest=Path(directory);dest.mkdir(parents=True,exist_ok=True)
    (dest/'result.json').write_text(json.dumps(result,indent=2,allow_nan=False)+'\n',encoding='utf8')
    (dest/'report.html').write_text(report_html(result),encoding='utf8')
    with (dest/'predictions.csv').open('w',newline='',encoding='utf8') as stream:
        writer=csv.writer(stream)
        from .data import FIELDS
        writer.writerow(['epoch_yr',*FIELDS])
        for epoch,row in zip(result['epochs_yr'],result['predictions']):writer.writerow([epoch,*row])
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    fig,axes=plt.subplots(1,2,figsize=(10,4),layout='constrained')
    times=np.array(result['epochs_yr']);track=np.array(result['center_track_mas'])
    err=np.array([[np.nan if v is None else v for v in row] for row in result['center_track_sigma_mas']])
    axes[0].errorbar(times,track[:,0],yerr=err[:,0],fmt='o',label='x center')
    axes[0].errorbar(times,track[:,1],yerr=err[:,1],fmt='s',label='y center')
    axes[0].set(xlabel='Observer epoch (year)',ylabel='Lens center (mas)');axes[0].legend(frameon=False)
    axes[1].plot(result['whitened_residuals'],'o',ms=3);axes[1].axhline(0,color='.6',lw=.8)
    axes[1].set(xlabel='Whitened component (full covariance)',ylabel='Residual / supplied noise')
    fig.suptitle('Blackhole Tracker · '+result['status'])
    fig.savefig(dest/'track_and_residuals.png',dpi=170);plt.close(fig)
    files=[p for p in dest.iterdir() if p.is_file() and p.name!='manifest.json']
    (dest/'manifest.json').write_text(json.dumps({'dataset_sha256':result['data_summary']['sha256'],
        'files':{p.name:hashlib.sha256(p.read_bytes()).hexdigest() for p in files}},indent=2)+'\n')
    return str(dest.resolve())
