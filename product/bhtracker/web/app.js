(() => {
  'use strict';
  const $ = id => document.getElementById(id);
  const state = { dataset: null, validated: false, validating: false, loading: false, result: null, examples: [], busy: false,
    jobId: null, polling: false, pollTimer: null, elapsedTimer: null, startedAt: null, toastTimer: null };
  const SVG = 'http://www.w3.org/2000/svg';
  const finite = value => typeof value === 'number' && Number.isFinite(value);
  const human = value => String(value ?? '').replaceAll('_', ' ').toLowerCase().replace(/^./, c => c.toUpperCase());
  const family = () => document.querySelector('input[name="family"]:checked').value;
  const configuration = () => ({ family: family(), nuisance: $('nuisance-select').value,
    starts: Math.max(1, Math.min(8, Math.round(Number($('fit-starts').value) || 2))) });

  function number(value, digits = 4) {
    if (!finite(value)) return 'Not available';
    if (value === 0) return '0';
    const absolute = Math.abs(value);
    if (absolute >= 1e6 || absolute < 1e-3) return value.toExponential(digits - 1);
    return new Intl.NumberFormat(undefined, { maximumSignificantDigits: digits }).format(value);
  }
  function vector(values) {
    if (!Array.isArray(values) || values.length < 2) return 'Not available';
    return `x ${number(values[0])} · y ${number(values[1])}`;
  }
  function uncertainty(value, units) {
    if (Array.isArray(value)) return `1σ: x ${finite(value[0]) ? '± ' + number(value[0]) : 'unresolved'} · y ${finite(value[1]) ? '± ' + number(value[1]) : 'unresolved'} ${units}`;
    return finite(value) ? `± ${number(value)} ${units} · local 1σ` : 'Uncertainty not identified';
  }
  function showError(message) { $('global-message').textContent = message; $('global-message').hidden = false; }
  function clearError() { $('global-message').hidden = true; $('global-message').textContent = ''; }
  function toast(message) {
    clearTimeout(state.toastTimer); $('toast').textContent = message; $('toast').hidden = false;
    state.toastTimer = setTimeout(() => { $('toast').hidden = true; }, 4500);
  }
  async function api(path, body) {
    const controller = new AbortController();
    const timer = setTimeout(() => controller.abort(), 30000);
    try {
      const response = await fetch(path, { method: body === undefined ? 'GET' : 'POST',
        headers: body === undefined ? {} : { 'Content-Type': 'application/json' },
        body: body === undefined ? undefined : JSON.stringify(body), signal: controller.signal,
        cache: 'no-store' });
      let payload;
      try { payload = await response.json(); }
      catch { throw new Error(`The local server returned an unreadable response (${response.status}).`); }
      if (!response.ok || payload?.error) throw new Error(payload?.error || `Request failed (${response.status}).`);
      return payload;
    } catch (error) {
      if (error.name === 'AbortError') throw new Error('The local server did not respond within 30 seconds.');
      throw error;
    } finally { clearTimeout(timer); }
  }
  function syncControls() {
    const hasData = !!state.dataset;
    const locked = state.busy || state.validating || state.loading;
    $('load-example').disabled = locked || !state.examples.length;
    $('example-select').disabled = locked || !state.examples.length;
    $('file-input').disabled = locked;
    $('validate-data').disabled = locked || !hasData;
    $('run-fit').disabled = locked || !state.validated;
    $('clear-data').disabled = locked;
    $('source-json').disabled = locked || !hasData;
    $('copy-input').disabled = !hasData;
    $('download-input').disabled = !hasData;
    document.querySelectorAll('input[name="family"], #nuisance-select, #fit-starts').forEach(el => { el.disabled = locked; });
    $('run-fit').querySelector('span').textContent = state.busy ? 'Inference running…' : 'Run inference';
    $('r1-note').hidden = family() !== 'r1';
    $('config-json').textContent = JSON.stringify(configuration(), null, 2);
  }
  function runMessage(message, busy = false) {
    $('run-state').textContent = message; $('run-state').classList.toggle('busy', busy);
  }
  function validationMessage(message, status = '') {
    const element = $('validation-state');
    element.className = `validation-state ${status}`;
    element.querySelector('.status-mark').textContent = status === 'success' ? '✓' : status === 'error' ? '!' : '○';
    element.querySelector('span:last-child').textContent = message;
    $('data-badge').className = `badge ${status || 'neutral'}`;
    $('data-badge').textContent = status === 'success' ? 'Validated' : status === 'error' ? 'Not valid' : state.dataset ? 'Needs validation' : 'No data';
  }
  function clearResults() {
    state.result = null;
    $('result-content').hidden = true; $('empty-results').hidden = false;
    $('raw-result-details').hidden = true; $('raw-result').textContent = '';
    $('empty-results').querySelector('h2').textContent = state.dataset ? 'Ready for inference.' : 'Start with the ordinary lens.';
    $('empty-results').querySelector('.empty-footnote').textContent = state.dataset ? 'No completed fit for this input and configuration.' : 'No fit has been run. No example results are displayed.';
  }
  function setDataset(dataset, name) {
    if (!dataset || typeof dataset !== 'object' || Array.isArray(dataset)) throw new Error('The dataset must be a JSON object.');
    state.dataset = dataset; state.validated = false;
    clearResults(); clearError();
    $('source-json').value = JSON.stringify(dataset, null, 2);
    $('loaded-name').textContent = name || dataset.title || dataset.name || 'Imported dataset';
    $('loaded-file').hidden = false; $('dataset-summary').hidden = true;
    validationMessage('Input loaded. Validate its measurements, units and covariance.');
    runMessage('Validate the dataset to enable the fit.'); syncControls();
  }
  function renderSummary(summary) {
    const container = $('dataset-summary'); container.replaceChildren();
    const entries = [
      ['Epochs / rows', summary?.rows ?? summary?.n_epochs ?? state.dataset?.epochs_yr?.length],
      ['Observer span', finite(summary?.span_observer_years) ? `${number(summary.span_observer_years)} yr` : null],
      ['Provenance', summary?.provenance?.kind ? human(summary.provenance.kind) : null],
      ['Covariance', summary?.covariance ? 'Validated full matrix' : null]
    ];
    for (const [label, value] of entries) {
      if (value === null || value === undefined) continue;
      const div = document.createElement('div'), dt = document.createElement('dt'), dd = document.createElement('dd');
      dt.textContent = label; dd.textContent = value; div.append(dt, dd); container.append(div);
    }
    container.hidden = !container.children.length;
  }
  async function validateDataset() {
    if (!state.dataset || state.busy || state.validating || state.loading) return;
    clearError(); state.validated = false; state.validating = true; syncControls();
    validationMessage('Checking the dataset…');
    try {
      const parsed = JSON.parse($('source-json').value);
      if (!parsed || typeof parsed !== 'object' || Array.isArray(parsed)) throw new Error('The dataset must be a JSON object.');
      state.dataset = parsed;
      const response = await api('/api/validate', { dataset: parsed });
      if (!response.valid) throw new Error(response.message || 'The dataset did not pass validation.');
      state.validated = true;
      validationMessage('Measurements and covariance passed validation.', 'success');
      renderSummary(response.summary); runMessage('Ready. Choose a model and run inference.');
    } catch (error) {
      validationMessage(error instanceof SyntaxError ? `Invalid JSON: ${error.message}` : error.message, 'error');
      $('dataset-summary').hidden = true; runMessage('Resolve the input error before fitting.');
    } finally { state.validating = false; syncControls(); }
  }
  async function loadExamples() {
    try {
      const examples = await api('/api/examples');
      if (!Array.isArray(examples)) throw new Error('The examples endpoint returned an invalid list.');
      state.examples = examples;
      $('example-select').replaceChildren();
      if (!examples.length) {
        const option = new Option('No built-in examples available', ''); $('example-select').append(option);
      } else {
        for (const example of examples) $('example-select').append(new Option(example.title || example.id, example.id));
        updateExampleDescription();
      }
    } catch (error) {
      $('example-select').replaceChildren(new Option('Examples unavailable', ''));
      $('example-description').textContent = `Could not load built-in examples: ${error.message}`;
    } finally { syncControls(); }
  }
  function updateExampleDescription() {
    const example = state.examples.find(item => String(item.id) === $('example-select').value);
    $('example-description').textContent = example?.description || '';
  }
  async function loadSelectedExample() {
    const id = $('example-select').value; if (!id || state.busy || state.loading) return;
    state.loading = true; syncControls(); clearError();
    try {
      const dataset = await api(`/api/examples/${encodeURIComponent(id)}`);
      const example = state.examples.find(item => String(item.id) === id);
      setDataset(dataset, example?.title || id);
    } catch (error) { showError(`Could not load the example: ${error.message}`); }
    finally { state.loading = false; syncControls(); }
  }
  function stopElapsed() { clearInterval(state.elapsedTimer); state.elapsedTimer = null; }
  function finishBusy() {
    state.busy = false; state.polling = false;
    stopElapsed(); clearTimeout(state.pollTimer); $('retry-job').hidden = true; syncControls();
  }
  async function submitFit() {
    if (!state.dataset || !state.validated || state.busy) return;
    clearError(); clearResults();
    state.busy = true; state.jobId = null; state.startedAt = Date.now(); syncControls();
    $('empty-results').querySelector('h2').textContent = 'Inference is running.';
    $('empty-results').querySelector('.empty-footnote').textContent = 'The local server is fitting your dataset. Results appear only when the job completes.';
    runMessage('Submitting the fit…', true);
    state.elapsedTimer = setInterval(() => {
      if (state.polling) runMessage(`Fitting ${family() === 'r1' ? 'the conditional R1 response' : 'the GR compact lens'} · ${Math.floor((Date.now() - state.startedAt) / 1000)} s elapsed`, true);
    }, 1000);
    try {
      const response = await api('/api/fit', { dataset: state.dataset, ...configuration() });
      if (!response.job_id) throw new Error('The server did not return a fit job identifier.');
      state.jobId = response.job_id; state.polling = true; pollJob();
    } catch (error) {
      finishBusy(); showError(`The fit could not start: ${error.message}`);
      runMessage('Fit not started. Review the error and retry.');
      $('empty-results').querySelector('h2').textContent = 'No completed fit.';
      $('empty-results').querySelector('.empty-footnote').textContent = 'The server did not return a completed inference result.';
    }
  }
  async function pollJob() {
    if (!state.jobId) return;
    state.polling = true; $('retry-job').hidden = true;
    try {
      const job = await api(`/api/jobs/${encodeURIComponent(state.jobId)}`);
      if (job.status === 'complete') {
        if (!job.result || typeof job.result !== 'object') throw new Error('The completed job has no result object.');
        state.result = job.result; finishBusy(); renderResult();
        runMessage(`Complete${finite(job.result.runtime_seconds) ? ` · ${number(job.result.runtime_seconds, 3)} s` : ''}. Inspect fit status and admission before interpreting estimates.`);
        toast('Fit complete. Review the diagnostics.');
      } else if (job.status === 'failed') {
        finishBusy(); showError(`The fit failed: ${job.error || 'No error detail was provided.'}`);
        runMessage('Fit failed. No inference result is displayed.');
        $('empty-results').querySelector('h2').textContent = 'No completed fit.';
        $('empty-results').querySelector('.empty-footnote').textContent = 'The failed job returned no result. Review the message before retrying.';
      } else if (job.status === 'running' || job.status === 'queued') {
        runMessage(`Fit running · ${Math.floor((Date.now() - state.startedAt) / 1000)} s elapsed`, true);
        state.pollTimer = setTimeout(pollJob, 1000);
      } else {
        throw new Error(`Unknown job status: ${String(job.status)}.`);
      }
    } catch (error) {
      state.polling = false; stopElapsed();
      showError(`Could not read the job status: ${error.message} The job may still be running on the local server. Use “Retry job status” to reconnect without starting another fit.`);
      runMessage('Status connection interrupted. The fit has not been declared failed.');
      $('retry-job').hidden = false;
    }
  }
  function renderResult() {
    const result = state.result; if (!result) return;
    clearError(); $('empty-results').hidden = true; $('result-content').hidden = false;
    const estimate = result.estimates || {}, sigma = result.uncertainty || {};
    $('mass-value').textContent = number(estimate.mass_msun);
    $('mass-uncertainty').textContent = uncertainty(sigma.mass_msun, 'M☉');
    $('center-value').textContent = vector(estimate.center_mas);
    $('center-uncertainty').textContent = uncertainty(sigma.center_mas, 'mas');
    $('motion-value').textContent = vector(estimate.proper_motion_mas_yr);
    $('motion-uncertainty').textContent = uncertainty(sigma.proper_motion_mas_yr, 'mas/yr');
    $('quality-value').textContent = finite(result.chi2) && finite(result.dof) ? `${number(result.chi2)} / ${number(result.dof)}` : 'Not available';
    const pValue = finite(result.goodness_of_fit_p) ? `p = ${number(result.goodness_of_fit_p, 3)}` : 'Goodness-of-fit probability unavailable';
    $('quality-detail').textContent = `${pValue}${result.optimizer_converged === false ? ' · optimizer not converged' : ''}`;
    const statuses = {
      TRACKED: ['success', 'The fitted motion meets the declared tracking criterion within this model and calibration.'],
      PARTIALLY_IDENTIFIED: ['warning', 'Some parameters are constrained, but the declared motion criterion is not met. Review the unresolved directions.'],
      INSUFFICIENT_INFORMATION: ['warning', 'The data do not identify the requested track. Point estimates are shown for audit; null uncertainties are unresolved, not zero.'],
      FIT_REJECTED: ['error', 'The fit did not pass its admission checks. Displayed estimates are diagnostic and do not establish a valid track.']
    };
    const [tone, description] = statuses[result.status] || ['neutral', 'Review the complete diagnostics and model assumptions.'];
    $('result-status-label').className = `badge ${tone}`;
    $('result-status-label').textContent = human(result.status || 'Status not provided');
    $('result-status-description').textContent = description;
    renderModelContext(result.model);
    renderEnvelope(result.conditional_envelope); renderDiagnostics(result.diagnostics);
    $('raw-result-details').hidden = false; $('raw-result').textContent = JSON.stringify(result, null, 2);
    renderCharts();
  }
  function renderModelContext(model) {
    const container = $('model-context'); container.replaceChildren();
    if (!model || model.family !== 'r1') { container.hidden = true; return; }
    const treatments = {
      fixed_regular_body_proxy: 'Fixed regular-body flux proxy',
      fitted_conditional_flux: 'Fitted conditional scalar flux'
    };
    const label = document.createElement('span'); label.className = 'model-context-label'; label.textContent = 'FLUX TREATMENT';
    const value = document.createElement('strong');
    value.textContent = treatments[model.flux_policy] || (model.flux_policy ? human(model.flux_policy) : 'Not specified in this result');
    container.append(label, value); container.hidden = false;
  }
  function renderEnvelope(envelope) {
    const status = envelope?.status || 'NOT_REQUESTED';
    $('envelope-status').textContent = human(status);
    const container = $('envelope-detail'); container.replaceChildren();
    if (!envelope || status === 'NOT_REQUESTED') {
      const p = document.createElement('p'); p.className = 'field-help';
      p.textContent = 'No conditional envelope was requested for this fit.'; container.append(p); return;
    }
    const description = document.createElement('p'); description.className = 'field-help';
    description.textContent = status === 'NOT_ADMITTED'
      ? 'This response is outside the model’s admitted approximation domain for the supplied configuration.'
      : 'Sensitivity is assessed only for boundary states inside the declared model-admission bound.';
    container.append(description);
    if (status === 'NOT_ADMITTED' && Array.isArray(envelope.reasons)) {
      for (const reason of envelope.reasons) {
        const p = document.createElement('p'); p.className = 'envelope-explanation';
        p.textContent = typeof reason === 'string' ? reason : reason.message || JSON.stringify(reason);
        container.append(p);
      }
    }
    const cards = document.createElement('div'); cards.className = 'envelope-measures';
    function card(label, value, detail, wide = false) {
      const article = document.createElement('div'); article.className = `envelope-measure${wide ? ' wide' : ''}`;
      const heading = document.createElement('div'), number = document.createElement('strong'), explanation = document.createElement('p');
      heading.className = 'envelope-measure-label'; heading.textContent = label;
      number.textContent = value; explanation.textContent = detail;
      article.append(heading, number, explanation); cards.append(article);
    }
    if (finite(envelope.boundary_amplitude_cap)) card('Admitted boundary norm', `≤ ${number(envelope.boundary_amplitude_cap, 3)}`,
      'Sufficient bound on the normalized scalar state. This is not a predicted astrophysical amplitude.');
    if (finite(envelope.maximum_admitted_local_snr)) card('Largest admitted local SNR', number(envelope.maximum_admitted_local_snr, 3),
      'Largest local signal-to-noise ratio after source and lens freedoms are profiled.');
    if (envelope.all_admitted_states_below_unit_snr === true) {
      card('Do the data resolve admitted states?', 'No',
        'Every state inside the stated norm bound remains below unit local SNR for this covariance.', true);
    } else if (status === 'UNRESOLVED') {
      card('Do the data resolve admitted states?', 'Not established',
        'The analysis does not resolve all conditional envelope directions within the declared model.', true);
    } else if (status === 'INFORMATION_AVAILABLE') {
      card('Do the data resolve admitted states?', 'Local sensitivity available',
        'The model retains conditional response information. This is a sensitivity result, not a deformation detection.', true);
    }
    if (cards.children.length) container.append(cards);
    if (envelope.unconstrained_estimate_inside_admitted_ball === false) {
      const note = document.createElement('div'); note.className = 'envelope-domain-note';
      const title = document.createElement('strong'), body = document.createElement('p');
      title.textContent = 'Unconstrained estimate outside the admitted domain';
      body.textContent = 'That unrestricted fit value is not a physical boundary-state measurement. It must not be interpreted as a detected deformation.';
      note.append(title, body); container.append(note);
    }
    const technical = document.createElement('details'); technical.className = 'envelope-technical';
    const summary = document.createElement('summary'); summary.textContent = 'Technical admission details';
    const pre = document.createElement('pre'); pre.className = 'admission-details'; pre.textContent = JSON.stringify(envelope, null, 2);
    technical.append(summary, pre); container.append(technical);
  }
  function renderDiagnostics(diagnostics) {
    const container = $('diagnostics-list'); container.replaceChildren();
    const list = Array.isArray(diagnostics) ? diagnostics : [];
    if (!list.length) {
      const p = document.createElement('p'); p.className = 'field-help'; p.textContent = 'The server supplied no additional diagnostics.'; container.append(p); return;
    }
    for (const item of list) {
      const level = ['warning', 'error'].includes(item.level) ? item.level : 'info';
      const row = document.createElement('div'); row.className = `diagnostic ${level}`;
      const icon = document.createElement('span'); icon.className = 'diagnostic-icon'; icon.textContent = level === 'info' ? 'i' : '!'; icon.setAttribute('aria-label', level);
      const body = document.createElement('div'), code = document.createElement('div'), message = document.createElement('p');
      code.className = 'diagnostic-code'; code.textContent = item.code || level.toUpperCase(); message.textContent = item.message || '';
      body.append(code, message); row.append(icon, body); container.append(row);
    }
  }

  function svgElement(tag, attributes = {}, text) {
    const element = document.createElementNS(SVG, tag);
    Object.entries(attributes).forEach(([key, value]) => element.setAttribute(key, value));
    if (text !== undefined) element.textContent = text;
    return element;
  }
  function extent(values) {
    const valid = values.filter(finite); if (!valid.length) return [-1, 1];
    let lo = Math.min(...valid), hi = Math.max(...valid);
    const span = hi - lo || Math.max(Math.abs(lo) * .15, .01);
    return [lo - .12 * span, hi + .12 * span];
  }
  function tick(value) {
    if (value === 0 || Math.abs(value) < 1e-13) return '0';
    const abs = Math.abs(value);
    return abs >= 1e4 || abs < .001 ? value.toExponential(1) : String(Number(value.toPrecision(4)));
  }
  function plotFrame(container, title, xExtent, yExtent, xLabel, yLabel, equalAspect = false, height = 300) {
    container.replaceChildren();
    const width = Math.max(340, container.clientWidth || 600);
    const margin = { l: 66, r: 20, t: 18, b: 47 }, pw = width - margin.l - margin.r, ph = height - margin.t - margin.b;
    let [xmin, xmax] = xExtent, [ymin, ymax] = yExtent;
    if (equalAspect) {
      const unit = Math.max((xmax - xmin) / pw, (ymax - ymin) / ph), cx = (xmin + xmax) / 2, cy = (ymin + ymax) / 2;
      xmin = cx - unit * pw / 2; xmax = cx + unit * pw / 2; ymin = cy - unit * ph / 2; ymax = cy + unit * ph / 2;
    }
    const x = value => margin.l + (value - xmin) / (xmax - xmin) * pw;
    const y = value => margin.t + ph - (value - ymin) / (ymax - ymin) * ph;
    const svg = svgElement('svg', { viewBox: `0 0 ${width} ${height}`, role: 'img', 'aria-label': title });
    svg.append(svgElement('title', {}, title));
    const count = width < 470 ? 4 : 5;
    for (let i = 0; i <= count; i++) {
      const xv = xmin + (xmax - xmin) * i / count, yv = ymin + (ymax - ymin) * i / count;
      svg.append(svgElement('line', { x1: x(xv), x2: x(xv), y1: margin.t, y2: margin.t + ph, stroke: '#e4ebee', 'stroke-width': 1 }));
      svg.append(svgElement('line', { x1: margin.l, x2: margin.l + pw, y1: y(yv), y2: y(yv), stroke: '#e4ebee', 'stroke-width': 1 }));
      svg.append(svgElement('text', { x: x(xv), y: height - 27, 'text-anchor': 'middle', fill: '#627987', 'font-family': 'Consolas,monospace', 'font-size': 10 }, tick(xv)));
      svg.append(svgElement('text', { x: margin.l - 9, y: y(yv) + 3, 'text-anchor': 'end', fill: '#627987', 'font-family': 'Consolas,monospace', 'font-size': 10 }, tick(yv)));
    }
    svg.append(svgElement('text', { x: margin.l + pw / 2, y: height - 5, 'text-anchor': 'middle', fill: '#526c7b', 'font-size': 11 }, xLabel));
    svg.append(svgElement('text', { transform: `translate(13 ${margin.t + ph / 2}) rotate(-90)`, 'text-anchor': 'middle', fill: '#526c7b', 'font-size': 11 }, yLabel));
    container.append(svg);
    return { svg, x, y, width, height };
  }
  function line(frame, points, color = '#007f93', dash = '') {
    const valid = points.filter(p => p.every(finite)); if (!valid.length) return;
    const path = valid.map((point, i) => `${i ? 'L' : 'M'}${frame.x(point[0])},${frame.y(point[1])}`).join(' ');
    const attributes = { d: path, fill: 'none', stroke: color, 'stroke-width': 2, 'stroke-linecap': 'round', 'stroke-linejoin': 'round' };
    if (dash) attributes['stroke-dasharray'] = dash;
    frame.svg.append(svgElement('path', attributes));
  }
  function marker(frame, xValue, yValue, options = {}) {
    if (!finite(xValue) || !finite(yValue)) return;
    const point = svgElement('circle', { cx: frame.x(xValue), cy: frame.y(yValue), r: options.radius || 3.2,
      fill: options.hollow ? '#fff' : '#172d3c', stroke: options.color || '#172d3c', 'stroke-width': 1.2 });
    if (options.label) point.append(svgElement('title', {}, options.label));
    frame.svg.append(point);
  }
  function errorBar(frame, xValue, yValue, sx, sy, color = '#849ca9') {
    const attributes = { stroke: color, 'stroke-width': 1, opacity: .65 };
    const x = frame.x(xValue), y = frame.y(yValue);
    if (finite(sx) && sx >= 0) {
      const left = frame.x(xValue - sx), right = frame.x(xValue + sx);
      frame.svg.append(svgElement('line', { ...attributes, x1: left, x2: right, y1: y, y2: y }));
      frame.svg.append(svgElement('line', { ...attributes, x1: left, x2: left, y1: y - 2, y2: y + 2 }));
      frame.svg.append(svgElement('line', { ...attributes, x1: right, x2: right, y1: y - 2, y2: y + 2 }));
    }
    if (finite(sy) && sy >= 0) {
      const top = frame.y(yValue + sy), bottom = frame.y(yValue - sy);
      frame.svg.append(svgElement('line', { ...attributes, x1: x, x2: x, y1: top, y2: bottom }));
      frame.svg.append(svgElement('line', { ...attributes, x1: x - 2, x2: x + 2, y1: top, y2: top }));
      frame.svg.append(svgElement('line', { ...attributes, x1: x - 2, x2: x + 2, y1: bottom, y2: bottom }));
    }
  }
  function chartEmpty(container, message) { container.replaceChildren(); const p = document.createElement('div'); p.className = 'chart-empty'; p.textContent = message; container.append(p); }
  function rows() {
    const result = state.result || {};
    return { observed: Array.isArray(result.observations) ? result.observations : [],
      predicted: Array.isArray(result.predictions) ? result.predictions : [],
      epochs: Array.isArray(result.epochs_yr) ? result.epochs_yr : [] };
  }
  function covariance() { return state.result?.covariance || state.dataset?.covariance || []; }
  function blockCovariance(index) {
    const full = covariance(), offset = index * 6;
    if (!Array.isArray(full) || !Array.isArray(full[offset]) || full[offset].length < offset + 6) return null;
    return Array.from({ length: 6 }, (_, i) => full[offset + i]?.slice(offset, offset + 6));
  }
  function variance(row, cov) {
    if (!cov || cov.some(r => !Array.isArray(r) || r.length < 6)) return null;
    let value = 0;
    for (let i = 0; i < 6; i++) for (let j = 0; j < 6; j++) {
      if (!finite(cov[i][j])) return null; value += row[i] * cov[i][j] * row[j];
    }
    return value >= -1e-16 ? Math.sqrt(Math.max(0, value)) : null;
  }
  function pairCenter(row, index) {
    if (!Array.isArray(row) || row.slice(0, 5).some(value => !finite(value))) return null;
    const k = Math.exp(row[4] / 2); if (!finite(k)) return null;
    const wp = 1 / (1 + k), wm = k / (1 + k);
    const point = [wp * row[0] + wm * row[2], wp * row[1] + wm * row[3]];
    const jx = [wp, 0, wm, 0, .5 * k * (row[2] - row[0]) / (1 + k) ** 2, 0];
    const jy = [0, wp, 0, wm, .5 * k * (row[3] - row[1]) / (1 + k) ** 2, 0];
    const cov = blockCovariance(index);
    return { point, sigma: [variance(jx, cov), variance(jy, cov)], index };
  }
  function renderTrack() {
    const container = $('track-chart'), { observed, predicted, epochs } = rows();
    const track = state.result?.center_track_mas;
    const trackSigma = state.result?.center_track_sigma_mas || [];
    const legend = container.closest('.chart-panel').querySelector('.plot-legend');
    function setLegend(observedLabel, fittedLabel) {
      const first = document.createElement('span'), second = document.createElement('span');
      const dot = document.createElement('i'), stroke = document.createElement('i');
      dot.className = 'legend-dot'; stroke.className = 'legend-line';
      first.append(dot, document.createTextNode(observedLabel)); second.append(stroke, document.createTextNode(fittedLabel));
      legend.replaceChildren(first, second);
    }
    if (Array.isArray(track) && track.length && track.every(p => Array.isArray(p) && p.length >= 2 && p.slice(0, 2).every(finite))) {
      setLegend('GR pair diagnostic', 'Fitted center');
      const measured = observed.map(pairCenter).filter(Boolean);
      const xs = [], ys = [];
      track.forEach((p, i) => { xs.push(p[0]); ys.push(p[1]); if (finite(trackSigma[i]?.[0])) xs.push(p[0] - trackSigma[i][0], p[0] + trackSigma[i][0]); if (finite(trackSigma[i]?.[1])) ys.push(p[1] - trackSigma[i][1], p[1] + trackSigma[i][1]); });
      measured.forEach(({ point: p, sigma: sig }) => { xs.push(p[0]); ys.push(p[1]); if (finite(sig[0])) xs.push(p[0] - sig[0], p[0] + sig[0]); if (finite(sig[1])) ys.push(p[1] - sig[1], p[1] + sig[1]); });
      const frame = plotFrame(container, 'Lens center trajectory with marginal uncertainties', extent(xs), extent(ys), 'x [mas]', 'y [mas]', true, 320);
      track.forEach((p, i) => errorBar(frame, p[0], p[1], trackSigma[i]?.[0], trackSigma[i]?.[1], '#3aa8b6'));
      line(frame, track);
      measured.forEach(({ point: p, sigma: sig, index }) => { errorBar(frame, p[0], p[1], sig[0], sig[1]); marker(frame, p[0], p[1], { label: `Pair-derived GR center, epoch ${epochs[index] ?? index}: (${number(p[0])}, ${number(p[1])}) mas` }); });
      if (track.length > 1) {
        [0, track.length - 1].forEach(i => { const p = track[i]; frame.svg.append(svgElement('text', { x: frame.x(p[0]) + 7, y: frame.y(p[1]) - 8, fill: '#007f93', 'font-size': 9 }, `t=${number(epochs[i] ?? i, 3)}`)); });
      }
      $('track-title').textContent = 'Lens-center trajectory';
      $('track-caption').textContent = 'Dots: centers inferred from the observed pairs using the isolated GR identity, a diagnostic rather than an independent observable. Cyan: fitted center track. Bars show marginal 1σ uncertainties where identified; the fit retains the full covariance.';
      return;
    }
    setLegend('Observed images', 'Fitted images');
    if (!observed.length && !predicted.length) { chartEmpty(container, 'The server did not supply astrometric observations or predictions.'); return; }
    const xs = [], ys = [];
    [...observed, ...predicted].forEach(row => { if (Array.isArray(row)) { xs.push(row[0], row[2]); ys.push(row[1], row[3]); } });
    const frame = plotFrame(container, 'Observed and fitted primary image positions', extent(xs), extent(ys), 'x [mas]', 'y [mas]', true, 320);
    line(frame, predicted.map(row => [row[0], row[1]]));
    line(frame, predicted.map(row => [row[2], row[3]]), '#007f93', '5 4');
    observed.forEach((row, i) => {
      const cov = blockCovariance(i);
      [0, 2].forEach(offset => {
        const sx = finite(cov?.[offset]?.[offset]) ? Math.sqrt(Math.max(0, cov[offset][offset])) : null;
        const sy = finite(cov?.[offset + 1]?.[offset + 1]) ? Math.sqrt(Math.max(0, cov[offset + 1][offset + 1])) : null;
        errorBar(frame, row[offset], row[offset + 1], sx, sy);
        marker(frame, row[offset], row[offset + 1], { hollow: offset === 2, label: `${offset ? 'Inner' : 'Outer'} image, epoch ${epochs[i] ?? i}` });
      });
    });
    $('track-title').textContent = 'Observed and fitted image plane';
    $('track-caption').textContent = 'Filled markers and solid line: positive-parity image. Open markers and dashed line: negative-parity image. Coordinate bars show supplied marginal 1σ uncertainty where available.';
  }
  function renderSeries() {
    const container = $('series-chart'), { observed, predicted, epochs } = rows();
    if (!observed.length || !epochs.length) { chartEmpty(container, 'The server did not supply an epoch series.'); return; }
    const index = Number($('series-select').value);
    const labels = ['Positive-parity image x [mas]', 'Positive-parity image y [mas]', 'Negative-parity image x [mas]', 'Negative-parity image y [mas]', 'ln(flux ratio)', 'Time delay [s]'];
    const values = [];
    const errors = observed.map((row, i) => {
      values.push(row[index]); const cov = blockCovariance(i), value = cov?.[index]?.[index];
      const sigma = finite(value) && value >= 0 ? Math.sqrt(value) : null;
      if (finite(sigma)) values.push(row[index] - sigma, row[index] + sigma);
      if (finite(predicted[i]?.[index])) values.push(predicted[i][index]);
      return sigma;
    });
    const frame = plotFrame(container, `Observed and fitted ${labels[index]} over epochs`, extent(epochs), extent(values), 'Epoch [observer yr]', labels[index], false, 270);
    const ordering = epochs.map((time, i) => ({ time, i })).sort((a, b) => a.time - b.time);
    line(frame, ordering.filter(item => finite(predicted[item.i]?.[index])).map(item => [item.time, predicted[item.i][index]]));
    observed.forEach((row, i) => { errorBar(frame, epochs[i], row[index], null, errors[i]); marker(frame, epochs[i], row[index], { label: `Epoch ${epochs[i]}: observed ${number(row[index])}${finite(errors[i]) ? ' ± ' + number(errors[i]) : ''}; fit ${number(predicted[i]?.[index])}` }); });
    $('series-caption').textContent = errors.some(finite)
      ? 'Markers and bars: measured values with input marginal 1σ uncertainties. Cyan: fitted model. The inference uses the complete covariance, including correlations not shown by individual bars.'
      : 'Markers: supplied measurements. Cyan: fitted model. The server did not supply marginal uncertainties for this display.';
  }
  function renderCharts() { if (state.result) { renderTrack(); renderSeries(); } }
  function download(content, filename, type) {
    const url = URL.createObjectURL(new Blob([content], { type }));
    const anchor = document.createElement('a'); anchor.href = url; anchor.download = filename;
    document.body.append(anchor); anchor.click(); anchor.remove(); setTimeout(() => URL.revokeObjectURL(url), 1000);
  }

  $('example-select').addEventListener('change', updateExampleDescription);
  $('load-example').addEventListener('click', loadSelectedExample);
  $('file-input').addEventListener('change', async event => {
    const file = event.target.files?.[0]; if (!file || state.busy || state.loading) return;
    state.loading = true; syncControls();
    try {
      if (file.size > 12 * 1024 * 1024) throw new Error('The JSON file exceeds the 12 MB input limit.');
      const text = (await file.text()).replace(/^\uFEFF/, '');
      setDataset(JSON.parse(text), file.name);
    } catch (error) { showError(`Could not import the dataset: ${error.message}`); }
    finally { event.target.value = ''; state.loading = false; syncControls(); }
  });
  $('clear-data').addEventListener('click', () => {
    state.dataset = null; state.validated = false; clearResults(); clearError();
    $('source-json').value = ''; $('loaded-file').hidden = true; $('dataset-summary').hidden = true;
    validationMessage('Load an example or import a dataset to begin.');
    runMessage('Validate a dataset to enable the fit.');
    $('empty-results').querySelector('h2').textContent = 'Start with the ordinary lens.';
    $('empty-results').querySelector('.empty-footnote').textContent = 'No fit has been run. No example results are displayed.';
    syncControls();
  });
  $('source-json').addEventListener('input', () => {
    state.validated = false; clearResults(); $('dataset-summary').hidden = true;
    validationMessage('Input edited. Validate the updated JSON before fitting.'); runMessage('Revalidate the edited dataset.'); syncControls();
  });
  $('validate-data').addEventListener('click', validateDataset);
  document.querySelectorAll('input[name="family"], #nuisance-select, #fit-starts').forEach(element => element.addEventListener('change', () => {
    if (element.id === 'fit-starts') element.value = configuration().starts;
    clearResults(); syncControls();
    runMessage(state.validated ? 'Configuration changed. Ready to run inference.' : 'Validate a dataset to enable the fit.');
  }));
  $('run-fit').addEventListener('click', submitFit);
  $('retry-job').addEventListener('click', () => { clearError(); state.polling = true; runMessage('Reconnecting to the fit job…', true); pollJob(); });
  $('series-select').addEventListener('change', renderSeries);
  $('copy-input').addEventListener('click', async () => {
    try { await navigator.clipboard.writeText($('source-json').value); toast('Dataset JSON copied.'); }
    catch { toast('Clipboard access is unavailable. Use Download input instead.'); }
  });
  $('download-input').addEventListener('click', () => download($('source-json').value, 'bhtracker-dataset.json', 'application/json'));
  $('export-json').addEventListener('click', () => { if (state.result) download(JSON.stringify(state.result, null, 2), 'bhtracker-result.json', 'application/json'); });
  $('export-html').addEventListener('click', async () => {
    if (!state.result) return;
    $('export-html').disabled = true;
    try {
      const report = await api('/api/report', { result: state.result });
      if (typeof report.html !== 'string' || !report.html) throw new Error('The server did not supply report HTML.');
      download(report.html, 'bhtracker-report.html', 'text/html;charset=utf-8'); toast('HTML report exported.');
    } catch (error) { showError(`Could not export the report: ${error.message}`); }
    finally { $('export-html').disabled = false; }
  });
  let resizeTimer;
  window.addEventListener('resize', () => { clearTimeout(resizeTimer); resizeTimer = setTimeout(renderCharts, 150); });
  window.addEventListener('beforeunload', () => { clearTimeout(state.pollTimer); stopElapsed(); });
  syncControls(); loadExamples();
})();
