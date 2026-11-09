/* GreenLight · charts/charts.js
 * Dependency-free SVG charts: trajectory (observed + predictive bands +
 * previous-version ghost + events), parameter posterior bars, phenotype bar,
 * sparkline. All charts are responsive (viewBox) and theme-aware (CSS vars).
 */
(function (root) {
  const GL = root.GL = root.GL || {};
  const U = GL.util, C = GL.charts = {};
  const NS = 'http://www.w3.org/2000/svg';
  const el = (tag, attrs, parent) => { const e = document.createElementNS(NS, tag); for (const k in attrs) e.setAttribute(k, attrs[k]); if (parent) parent.appendChild(e); return e; };
  const t = (k, p) => GL.i18n.t(k, p);

  function niceTicks(lo, hi, n) {
    const span = hi - lo, step0 = Math.pow(10, Math.floor(Math.log10(span / n)));
    const err = span / n / step0; const step = err >= 7.5 ? 10 * step0 : err >= 3.5 ? 5 * step0 : err >= 1.5 ? 2 * step0 : step0;
    const out = []; for (let v = Math.ceil(lo / step) * step; v <= hi + 1e-9; v += step) out.push(+v.toFixed(10)); return out;
  }
  function logTicks(lo, hi) { const out = []; for (let e = Math.floor(Math.log10(lo)); e <= Math.ceil(Math.log10(hi)); e++) for (const m of [1, 2, 5]) { const v = m * Math.pow(10, e); if (v >= lo && v <= hi) out.push(v); } return out; }

  // opt: { grid, bands:{q05..q95}, prevBands, obs:[{day,value}], obsFuture:[{day,value}], asOfDay, variable:{label,unit,log,thresholds}, events:[{type,day_start,day_end}], xRange:[a,b], height, ldWindow:[-5,-3] }
  C.trajectory = function (container, opt) {
    container.innerHTML = '';
    const W = 760, H = opt.height || 300, m = { l: 52, r: 16, t: 14, b: 42 };
    const svg = el('svg', { viewBox: `0 0 ${W} ${H}`, class: 'gl-chart gl-traj', role: 'img' }, container);
    const x0 = opt.xRange ? opt.xRange[0] : -6, x1 = opt.xRange ? opt.xRange[1] : 120;
    const allVals = [];
    if (opt.bands) { allVals.push(...opt.bands.q05, ...opt.bands.q95); }
    if (opt.prevBands) { allVals.push(...opt.prevBands.q05, ...opt.prevBands.q95); }
    (opt.obs || []).forEach(o => allVals.push(o.value)); (opt.obsFuture || []).forEach(o => allVals.push(o.value));
    (opt.variable.thresholds || []).forEach(v => allVals.push(v));
    const log = !!opt.variable.log;
    const floor = log ? (opt.variable.id === 'plt' ? 1 : opt.variable.id === 'tumor_burden' ? 0.1 : 0.01) : 0;
    let lo = Math.min(...allVals.filter(v => isFinite(v))), hi = Math.max(...allVals.filter(v => isFinite(v)));
    if (!isFinite(lo)) { lo = 0; hi = 1; }
    if (log) { lo = Math.max(floor, lo * 0.8); hi = hi * 1.25; if (hi / lo < 4) hi = lo * 4; } else { lo = 0; hi = hi * 1.1 || 1; }
    const X = (d) => m.l + (d - x0) / (x1 - x0) * (W - m.l - m.r);
    const Y = log ? (v) => m.t + (1 - (Math.log(Math.max(v, lo)) - Math.log(lo)) / (Math.log(hi) - Math.log(lo))) * (H - m.t - m.b)
                  : (v) => m.t + (1 - (v - lo) / (hi - lo)) * (H - m.t - m.b);
    // LD shading
    const ld = opt.ldWindow || [-5, -3];
    el('rect', { x: X(ld[0]), y: m.t, width: X(ld[1] + 1) - X(ld[0]), height: H - m.t - m.b, class: 'gl-ld' }, svg);
    // gridlines + y axis
    const ticks = log ? logTicks(lo, hi) : niceTicks(lo, hi, 5);
    ticks.forEach(v => { el('line', { x1: m.l, x2: W - m.r, y1: Y(v), y2: Y(v), class: 'gl-grid' }, svg); const tx = el('text', { x: m.l - 6, y: Y(v) + 4, class: 'gl-tick', 'text-anchor': 'end' }, svg); tx.textContent = U.fmt(v, v >= 100 ? 0 : v >= 10 ? 0 : v >= 1 ? 1 : 2); });
    const xt = niceTicks(x0, x1, 8);
    xt.forEach(d => { el('line', { x1: X(d), x2: X(d), y1: H - m.b, y2: H - m.b + 4, class: 'gl-axis' }, svg); const tx = el('text', { x: X(d), y: H - m.b + 16, class: 'gl-tick', 'text-anchor': 'middle' }, svg); tx.textContent = d; });
    el('line', { x1: m.l, x2: W - m.r, y1: H - m.b, y2: H - m.b, class: 'gl-axis' }, svg);
    const xl = el('text', { x: (m.l + W - m.r) / 2, y: H - 4, class: 'gl-axis-label', 'text-anchor': 'middle' }, svg); xl.textContent = t('chart.xaxis');
    const yl = el('text', { x: 12, y: m.t + 4, class: 'gl-axis-label', transform: `rotate(-90 12 ${m.t + 4})`, 'text-anchor': 'end' }, svg); yl.textContent = opt.variable.label + (opt.variable.unit ? ' (' + opt.variable.unit + ')' : '');
    // thresholds
    (opt.variable.thresholds || []).forEach(v => { if (v > lo && v < hi) { el('line', { x1: m.l, x2: W - m.r, y1: Y(v), y2: Y(v), class: 'gl-thr' }, svg); const tx = el('text', { x: W - m.r - 2, y: Y(v) - 3, class: 'gl-thr-label', 'text-anchor': 'end' }, svg); tx.textContent = v; } });
    const band = (b, i0, i1, cls, qlo, qhi) => {
      let d = '';
      for (let i = i0; i <= i1; i++) d += (i === i0 ? 'M' : 'L') + X(opt.grid[i]).toFixed(1) + ',' + Y(b[qhi][i]).toFixed(1);
      for (let i = i1; i >= i0; i--) d += 'L' + X(opt.grid[i]).toFixed(1) + ',' + Y(b[qlo][i]).toFixed(1);
      el('path', { d: d + 'Z', class: cls }, svg);
    };
    const line = (vals, i0, i1, cls) => { let d = ''; for (let i = i0; i <= i1; i++) d += (i === i0 ? 'M' : 'L') + X(opt.grid[i]).toFixed(1) + ',' + Y(vals[i]).toFixed(1); el('path', { d, class: cls }, svg); };
    const grid = opt.grid || [];
    const iAs = Math.max(0, grid.findIndex(d => d >= opt.asOfDay));
    const n = grid.length - 1;
    if (opt.prevBands) { band(opt.prevBands, 0, n, 'gl-band-prev', 'q05', 'q95'); line(opt.prevBands.q50, 0, n, 'gl-line-prev'); }
    if (opt.bands) {
      band(opt.bands, 0, n, 'gl-band90', 'q05', 'q95'); band(opt.bands, 0, n, 'gl-band50', 'q25', 'q75');
      if (iAs > 0) line(opt.bands.q50, 0, iAs, 'gl-line-fit');
      line(opt.bands.q50, iAs, n, 'gl-line-forecast');
    }
    // events
    const evY = m.t + 2;
    (opt.events || []).forEach(e => {
      const cls = 'gl-ev gl-ev-' + e.type;
      const xa = X(e.day_start), xb = X(Math.max(e.day_end === null || e.day_end === undefined ? e.day_start : e.day_end, e.day_start) + 1);
      const r = el('rect', { x: xa, y: evY, width: Math.max(3, xb - xa), height: 6, rx: 1, class: cls }, svg);
      const title = el('title', {}, r); title.textContent = t('event.' + e.type) + ' ' + U.dayLabel(e.day_start) + (e.day_end > e.day_start ? '–' + U.dayLabel(e.day_end) : '');
    });
    // today line
    el('line', { x1: X(opt.asOfDay), x2: X(opt.asOfDay), y1: m.t, y2: H - m.b, class: 'gl-today' }, svg);
    const tl = el('text', { x: X(opt.asOfDay) + 4, y: m.t + 22, class: 'gl-today-label' }, svg); tl.textContent = t('chart.today') + ' ' + U.dayLabel(opt.asOfDay);
    // observations
    (opt.obsFuture || []).forEach(o => { if (o.day >= x0 && o.day <= x1) el('circle', { cx: X(o.day), cy: Y(Math.max(o.value, floor)), r: 2.6, class: 'gl-obs-future' }, svg); });
    (opt.obs || []).forEach(o => { if (o.day >= x0 && o.day <= x1) { const c = el('circle', { cx: X(o.day), cy: Y(Math.max(o.value, floor)), r: 3, class: 'gl-obs' }, svg); const ti = el('title', {}, c); ti.textContent = U.dayLabel(o.day) + ': ' + o.value + ' ' + (opt.variable.unit || ''); } });
    // hover readout
    const hover = el('g', { class: 'gl-hover', style: 'display:none' }, svg);
    const hl = el('line', { y1: m.t, y2: H - m.b, class: 'gl-hover-line' }, hover);
    const hb = el('rect', { rx: 4, class: 'gl-hover-box', width: 150, height: 36 }, hover);
    const ht1 = el('text', { class: 'gl-hover-text' }, hover), ht2 = el('text', { class: 'gl-hover-text' }, hover);
    svg.addEventListener('mousemove', (ev) => {
      const pt = svg.createSVGPoint(); pt.x = ev.clientX; pt.y = ev.clientY; const p = pt.matrixTransform(svg.getScreenCTM().inverse());
      if (p.x < m.l || p.x > W - m.r) { hover.style.display = 'none'; return; }
      const d = Math.round(x0 + (p.x - m.l) / (W - m.l - m.r) * (x1 - x0)); const i = grid.indexOf(d);
      if (i < 0 || !opt.bands) { hover.style.display = 'none'; return; }
      hover.style.display = '';
      hl.setAttribute('x1', X(d)); hl.setAttribute('x2', X(d));
      const bx = Math.min(X(d) + 8, W - m.r - 152);
      hb.setAttribute('x', bx); hb.setAttribute('y', m.t + 30);
      ht1.setAttribute('x', bx + 6); ht1.setAttribute('y', m.t + 44); ht1.textContent = U.dayLabel(d) + ' · ' + t('chart.median') + ' ' + U.fmt(opt.bands.q50[i]);
      ht2.setAttribute('x', bx + 6); ht2.setAttribute('y', m.t + 58); ht2.textContent = '90%: ' + U.fmt(opt.bands.q05[i]) + ' – ' + U.fmt(opt.bands.q95[i]);
    });
    svg.addEventListener('mouseleave', () => { hover.style.display = 'none'; });
    return svg;
  };

  // params: [{ id, label, median, lo, hi, prev:{median,lo,hi}|null, range:[lo,hi], log, unit }]
  C.paramBars = function (container, params) {
    container.innerHTML = '';
    const rowH = 34, W = 420, H = params.length * rowH + 8, m = { l: 150, r: 78 };
    const svg = el('svg', { viewBox: `0 0 ${W} ${H}`, class: 'gl-chart gl-params' }, container);
    params.forEach((p, i) => {
      const y = 8 + i * rowH + rowH / 2;
      const [a, b] = p.range; const tr = p.log ? (v) => Math.log(Math.max(v, a)) : (v) => v;
      const X = (v) => m.l + (tr(U.clamp(v, a, b)) - tr(a)) / (tr(b) - tr(a)) * (W - m.l - m.r);
      const lab = el('text', { x: m.l - 8, y: y + 4, class: 'gl-param-label', 'text-anchor': 'end' }, svg); lab.textContent = p.label;
      el('line', { x1: m.l, x2: W - m.r, y1: y, y2: y, class: 'gl-param-track' }, svg);
      if (p.prev) { el('rect', { x: X(p.prev.lo), y: y - 9, width: Math.max(2, X(p.prev.hi) - X(p.prev.lo)), height: 18, rx: 3, class: 'gl-param-prev' }, svg); el('line', { x1: X(p.prev.median), x2: X(p.prev.median), y1: y - 11, y2: y + 11, class: 'gl-param-prev-med' }, svg); }
      el('rect', { x: X(p.lo), y: y - 6, width: Math.max(2, X(p.hi) - X(p.lo)), height: 12, rx: 3, class: 'gl-param-ci' }, svg);
      el('line', { x1: X(p.median), x2: X(p.median), y1: y - 9, y2: y + 9, class: 'gl-param-med' }, svg);
      const v = el('text', { x: W - m.r + 6, y: y + 4, class: 'gl-param-value' }, svg); v.textContent = p.fmt ? p.fmt(p.median) : U.fmt(p.median);
      if (p.help) { const g = el('g', { class: 'help-svg', 'data-help': p.help, transform: `translate(${W - 10},${y})`, role: 'button', tabindex: '0' }, svg); el('circle', { r: 6.5 }, g); const q = el('text', { y: 3.2, 'text-anchor': 'middle' }, g); q.textContent = '?'; const ti = el('title', {}, g); ti.textContent = t('help.open'); }
    });
    return svg;
  };

  // stacked phenotype bar: {quick, second_wave, [delayed,] persistent}
  C.phenotypeBar = function (container, probs, opt) {
    opt = opt || {};
    container.innerHTML = '';
    const W = 300, H = opt.height || 18;
    const svg = el('svg', { viewBox: `0 0 ${W} ${H}`, class: 'gl-chart gl-pheno', preserveAspectRatio: 'none' }, container);
    let x = 0;
    [['quick', probs.quick], ['second_wave', probs.second_wave], ['delayed', probs.delayed], ['persistent', probs.persistent]].forEach(([k, p]) => {
      const w = W * (p || 0);
      if (w > 0) { const r = el('rect', { x, y: 0, width: w, height: H, class: 'gl-ph-' + k }, svg); const ti = el('title', {}, r); ti.textContent = t('pheno.' + k) + ' ' + U.pct(p); if (w > 36) { const tx = el('text', { x: x + w / 2, y: H / 2 + 4, class: 'gl-ph-text', 'text-anchor': 'middle' }, svg); tx.textContent = U.pct(p); } }
      x += w;
    });
    return svg;
  };

  C.sparkline = function (container, values, opt) {
    opt = opt || {}; container.innerHTML = '';
    const W = opt.width || 120, H = opt.height || 28;
    const svg = el('svg', { viewBox: `0 0 ${W} ${H}`, class: 'gl-chart gl-spark' }, container);
    const v = values.filter(x => isFinite(x)); if (!v.length) return svg;
    const lo = Math.min(...v), hi = Math.max(...v);
    const X = (i) => 2 + i / Math.max(1, values.length - 1) * (W - 4), Y = (x) => 2 + (1 - (x - lo) / ((hi - lo) || 1)) * (H - 4);
    let d = ''; values.forEach((x, i) => { d += (i === 0 ? 'M' : 'L') + X(i).toFixed(1) + ',' + Y(x).toFixed(1); });
    el('path', { d, class: 'gl-spark-line' }, svg);
    return svg;
  };
})(typeof window !== 'undefined' ? window : global);
