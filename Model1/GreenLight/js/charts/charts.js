/* GreenLight · charts/charts.js
 * Small SVG charts with fixed axes (the clinician view keeps every axis fixed per patient,
 * so scrubbing and playback never rescale): lineage forecast (band + median + observations),
 * prediction evolution (stacked P(phenotype) by as-of day), ferritin driver, heat map, scatter.
 */
(function (root) {
  const GL = root.GL = root.GL || {};
  const Ch = GL.charts = {};
  const f = (x) => (Math.round(x * 10) / 10).toString();
  Ch.scale = function (d0, d1, r0, r1, log) {
    if (log) { const l0 = Math.log(d0), l1 = Math.log(d1); return (v) => r0 + (Math.log(Math.max(v, d0 * 1e-3)) - l0) / (l1 - l0) * (r1 - r0); }
    return (v) => r0 + (v - d0) / (d1 - d0) * (r1 - r0);
  };
  Ch.logTicks = function (lo, hi) {
    const out = [];
    for (let e = Math.floor(Math.log10(lo)); e <= Math.ceil(Math.log10(hi)); e++) [1, 2, 5].forEach(m => { const v = m * Math.pow(10, e); if (v >= lo * 0.999 && v <= hi * 1.001) out.push(v); });
    return out.length > 6 ? out.filter((v, i) => /^1/.test(String(v)) || i % 2 === 0) : out;
  };
  const fmtTick = (v) => v >= 1000 ? (v / 1000) + 'k' : v >= 1 ? String(Math.round(v * 100) / 100) : String(v);
  const path = (xs, ys, X, Y) => { let s = ''; xs.forEach((x, i) => { if (ys[i] === null || ys[i] === undefined || !isFinite(ys[i])) return; s += (s ? 'L' : 'M') + f(X(x)) + ',' + f(Y(ys[i])); }); return s; };
  const area = (xs, lo, hi, X, Y) => {
    let s = '';
    xs.forEach((x, i) => { s += (i ? 'L' : 'M') + f(X(x)) + ',' + f(Y(hi[i])); });
    for (let i = xs.length - 1; i >= 0; i--) s += 'L' + f(X(xs[i])) + ',' + f(Y(lo[i]));
    return s + 'Z';
  };
  // piecewise time axis shared by the charts and the scrubber: days -10..30 take 64% of the width
  Ch.timeScale = function (r0, r1, d0, d1, split) {
    d0 = d0 === undefined ? -10 : d0; d1 = d1 === undefined ? 90 : d1; split = split === undefined ? 30 : split;
    const w = r1 - r0, k = 0.64;
    return (d) => d <= split ? r0 + (d - d0) / (split - d0) * w * k : r0 + w * k + (d - split) / (d1 - split) * w * (1 - k);
  };

  // lineage chart: o = { w, h, yDomain:[lo,hi], log, grid, band:{lo,med,hi}, obs, asOf, horizon, color, thr:[{v,label}],
  //                      whatif:{days,v}, truth:{days,v}, gcsf:[[s,e]], future:'hollow'|'hide' }
  Ch.lineage = function (o) {
    const m = { l: 34, r: 8, t: 6, b: 18 }, W = o.w, H = o.h;
    const X = Ch.timeScale(m.l, W - m.r), Y = Ch.scale(o.yDomain[0], o.yDomain[1], H - m.b, m.t, o.log);
    const yt = o.log ? Ch.logTicks(o.yDomain[0], o.yDomain[1]) : [0, 0.25, 0.5, 0.75, 1].map(q => o.yDomain[0] + q * (o.yDomain[1] - o.yDomain[0]));
    const xt = [-10, 0, 7, 14, 21, 30, 60, 90];
    let s = `<rect class="past" x="${m.l}" y="${m.t}" width="${f(Math.max(0, X(Math.min(o.asOf, 90)) - m.l))}" height="${H - m.t - m.b}"/>`;
    s += '<g class="grid">' + yt.map(v => `<line x1="${m.l}" x2="${W - m.r}" y1="${f(Y(v))}" y2="${f(Y(v))}"/>`).join('') + '</g>';
    s += '<g class="ax">' + yt.map(v => `<text x="${m.l - 4}" y="${f(Y(v) + 3)}" text-anchor="end">${fmtTick(v)}</text>`).join('') +
      xt.map(d => `<text x="${f(X(d))}" y="${H - 5}" text-anchor="${d === 90 ? 'end' : 'middle'}">${d === 0 ? 'D0' : (d > 0 ? '+' : '') + d}</text>`).join('') +
      `<line x1="${m.l}" x2="${W - m.r}" y1="${H - m.b}" y2="${H - m.b}"/></g>`;
    (o.thr || []).forEach(t => { if (t.v > o.yDomain[0] && t.v < o.yDomain[1]) s += `<line class="thr" x1="${m.l}" x2="${W - m.r}" y1="${f(Y(t.v))}" y2="${f(Y(t.v))}"/><text class="thr-lbl" x="${W - m.r - 2}" y="${f(Y(t.v) - 2)}" text-anchor="end">${t.label}</text>`; });
    (o.gcsf || []).forEach(w => { s += `<rect class="gcsf-bar" x="${f(X(w[0]))}" y="${H - m.b - 4}" width="${f(Math.max(2, X(w[1]) - X(w[0])))}" height="4"><title>G-CSF</title></rect>`; });
    if (o.band) {
      s += `<path class="band" fill="${o.color}" d="${area(o.grid, o.band.lo, o.band.hi, X, Y)}"/>`;
      s += `<path class="med" stroke="${o.color}" d="${path(o.grid, o.band.med, X, Y)}"/>`;
    }
    if (o.truth) s += `<path class="truth" d="${path(o.truth.days, o.truth.v, X, Y)}"/>`;
    if (o.whatif) s += `<path class="whatif" stroke="${o.color}" d="${path(o.whatif.days, o.whatif.v, X, Y)}"/>`;
    (o.obs || []).forEach(p => {
      const d = p[0], v = p[1], cens = p[2] === 1, keep = p[3] === 1;
      if (d < -10 || d > 90) return;
      const future = d > o.asOf;
      if (future && (o.future === 'hide' || d > o.horizon)) return;
      const cx = f(X(d)), cy = f(Y(Math.max(v, o.yDomain[0]))), cls = 'obs' + (future ? ' obs-future' : '');
      s += cens ? `<path class="${cls}" fill="${o.color}" stroke="${o.color}" d="M${cx - 3.5},${cy - 3}L${+cx + 3.5},${cy - 3}L${cx},${+cy + 3.5}Z"/>`
        : `<circle class="${cls}" cx="${cx}" cy="${cy}" r="${keep ? 2.6 : 2}" fill="${keep ? o.color : 'none'}" stroke="${o.color}"/>`;
    });
    s += `<line class="asof-line" x1="${f(X(o.asOf))}" x2="${f(X(o.asOf))}" y1="${m.t}" y2="${H - m.b}"/>`;
    return `<svg viewBox="0 0 ${W} ${H}" preserveAspectRatio="none" role="img">${s}</svg>`;
  };

  // prediction evolution: stacked P(quick, second wave, delayed) by as-of day
  Ch.evolution = function (o) {
    const m = { l: 34, r: 8, t: 6, b: 16 }, W = o.w, H = o.h;
    const X = Ch.timeScale(m.l, W - m.r, -10, 90), Y = Ch.scale(0, 1, H - m.b, m.t);
    const ev = o.evolution.filter(e => e[0] <= o.horizon);
    if (!ev.length) return `<svg viewBox="0 0 ${W} ${H}"></svg>`;
    const xs = ev.map(e => e[0]);
    const q = ev.map(e => e[1]), qs = ev.map(e => e[1] + e[2]);
    let s = '';
    const band = (lo, hi, cls) => `<path class="${cls}" d="${area(xs, lo, hi, X, Y)}" opacity=".85"/>`;
    s += band(xs.map(() => 0), q, 'ph-q-f') + band(q, qs, 'ph-s-f') + band(qs, xs.map(() => 1), 'ph-d-f');
    s = s.replace(/class="ph-q-f"/g, 'fill="var(--ph-quick)"').replace(/class="ph-s-f"/g, 'fill="var(--ph-sw)"').replace(/class="ph-d-f"/g, 'fill="var(--ph-del)"');
    s += '<g class="ax">' + [0, 0.5, 1].map(v => `<text x="${m.l - 4}" y="${f(Y(v) + 3)}" text-anchor="end">${Math.round(v * 100)}%</text>`).join('') +
      [-5, 0, 7, 14, 21, 30, 60, 90].map(d => `<text x="${f(X(d))}" y="${H - 3}" text-anchor="${d === 90 ? 'end' : 'middle'}">${d === 0 ? 'D0' : (d > 0 ? '+' : '') + d}</text>`).join('') + '</g>';
    if (o.label) s += `<text x="${W - m.r}" y="${m.t + 10}" text-anchor="end" class="thr-lbl" style="font-weight:700;fill:var(--text)">${o.label}</text>`;
    s += `<line class="asof-line" x1="${f(X(o.asOf))}" x2="${f(X(o.asOf))}" y1="${m.t}" y2="${H - m.b}" style="stroke-width:2;opacity:.9"/>`;
    return `<svg viewBox="0 0 ${W} ${H}" preserveAspectRatio="none" role="img">${s}</svg>`;
  };

  // ferritin driver strip: observed log(ferritin / LD value) and the as-of driver (observed + kernel)
  Ch.ferritin = function (o) {
    const m = { l: 34, r: 8, t: 4, b: 14 }, W = o.w, H = o.h;
    const X = Ch.timeScale(m.l, W - m.r, -10, 90), Y = Ch.scale(-0.2, o.yMax, H - m.b, m.t);
    let s = `<rect class="past" x="${m.l}" y="${m.t}" width="${f(Math.max(0, X(Math.min(o.asOf, 90)) - m.l))}" height="${H - m.t - m.b}"/>`;
    s += `<line class="thr" x1="${m.l}" x2="${W - m.r}" y1="${f(Y(0))}" y2="${f(Y(0))}"/>`;
    if (o.drv) { const xs = o.drv.v.map((_, i) => o.drv.t0 + i * o.drv.dt); const sel = xs.map(x => x >= -10 && x <= 90);
      s += `<path class="med" stroke="var(--amber)" d="${path(xs.filter((_, i) => sel[i]), o.drv.v.filter((_, i) => sel[i]), X, Y)}"/>`; }
    (o.ferritin || []).forEach(p => { if (p[0] < -10 || p[0] > 90 || p[0] > o.horizon) return; const fut = p[0] > o.asOf;
      s += `<circle cx="${f(X(p[0]))}" cy="${f(Y(Math.log(p[1] / o.base)))}" r="2.4" fill="${fut ? 'none' : 'var(--amber)'}" stroke="var(--amber)" opacity="${fut ? .5 : 1}"/>`; });
    if (o.cart) { const sel = o.cart.days.map((d, i) => d >= -10 && d <= 90); const mx = Math.max(...o.cart.v, 1e-9);
      const Yc = Ch.scale(0, mx, H - m.b, m.t);
      s += `<path class="med" stroke="var(--l-alc)" style="stroke-dasharray:4 3;stroke-width:1.4" d="${path(o.cart.days.filter((_, i) => sel[i]), o.cart.v.filter((_, i) => sel[i]), X, Yc)}"/>`; }
    s += `<g class="ax"><text x="${m.l - 4}" y="${f(Y(0) + 3)}" text-anchor="end">0</text><text x="${m.l - 4}" y="${m.t + 8}" text-anchor="end">${fmtTick(Math.round(o.yMax * 10) / 10)}</text></g>`;
    s += `<line class="asof-line" x1="${f(X(o.asOf))}" x2="${f(X(o.asOf))}" y1="${m.t}" y2="${H - m.b}"/>`;
    return `<svg viewBox="0 0 ${W} ${H}" preserveAspectRatio="none" role="img">${s}</svg>`;
  };

  // heat map: rows (id, label colour) x days with presence
  Ch.coverage = function (o) {
    const rows = o.rows, d0 = -30, d1 = 90, cw = 6, rh = 9, lw = 48, W = lw + (d1 - d0 + 1) * cw + 8, H = rows.length * rh + 22;
    let s = '';
    rows.forEach((r, i) => {
      s += `<text x="${lw - 14}" y="${i * rh + 8}" text-anchor="end" style="font-size:8px;fill:var(--muted)">${r.id}</text><rect x="${lw - 10}" y="${i * rh + 1}" width="7" height="7" rx="2" fill="${r.color}"/>`;
      const set = new Set(r.days);
      for (let d = d0; d <= d1; d++) if (set.has(d)) s += `<rect x="${lw + (d - d0) * cw}" y="${i * rh + 1}" width="${cw - 1}" height="${rh - 2}" fill="var(--l-anc)" opacity=".75"/>`;
    });
    [-30, -5, 0, 14, 30, 60, 90].forEach(d => { s += `<text x="${lw + (d - d0) * cw + 2}" y="${rows.length * rh + 14}" text-anchor="middle" style="font-size:9px;fill:var(--muted)">${d}</text>`; });
    s += `<line x1="${lw + (0 - d0) * cw}" x2="${lw + (0 - d0) * cw}" y1="0" y2="${rows.length * rh}" stroke="var(--red)" stroke-width="1"/>`;
    return `<svg class="heat" viewBox="0 0 ${W} ${H}" role="img">${s}</svg>`;
  };

  // scatter: points [{x, y, label}], log axes optional
  Ch.scatter = function (o) {
    const m = { l: 40, r: 10, t: 10, b: 30 }, W = o.w || 360, H = o.h || 260;
    const X = Ch.scale(o.dom[0], o.dom[1], m.l, W - m.r, o.log), Y = Ch.scale(o.dom[0], o.dom[1], H - m.b, m.t, o.log);
    const t = o.log ? Ch.logTicks(o.dom[0], o.dom[1]) : [o.dom[0], (o.dom[0] + o.dom[1]) / 2, o.dom[1]];
    let s = `<line x1="${X(o.dom[0])}" y1="${Y(o.dom[0])}" x2="${X(o.dom[1])}" y2="${Y(o.dom[1])}" stroke="var(--faint)" stroke-dasharray="4 3"/>`;
    s += '<g class="ax">' + t.map(v => `<text x="${f(X(v))}" y="${H - m.b + 14}" text-anchor="middle">${fmtTick(v)}</text><text x="${m.l - 4}" y="${f(Y(v) + 3)}" text-anchor="end">${fmtTick(v)}</text>`).join('') +
      `<text x="${(m.l + W - m.r) / 2}" y="${H - 4}" text-anchor="middle">${o.xLabel || ''}</text><text x="10" y="${(m.t + H - m.b) / 2}" text-anchor="middle" transform="rotate(-90 10 ${(m.t + H - m.b) / 2})">${o.yLabel || ''}</text></g>`;
    o.pts.forEach(p => { if (!(p.x > 0 || !o.log) || !(p.y > 0 || !o.log)) return; s += `<circle cx="${f(X(p.x))}" cy="${f(Y(p.y))}" r="4" fill="${p.color || 'var(--dp)'}" opacity=".85"><title>${p.label}</title></circle>`; });
    return `<svg viewBox="0 0 ${W} ${H}" style="width:100%;height:auto" role="img">${s}</svg>`;
  };
})(typeof window !== 'undefined' ? window : global);
