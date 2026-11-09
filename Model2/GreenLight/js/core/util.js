/* GreenLight · core/util.js
 * Shared helpers: namespace, seeded RNG, stats, formatting.
 * Works in browser (window.GL) and node (global.GL) without a bundler.
 */
(function (root) {
  const GL = root.GL = root.GL || {};
  const U = GL.util = {};

  // ---- deterministic RNG (mulberry32) ------------------------------------
  U.rng = function (seed) {
    let a = (seed >>> 0) || 1;
    const next = function () {
      a |= 0; a = a + 0x6D2B79F5 | 0;
      let t = Math.imul(a ^ a >>> 15, 1 | a);
      t = t + Math.imul(t ^ t >>> 7, 61 | t) ^ t;
      return ((t ^ t >>> 14) >>> 0) / 4294967296;
    };
    const r = { next };
    r.uniform = (lo, hi) => lo + (hi - lo) * next();
    r.normal = function () {
      let u = 0, v = 0;
      while (u === 0) u = next();
      while (v === 0) v = next();
      return Math.sqrt(-2 * Math.log(u)) * Math.cos(2 * Math.PI * v);
    };
    r.lognormal = (medianV, sdlog) => medianV * Math.exp(sdlog * r.normal());
    r.pick = (arr) => arr[Math.floor(next() * arr.length)];
    return r;
  };
  U.hashString = function (s) {
    let h = 2166136261;
    for (let i = 0; i < s.length; i++) { h ^= s.charCodeAt(i); h = Math.imul(h, 16777619); }
    return h >>> 0;
  };

  // ---- stats ---------------------------------------------------------------
  U.clamp = (x, lo, hi) => Math.min(hi, Math.max(lo, x));
  U.sum = (a) => a.reduce((s, x) => s + x, 0);
  U.mean = (a) => a.length ? U.sum(a) / a.length : NaN;
  U.median = (a) => U.quantile(a, 0.5);
  U.quantile = function (arr, q) {
    if (!arr.length) return NaN;
    const a = arr.slice().sort((x, y) => x - y);
    const pos = (a.length - 1) * q, lo = Math.floor(pos), hi = Math.ceil(pos);
    return a[lo] + (a[hi] - a[lo]) * (pos - lo);
  };
  // weighted quantile: values v[], weights w[] (normalised or not)
  U.wquantile = function (v, w, q) {
    const idx = v.map((_, i) => i).sort((a, b) => v[a] - v[b]);
    const tot = U.sum(w);
    let acc = 0;
    for (const i of idx) { acc += w[i] / tot; if (acc >= q) return v[i]; }
    return v[idx[idx.length - 1]];
  };
  U.wmean = function (v, w) {
    const tot = U.sum(w); let s = 0;
    for (let i = 0; i < v.length; i++) s += v[i] * w[i];
    return s / tot;
  };
  U.normalise = function (w) { const t = U.sum(w); return w.map(x => x / t); };
  U.ess = function (w) { const n = U.normalise(w); return 1 / U.sum(n.map(x => x * x)); };
  U.logsumexp = function (a) { const m = Math.max(...a); return m + Math.log(U.sum(a.map(x => Math.exp(x - m)))); };
  U.sigmoid = (x) => 1 / (1 + Math.exp(-x));

  // ---- formatting ----------------------------------------------------------
  U.fmt = function (x, d) {
    if (x === null || x === undefined || Number.isNaN(x)) return '–';
    if (d === undefined) d = Math.abs(x) >= 100 ? 0 : Math.abs(x) >= 10 ? 1 : 2;
    return Number(x).toFixed(d);
  };
  U.pct = (p, d) => (p === null || p === undefined || Number.isNaN(p)) ? '–' : (100 * p).toFixed(d === undefined ? 0 : d) + '%';
  U.dayLabel = (d) => (d >= 0 ? 'D+' + d : 'D' + d);
  U.deepClone = (o) => JSON.parse(JSON.stringify(o));
  U.uid = (function () { let n = 0; return (p) => (p || 'id') + '_' + Date.now().toString(36) + '_' + (n++).toString(36); })();
  U.interp = function (xs, ys, x) {
    if (x <= xs[0]) return ys[0];
    if (x >= xs[xs.length - 1]) return ys[ys.length - 1];
    let lo = 0, hi = xs.length - 1;
    while (hi - lo > 1) { const m = (lo + hi) >> 1; if (xs[m] <= x) lo = m; else hi = m; }
    const t = (x - xs[lo]) / (xs[hi] - xs[lo]);
    return ys[lo] + t * (ys[hi] - ys[lo]);
  };
  U.escapeHtml = (s) => String(s).replace(/[&<>"']/g, c => ({ '&': '&amp;', '<': '&lt;', '>': '&gt;', '"': '&quot;', "'": '&#39;' }[c]));
})(typeof window !== 'undefined' ? window : global);
