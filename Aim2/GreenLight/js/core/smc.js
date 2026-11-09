/* GreenLight · core/smc.js
 * Shared particle-posterior helpers used by local modules:
 * tempering to a target effective sample size, systematic resampling,
 * jitter (rejuvenation), weighted summaries and quantile bands.
 */
(function (root) {
  const GL = root.GL = root.GL || {};
  const U = GL.util;
  const SMC = GL.smc = {};

  SMC.temper = function (logw, target) {
    const w = (b) => { const m = Math.max(...logw); return logw.map(x => Math.exp(b * (x - m))); };
    if (U.ess(w(1)) >= target) return { beta: 1, w: w(1) };
    let lo = 0, hi = 1;
    for (let k = 0; k < 30; k++) { const mid = (lo + hi) / 2; if (U.ess(w(mid)) >= target) lo = mid; else hi = mid; }
    return { beta: lo, w: w(lo) };
  };
  SMC.resample = function (parts, w, m, seed) {
    const r = U.rng(seed), cw = [], tot = U.sum(w); let acc = 0;
    for (const x of w) { acc += x / tot; cw.push(acc); }
    const out = []; let j = 0; const u0 = r.next() / m;
    for (let i = 0; i < m; i++) { const u = u0 + i / m; while (cw[j] < u && j < cw.length - 1) j++; out.push(parts[j]); }
    return out;
  };
  // jitter: spec = { key: { log: bool, sd: number, lo, hi } }  (sd on log scale when log)
  SMC.jitter = function (parts, spec, seed) {
    const r = U.rng(seed);
    return parts.map(p => {
      const q = { ...p };
      for (const k in spec) {
        const s = spec[k];
        if (s.log) q[k] = U.clamp(p[k] * Math.exp(s.sd * r.normal()), s.lo, s.hi);
        else q[k] = U.clamp(p[k] + s.sd * r.normal(), s.lo, s.hi);
      }
      return q;
    });
  };
  SMC.summarise = function (parts, w, key) {
    const v = parts.map(p => p[key]);
    return { median: U.wquantile(v, w, 0.5), lo: U.wquantile(v, w, 0.05), hi: U.wquantile(v, w, 0.95), mean: U.wmean(v, w) };
  };
  // sort once, read five quantiles
  function pushQ(out, col) {
    col.sort((a, b) => a - b);
    const n = col.length, at = (q) => { const pos = (n - 1) * q, lo = Math.floor(pos), hi = Math.ceil(pos); return col[lo] + (col[hi] - col[lo]) * (pos - lo); };
    out.q05.push(at(0.05)); out.q25.push(at(0.25)); out.q50.push(at(0.5)); out.q75.push(at(0.75)); out.q95.push(at(0.95));
  }
  SMC.bands = function (grid, arrs) {
    const out = { q05: [], q25: [], q50: [], q75: [], q95: [] };
    for (let i = 0; i < grid.length; i++) pushQ(out, arrs.map(a => a[i]));
    return out;
  };
  // predictive bands: particle spread convolved with log-normal measurement noise (sd on log scale)
  const Z = [-1.645, -0.675, 0, 0.675, 1.645];
  SMC.predictiveBands = function (grid, arrs, sd) {
    if (!sd) return SMC.bands(grid, arrs);
    const mult = Z.map(z => Math.exp(sd * z));
    const out = { q05: [], q25: [], q50: [], q75: [], q95: [] };
    for (let i = 0; i < grid.length; i++) {
      const col = [];
      for (const a of arrs) for (const m of mult) col.push(a[i] * m);
      pushQ(out, col);
    }
    return out;
  };
  SMC.val = (a) => ({ value: U.median(a), lo: U.quantile(a, 0.1), hi: U.quantile(a, 0.9), mean: U.mean(a) });

  // Generic two-stage fit: prior -> weight -> resample -> jitter -> re-weight.
  // simulateLL(particle) -> {ll, n}
  SMC.fit = function (prior, simulateLL, jitterSpec, seed, opts) {
    opts = opts || {};
    const n = prior.length, target = Math.max(opts.minEss || 40, (opts.essFrac || 0.1) * n);
    const ll1 = new Array(n); let nObs = 0;
    for (let i = 0; i < n; i++) { const r = simulateLL(prior[i]); ll1[i] = r.ll; nObs = r.n; }
    if (nObs === 0) return { particles: prior, weights: new Array(n).fill(1 / n), ess: n, beta: 1, nObs: 0, logEvidence: 0, stages: 1 };
    const t1 = SMC.temper(ll1, target);
    let parts = SMC.resample(prior, t1.w, n, seed + 101);
    parts = SMC.jitter(parts, jitterSpec, seed + 202);
    const ll2 = new Array(n);
    for (let i = 0; i < n; i++) ll2[i] = simulateLL(parts[i]).ll;
    const t2 = SMC.temper(ll2, target);
    const w = U.normalise(t2.w);
    return { particles: parts, weights: w, ess: U.ess(w), beta: t2.beta, nObs, logEvidence: U.logsumexp(ll1) - Math.log(n), stages: 2 };
  };
})(typeof window !== 'undefined' ? window : global);
