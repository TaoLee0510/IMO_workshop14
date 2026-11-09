/* GreenLight · tools/test_model3_ode.js
 * Consistency test: the JavaScript port of Model3 (js/engine/model3_core.js) against
 * reference trajectories from the R/C++ engine (tools/ref_model3.json, written by step G,
 * Model3/Code/R/an_G_export.R): daily grid, digital-patient designs incl. bridging, G-CSF,
 * TPO-RA and a CD34+ boost; the Model2 cases plus cases that exercise the Model3 terms (neutrophil
 * chains n_tr_n / nk_n 8/6, 5/1 and the Model2 chain 3/2; P1 th_cp > 0; P2 d_zpl != 0).
 * Pass: (1) C.buildPar reproduces R build_par() (incl. ktr_n = n_tr_n / MTT_N and
 * e_z_pl = plogis(qlogis(e_z) + d_zpl)) to 1e-12 relative; (2) max relative difference < 1e-8 on
 * ANC, PLT, AMC, ALC (relative to max(|ref|, 1e-3)); S, clone fraction, CAR T and Z are reported.
 * Usage: node tools/test_model3_ode.js
 */
const path = require('path'), fs = require('fs');
require(path.join(__dirname, '..', 'js', 'engine', 'model3_core.js'));
const C = global.GL.model3;
const refs = JSON.parse(fs.readFileSync(path.join(__dirname, 'ref_model3.json'), 'utf8'));
const TOL = 1e-8, TOL_PAR = 1e-12;
const one = (x) => Array.isArray(x) ? x[0] : x;
let worst = 0, worstPar = 0, fail = false;
for (const r of refs) {
  const par = Float64Array.from(C.PAR.map(k => one(r.par[k])));
  if (C.PAR.some(k => r.par[k] === undefined) || Object.keys(r.par).length !== C.PAR.length) { console.log('FAIL parameter layout differs from R build_par: ' + r.name); fail = true; }
  const s = {}; for (const k in r.struct) s[k] = one(r.struct[k]);
  // (1) build_par port: individual values taken from the R vector, struct and model id as exported
  const ind = {}; ['r_s', 'f0', 'th_as', 'th_ap', 'th_bs', 'rho', 'a0'].forEach(k => { ind[k] = one(r.par[k]); });
  const bp = C.buildPar(ind, s, C.spec(r.model));
  let dp = 0; C.PAR.forEach((k, i) => { dp = Math.max(dp, Math.abs(bp[i] - par[i]) / Math.max(Math.abs(par[i]), 1e-12)); });
  worstPar = Math.max(worstPar, dp); if (!(dp < TOL_PAR)) fail = true;
  const des = Object.assign({}, r.design, {
    gcsf: (r.design.gcsf || []).map(w => [w[0], w[1]]), tpo: (r.design.tpo || []).map(w => [w[0], w[1]]),
    cd34: (r.design.cd34 || []).map(c => ({ day: [].concat(c.day)[0], add: [].concat(c.add)[0] }))
  });
  const t0 = Date.now();
  const sim = C.simulate(par, { design: des, struct: s, drv: { t0: r.drv.t0, dt: r.drv.dt, v: Float64Array.from(r.drv.v) }, tEnd: Math.max(...r.days) });
  const ms = Date.now() - t0;
  const map = { anc_rel: 'anc_rel', plt_rel: 'plt_rel', amc_rel: 'amc_rel', alc_rel: 'alc_rel', S: 'S', f_clone: 'f_clone', cart: 'cart', z: 'z' };
  const row = {};
  for (const [v, w] of Object.entries(map)) {
    let m = 0;
    r.days.forEach((d, i) => {
      const j = sim ? sim.days.indexOf(d) : -1; if (j < 0) { m = Infinity; return; }
      const a = sim[w][j], b = r.out[v][i];
      m = Math.max(m, Math.abs(a - b) / Math.max(Math.abs(b), 1e-3));
    });
    row[v] = m;
    if (['anc_rel', 'plt_rel', 'amc_rel', 'alc_rel'].includes(v)) { worst = Math.max(worst, m); if (!(m < TOL)) fail = true; }
  }
  console.log(r.name.padEnd(46), `N ${one(r.par.n_tr_n)}/${one(r.par.nk_n)} par ${dp.toExponential(0)}`.padEnd(16),
    Object.entries(row).map(([k, x]) => k + ' ' + x.toExponential(1)).join('  '), ` (${ms} ms, ${r.days.length} days)`);
}
console.log((fail ? 'FAIL' : 'PASS') + `: ${refs.length} cases; build_par max relative difference ${worstPar.toExponential(2)} (tolerance ${TOL_PAR}); ` +
  'max relative difference on ANC/PLT/AMC/ALC = ' + worst.toExponential(2) + ' (tolerance ' + TOL + ')');
process.exit(fail ? 1 : 0);
