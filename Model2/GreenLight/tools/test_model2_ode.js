/* GreenLight · tools/test_model2_ode.js
 * Consistency test: the JavaScript port of Model2 (js/engine/model2_core.js) against
 * reference trajectories from the R/C++ engine (tools/ref_model2.json, written by step G,
 * Model2/Code/R/an_G_export.R): daily grid, real patient designs incl. bridging, G-CSF,
 * TPO-RA and a CD34+ boost. Pass: max relative difference < 1e-8 on ANC, PLT, AMC, ALC
 * (relative to max(|ref|, 1e-3)); S, clone fraction, CAR T and Z are reported.
 * Usage: node tools/test_model2_ode.js
 */
const path = require('path'), fs = require('fs');
require(path.join(__dirname, '..', 'js', 'engine', 'model2_core.js'));
const C = global.GL.model2;
const refs = JSON.parse(fs.readFileSync(path.join(__dirname, 'ref_model2.json'), 'utf8'));
const TOL = 1e-8;
let worst = 0, fail = false;
for (const r of refs) {
  const par = Float64Array.from(C.PAR.map(k => Array.isArray(r.par[k]) ? r.par[k][0] : r.par[k]));
  const s = {}; for (const k in r.struct) s[k] = Array.isArray(r.struct[k]) ? r.struct[k][0] : r.struct[k];
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
      const j = sim.days.indexOf(d); if (j < 0) { m = Infinity; return; }
      const a = sim[w][j], b = r.out[v][i];
      m = Math.max(m, Math.abs(a - b) / Math.max(Math.abs(b), 1e-3));
    });
    row[v] = m;
    if (['anc_rel', 'plt_rel', 'amc_rel', 'alc_rel'].includes(v)) { worst = Math.max(worst, m); if (!(m < TOL)) fail = true; }
  }
  console.log(r.name.padEnd(38), Object.entries(row).map(([k, x]) => k + ' ' + x.toExponential(1)).join('  '), ` (${ms} ms, ${r.days.length} days)`);
}
console.log((fail ? 'FAIL' : 'PASS') + ': max relative difference on ANC/PLT/AMC/ALC = ' + worst.toExponential(2) + ' (tolerance ' + TOL + ')');
process.exit(fail ? 1 : 0);
