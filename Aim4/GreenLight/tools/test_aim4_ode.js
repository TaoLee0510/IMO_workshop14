/* GreenLight · tools/test_aim4_ode.js
 * Consistency test: the JavaScript port of the Aim 4 parent model
 * (js/modules/aim4_core.js) against reference trajectories from the R/C++
 * integrator (tools/ref_trajectories.json, written by step G or
 * tools/build_aim4_data.R). Pass: max relative difference < 1% on ANC, PLT, AMC, ALC and on the mean
 * time and number of divisions since MPP entry of each (track)
 * (relative to max(|ref|, 1e-3)).
 * Usage: node tools/test_aim4_ode.js
 */
const path = require('path'), fs = require('fs');
require(path.join(__dirname, '..', 'js', 'modules', 'aim4_core.js'));
const C = global.GL.aim4;
const refs = JSON.parse(fs.readFileSync(path.join(__dirname, 'ref_trajectories.json'), 'utf8'));
let worst = 0, fail = false;
for (const r of refs) {
  const par = Float64Array.from(C.PAR.map(k => r.par[k]));
  const des = Object.assign({}, r.design, {
    gcsf: (r.design.gcsf || []).map(w => [w[0], w[1]]), tpo: (r.design.tpo || []).map(w => [w[0], w[1]]),
    cd34: (r.design.cd34 || []).map(c => ({ day: c.day, add: c.add }))
  });
  const t0 = Date.now();
  const sim = C.simulate(par, { design: des, struct: r.struct, drv: { t0: r.drv.t0, dt: r.drv.dt, v: Float64Array.from(r.drv.v) }, tEnd: Math.max(...r.days), track: true });
  const ms = Date.now() - t0;
  const row = {};
  for (const v of ['anc_rel', 'plt_rel', 'amc_rel', 'alc_rel', 'S', 'f_clone', 'age_anc', 'age_plt', 'age_amc', 'age_alc', 'gen_anc', 'gen_plt', 'gen_amc', 'gen_alc']) {
    if (!r.out[v]) continue;
    let m = 0;
    r.days.forEach((d, i) => {
      const j = sim.days.indexOf(d); if (j < 0) return;
      const a = sim[v][j], b = r.out[v][i];
      m = Math.max(m, Math.abs(a - b) / Math.max(Math.abs(b), 1e-3));
    });
    row[v] = m;
    if (v !== 'S' && v !== 'f_clone') { worst = Math.max(worst, m); if (m >= 0.01) fail = true; }
  }
  console.log(r.name.padEnd(38), Object.entries(row).map(([k, x]) => k + ' ' + (100 * x).toExponential(1) + '%').join(' '), ` (${ms} ms)`);
}
console.log((fail ? 'FAIL' : 'PASS') + ': max relative difference on ANC/PLT/AMC/ALC and their age / generation since MPP entry = ' + (100 * worst).toExponential(2) + '%');
process.exit(fail ? 1 : 0);
