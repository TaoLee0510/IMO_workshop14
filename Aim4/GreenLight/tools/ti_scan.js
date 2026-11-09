const GL = require('./load_core'); const S = GL.schema, TI = GL.modules.tumor_immune, U = GL.util;
const snap = { patient_id: 'X', asOfDay: -5, history: {}, mutations: { tested: false, variants: [] }, observations: [], events: [], ld: null, bridging: 'none' };
const rows = [];
for (const Emax of [1.5, 6, 20]) for (const T0 of [15, 80, 400]) {
  const r = TI.simulate({ rho: 1.3, Emax, r: 0.03, T0, e0: 1 }, snap, S.defaultPlan(), {}, 120);
  const pk = Math.max(...r.E), pd = r.days[r.E.indexOf(pk)];
  rows.push(`Emax=${Emax} T0=${T0}: pk=${pk.toFixed(1)}@d${pd} E28=${r.E[r.days.indexOf(28)].toFixed(2)} T30=${(100*r.T[r.days.indexOf(30)]/T0).toFixed(0)}% T90=${(100*r.T[r.days.indexOf(90)]/T0).toFixed(0)}%`);
}
console.log(rows.join('\n'));
const prior = TI.priorParticles(snap, 400, 1);
let resp = 0, cr = 0; const peaks = [], pds = [];
for (const p of prior) { const r = TI.simulate(p, snap, S.defaultPlan(), {}, 120); const pk = Math.max(...r.E); peaks.push(pk); pds.push(r.days[r.E.indexOf(pk)]); if (r.T[r.days.indexOf(90)] / r.T0 < TI.POP.respThr) resp++; if (r.T[r.days.indexOf(90)] < TI.POP.crThr) cr++; }
console.log('prior P(resp)=' + (resp / 400).toFixed(2), 'P(cr)=' + (cr / 400).toFixed(2), 'peak median=' + U.median(peaks).toFixed(1), '[', U.quantile(peaks, 0.1).toFixed(1), U.quantile(peaks, 0.9).toFixed(1), '] peak day median=' + U.median(pds), '[', U.quantile(pds, 0.1), U.quantile(pds, 0.9), ']');
