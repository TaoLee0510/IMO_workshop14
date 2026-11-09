const GL = require('./load_core'); const S = GL.schema, TI = GL.modules.tumor_immune, U = GL.util;
const snap = { patient_id: 'X', asOfDay: -5, history: {}, mutations: { tested: false, variants: [] }, observations: [], events: [], ld: null, bridging: 'none' };
const POP = TI.POP;
for (const rho of [1.6, 2.0]) for (const tauE of [7, 9]) for (const KT of [3, 6]) for (const mu of [0.12, 0.15]) {
  Object.assign(POP, { tauE, KT, mu });
  const typ = TI.simulate({ rho, Emax: 6, r: 0.03, T0: 80, e0: 1 }, snap, S.defaultPlan(), {}, 120);
  const pk = Math.max(...typ.E), pd = typ.days[typ.E.indexOf(pk)];
  const small = TI.simulate({ rho, Emax: 6, r: 0.03, T0: 15, e0: 1 }, snap, S.defaultPlan(), {}, 120);
  const pk2 = Math.max(...small.E), pd2 = small.days[small.E.indexOf(pk2)];
  const prior = TI.priorParticles(snap, 400, 1);
  let resp = 0; const peaks = [], pds = [];
  for (const p of prior) { p.rho = p.rho / 1.3 * rho; const r = TI.simulate(p, snap, S.defaultPlan(), {}, 120); const pkk = Math.max(...r.E); peaks.push(pkk); pds.push(r.days[r.E.indexOf(pkk)]); if (r.T[r.days.indexOf(90)] / r.T0 < POP.respThr) resp++; }
  console.log(`rho=${rho} tauE=${tauE} KT=${KT} mu=${mu}: typical pk=${pk.toFixed(1)}@d${pd} E28=${typ.E[typ.days.indexOf(28)].toFixed(2)} T90=${(100*typ.T[typ.days.indexOf(90)]/80).toFixed(0)}% | small-tumour pk=${pk2.toFixed(1)}@d${pd2} T90=${(100*small.T[small.days.indexOf(90)]/15).toFixed(0)}% | prior P(resp)=${(resp/400).toFixed(2)} peak med=${U.median(peaks).toFixed(1)} day med=${U.median(pds)} [${U.quantile(pds,0.1)}-${U.quantile(pds,0.9)}]`);
}
