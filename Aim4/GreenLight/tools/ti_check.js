const GL = require('./load_core'); const S = GL.schema, TI = GL.modules.tumor_immune;
const snap = { patient_id: 'X', asOfDay: -5, history: {}, mutations: { tested: false, variants: [] }, observations: [], events: [], ld: null, bridging: 'none' };
for (const rho of [0.45, 0.68, 1.0]) for (const T0 of [20, 80, 400]) {
  const r = TI.simulate({ rho, r: 0.03, T0, e0: 1 }, snap, S.defaultPlan(), {}, 120);
  const pk = Math.max(...r.E), pd = r.days[r.E.indexOf(pk)];
  console.log(`rho=${rho} T0=${T0}: peakE=${pk.toFixed(2)} at d${pd}  E(d7)=${r.E[r.days.indexOf(7)].toFixed(2)} E(d14)=${r.E[r.days.indexOf(14)].toFixed(2)} E(d28)=${r.E[r.days.indexOf(28)].toFixed(2)}  T30=${(100*r.T[r.days.indexOf(30)]/T0).toFixed(1)}% T90=${(100*r.T[r.days.indexOf(90)]/T0).toFixed(1)}%`);
}
