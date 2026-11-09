const GL = require('./load_core');
const S = GL.schema, M = GL.modules.marrow;
const snap = { patient_id: 'X', asOfDay: -5, history: {}, costim: 'CD28', mutations: { tested: false, variants: [] }, observations: [], events: [], ld: null };
const plan = S.defaultPlan();
plan.stages.cytopenia.gcsf = process.argv[2] || 'none';
const fmt = (x) => x.toFixed(2).padStart(5);
const days = [-5, 0, 3, 7, 10, 14, 21, 28, 35, 42, 50, 60, 75, 90];
console.log('plan gcsf=', plan.stages.cytopenia.gcsf);
for (const f0 of [0, 0.3]) for (const R of [0.9, 0.5, 0.25]) for (const g of [30, 120, 300]) for (const q of [1/200, 1/80, 1/35]) {
  const p = { R, g, q, f0, tau: 1, N0: 3.6, PLT0: 230 };
  const r = M.simulate(p, snap, plan, {}, 120);
  const ph = S.classifyPhenotype(r.days, r.anc);
  const idx = days.map(d => r.days.indexOf(d));
  console.log(`f0=${f0} R=${R} g=${g.toString().padStart(3)} q=1/${(1/q).toFixed(0).padStart(3)} ${ph.padEnd(12)} ANC:` + idx.map(i => fmt(r.anc[i])).join('') + '  CRPmax=' + Math.max(...r.crp).toFixed(0) + ' PLT60=' + r.plt[r.days.indexOf(60)].toFixed(0) + ' gcsf=' + r.pol.gcsfCourses);
}
