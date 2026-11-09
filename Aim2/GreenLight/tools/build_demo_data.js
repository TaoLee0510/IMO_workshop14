/* Build GreenLight demo data from Aim2/TestData (synthetic cohort, same contract as real Data/).
 * Usage: node tools/build_demo_data.js [path/to/TestData]
 * Writes data/demo_patients.js (browser, file:// friendly) and data/demo_patients.json.
 * Tumour–immune observations (car_copies, tumor_burden) do not exist in TestData; they are
 * generated here from the placeholder module with per-patient hidden parameters and are
 * labelled synthetic in the data file.
 */
const fs = require('fs'), path = require('path');
const GL = require('./load_core');
const S = GL.schema, TI = GL.modules.tumor_immune, U = GL.util;
const src = process.argv[2] || path.join(__dirname, '..', '..', 'TestData');
const csv = (f) => {
  const lines = fs.readFileSync(path.join(src, f), 'utf8').trim().split(/\r?\n/);
  const h = lines[0].split(',');
  return lines.slice(1).map(l => { const c = l.split(','); const o = {}; h.forEach((k, i) => { o[k] = c[i] === undefined ? '' : c[i]; }); return o; });
};
const num = (x) => (x === '' || x === undefined) ? null : Number(x);
const patients = csv('patients.csv'), labs = csv('labs_long.csv'), inflam = csv('inflam_long.csv'), events = csv('events.csv'), tested = csv('mutation_tested.csv'), muts = csv('mutations_pre.csv');

// status 'active': in treatment, data arrive batch by batch (dynamic prediction demo)
// status 'completed': full record available (training / validation cohort), viewable at any time point
const DEMO = [
  { id: 'P148', status: 'active', initial: 0, tag: 'persistent_tp53', ti: { rho: 2.0, Emax: 3.0, r: 0.04, T0: 90, e0: 1 } },
  { id: 'P083', status: 'active', initial: 2, tag: 'second_wave_tp53_low', ti: { rho: 1.9, Emax: 1.8, r: 0.05, T0: 200, e0: 1 } },
  { id: 'P002', status: 'active', initial: 3, tag: 'quick', ti: { rho: 2.0, Emax: 5.0, r: 0.03, T0: 40, e0: 1 } },
  { id: 'P001', status: 'active', initial: 4, tag: 'second_wave_crs4', ti: { rho: 2.2, Emax: 9.0, r: 0.03, T0: 150, e0: 1 } },
  { id: 'P004', status: 'active', initial: 5, tag: 'quick_reduced_ld', ti: { rho: 2.0, Emax: 4.0, r: 0.03, T0: 30, e0: 1 } },
  { id: 'P010', status: 'active', initial: 6, tag: 'second_wave_crs4_reduced_ld', ti: { rho: 2.3, Emax: 12.0, r: 0.03, T0: 300, e0: 1 } },
  { id: 'P026', status: 'active', initial: 7, tag: 'persistent_crs4_cd34', ti: { rho: 2.1, Emax: 8.0, r: 0.03, T0: 120, e0: 1 } },
  { id: 'P114', status: 'active', initial: 7, tag: 'slow_dnmt3a', ti: { rho: 2.0, Emax: 4.5, r: 0.03, T0: 25, e0: 1 } },
  { id: 'P005', status: 'completed', initial: 8, tag: 'quick_relapse', ti: { rho: 2.0, Emax: 2.0, r: 0.06, T0: 120, e0: 1 } },
  { id: 'P011', status: 'completed', initial: 8, tag: 'quick', ti: { rho: 2.0, Emax: 4.0, r: 0.03, T0: 20, e0: 1 } },
  { id: 'P014', status: 'completed', initial: 8, tag: 'quick', ti: { rho: 1.9, Emax: 5.0, r: 0.03, T0: 35, e0: 1 } },
  { id: 'P017', status: 'completed', initial: 8, tag: 'quick', ti: { rho: 2.1, Emax: 6.0, r: 0.03, T0: 60, e0: 1 } },
  { id: 'P003', status: 'completed', initial: 8, tag: 'second_wave_early_relapse', ti: { rho: 1.8, Emax: 1.6, r: 0.07, T0: 220, e0: 1 } },
  { id: 'P006', status: 'completed', initial: 8, tag: 'second_wave', ti: { rho: 2.0, Emax: 5.0, r: 0.03, T0: 70, e0: 1 } },
  { id: 'P009', status: 'completed', initial: 8, tag: 'second_wave', ti: { rho: 2.1, Emax: 7.0, r: 0.03, T0: 90, e0: 1 } },
  { id: 'P012', status: 'completed', initial: 8, tag: 'second_wave_partial', ti: { rho: 1.9, Emax: 3.0, r: 0.05, T0: 150, e0: 1 } },
  { id: 'P016', status: 'completed', initial: 8, tag: 'second_wave', ti: { rho: 2.2, Emax: 8.0, r: 0.03, T0: 110, e0: 1 } },
  { id: 'P039', status: 'completed', initial: 8, tag: 'persistent_cd34', ti: { rho: 2.3, Emax: 10.0, r: 0.03, T0: 200, e0: 1 } },
  { id: 'P044', status: 'completed', initial: 8, tag: 'persistent', ti: { rho: 2.2, Emax: 9.0, r: 0.03, T0: 160, e0: 1 } },
  { id: 'P062', status: 'completed', initial: 8, tag: 'persistent_ppm1d', ti: { rho: 2.0, Emax: 6.0, r: 0.03, T0: 80, e0: 1 } },
  { id: 'P202', status: 'completed', initial: 8, tag: 'persistent_tet2_tp53', ti: { rho: 1.9, Emax: 2.5, r: 0.05, T0: 130, e0: 1 } }
];
const CHECKPOINTS = [-5, 0, 7, 14, 21, 28, 42, 60, 90];
const out = [];
for (const d of DEMO) {
  const p = patients.find(x => x.patient_id === d.id);
  if (!p) throw new Error('missing ' + d.id);
  const obs = [];
  for (const r of labs.filter(x => x.patient_id === d.id)) for (const v of ['anc', 'plt', 'wbc', 'hgb', 'alc']) if (r[v] !== '') obs.push({ day: Number(r.day), variable: v, value: Number(r[v]), source: 'cbc' });
  for (const r of inflam.filter(x => x.patient_id === d.id)) for (const v of ['crp', 'ferritin']) if (r[v] !== '') obs.push({ day: Number(r.day), variable: v, value: Number(r[v]), source: 'chem' });
  const ev = events.filter(x => x.patient_id === d.id).map(e => ({ type: e.event, day_start: Number(e.start_day), day_end: Number(e.end_day), dose: num(e.dose) }));
  const t = tested.find(x => x.patient_id === d.id);
  const variants = muts.filter(x => x.patient_id === d.id).map(m => ({ gene: m.gene, protein_change: m.protein_change, vaf: Number(m.vaf), depth: Number(m.depth) }));
  const pt = {
    patient_id: d.id, label: d.id,
    demographics: { age: Number(p.age), sex: p.sex },
    product: p.product, costim: p.costim,
    history: { prior_lines: Number(p.prior_lines), prior_asct: Number(p.prior_asct), prior_bendamustine: Number(p.prior_bendamustine), months_since_last_chemo: Number(p.months_since_last_chemo), marrow_involvement: Number(p.marrow_involvement) },
    bridging: p.bridging, ld: { flu_total: Number(p.ld_flu_total), cy_total: Number(p.ld_cy_total), dose_reduced: Number(p.ld_dose_reduced) }, car_dose_frac: 1,
    crs: p.crs_grade !== '' && Number(p.crs_grade) > 0 ? { grade: Number(p.crs_grade), onset_day: num(p.crs_onset_day), end_day: num(p.crs_end_day) } : null,
    mutations: { tested: !!t, sample_day: t ? Number(t.sample_day) : null, sample_timepoint: t ? t.sample_timepoint : null, panel: t ? t.panel_version : null, variants },
    observations: obs, events: ev,
    recorded: { relapse_day: num(p.relapse_day), death_day: num(p.death_day), last_followup_day: num(p.last_followup_day), phenotype_clinical: p.phenotype_clinical || null },
    demo: { status: d.status, checkpoints: CHECKPOINTS, initialCheckpoint: d.initial, tag: d.tag, synthetic_variables: ['car_copies', 'tumor_burden'] }
  };
  // synthetic tumour–immune observations from the placeholder model
  const snap = S.snapshot(pt, 200);
  const rec = TI.simulate(d.ti, snap, S.defaultPlan(), {}, 120);
  const rng = U.rng(U.hashString(d.id) + 99);
  pt.observations.push({ day: -7, variable: 'tumor_burden', value: +(d.ti.T0 * Math.exp(0.1 * rng.normal())).toFixed(1), source: 'pet_synthetic' });
  for (const day of [3, 7, 10, 14, 21, 28, 42, 60]) { const i = rec.days.indexOf(day); pt.observations.push({ day, variable: 'car_copies', value: +(Math.max(0.005, rec.E[i] * Math.exp(0.25 * rng.normal()))).toFixed(3), source: 'qpcr_synthetic' }); }
  for (const day of [30, 90]) { const i = rec.days.indexOf(day); pt.observations.push({ day, variable: 'tumor_burden', value: +(Math.max(0.1, rec.T[i] * Math.exp(0.2 * rng.normal()))).toFixed(1), source: 'pet_synthetic' }); }
  pt.observations.sort((a, b) => a.day - b.day || a.variable.localeCompare(b.variable));
  out.push(pt);
  console.log(d.id, d.status, p.phenotype_clinical, 'obs', obs.length, 'events', ev.length, 'variants', variants.length, 'TI resp:', (rec.T[rec.days.indexOf(90)] / d.ti.T0 < 0.1), 'peakE', Math.max(...rec.E).toFixed(1));
}
const payload = { generated: new Date().toISOString(), source: 'Aim2/TestData (synthetic cohort)', patients: out };
fs.writeFileSync(path.join(__dirname, '..', 'data', 'demo_patients.json'), JSON.stringify(payload));
fs.writeFileSync(path.join(__dirname, '..', 'data', 'demo_patients.js'), '/* generated by tools/build_demo_data.js — do not edit */\nwindow.GL_DEMO_DATA = ' + JSON.stringify(payload) + ';\n');
console.log('wrote data/demo_patients.js (' + out.length + ' patients)');
