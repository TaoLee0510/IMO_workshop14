/* Build GreenLight's demo cohort from the Aim 4 parent model (v1.2): fully synthetic patients
 * (no value comes from a real record). Each patient has hand-set "true" parameters chosen to show
 * one Aim 4 mechanism; counts are simulated with js/modules/aim4_core.js (the same equations as
 * Aim4/Code/src/aim4_model.cpp) with the population constants of data/aim4_pop.js, plus
 * measurement noise. Ferritin follows a Gamma-shaped course; CRP is derived from it. Tumour–immune
 * observations (car_copies, tumor_burden) come from the placeholder tumour–immune module.
 * Usage: node tools/build_demo_data.js
 * Writes data/demo_patients.js (browser, file:// friendly) and data/demo_patients.json.
 */
const fs = require('fs'), path = require('path');
const GL = require('./load_core');
const S = GL.schema, TI = GL.modules.tumor_immune, U = GL.util, C = GL.aim4;
const pop = global.GL_AIM4_POP;
if (!pop) throw new Error('data/aim4_pop.js not found');
const sp = C.spec(pop.model), st = pop.struct;

const B = { r_s: 0.08, th_ap: 2.3, th_as: 0.12, th_bs: 0.02, f0: 0.01 };   // typical patient
// status 'active': in treatment, data arrive checkpoint by checkpoint; 'completed': full record
// truth: individual parameters; fer: ferritin course (A = peak log-ratio, tp = peak day); base: set points
const DEMO = [
  { id: 'D01', status: 'active', initial: 0, tag: 'b_second_wave', age: 63, sex: 'M', product: 'A', costim: 'CD28', lines: 3, asct: 1, bridging: 'none', egfr: 82,
    truth: { ...B, th_bs: 0.06 }, fer: { A: 2.6, tp: 9, base: 900 }, base: { anc: 3.4, plt: 190, amc: 0.55, alc: 1.1, hgb: 12.1 }, crs: 3, ti: { rho: 2.2, Emax: 9.0, r: 0.03, T0: 150, e0: 1 } },
  { id: 'D02', status: 'active', initial: 2, tag: 'a_second_wave', age: 71, sex: 'F', product: 'B', costim: '4-1BB', lines: 4, asct: 1, bridging: 'chemo', egfr: 64,
    truth: { ...B, r_s: 0.035, th_as: 1.3, th_bs: 0.005 }, fer: { A: 0.6, tp: 7, base: 450 }, base: { anc: 4.2, plt: 160, amc: 0.6, alc: 0.8, hgb: 11.2 }, crs: 1, ti: { rho: 1.9, Emax: 1.8, r: 0.05, T0: 200, e0: 1 } },
  { id: 'D03', status: 'active', initial: 3, tag: 'quick', age: 52, sex: 'M', product: 'A', costim: 'CD28', lines: 2, asct: 0, bridging: 'none', egfr: 98,
    truth: { ...B }, fer: { A: 0.8, tp: 7, base: 380 }, base: { anc: 3.6, plt: 240, amc: 0.5, alc: 1.4, hgb: 13.1 }, crs: 1, ti: { rho: 2.0, Emax: 5.0, r: 0.03, T0: 40, e0: 1 } },
  { id: 'D04', status: 'active', initial: 4, tag: 'clone_ddr', age: 74, sex: 'M', product: 'B', costim: '4-1BB', lines: 4, asct: 1, bridging: 'chemo', egfr: 70,
    truth: { ...B, f0: 0.4, th_as: 0.6, r_s: 0.05, th_bs: 0.03 }, fer: { A: 1.2, tp: 8, base: 620 }, base: { anc: 3.2, plt: 150, amc: 0.5, alc: 0.9, hgb: 11.0 }, crs: 2,
    mut: [{ gene: 'TP53', protein_change: 'p.R248Q' }, { gene: 'PPM1D', protein_change: 'p.R458*', share: 0.15 }], ti: { rho: 2.0, Emax: 4.5, r: 0.03, T0: 120, e0: 1 } },
  { id: 'D05', status: 'active', initial: 5, tag: 'quick_reduced_ld', age: 66, sex: 'F', product: 'C', costim: '4-1BB', lines: 3, asct: 0, bridging: 'radiation', egfr: 55, ld: 0.75,
    truth: { ...B, r_s: 0.07, th_as: 0.5 }, fer: { A: 0.8, tp: 8, base: 520 }, base: { anc: 3.8, plt: 210, amc: 0.45, alc: 1.0, hgb: 12.4 }, crs: 1, ti: { rho: 2.0, Emax: 4.0, r: 0.03, T0: 30, e0: 1 } },
  { id: 'D06', status: 'active', initial: 6, tag: 'persistent_chronic_inflammation', age: 58, sex: 'M', product: 'A', costim: 'CD28', lines: 3, asct: 1, bridging: 'none', egfr: 88,
    truth: { ...B, r_s: 0.03, th_as: 0.5, th_bs: 0.08, th_ap: 3 }, fer: { A: 3.2, tp: 28, base: 1100 }, base: { anc: 3.0, plt: 170, amc: 0.5, alc: 1.2, hgb: 11.8 }, crs: 3, ti: { rho: 2.3, Emax: 12.0, r: 0.03, T0: 300, e0: 1 } },
  { id: 'D07', status: 'active', initial: 7, tag: 'delayed_low_reserve', age: 76, sex: 'F', product: 'B', costim: '4-1BB', lines: 5, asct: 1, bridging: 'chemo', egfr: 48,
    truth: { ...B, r_s: 0.02, th_as: 2.5 }, fer: { A: 1.0, tp: 8, base: 700 }, base: { anc: 2.8, plt: 130, amc: 0.45, alc: 0.7, hgb: 10.6 }, crs: 2, gcsf: true, ti: { rho: 2.1, Emax: 8.0, r: 0.03, T0: 120, e0: 1 } },
  { id: 'D08', status: 'active', initial: 7, tag: 'quick_high_ferritin', age: 49, sex: 'F', product: 'A', costim: 'CD28', lines: 2, asct: 0, bridging: 'none', egfr: 105,
    truth: { ...B, r_s: 0.12 }, fer: { A: 1.8, tp: 8, base: 300 }, base: { anc: 4.0, plt: 260, amc: 0.6, alc: 1.6, hgb: 13.4 }, crs: 2, ti: { rho: 2.0, Emax: 4.5, r: 0.03, T0: 25, e0: 1 } },
  { id: 'D09', status: 'completed', initial: 8, tag: 'quick', age: 61, sex: 'M', product: 'A', costim: 'CD28', lines: 2, asct: 1, bridging: 'none', egfr: 91,
    truth: { ...B }, fer: { A: 0.7, tp: 6, base: 420 }, base: { anc: 3.8, plt: 220, amc: 0.55, alc: 1.3, hgb: 13.0 }, crs: 1, ti: { rho: 2.0, Emax: 4.0, r: 0.03, T0: 20, e0: 1 } },
  { id: 'D10', status: 'completed', initial: 8, tag: 'quick', age: 45, sex: 'F', product: 'C', costim: '4-1BB', lines: 2, asct: 0, bridging: 'radiation', egfr: 112,
    truth: { ...B, r_s: 0.1 }, fer: { A: 1.1, tp: 8, base: 250 }, base: { anc: 3.3, plt: 250, amc: 0.5, alc: 1.5, hgb: 12.8 }, crs: 1, ti: { rho: 1.9, Emax: 5.0, r: 0.03, T0: 35, e0: 1 } },
  { id: 'D11', status: 'completed', initial: 8, tag: 'a_second_wave', age: 68, sex: 'M', product: 'B', costim: '4-1BB', lines: 4, asct: 1, bridging: 'chemo', egfr: 74,
    truth: { ...B, r_s: 0.035, th_as: 1.3, th_bs: 0.005 }, fer: { A: 0.5, tp: 7, base: 560 }, base: { anc: 4.4, plt: 175, amc: 0.6, alc: 0.9, hgb: 11.6 }, crs: 1, ti: { rho: 2.1, Emax: 6.0, r: 0.03, T0: 60, e0: 1 } },
  { id: 'D12', status: 'completed', initial: 8, tag: 'b_second_wave', age: 57, sex: 'F', product: 'A', costim: 'CD28', lines: 3, asct: 0, bridging: 'none', egfr: 95,
    truth: { ...B, th_bs: 0.08, r_s: 0.06 }, fer: { A: 2.8, tp: 10, base: 820 }, base: { anc: 3.5, plt: 205, amc: 0.5, alc: 1.2, hgb: 12.2 }, crs: 3, ti: { rho: 2.2, Emax: 8.0, r: 0.03, T0: 110, e0: 1 } },
  { id: 'D13', status: 'completed', initial: 8, tag: 'ab_delayed', age: 70, sex: 'M', product: 'B', costim: '4-1BB', lines: 4, asct: 1, bridging: 'chemo', egfr: 66,
    truth: { ...B, r_s: 0.04, th_as: 1.2, th_bs: 0.06 }, fer: { A: 2.2, tp: 9, base: 760 }, base: { anc: 3.1, plt: 160, amc: 0.5, alc: 0.8, hgb: 11.1 }, crs: 2, gcsf: true, ti: { rho: 2.0, Emax: 7.0, r: 0.03, T0: 90, e0: 1 } },
  { id: 'D14', status: 'completed', initial: 8, tag: 'delayed_low_reserve', age: 73, sex: 'F', product: 'A', costim: 'CD28', lines: 5, asct: 1, bridging: 'chemo', egfr: 52,
    truth: { ...B, r_s: 0.015, th_as: 3.0, th_ap: 3 }, fer: { A: 1.2, tp: 8, base: 680 }, base: { anc: 2.6, plt: 120, amc: 0.4, alc: 0.6, hgb: 10.4 }, crs: 2, gcsf: true, tpo: 28, ti: { rho: 2.0, Emax: 6.0, r: 0.03, T0: 80, e0: 1 } },
  { id: 'D15', status: 'completed', initial: 8, tag: 'persistent_chronic_inflammation_cd34', age: 64, sex: 'M', product: 'A', costim: 'CD28', lines: 3, asct: 1, bridging: 'none', egfr: 84,
    truth: { ...B, r_s: 0.01, th_as: 2.0, th_bs: 0.08, th_ap: 3 }, fer: { A: 3.5, tp: 30, base: 1300 }, base: { anc: 3.2, plt: 180, amc: 0.5, alc: 1.0, hgb: 11.5 }, crs: 3, cd34: 60, ti: { rho: 2.3, Emax: 10.0, r: 0.03, T0: 200, e0: 1 } },
  { id: 'D16', status: 'completed', initial: 8, tag: 'clone_ddr_delayed', age: 75, sex: 'F', product: 'B', costim: '4-1BB', lines: 4, asct: 1, bridging: 'chemo', egfr: 60,
    truth: { ...B, f0: 0.55, r_s: 0.02, th_as: 2.5, th_bs: 0.06, th_ap: 3 }, fer: { A: 2.0, tp: 9, base: 900 }, base: { anc: 3.0, plt: 150, amc: 0.45, alc: 0.7, hgb: 10.9 }, crs: 2,
    mut: [{ gene: 'PPM1D', protein_change: 'p.R458*' }, { gene: 'TP53', protein_change: 'p.Y220C', share: 0.3 }], ti: { rho: 2.0, Emax: 2.5, r: 0.05, T0: 130, e0: 1 } }
];
const CHECKPOINTS = [-5, 0, 7, 14, 21, 28, 42, 60, 90];
const CBC_DAYS = [-30, -21, -14, -10].concat(Array.from({ length: 29 }, (_, i) => i - 7), Array.from({ length: 10 }, (_, i) => 23 + 2 * i), Array.from({ length: 10 }, (_, i) => 45 + 5 * i));
const CHEM_DAYS = [-7, -5, -3].concat(Array.from({ length: 15 }, (_, i) => i), Array.from({ length: 8 }, (_, i) => 16 + 2 * i), [35, 42, 49, 56, 63, 70, 77, 84, 90]);
const SD = { anc: 0.15, plt: 0.12, amc: 0.2, alc: 0.2, wbc: 0.12, hgb: 0.04, ferritin: 0.12, crp: 0.2 };
const FLOOR = { anc: 0.1, amc: 0.01, alc: 0.02 };
const T0 = -31, H = (pop.kernel_pop && pop.kernel_pop.h) || 2;

function simulate(d, gcsf, tpo, cd34) {
  const design = { t0: T0, egfr: d.egfr, bridging: d.bridging === 'chemo' ? 'cy' : 'none', apheresis_day: -28,
    ld_scale: d.ld || 1, gcsf, tpo, cd34 };
  const drv = C.buildDriver({ t0: T0, tEnd: 130, smooth: null, upTo: 0, kernel: { A: d.fer.A, tp: d.fer.tp }, h: H });
  const sim = C.simulate(C.buildPar(d.truth, st, sp), { design, struct: st, drv, tEnd: 90 });
  if (!sim) throw new Error('simulation failed for ' + d.id);
  return { sim, drv };
}
const at = (sim, key, day) => sim[key][sim.days.indexOf(day)];

const out = [];
for (const d of DEMO) {
  const rng = U.rng(U.hashString(d.id) + 7);
  // supportive care: reactive G-CSF (5 days when the ANC falls below 1.0 after D+14), TPO-RA, CD34+ boost
  let gcsf = [];
  const tpo = d.tpo ? [[d.tpo, d.tpo + 60]] : [];
  const cd34 = d.cd34 ? [{ day: d.cd34, add: pop.cd34_add || 0.3 }] : [];
  let { sim, drv } = simulate(d, gcsf, tpo, cd34);
  if (d.gcsf) {
    const low = sim.days.find(day => day >= 15 && day <= 30 && at(sim, 'anc_rel', day) * d.base.anc < 1.0);
    if (low !== undefined) { gcsf = [[low, low + 5]]; ({ sim, drv } = simulate(d, gcsf, tpo, cd34)); }
  }
  const obs = [];
  const noisy = (v, x) => x * Math.exp(SD[v] * rng.normal());
  const push = (day, v, x, src) => {
    let val = noisy(v, x), cens = false;
    if (FLOOR[v] && val < FLOOR[v]) { val = FLOOR[v]; cens = true; }
    const o = { day, variable: v, value: +val.toPrecision(3), source: src };
    if (cens) o.cens = true;
    obs.push(o);
  };
  for (const day of CBC_DAYS) {
    const anc = at(sim, 'anc_rel', day) * d.base.anc, plt = at(sim, 'plt_rel', day) * d.base.plt;
    const amc = at(sim, 'amc_rel', day) * d.base.amc, alc = at(sim, 'alc_rel', day) * d.base.alc;
    push(day, 'anc', anc, 'cbc'); push(day, 'plt', Math.max(plt, 1), 'cbc'); push(day, 'amc', amc, 'cbc'); push(day, 'alc', alc, 'cbc');
    push(day, 'wbc', anc + amc + alc + 0.15, 'cbc');
    push(day, 'hgb', d.base.hgb * (0.82 + 0.18 * Math.min(1, at(sim, 'plt_rel', day))), 'cbc');
  }
  const fer0 = d.fer.base, Dk = (t) => C.kernel(t, d.fer.A, d.fer.tp, H);
  for (const day of CHEM_DAYS) {
    push(day, 'ferritin', fer0 * Math.exp(Dk(day)), 'chem');
    push(day, 'crp', 3 + 40 * Dk(day + 3), 'chem');
  }
  // events
  const ev = [];
  if (d.crs >= 2) ev.push({ type: 'tocilizumab', day_start: 3, day_end: 3, dose: 8 });
  if (d.crs >= 3) ev.push({ type: 'steroid', day_start: 4, day_end: 8, dose: 10 });
  gcsf.forEach(w => ev.push({ type: 'gcsf', day_start: w[0], day_end: w[1] - 1, dose: 5 }));
  tpo.forEach(w => ev.push({ type: 'tpo_ra', day_start: w[0], day_end: w[0] + 28, dose: 75 }));
  cd34.forEach(c => ev.push({ type: 'cd34_boost', day_start: c.day, day_end: c.day, dose: 2 }));
  for (const day of CBC_DAYS) if (day >= 3 && at(sim, 'plt_rel', day) * d.base.plt < 15) { ev.push({ type: 'plt_tx', day_start: day, day_end: day, dose: 1 }); }
  // clonal haematopoiesis panel: VAF = f0 / 2 (heterozygous) for the leading variant
  const variants = (d.mut || []).map((m, i) => ({ gene: m.gene, protein_change: m.protein_change,
    vaf: +((i === 0 ? d.truth.f0 / 2 : d.truth.f0 / 2 * (m.share || 0.2)) * Math.exp(0.08 * rng.normal())).toFixed(4), depth: 1500 + Math.round(1500 * rng.uniform(0, 1)) }));
  const series = obs.filter(o => o.variable === 'anc').map(o => ({ day: o.day, value: o.value, cens: !!o.cens }));
  const pheno = C.phenotype(series, true);
  const pt = {
    patient_id: d.id, label: d.id,
    demographics: { age: d.age, sex: d.sex },
    product: d.product, costim: d.costim,
    history: { prior_lines: d.lines, prior_asct: d.asct, prior_bendamustine: d.lines >= 4 ? 1 : 0, months_since_last_chemo: +(0.8 + 4 * rng.uniform(0, 1)).toFixed(1), marrow_involvement: 0 },
    bridging: d.bridging, ld: { flu_total: Math.round(90 * (d.ld || 1)), cy_total: Math.round(1500 * (d.ld || 1)), dose_reduced: d.ld ? 1 : 0 }, car_dose_frac: 1,
    crs: d.crs > 0 ? { grade: d.crs, onset_day: 2, end_day: d.crs >= 3 ? 12 : 7 } : null,
    mutations: { tested: true, sample_day: -32, sample_timepoint: 'apheresis', panel: 'Myeloid panel (synthetic)', variants },
    aim4: { egfr: d.egfr, bridging_class: d.bridging === 'chemo' ? 'cy' : 'none', apheresis_day: -28 },
    observations: obs, events: ev,
    recorded: { relapse_day: null, death_day: null, last_followup_day: 90, phenotype_clinical: pheno },
    demo: { status: d.status, checkpoints: CHECKPOINTS, initialCheckpoint: d.initial, tag: d.tag, synthetic_variables: 'all' }
  };
  // synthetic tumour–immune observations from the placeholder model
  const snap = S.snapshot(pt, 200);
  const rec = TI.simulate(d.ti, snap, S.defaultPlan(), {}, 120);
  pt.observations.push({ day: -7, variable: 'tumor_burden', value: +(d.ti.T0 * Math.exp(0.1 * rng.normal())).toFixed(1), source: 'pet_synthetic' });
  for (const day of [3, 7, 10, 14, 21, 28, 42, 60]) { const i = rec.days.indexOf(day); pt.observations.push({ day, variable: 'car_copies', value: +(Math.max(0.005, rec.E[i] * Math.exp(0.25 * rng.normal()))).toFixed(3), source: 'qpcr_synthetic' }); }
  for (const day of [30, 90]) { const i = rec.days.indexOf(day); pt.observations.push({ day, variable: 'tumor_burden', value: +(Math.max(0.1, rec.T[i] * Math.exp(0.2 * rng.normal()))).toFixed(1), source: 'pet_synthetic' }); }
  pt.observations.sort((a, b) => a.day - b.day || a.variable.localeCompare(b.variable));
  out.push(pt);
  console.log(d.id.padEnd(4), d.status.padEnd(9), d.tag.padEnd(38), 'phenotype', String(pheno).padEnd(12), 'obs', obs.length, 'events', ev.length, 'gcsf', JSON.stringify(gcsf));
}
const payload = { generated: new Date().toISOString(), source: 'Aim 4 parent model v1.2 (synthetic demo cohort; tools/build_demo_data.js)', patients: out };
fs.writeFileSync(path.join(__dirname, '..', 'data', 'demo_patients.json'), JSON.stringify(payload));
fs.writeFileSync(path.join(__dirname, '..', 'data', 'demo_patients.js'), '/* generated by tools/build_demo_data.js — do not edit */\nwindow.GL_DEMO_DATA = ' + JSON.stringify(payload) + ';\n');
console.log('wrote data/demo_patients.js (' + out.length + ' patients)');
