/* GreenLight · core/plans.js
 * Candidate treatment-plan library, plan description helpers and the
 * multi-objective utility that ranks plans across modules.
 * Doses/thresholds are illustrative — not clinical guidance.
 */
(function (root) {
  const GL = root.GL = root.GL || {};
  const S = GL.schema, U = GL.util;
  const P = GL.plans = {};

  function make(id, nameKey, mod) {
    const p = S.defaultPlan(); p.id = id; p.nameKey = nameKey; p.custom = false;
    mod(p.stages); return p;
  }
  P.library = function () {
    return [
      make('soc', 'plan.soc', () => {}),
      make('soc_proph_gcsf', 'plan.soc_proph_gcsf', s => { s.cytopenia.gcsf = 'prophylactic'; s.cytopenia.gcsf_start_day = 5; }),
      make('reduced_ld', 'plan.reduced_ld', s => { s.ld.flu_frac = 0.75; s.ld.cy_frac = 0.75; }),
      make('soc_cd34', 'plan.soc_cd34', s => { s.cytopenia.cd34_backup = true; s.cytopenia.cd34_boost = 'if_persistent'; s.cytopenia.cd34_boost_day = 45; }),
      make('reduced_ld_cd34', 'plan.reduced_ld_cd34', s => { s.ld.flu_frac = 0.75; s.ld.cy_frac = 0.75; s.cytopenia.cd34_backup = true; s.cytopenia.cd34_boost = 'if_persistent'; s.cytopenia.cd34_boost_day = 45; }),
      make('early_antiinflam', 'plan.early_antiinflam', s => { s.crs.toci_grade = 1; s.crs.anakinra_prophylaxis = true; }),
      make('early_antiinflam_cd34', 'plan.early_antiinflam_cd34', s => { s.crs.toci_grade = 1; s.crs.anakinra_prophylaxis = true; s.cytopenia.cd34_backup = true; s.cytopenia.cd34_boost = 'if_persistent'; s.cytopenia.cd34_boost_day = 45; s.cytopenia.tpo_ra = true; }),
      make('soc_tpo', 'plan.soc_tpo', s => { s.cytopenia.tpo_ra = true; s.cytopenia.tpo_plt_thr = 50; }),
      make('intensified', 'plan.intensified', s => { s.infusion.car_dose_frac = 1.25; s.cytopenia.gcsf = 'prophylactic'; s.cytopenia.gcsf_start_day = 5; })
    ];
  };

  // plan parameter catalogue (drives the Plan Builder UI and the diff view)
  P.fields = [
    { stage: 'bridging', key: 'type', labelKey: 'field.bridging.type', type: 'select', options: ['as_recorded', 'none', 'chemo', 'radiation'] },
    { stage: 'ld', key: 'flu_frac', labelKey: 'field.ld.flu_frac', unit: '% of 30 mg/m² ×3', type: 'range', min: 0.5, max: 1.25, step: 0.05, fmt: v => Math.round(v * 100) + '%' },
    { stage: 'ld', key: 'cy_frac', labelKey: 'field.ld.cy_frac', unit: '% of 500 mg/m² ×3', type: 'range', min: 0.5, max: 1.25, step: 0.05, fmt: v => Math.round(v * 100) + '%' },
    { stage: 'infusion', key: 'car_dose_frac', labelKey: 'field.infusion.car_dose_frac', unit: '% of label dose', type: 'range', min: 0.5, max: 1.5, step: 0.05, fmt: v => Math.round(v * 100) + '%' },
    { stage: 'crs', key: 'toci_grade', labelKey: 'field.crs.toci_grade', unit: 'CRS grade', type: 'select', options: [1, 2, 3] },
    { stage: 'crs', key: 'steroid_grade', labelKey: 'field.crs.steroid_grade', unit: 'CRS grade', type: 'select', options: [2, 3] },
    { stage: 'crs', key: 'anakinra_prophylaxis', labelKey: 'field.crs.anakinra_prophylaxis', type: 'bool' },
    { stage: 'cytopenia', key: 'gcsf', labelKey: 'field.cytopenia.gcsf', type: 'select', options: ['none', 'reactive', 'prophylactic'] },
    { stage: 'cytopenia', key: 'gcsf_start_day', labelKey: 'field.cytopenia.gcsf_start_day', unit: 'day', type: 'range', min: 3, max: 28, step: 1, fmt: v => 'D+' + v },
    { stage: 'cytopenia', key: 'gcsf_anc_thr', labelKey: 'field.cytopenia.gcsf_anc_thr', unit: '×10⁹/L', type: 'range', min: 0.2, max: 1.0, step: 0.1, fmt: v => v.toFixed(1) },
    { stage: 'cytopenia', key: 'cd34_backup', labelKey: 'field.cytopenia.cd34_backup', type: 'bool' },
    { stage: 'cytopenia', key: 'cd34_boost', labelKey: 'field.cytopenia.cd34_boost', type: 'select', options: ['none', 'if_persistent'] },
    { stage: 'cytopenia', key: 'cd34_boost_day', labelKey: 'field.cytopenia.cd34_boost_day', unit: 'day', type: 'range', min: 30, max: 90, step: 5, fmt: v => 'D+' + v },
    { stage: 'cytopenia', key: 'tpo_ra', labelKey: 'field.cytopenia.tpo_ra', type: 'bool' },
    { stage: 'cytopenia', key: 'plt_tx_thr', labelKey: 'field.cytopenia.plt_tx_thr', unit: '×10⁹/L', type: 'range', min: 10, max: 50, step: 5, fmt: v => String(v) },
    { stage: 'followup', key: 'cbc_interval_days', labelKey: 'field.followup.cbc_interval_days', unit: 'days', type: 'range', min: 1, max: 7, step: 1, fmt: v => v + ' d' }
  ];

  P.diff = function (a, b) {
    const out = [];
    for (const f of P.fields) {
      const va = a.stages[f.stage][f.key], vb = b.stages[f.stage][f.key];
      if (va !== vb) out.push({ stage: f.stage, key: f.key, labelKey: f.labelKey, from: va, to: vb, field: f });
    }
    return out;
  };
  P.formatValue = function (field, v) {
    if (field.type === 'bool') return v ? '✓' : '—';
    if (field.fmt) return field.fmt(Number(v));
    const k = 'opt.' + field.key + '.' + v;
    return GL.i18n && GL.i18n.has(k) ? GL.i18n.t(k) : String(v);
  };
  // plans that differ only in stages already administered are the same plan from today on
  P.effectiveKey = function (plan, asOfDay) {
    const st = {};
    for (const stg of S.stageOrder) st[stg] = S.stageLocked(stg, asOfDay) ? 'locked' : plan.stages[stg];
    return JSON.stringify(st);
  };

  // ---- utility ------------------------------------------------------------------
  // score = w_eff·P(response)
  //       − w_rec·[1.0·P(persistent) + 0.8·P(delayed) + 0.6·P(second wave)]   (recovery-phenotype risk;
  //                                                              P(delayed) only from the Aim 4 module)
  //       − w_tox·[0.5·P(CRS≥3) + 0.3·min(1, sevDays/30) + 0.25·P(PLT<20)]
  //       − burden
  // The recovery term is separate so that "avoid persistent cytopenia and the second wave" is an
  // explicit, physician-adjustable priority rather than one entry in a toxicity sum.
  P.defaultWeights = { eff: 1.0, rec: 1.2, tox: 1.0 };
  P.recoveryWeights = { persistent: 1.0, delayed: 0.8, second_wave: 0.6 };
  P.normWeights = (w) => Object.assign({}, P.defaultWeights, w || {});
  // outcomes of the marrow model (Aim 4 parent model, module 'marrow4')
  P.marrowOut = (outcomes) => (outcomes && outcomes.marrow4) || {};
  P.phenoKeys = ['quick', 'second_wave', 'delayed', 'persistent'];
  // {class: probability} for the classes the module provides (quick, second wave, delayed, persistent)
  P.phenoProbs = (m) => { const o = {}; P.phenoKeys.forEach(k => { if (m && m['p_' + k]) o[k] = m['p_' + k].value; }); return o; };
  // outcomes: { marrow4: {...}, tumor_immune: {...} } (whatever modules are present)
  P.utility = function (outcomes, plan, weights) {
    weights = P.normWeights(weights);
    const m = P.marrowOut(outcomes), t = outcomes.tumor_immune || {};
    const g = (o, k) => (o && o[k] && typeof o[k].value === 'number') ? o[k].value : null;
    const comp = {};
    comp.efficacy = g(t, 'p_response') !== null ? g(t, 'p_response') : null;
    const rec = [];
    if (g(m, 'p_persistent') !== null) rec.push(P.recoveryWeights.persistent * g(m, 'p_persistent'));
    if (g(m, 'p_delayed') !== null) rec.push(P.recoveryWeights.delayed * g(m, 'p_delayed'));
    if (g(m, 'p_second_wave') !== null) rec.push(P.recoveryWeights.second_wave * g(m, 'p_second_wave'));
    comp.recovery = rec.length ? U.sum(rec) : null;
    const tox = [];
    if (g(m, 'p_crs3') !== null) tox.push(0.5 * g(m, 'p_crs3'));
    if (g(m, 'sev_neut_days') !== null) tox.push(0.3 * Math.min(1, g(m, 'sev_neut_days') / 30));
    if (g(m, 'p_plt20') !== null) tox.push(0.25 * g(m, 'p_plt20'));
    comp.toxicity = tox.length ? U.sum(tox) : null;
    let burden = 0;
    const s = plan.stages;
    if (g(m, 'gcsf_courses') !== null) burden += 0.03 * g(m, 'gcsf_courses');
    if (g(m, 'plt_tx') !== null) burden += 0.02 * g(m, 'plt_tx');
    if (s.cytopenia.cd34_backup) burden += 0.05;
    if (s.crs.anakinra_prophylaxis) burden += 0.03;
    if (s.cytopenia.tpo_ra) burden += 0.02;
    if (s.cytopenia.gcsf === 'prophylactic') burden += 0.02;
    if (s.infusion.car_dose_frac > 1) burden += 0.04;
    comp.burden = burden;
    const score = weights.eff * (comp.efficacy || 0) - weights.rec * (comp.recovery || 0) - weights.tox * (comp.toxicity || 0) - burden;
    return { score, components: comp };
  };
  // Which candidate is best on each recovery pattern. list: ranked [{ id, plan, outcomes, rank }].
  // Ties go to the better-ranked plan (the list is in rank order).
  P.phenoBest = function (list) {
    const g = (x, k) => { const m = P.marrowOut(x.outcomes); return m[k] ? m[k].value : null; };
    const pick = (k, better) => { let b = null; for (const x of list) { const v = g(x, k); if (v === null) continue; if (!b || better(v, b.value)) b = { id: x.id, plan: x.plan, rank: x.rank, value: v }; } return b; };
    return { quick: pick('p_quick', (a, b) => a > b), second_wave: pick('p_second_wave', (a, b) => a < b), delayed: pick('p_delayed', (a, b) => a < b), persistent: pick('p_persistent', (a, b) => a < b) };
  };
})(typeof window !== 'undefined' ? window : global);
