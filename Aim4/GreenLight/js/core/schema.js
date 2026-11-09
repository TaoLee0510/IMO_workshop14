/* GreenLight · core/schema.js
 * Data contract (patient snapshot, long-format observations, events),
 * treatment-plan contract (staged intervention timeline) and the module
 * interface contract. Every adapter and module speaks these shapes.
 * See docs/DATA_API.md and docs/MODULE_API.md.
 */
(function (root) {
  const GL = root.GL = root.GL || {};
  const S = GL.schema = {};

  // ---- observed variables (long format: {day, variable, value}) -----------
  // Units are display units; `log` tells charts/likelihoods to work on log scale.
  S.variables = {
    anc:          { labelKey: 'var.anc',          unit: '×10⁹/L', log: true,  thresholds: [0.5, 1.0], group: 'blood' },
    plt:          { labelKey: 'var.plt',          unit: '×10⁹/L', log: true,  thresholds: [20, 50],   group: 'blood' },
    hgb:          { labelKey: 'var.hgb',          unit: 'g/dL',   log: false, thresholds: [8],        group: 'blood' },
    wbc:          { labelKey: 'var.wbc',          unit: '×10⁹/L', log: true,  thresholds: [],         group: 'blood' },
    alc:          { labelKey: 'var.alc',          unit: '×10⁹/L', log: true,  thresholds: [],         group: 'blood' },
    amc:          { labelKey: 'var.amc',          unit: '×10⁹/L', log: true,  thresholds: [0.1],      group: 'blood' },
    crp:          { labelKey: 'var.crp',          unit: 'mg/L',   log: true,  thresholds: [50],       group: 'inflammation' },
    ferritin:     { labelKey: 'var.ferritin',     unit: 'ng/mL',  log: true,  thresholds: [],         group: 'inflammation' },
    car_copies:   { labelKey: 'var.car_copies',   unit: '×10³ copies/µg', log: true, thresholds: [], group: 'immune' },
    tumor_burden: { labelKey: 'var.tumor_burden', unit: 'mL (MTV)', log: true, thresholds: [],       group: 'tumor' },
    ldh:          { labelKey: 'var.ldh',          unit: 'U/L',    log: true,  thresholds: [250],      group: 'tumor' }
  };

  S.eventTypes = ['tocilizumab', 'steroid', 'anakinra', 'gcsf', 'tpo_ra', 'plt_tx', 'rbc_tx', 'cd34_boost', 'cd34_backup', 'bridging', 'lymphodepletion', 'infusion'];

  // ---- patient snapshot ----------------------------------------------------
  // What a module sees. Built by GL.schema.snapshot(patient, asOfDay): only
  // observations/events that were knowable at asOfDay are included.
  S.snapshot = function (patient, asOfDay) {
    const obs = patient.observations.filter(o => o.day <= asOfDay);
    const ev = patient.events.filter(e => e.day_start <= asOfDay).map(e => ({
      ...e, day_end: Math.min(e.day_end === undefined || e.day_end === null ? e.day_start : e.day_end, Math.max(asOfDay, e.day_start))
    }));
    const crsKnown = patient.crs && patient.crs.onset_day !== null && patient.crs.onset_day !== undefined && patient.crs.onset_day <= asOfDay;
    return {
      patient_id: patient.patient_id,
      asOfDay,
      demographics: patient.demographics,
      product: patient.product, costim: patient.costim,
      history: patient.history,            // prior lines, ASCT, bendamustine, months since chemo, marrow involvement
      bridging: patient.bridging,
      ld: patient.ld,                      // recorded LD doses (only meaningful once asOfDay >= -5)
      crs: crsKnown ? patient.crs : null,
      mutations: patient.mutations,        // pre-treatment panel (known from the start)
      aim4: patient.aim4,                  // Aim 4 design inputs (eGFR, bridging class, set points, covariates) when known
      observations: obs,
      events: ev
    };
  };

  S.obsByVar = function (snapshot, variable) {
    return snapshot.observations.filter(o => o.variable === variable).sort((a, b) => a.day - b.day);
  };

  // ---- treatment plan contract --------------------------------------------
  // A plan is a staged intervention timeline. Modules read only the keys they
  // declare in consumes.planKeys. Stages before the snapshot day are "locked"
  // (already administered) and are overridden by the recorded history.
  S.stageOrder = ['bridging', 'ld', 'infusion', 'crs', 'cytopenia', 'followup'];
  S.stageWindows = { bridging: [-21, -7], ld: [-5, -3], infusion: [0, 0], crs: [0, 28], cytopenia: [10, 90], followup: [0, 120] };
  S.stageLocked = function (stage, asOfDay) {
    const w = S.stageWindows[stage];
    if (stage === 'crs' || stage === 'cytopenia' || stage === 'followup') return false; // policies keep applying to the future
    return asOfDay >= w[0];
  };
  S.defaultPlan = function () {
    return {
      id: 'soc', nameKey: 'plan.soc', custom: false,
      stages: {
        bridging:  { type: 'as_recorded' },               // none | chemo | radiation | as_recorded
        ld:        { flu_frac: 1.0, cy_frac: 1.0 },       // fraction of standard Flu 30 mg/m² ×3, Cy 500 mg/m² ×3
        infusion:  { car_dose_frac: 1.0 },                // fraction of label dose
        crs:       { toci_grade: 2, steroid_grade: 3, anakinra_prophylaxis: false },
        cytopenia: { gcsf: 'reactive', gcsf_start_day: 14, gcsf_anc_thr: 0.5, gcsf_days: 5,
                     cd34_backup: false, cd34_boost: 'none', cd34_boost_day: 60,
                     tpo_ra: false, tpo_plt_thr: 50, plt_tx_thr: 20 },
        followup:  { cbc_interval_days: 2 }
      }
    };
  };
  S.validatePlan = function (plan) {
    const errs = [];
    if (!plan || !plan.stages) return ['plan.stages missing'];
    for (const st of S.stageOrder) if (!plan.stages[st]) errs.push('stage missing: ' + st);
    const c = plan.stages.cytopenia || {};
    if (c.gcsf && !['none', 'reactive', 'prophylactic'].includes(c.gcsf)) errs.push('cytopenia.gcsf invalid');
    const ld = plan.stages.ld || {};
    if (ld.flu_frac < 0 || ld.flu_frac > 1.5) errs.push('ld.flu_frac out of range');
    if (ld.cy_frac < 0 || ld.cy_frac > 1.5) errs.push('ld.cy_frac out of range');
    return errs;
  };

  // ---- module interface contract ------------------------------------------
  // A module object must provide:
  //   meta: { id, nameKey, owner, version, kind: 'local'|'remote', placeholder?: bool }
  //   consumes: { variables: [..], planKeys: [..], signals: [..] }
  //   provides: { signals: [..], params: [{id,labelKey,unit?,log?}], outcomes: [{id,labelKey,kind:'prob'|'value',direction:'higher'|'lower',unit?}] }
  //   async fit(snapshot, ctx)            -> state (serialisable)
  //   async forecast(state, plan, ctx)    -> { grid, variables:{v:{q05,q25,q50,q75,q95}}, outcomes:{id:{value,lo,hi}}, signals:{id:{grid,values}}, simEvents:[] , extra:{} }
  // ctx = { asOfDay, horizon, nParticles, seed, signals:{id:{grid,values}}, snapshot }
  S.validateModule = function (m) {
    const errs = [];
    if (!m || !m.meta || !m.meta.id) errs.push('meta.id missing');
    if (typeof m.fit !== 'function') errs.push('fit() missing');
    if (typeof m.forecast !== 'function') errs.push('forecast() missing');
    if (!m.consumes) errs.push('consumes missing');
    if (!m.provides) errs.push('provides missing');
    return errs;
  };

  // ---- phenotype rule (mirrors Aim2 Code/config.R thresholds) --------------
  S.phenotypeRule = { recover_thr: 1.0, recover_run: 3, early_window: [0, 28], dip_window: [21, 60], dip_thr: 1.0, dip_run: 2, persist_window: [45, 75], persist_thr: 0.5 };
  // daily grid days[] with anc[] -> 'quick' | 'second_wave' | 'persistent'
  S.classifyPhenotype = function (days, anc, rule) {
    rule = rule || S.phenotypeRule;
    const inWin = (d, w) => d >= w[0] && d <= w[1];
    const pw = [];
    for (let i = 0; i < days.length; i++) if (inWin(days[i], rule.persist_window)) pw.push(anc[i]);
    if (pw.length && GL.util.median(pw) < rule.persist_thr) return 'persistent';
    // recovery: `recover_run` consecutive days >= thr inside early window
    let run = 0, recoveredAt = null;
    for (let i = 0; i < days.length; i++) {
      if (!inWin(days[i], rule.early_window)) { continue; }
      run = anc[i] >= rule.recover_thr ? run + 1 : 0;
      if (run >= rule.recover_run) { recoveredAt = days[i]; break; }
    }
    if (recoveredAt === null) return 'quick'; // slow but monotone: counted with quick (no second wave)
    run = 0;
    for (let i = 0; i < days.length; i++) {
      if (days[i] <= recoveredAt || !inWin(days[i], rule.dip_window)) continue;
      run = anc[i] < rule.dip_thr ? run + 1 : 0;
      if (run >= rule.dip_run) return 'second_wave';
    }
    return 'quick';
  };
})(typeof window !== 'undefined' ? window : global);
