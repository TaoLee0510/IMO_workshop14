/* GreenLight · modules/marrow4.js
 * Aim 4 marrow module — "who & when": the parent model of Aim 4, v1.2 (stem-cell
 * pool with a two-clone layer -> MPP -> CLP / CMP -> GMP (neutrophils, monocytes)
 * / MEP (platelets), CLP -> lymphocytes; chemotherapy exposure a, CAR T-phase
 * driver b from ferritin)
 * running in the browser through js/modules/aim4_core.js, with a particle
 * posterior over the five patient parameters
 *   r_s (stem-cell regrowth = reserve), th_ap (chemotherapy hit on progenitors),
 *   th_as (chemotherapy hit on stem cells), th_bs (CAR T-phase hit on stem cells),
 *   f0 (clone fraction)
 * Population constants come from data/aim4_pop.js (window.GL_AIM4_POP): prior
 * centres until step G of Aim4/Code exports the fitted population.
 * Outcomes: recovery phenotype (quick / second wave / delayed / persistent) and
 * two counterfactual readouts at the stem-cell level - when (a / b / a+b) and
 * why (reserve / clone) - as in Aim4/Code/R/attribution.R; the mean time and number
 * of divisions since MPP entry of the circulating cells (model moments); forecasts of
 * ANC, PLT, AMC, ALC (lymphoid branch, not fitted to ALC) and ferritin.
 */
(function (root) {
  const GL = root.GL = root.GL || {};
  const U = GL.util, S = GL.schema, C = GL.aim4;
  GL.modules = GL.modules || {};

  const POP = () => root.GL_AIM4_POP;
  const spec = () => C.spec(POP().model);
  const SIGMA_FALLBACK = { anc: 0.3, plt: 0.25, amc: 0.35 };
  const T_END = 120, GRID_FROM = -10;

  // ---- normal log-CDF (bounded observations) ---------------------------------------------
  // log erfc(x) for x >= 0 (Numerical Recipes erfcc, in log form: no underflow)
  function logErfc(x) {
    const t = 1 / (1 + 0.5 * x);
    return Math.log(t) - x * x - 1.26551223 + t * (1.00002368 + t * (0.37409196 + t * (0.09678418 + t * (-0.18628806 + t * (0.27886807 + t * (-1.13520398 + t * (1.48851587 + t * (-0.82215223 + t * 0.17087277))))))));
  }
  function logPhi(z) {
    if (z < 0) return Math.log(0.5) + logErfc(-z / Math.SQRT2);
    return Math.log(1 - 0.5 * Math.exp(logErfc(z / Math.SQRT2)));
  }

  // ---- patient information ----------------------------------------------------------------
  // Real-cohort records carry aim4 = {egfr, bridging_class, apheresis_day, t0, base, ferritin_base, cov};
  // other records get defaults from what the snapshot shows.
  function patientInfo(snapshot) {
    const a = snapshot.aim4 || {};
    const pop = POP();
    const obs = snapshot.observations.filter(o => !o.excluded && o.value > 0);
    const by = (v) => obs.filter(o => o.variable === v).sort((x, y) => x.day - y.day);
    const setpoint = (v, dflt) => {
      if (a.base && a.base[v] > 0) return a.base[v];
      const ok = by(v).filter(o => !o.cens && !o.upper);
      const pre = ok.filter(o => o.day <= -6);
      if (pre.length) return U.median(pre.map(o => o.value));
      const early = ok.filter(o => o.day <= -3);
      if (early.length) return U.median(early.slice(0, 2).map(o => o.value));
      return dflt;
    };
    const fer = by('ferritin');
    let ferBase = a.ferritin_base;
    if (!(ferBase > 0)) { const pre = fer.filter(o => o.day <= -1); ferBase = pre.length ? U.median(pre.map(o => o.value)) : fer.length ? fer[0].value : 500; }
    const firstDay = obs.length ? Math.min(...obs.map(o => o.day)) : -31;
    const h = snapshot.history || {}, dem = snapshot.demographics || {}, cc = pop.covariates;
    const cov = a.cov || {
      age_c: ((dem.age || cc.age_ref) - cc.age_ref) / cc.age_scale,
      prior_lines_c: (h.prior_lines === undefined ? cc.prior_lines_median : h.prior_lines) - cc.prior_lines_median,
      prior_asct: h.prior_asct ? 1 : 0, log_rdw_c: 0, mcv_c: 0
    };
    return {
      base: { anc: setpoint('anc', 3.5), plt: setpoint('plt', 200), amc: setpoint('amc', 0.55), alc: setpoint('alc', 1.0) },
      ferBase, ferObs: fer,
      bridging: a.bridging_class || (snapshot.bridging === 'chemo' ? 'cy' : 'none'),
      egfr: a.egfr > 0 ? a.egfr : 90,
      apheresis: a.apheresis_day === undefined || a.apheresis_day === null ? -27 : a.apheresis_day,
      t0: a.t0 !== undefined && a.t0 !== null ? a.t0 : Math.min(-31, Math.floor(firstDay) - 1),
      cov, rEngine: a.r_engine || null
    };
  }

  // inputs: recorded LD doses once LD has started, plan doses before; recorded supportive care
  function ldScale(snapshot, plan) {
    if (snapshot.asOfDay >= -5 && snapshot.ld) {
      const l = snapshot.ld;
      if (l.flu_total !== undefined && l.cy_total !== undefined) return 0.5 * l.flu_total / 90 + 0.5 * l.cy_total / 1500;
      return l.dose_reduced ? 0.75 : 1;
    }
    return 0.5 * plan.stages.ld.flu_frac + 0.5 * plan.stages.ld.cy_frac;
  }
  function bridgingClass(info, plan, asOf) {
    const ty = plan.stages.bridging.type;
    if (ty === 'as_recorded' || S.stageLocked('bridging', asOf)) return info.bridging;
    return ty === 'chemo' ? (info.bridging === 'none' ? 'cy' : info.bridging) : 'none';
  }
  function designFor(info, snapshot, plan, asOf) {
    const gcsf = [], tpo = [], cd34 = [];
    for (const e of snapshot.events) {
      const s = Math.round(e.day_start), en = Math.round(e.day_end === undefined || e.day_end === null ? e.day_start : e.day_end);
      if (e.type === 'gcsf') gcsf.push([s, en + 1]);
      if (e.type === 'tpo_ra') tpo.push([s, Math.max(en, s + 60)]);
      if (e.type === 'cd34_boost') cd34.push({ day: s, add: POP().cd34_add || 0.3 });
    }
    return { t0: info.t0, egfr: info.egfr, bridging: bridgingClass(info, plan, asOf), apheresis_day: info.apheresis,
      ld_scale: ldScale(snapshot, plan), gcsf, tpo, cd34 };
  }

  // anti-inflammatory policies scale the CAR T-phase driver after the decision day
  // (illustrative factors; the standard of care - tocilizumab at grade 2, steroids at grade 3 - is 1)
  function crsFactor(crs) {
    let f = 1;
    if (crs.toci_grade === 1) f *= 0.8; else if (crs.toci_grade === 3) f *= 1.1;
    if (crs.steroid_grade === 2) f *= 0.9;
    if (crs.anakinra_prophylaxis) f *= 0.6;
    return f;
  }

  // ---- CAR T-phase driver ------------------------------------------------------------------
  function kernelDraw(fitted, rng, jitter) {
    const kp = POP().kernel_pop;
    if (fitted) return { A: fitted.A * Math.exp(jitter * rng.normal()), tp: U.clamp(fitted.tp * Math.exp(jitter * rng.normal()), 2, 25) };
    const z1 = rng.normal(), z2 = rng.normal();
    const c = kp.cov, l11 = Math.sqrt(c[0][0]), l21 = c[1][0] / l11, l22 = Math.sqrt(Math.max(c[1][1] - l21 * l21, 1e-6));
    return { A: Math.exp(kp.mu[0] + l11 * z1), tp: U.clamp(Math.exp(kp.mu[1] + l21 * z1 + l22 * z2), 2, 25) };
  }
  function driverFor(info, asOf, kernel) {
    const sm = C.smoothFerritin(info.ferObs, info.ferBase, asOf);
    return C.buildDriver({ t0: info.t0, tEnd: 130, smooth: sm, upTo: asOf, kernel, h: POP().kernel_pop.h });
  }
  function kernelMedian() { const kp = POP().kernel_pop; return { A: Math.exp(kp.mu[0]), tp: Math.exp(kp.mu[1]) }; }

  // ---- prior particles -----------------------------------------------------------------------
  const NAT = { r_s: Math.exp, th_ap: Math.exp, th_as: Math.exp, th_bs: Math.exp, f0: (x) => 1 / (1 + Math.exp(-x)) };
  // f0: population prior, or 2 x the largest pre-treatment VAF when a CH panel is on record
  // (demo records; the real cohort has no sequencing)
  function priorParticles(info, n, seed, mutations) {
    const pop = POP(), sp = spec(), r = U.rng(seed);
    const mut = mutations || {}, vafs = (mut.variants || []).map(v => v.vaf).filter(x => x > 0);
    const parts = [];
    for (let i = 0; i < n; i++) {
      const p = { r_s: 0, f0: 0, th_as: 0, th_ap: 0, th_bs: 0 };
      for (const k of sp.re) {
        const b = (pop.beta && pop.beta[k]) || {};
        let m = pop.mu[k];
        for (const c in b) m += b[c] * (info.cov[c] || 0);
        p[k] = NAT[k](m + pop.omega[k] * r.normal());
      }
      if (sp.re.includes('f0') && vafs.length) p.f0 = U.clamp(2 * Math.max(...vafs) * Math.exp(0.35 * r.normal()), 0.005, 0.9);
      else if (sp.re.includes('f0') && mut.tested) p.f0 = U.clamp(p.f0, 1e-4, 0.02);
      parts.push(p);
    }
    return parts;
  }
  function jitterSpec() {
    const sp = spec(), js = {};
    const all = { r_s: { log: true, sd: 0.15, lo: 0.002, hi: 2 }, th_ap: { log: true, sd: 0.15, lo: 0.01, hi: 40 },
      th_as: { log: true, sd: 0.2, lo: 1e-4, hi: 20 }, th_bs: { log: true, sd: 0.2, lo: 1e-5, hi: 5 }, f0: { log: true, sd: 0.25, lo: 1e-4, hi: 0.95 } };
    sp.re.forEach(k => { js[k] = all[k]; });
    return js;
  }

  // ---- likelihood -------------------------------------------------------------------------------
  const LIK_VARS = ['anc', 'plt', 'amc'];
  function loglik(sim, info, obs, sigma) {
    let ll = 0, n = 0;
    const t0 = sim.days[0];
    for (const o of obs) {
      const i = Math.round(o.day) - t0; if (i < 0 || i >= sim.days.length) continue;
      const rel = o.variable === 'anc' ? sim.anc_rel[i] : o.variable === 'plt' ? sim.plt_rel[i] : sim.amc_rel[i];
      const pred = Math.max(rel * info.base[o.variable], 1e-8), s = sigma[o.variable];
      const z = (Math.log(o.value) - Math.log(pred)) / s;
      ll += (o.cens || o.upper) ? logPhi(z) : -0.5 * z * z - Math.log(s);
      n++;
    }
    return { ll, n };
  }

  // ---- forecasting helpers --------------------------------------------------------------------
  // supportive-care policies of the plan, applied online to days after asOf
  function makeHook(plan, info, asOf, pol, struct) {
    const cy = plan.stages.cytopenia, add = POP().cd34_add || 0.3;
    return function (day, y, u) {
      if (day <= asOf) return;
      const anc = y[C.IDX.BN] * info.base.anc, plt = y[C.IDX.BP] * info.base.plt;
      if (day >= 0) {
        if (cy.gcsf === 'prophylactic' && day >= cy.gcsf_start_day && day <= 60 && day >= pol.gcsfUntil && anc < 1.0 && pol.gcsfCourses < 6) { pol.gcsfUntil = day + cy.gcsf_days; pol.gcsfCourses++; }
        if (cy.gcsf === 'reactive' && day >= cy.gcsf_start_day && day <= 90 && day >= pol.gcsfUntil && anc < cy.gcsf_anc_thr && pol.gcsfCourses < 5) { pol.gcsfUntil = day + cy.gcsf_days; pol.gcsfCourses++; }
        if (cy.tpo_ra && day >= 21 && plt < cy.tpo_plt_thr) pol.tpoFrom = Math.min(pol.tpoFrom, day);
        if (day >= 3 && plt < cy.plt_tx_thr && day - pol.lastPltTx >= 3) { pol.lastPltTx = day; pol.pltTx++; }
        if (cy.cd34_backup && cy.cd34_boost === 'if_persistent' && day === cy.cd34_boost_day && anc < 0.5 && !pol.boosted) { y[0] += add; pol.boosted = true; pol.boosts++; }
      }
      if (day < pol.gcsfUntil) u.gcsf = 1;
      if (day >= pol.tpoFrom) u.tpo = 1;
    };
  }

  function observedAnc(snapshot) {
    return snapshot.observations.filter(o => o.variable === 'anc' && !o.excluded).map(o => ({ day: o.day, value: o.value, cens: !!o.cens }));
  }

  function trapz(x, y, from) { let s = 0; for (let i = 1; i < x.length; i++) if (x[i - 1] >= from) s += (x[i] - x[i - 1]) * (y[i] + y[i - 1]) / 2; return s; }

  // counterfactual readouts on latent trajectories (no supportive care), as attribution.R
  function attribution(state, snapshot, info) {
    const pop = POP(), sp = spec(), st = pop.struct, asOf = state.asOfDay;
    const m = Math.min(40, state.particles.length);
    const parts = GL.smc.resample(state.particles, state.weights, m, state.seed + 31);
    const rng = U.rng(state.seed + 37);
    const design = designFor(info, snapshot, S.defaultPlan(), asOf);
    const rMed = Math.exp(pop.mu.r_s);
    const out = { when: {}, who: {}, f_a: [], f_b: [], f30: [], S30: [], n: 0, age: {}, gen: {}, age0: {}, gen0: {} };
    const LINES = ['anc', 'plt', 'amc', 'alc'];
    LINES.forEach(v => { out.age[v] = []; out.gen[v] = []; });
    const dayNow = Math.max(-25, Math.min(90, Math.round(asOf)));
    const cls = (p, drv, track) => {
      const sim = C.simulate(C.buildPar(p, st, sp), { design, struct: st, drv, tEnd: 90, track });
      if (!sim) return null;
      const series = [];
      for (let i = 0; i < sim.days.length; i++) if (sim.days[i] >= 0) series.push({ day: sim.days[i], value: sim.anc_rel[i] * info.base.anc });
      return { ph: C.phenotype(series, true), sim };
    };
    for (const p of parts) {
      const drv = driverFor(info, asOf, kernelDraw(state.kernel, rng, 0.15));
      const act = cls(p, drv, true); if (!act) continue;
      const noA = cls({ ...p, th_as: 0 }, drv).ph, noBS = cls({ ...p, th_bs: 0 }, drv).ph, noABS = cls({ ...p, th_as: 0, th_bs: 0 }, drv).ph;
      const noC = p.f0 > 0 ? cls({ ...p, f0: 0 }, drv).ph : act.ph, res = cls({ ...p, r_s: rMed }, drv).ph;
      const w = C.whenClass(act.ph, noA, noBS, noABS), h = C.whoClass(act.ph, noC, res) || 'none';
      out.when[w] = (out.when[w] || 0) + 1; out.who[h] = (out.who[h] || 0) + 1;
      const s = act.sim, i30 = s.days.indexOf(30);
      out.f_a.push(1 - Math.exp(-(sp.re.includes('th_as') ? p.th_as : 0) * trapz(s.days, s.C, -Infinity)));
      out.f_b.push(1 - Math.exp(-(sp.re.includes('th_bs') ? p.th_bs : 0) * trapz(s.days, s.I, 0)));
      out.f30.push(s.f_clone[i30]); out.S30.push(s.S[i30]);
      // mean time / divisions since MPP entry of the circulating cells on the viewing day, and pre-treatment
      const iNow = s.days.indexOf(dayNow);
      LINES.forEach(v => { out.age[v].push(s['age_' + v][iNow]); out.gen[v].push(s['gen_' + v][iNow]); out.age0[v] = s['age_' + v][0]; out.gen0[v] = s['gen_' + v][0]; });
      out.n++;
    }
    const pr = (o, k) => out.n ? (o[k] || 0) / out.n : null;
    return {
      n: out.n,
      when: { none: pr(out.when, 'none'), a: pr(out.when, 'a'), b: pr(out.when, 'b'), ab: (pr(out.when, 'ab_joint') || 0) + (pr(out.when, 'ab_redundant') || 0), other: pr(out.when, 'other') },
      who: { clone: pr(out.who, 'clone'), reserve: pr(out.who, 'reserve'), both: pr(out.who, 'both'), other: pr(out.who, 'other'), none: pr(out.who, 'none') },
      f_a: GL.smc.val(out.f_a), f_b: GL.smc.val(out.f_b), f30: GL.smc.val(out.f30), S30: GL.smc.val(out.S30),
      sinceMpp: { day: dayNow, lines: LINES.map(v => ({ id: v, age: GL.smc.val(out.age[v]), gen: GL.smc.val(out.gen[v]), age0: out.age0[v], gen0: out.gen0[v] })) }
    };
  }

  // ---- module -----------------------------------------------------------------------------------
  const M = {
    meta: { id: 'marrow4', nameKey: 'mod.marrow4.name', descKey: 'mod.marrow4.desc', owner: 'Aim 4 · Team Green', version: '1.2', kind: 'local', placeholder: false, engine: 'JS port of the Aim 4 C++ model v1.2 (RK45, 24 states + age/generation moments) · particle posterior' },
    consumes: { variables: ['anc', 'plt', 'amc', 'ferritin'], planKeys: ['bridging.type', 'ld.flu_frac', 'ld.cy_frac', 'crs.toci_grade', 'crs.steroid_grade', 'crs.anakinra_prophylaxis', 'cytopenia.*'], signals: [] },
    provides: {
      signals: ['inflammation'],
      variables: ['anc', 'plt', 'amc', 'alc', 'ferritin'],
      params: [
        { id: 'r_s', labelKey: 'param.r_s', unit: '/day', log: true, lo: 0.005, hi: 0.6, direction: 'higher', descKey: 'help.param.r_s' },
        { id: 'th_ap', labelKey: 'param.th_ap', unit: '/exposure', log: true, lo: 0.05, hi: 20, direction: 'lower', descKey: 'help.param.th_ap' },
        { id: 'th_as', labelKey: 'param.th_as', unit: '/exposure', log: true, lo: 0.002, hi: 5, direction: 'lower', descKey: 'help.param.th_as' },
        { id: 'th_bs', labelKey: 'param.th_bs', unit: '/(log-ratio·day)', log: true, lo: 0.0002, hi: 0.5, direction: 'lower', descKey: 'help.param.th_bs' },
        { id: 'f0', labelKey: 'param.f0', unit: 'fraction', lo: 0, hi: 0.8, direction: 'lower', descKey: 'help.param.f0' }
      ],
      outcomes: [
        { id: 'p_quick', labelKey: 'out.p_quick', kind: 'prob', direction: 'higher' },
        { id: 'p_second_wave', labelKey: 'out.p_second_wave', kind: 'prob', direction: 'lower' },
        { id: 'p_delayed', labelKey: 'out.p_delayed', kind: 'prob', direction: 'lower' },
        { id: 'p_persistent', labelKey: 'out.p_persistent', kind: 'prob', direction: 'lower' },
        { id: 'p_a_driven', labelKey: 'out.p_a_driven', kind: 'prob', direction: 'lower' },
        { id: 'p_b_driven', labelKey: 'out.p_b_driven', kind: 'prob', direction: 'lower' },
        { id: 'p_ab_driven', labelKey: 'out.p_ab_driven', kind: 'prob', direction: 'lower' },
        { id: 'p_clone', labelKey: 'out.p_clone', kind: 'prob', direction: 'lower' },
        { id: 'p_reserve', labelKey: 'out.p_reserve', kind: 'prob', direction: 'lower' },
        { id: 'f_a', labelKey: 'out.f_a', kind: 'value', direction: 'lower', unit: 'fraction' },
        { id: 'f_b', labelKey: 'out.f_b', kind: 'value', direction: 'lower', unit: 'fraction' },
        { id: 'f_clone30', labelKey: 'out.f_clone30', kind: 'value', direction: 'lower', unit: 'fraction' },
        { id: 'age_anc', labelKey: 'out.age_anc', kind: 'value', direction: 'neutral', unit: 'd' },
        { id: 'gen_anc', labelKey: 'out.gen_anc', kind: 'value', direction: 'neutral', unit: 'n' },
        { id: 'sev_neut_days', labelKey: 'out.sev_neut_days', kind: 'value', direction: 'lower', unit: 'd' },
        { id: 'anc_nadir_late', labelKey: 'out.anc_nadir_late', kind: 'value', direction: 'higher', unit: '×10⁹/L' },
        { id: 'p_plt20', labelKey: 'out.p_plt20', kind: 'prob', direction: 'lower' },
        { id: 'gcsf_courses', labelKey: 'out.gcsf_courses', kind: 'value', direction: 'lower', unit: 'n' },
        { id: 'plt_tx', labelKey: 'out.plt_tx', kind: 'value', direction: 'lower', unit: 'n' }
      ]
    },
    patientInfo, logPhi,

    async fit(snapshot, ctx) {
      const pop = POP(); if (!pop) throw new Error('data/aim4_pop.js not loaded');
      const n = Math.min((ctx && ctx.nParticles) || 400, 400);
      const seed = (ctx && ctx.seed) || U.hashString(snapshot.patient_id + ':aim4');
      const info = patientInfo(snapshot), sp = spec(), st = pop.struct, asOf = snapshot.asOfDay;
      const plan = S.defaultPlan();
      const design = designFor(info, snapshot, plan, asOf);
      const kernel = C.fitKernel(info.ferObs, info.ferBase, asOf, pop.kernel_pop.h);
      const drv = driverFor(info, asOf, kernel || kernelMedian());
      const tEnd = Math.max(Math.ceil(asOf), info.t0 + 2);
      const obs = snapshot.observations.filter(o => LIK_VARS.includes(o.variable) && !o.excluded && o.value > 0 && o.day > info.t0 && o.day <= tEnd);
      const sigma = Object.assign({}, SIGMA_FALLBACK, pop.sigma);
      const simLL = (p) => {
        const sim = C.simulate(C.buildPar(p, st, sp), { design, struct: st, drv, tEnd });
        return sim ? loglik(sim, info, obs, sigma) : { ll: -1e10, n: obs.length };
      };
      const fit = GL.smc.fit(priorParticles(info, n, seed, snapshot.mutations), simLL, jitterSpec(), seed);
      const params = {};
      for (const k of ['r_s', 'th_ap', 'th_as', 'th_bs', 'f0']) params[k] = GL.smc.summarise(fit.particles, fit.weights, k);
      return { moduleId: 'marrow4', model: pop.model, popSource: pop.source, asOfDay: asOf, seed, n, nObs: fit.nObs, beta: fit.beta, ess: fit.ess,
        particles: fit.particles, weights: fit.weights, params, logEvidence: fit.logEvidence, kernel, base: info.base };
    },

    async forecast(state, plan, ctx) {
      const pop = POP(), sp = spec(), st = pop.struct;
      const m = Math.min((ctx && ctx.nForecast) || 60, 60), horizon = Math.min((ctx && ctx.horizon) || T_END, T_END);
      const snapshot = ctx.snapshot, asOf = state.asOfDay;
      const info = patientInfo(snapshot);
      const parts = GL.smc.resample(state.particles, state.weights, m, state.seed + 7);
      const rng = U.rng(state.seed + 11);
      const design = designFor(info, snapshot, plan, asOf);
      const dScale = { after: Math.max(asOf, 0), factor: crsFactor(plan.stages.crs) };
      const obsAnc = observedAnc(snapshot);
      const recs = [], ph = { quick: 0, second_wave: 0, delayed: 0, persistent: 0 };
      const sev = [], nadir = [], plt20 = [], gcsf = [], ptx = [];
      for (const p of parts) {
        const drv = driverFor(info, asOf, kernelDraw(state.kernel, rng, state.kernel ? 0.15 : 0));
        const pol = { gcsfUntil: -1e9, gcsfCourses: 0, tpoFrom: 1e9, lastPltTx: -1e9, pltTx: 0, boosted: false, boosts: 0 };
        const sim = C.simulate(C.buildPar(p, st, sp), { design, struct: st, drv, tEnd: horizon, dScale, hook: makeHook(plan, info, asOf, pol, st) });
        if (!sim) continue;
        const r = { days: [], anc: [], plt: [], amc: [], alc: [], fer: [], D: [] };
        for (let i = 0; i < sim.days.length; i++) {
          const d = sim.days[i]; if (d < GRID_FROM) continue;
          r.days.push(d); r.anc.push(sim.anc_rel[i] * info.base.anc); r.plt.push(sim.plt_rel[i] * info.base.plt);
          r.amc.push(sim.amc_rel[i] * info.base.amc); r.alc.push(sim.alc_rel[i] * info.base.alc);
          r.D.push(sim.D[i]); r.fer.push(info.ferBase * Math.exp(Math.max(sim.D[i], 0)));
        }
        // phenotype: observed counts up to asOf, model afterwards (same rule as the R pipeline)
        const series = obsAnc.slice();
        for (let i = 0; i < r.days.length; i++) if (r.days[i] > asOf) series.push({ day: r.days[i], value: r.anc[i] });
        const cls = C.phenotype(series, true) || 'quick';
        ph[cls]++;
        let s = 0, nd = Infinity, p20 = 0;
        for (let i = 0; i < r.days.length; i++) {
          const d = r.days[i];
          if (d >= 0 && d <= 90 && r.anc[i] < 0.5) s++;
          if (d >= 14 && d <= 60 && r.anc[i] < nd) nd = r.anc[i];
          if (d >= 10 && d <= 90 && r.plt[i] < 20) p20 = 1;
        }
        sev.push(s); nadir.push(nd); plt20.push(p20); gcsf.push(pol.gcsfCourses); ptx.push(pol.pltTx);
        recs.push(r);
      }
      if (!recs.length) throw new Error('marrow4: all simulations failed');
      const k = recs.length, grid = recs[0].days;
      const sg = Object.assign({}, SIGMA_FALLBACK, pop.sigma);
      const q = (key, sd) => GL.smc.predictiveBands(grid, recs.map(r => r[key]), sd);
      // ALC: lymphoid branch (not fitted to ALC; observed ALC after infusion includes CAR T cells)
      const variables = { anc: q('anc', sg.anc), plt: q('plt', sg.plt), amc: q('amc', sg.amc), alc: q('alc', 0.35), ferritin: q('fer', 0.25) };
      if (!state.attr) state.attr = attribution(state, snapshot, info);
      const at = state.attr, val = GL.smc.val;
      const outcomes = {
        p_quick: { value: ph.quick / k }, p_second_wave: { value: ph.second_wave / k }, p_delayed: { value: ph.delayed / k }, p_persistent: { value: ph.persistent / k },
        p_a_driven: { value: at.when.a }, p_b_driven: { value: at.when.b }, p_ab_driven: { value: at.when.ab },
        p_clone: { value: (at.who.clone || 0) + (at.who.both || 0) }, p_reserve: { value: (at.who.reserve || 0) + (at.who.both || 0) },
        f_a: at.f_a, f_b: at.f_b, f_clone30: at.f30,
        age_anc: at.sinceMpp.lines[0].age, gen_anc: at.sinceMpp.lines[0].gen,
        sev_neut_days: val(sev), anc_nadir_late: val(nadir), p_plt20: { value: U.mean(plt20) }, gcsf_courses: val(gcsf), plt_tx: val(ptx)
      };
      const dMed = grid.map((_, i) => U.median(recs.map(r => r.D[i])));
      return { grid, variables, outcomes, signals: { inflammation: { grid, values: dMed } },
        extra: { attribution: at, model: pop.model, popSource: pop.source, rEngine: info.rEngine } };
    }
  };

  GL.modules.marrow4 = M;
  if (GL.registry) GL.registry.register(M);
})(typeof window !== 'undefined' ? window : global);
