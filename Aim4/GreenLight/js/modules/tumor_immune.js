/* GreenLight · modules/tumor_immune.js
 * PLACEHOLDER tumor–immune module. Stands in for the ODE model being built
 * by the other Aim (tumor–CAR T interplay). Variable names and equations are
 * illustrative and will be replaced when that model's contract arrives:
 *   E  CAR T effector cells (×10³ copies/µg)
 *        dE/dt = rho·e^(−t/τe)·E·T/(T+KT)·(1 − E/Emax) − mu·E   (antigen-driven, capacity-limited, waning)
 *   T  tumour burden (mL metabolic tumour volume)
 *        dT/dt = r·T·(1 − T/Tmax) − kappa·T·E/(Es + E)           (saturating kill)
 * Patient parameters (particle posterior): rho (expansion rate), Emax (expansion capacity), r (tumour growth), T0, e0 (dose scale).
 * Provides signal car_expansion(t) = E(t)/E_ref consumed by the marrow module.
 */
(function (root) {
  const GL = root.GL = root.GL || {};
  const U = GL.util, S = GL.schema;
  GL.modules = GL.modules || {};

  const POP = { KT: 3, mu: 0.12, tauE: 9, Tmax: 2000, kappa: 0.3, Es: 2, Eref: 5, E0base: 0.02, ldBoost: 0.6, bridgingReduce: { none: 1, as_recorded: 1, chemo: 0.6, radiation: 0.5 }, sterRho: 0.8, respThr: 0.1, crThr: 1.0 };

  function simulate(p, snapshot, plan, ctx, tEnd) {
    const asOf = snapshot.asOfDay;
    // locked stages use the record, open stages use the plan
    const ldOpen = asOf <= -5, infOpen = asOf < 0, brOpen = asOf < -6;
    let ldFrac = 1;
    if (!ldOpen && snapshot.ld) ldFrac = snapshot.ld.flu_total !== undefined ? 0.5 * (snapshot.ld.flu_total / 90 + snapshot.ld.cy_total / 1500) : (snapshot.ld.dose_reduced ? 0.75 : 1);
    else ldFrac = 0.5 * (plan.stages.ld.flu_frac + plan.stages.ld.cy_frac);
    const doseFrac = infOpen ? plan.stages.infusion.car_dose_frac : (snapshot.car_dose_frac || 1);
    const brType = brOpen ? plan.stages.bridging.type : (snapshot.bridging || 'none');
    const brRed = POP.bridgingReduce[brType === 'as_recorded' ? (snapshot.bridging || 'none') : brType] || 1;
    const rho = p.rho * (1 - POP.ldBoost + POP.ldBoost * ldFrac);
    const E0 = POP.E0base * doseFrac * p.e0;
    const T0 = p.T0 * (brOpen ? brRed : 1); // recorded bridging already reflected in baseline measurement
    // steroid use (recorded) dampens expansion
    const ster = snapshot.events.filter(e => e.type === 'steroid').map(e => [e.day_start, (e.day_end || e.day_start) + 1]);
    const y0 = new Float64Array([0, T0]);
    const rec = { days: [], E: [], T: [] };
    const rhs = function (t, y, dy) {
      const onSter = ster.some(w => t >= w[0] && t < w[1]);
      const rh = t >= 0 ? rho * Math.exp(-t / POP.tauE) * (onSter ? POP.sterRho : 1) : 0;
      const E = y[0], T = y[1];
      dy[0] = rh * E * T / (T + POP.KT) * (1 - E / p.Emax) - POP.mu * E;
      dy[1] = p.r * T * (1 - T / POP.Tmax) - POP.kappa * T * E / (POP.Es + E);
    };
    const onStep = function (t, y) {
      if (Math.abs(t - Math.round(t)) < 1e-9) {
        const d = Math.round(t);
        if (d === 0) y[0] = E0; // infusion
        rec.days.push(d); rec.E.push(y[0]); rec.T.push(y[1]);
      }
    };
    GL.ode.integrate({ t0: -6, t1: tEnd, dt: 0.25, y0, rhs, onStep });
    rec.T0 = T0;
    return rec;
  }

  function priorParticles(snapshot, n, seed) {
    const r = U.rng(seed + 11);
    const tb = S.obsByVar(snapshot, 'tumor_burden').filter(o => o.day <= 0);
    const T0med = tb.length ? tb[tb.length - 1].value : 60;
    const parts = [];
    for (let i = 0; i < n; i++) {
      parts.push({ rho: r.lognormal(2.0, 0.2), Emax: r.lognormal(6, 0.8), r: r.lognormal(0.03, 0.6), T0: tb.length ? T0med * Math.exp(0.15 * r.normal()) : r.lognormal(60, 0.9), e0: r.lognormal(1, 0.5) });
    }
    return parts;
  }
  const SIG = { car_copies: 0.6, tumor_burden: 0.45 };
  function loglik(rec, snapshot) {
    let ll = 0, n = 0; const idx = {}; rec.days.forEach((d, i) => { idx[d] = i; });
    for (const o of snapshot.observations) {
      if (!(o.variable in SIG) || o.day < 0 && o.variable === 'car_copies') continue;
      const i = idx[o.day]; if (i === undefined) continue;
      const pred = o.variable === 'car_copies' ? rec.E[i] : rec.T[i];
      const z = (Math.log(o.value + 0.05) - Math.log(pred + 0.05)) / SIG[o.variable];
      ll += -0.5 * z * z; n++;
    }
    return { ll, n };
  }
  const JITTER = { rho: { log: true, sd: 0.1, lo: 0.1, hi: 3 }, Emax: { log: true, sd: 0.2, lo: 0.1, hi: 200 }, r: { log: true, sd: 0.2, lo: 0.002, hi: 0.4 }, T0: { log: true, sd: 0.1, lo: 0.5, hi: 1900 }, e0: { log: true, sd: 0.15, lo: 0.1, hi: 10 } };
  const M = {
    meta: { id: 'tumor_immune', nameKey: 'mod.ti.name', descKey: 'mod.ti.desc', owner: 'Tumour–immune Aim (TBD)', version: '0.1-placeholder', kind: 'local', placeholder: true, engine: 'JS RK4 placeholder · particle posterior' },
    consumes: { variables: ['car_copies', 'tumor_burden'], planKeys: ['bridging.type', 'ld.flu_frac', 'ld.cy_frac', 'infusion.car_dose_frac'], signals: [] },
    provides: {
      signals: ['car_expansion'],
      params: [
        { id: 'rho', labelKey: 'param.rho', unit: 'd⁻¹', log: true, lo: 0.2, hi: 3, direction: 'higher', descKey: 'param.rho.desc' },
        { id: 'Emax', labelKey: 'param.Emax', unit: '×10³ copies/µg', log: true, lo: 0.1, hi: 200, direction: 'higher', descKey: 'param.Emax.desc' },
        { id: 'r', labelKey: 'param.r', unit: 'd⁻¹', log: true, lo: 0.003, hi: 0.3, direction: 'lower', descKey: 'param.r.desc' },
        { id: 'T0', labelKey: 'param.T0', unit: 'mL', log: true, lo: 1, hi: 1500, direction: 'lower', descKey: 'param.T0.desc' }
      ],
      outcomes: [
        { id: 'p_response', labelKey: 'out.p_response', kind: 'prob', direction: 'higher' },
        { id: 'p_cr', labelKey: 'out.p_cr', kind: 'prob', direction: 'higher' },
        { id: 'peak_expansion', labelKey: 'out.peak_expansion', kind: 'value', direction: 'higher', unit: '×10³ copies/µg' },
        { id: 'tumor_d90', labelKey: 'out.tumor_d90', kind: 'value', direction: 'lower', unit: '% of baseline' }
      ],
      variables: ['car_copies', 'tumor_burden']
    },
    POP, simulate, priorParticles,
    async fit(snapshot, ctx) {
      const n = (ctx && ctx.nParticles) || 400;
      const seed = (ctx && ctx.seed) || U.hashString(snapshot.patient_id);
      const prior = priorParticles(snapshot, n, seed);
      const plan = S.defaultPlan();
      const tEnd = Math.max(snapshot.asOfDay, -5);
      const simLL = (p) => loglik(simulate(p, snapshot, plan, ctx || {}, tEnd), snapshot);
      const fit = GL.smc.fit(prior, simLL, JITTER, seed);
      const params = {}; for (const k of ['rho', 'Emax', 'r', 'T0', 'e0']) params[k] = GL.smc.summarise(fit.particles, fit.weights, k);
      return { moduleId: 'tumor_immune', asOfDay: snapshot.asOfDay, seed, n, nObs: fit.nObs, beta: fit.beta, ess: fit.ess, particles: fit.particles, weights: fit.weights, params, logEvidence: fit.logEvidence };
    },
    async forecast(state, plan, ctx) {
      const m = (ctx && ctx.nForecast) || 80, horizon = (ctx && ctx.horizon) || 120;
      const snapshot = ctx.snapshot;
      const parts = GL.smc.resample(state.particles, state.weights, m, state.seed + 13);
      const recs = parts.map(p => simulate(p, snapshot, plan, ctx, horizon));
      const grid = recs[0].days;
      const q = (arrs, sd) => GL.smc.predictiveBands(grid, arrs, sd);
      const variables = { car_copies: q(recs.map(r => r.E), 0.35), tumor_burden: q(recs.map(r => r.T), 0.3) };
      const i90 = grid.indexOf(90), i30 = grid.indexOf(30);
      const resp = [], cr = [], peak = [], t90 = [];
      for (const r of recs) {
        const frac90 = r.T[i90] / r.T0;
        resp.push(frac90 < POP.respThr ? 1 : 0); cr.push(r.T[i90] < POP.crThr ? 1 : 0);
        peak.push(Math.max(...r.E)); t90.push(100 * frac90);
      }
      const val = GL.smc.val;
      const outcomes = { p_response: { value: U.mean(resp) }, p_cr: { value: U.mean(cr) }, peak_expansion: val(peak), tumor_d90: val(t90) };
      // signal: posterior-median expansion normalised to the reference peak
      const signals = { car_expansion: { grid, values: variables.car_copies.q50.map(v => v / POP.Eref) } };
      return { grid, variables, outcomes, signals, extra: {} };
    }
  };
  GL.modules.tumor_immune = M;
  if (GL.registry) GL.registry.register(M);
})(typeof window !== 'undefined' ? window : global);
