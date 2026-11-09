/* GreenLight · core/orchestrator.js
 * Runs the registered modules in dependency order, passes shared signals
 * downstream, evaluates candidate plans, ranks them, and packages everything
 * into an immutable TwinVersion. Also serves what-if evaluation.
 */
(function (root) {
  const GL = root.GL = root.GL || {};
  const S = GL.schema, U = GL.util, R = GL.registry, P = GL.plans;
  const O = GL.orchestrator = {};

  O.defaults = { horizon: 120, nParticles: 500, nForecast: 80 };

  const round = (x) => (typeof x === 'number' && isFinite(x)) ? Number(x.toPrecision(4)) : x;
  function compactForecast(fc) {
    const vars = {};
    for (const v in fc.variables) { vars[v] = {}; for (const q in fc.variables[v]) vars[v][q] = fc.variables[v][q].map(round); }
    const outs = {};
    for (const k in fc.outcomes) { outs[k] = {}; for (const q in fc.outcomes[k]) outs[k][q] = round(fc.outcomes[k][q]); }
    return { grid: fc.grid, variables: vars, outcomes: outs };
  }

  // run the module chain for one plan given fitted states
  async function chainForecast(mods, states, plan, ctxBase) {
    const signals = {}, forecasts = {}, outcomes = {};
    for (const mod of mods) {
      const fc = await mod.forecast(states[mod.meta.id], plan, { ...ctxBase, signals: { ...signals } });
      Object.assign(signals, fc.signals || {});
      forecasts[mod.meta.id] = fc; outcomes[mod.meta.id] = fc.outcomes;
    }
    return { forecasts, outcomes, signals };
  }

  O.buildVersion = async function (patient, asOfDay, opts) {
    opts = opts || {};
    const t0 = Date.now();
    const snapshot = S.snapshot(patient, asOfDay);
    const mods = R.order();
    const ctxBase = { asOfDay, horizon: opts.horizon || O.defaults.horizon, nParticles: opts.nParticles || O.defaults.nParticles, nForecast: opts.nForecast || O.defaults.nForecast, snapshot };
    const basePlan = opts.activePlan || S.defaultPlan();
    // 1. fit chain (upstream signals from the base plan)
    const states = {}, moduleInfo = {};
    let signals = {};
    for (const mod of mods) {
      const tf = Date.now();
      const state = await mod.fit(snapshot, { ...ctxBase, signals: { ...signals } });
      const fc = await mod.forecast(state, basePlan, { ...ctxBase, signals: { ...signals } });
      Object.assign(signals, fc.signals || {});
      states[mod.meta.id] = state;
      moduleInfo[mod.meta.id] = { meta: mod.meta, coverage: R.coverage(mod, snapshot), nObs: state.nObs, ess: round(state.ess), beta: round(state.beta), logEvidence: round(state.logEvidence), fitMs: Date.now() - tf, params: state.params };
    }
    // 2. candidate plans
    const seen = new Set();
    const plans = P.library().concat(opts.extraPlans || []).filter(pl => { const k = P.effectiveKey(pl, asOfDay); if (seen.has(k)) return false; seen.add(k); return true; });
    const weights = opts.weights || P.defaultWeights;
    const planResults = {};
    for (const plan of plans) {
      const r = await chainForecast(mods, states, plan, ctxBase);
      const u = P.utility(r.outcomes, plan, weights);
      const fcs = {}; for (const id in r.forecasts) fcs[id] = compactForecast(r.forecasts[id]);
      planResults[plan.id] = { plan, forecasts: fcs, outcomes: r.outcomes, score: round(u.score), components: u.components };
    }
    const ranking = Object.keys(planResults).sort((a, b) => planResults[b].score - planResults[a].score);
    // snapshot summary for the UI / diff
    const nObs = {}; snapshot.observations.forEach(o => { nObs[o.variable] = (nObs[o.variable] || 0) + 1; });
    const lastDay = snapshot.observations.length ? Math.max(...snapshot.observations.map(o => o.day)) : null;
    return {
      id: U.uid('v'), patient_id: patient.patient_id, asOfDay, createdAt: new Date().toISOString(), computeMs: Date.now() - t0,
      dataSummary: { nObs, lastDay, nEvents: snapshot.events.length, crsKnown: !!snapshot.crs, mutations: snapshot.mutations },
      modules: moduleInfo, states, plans: planResults, ranking, activePlanId: basePlan.id, weights, signals: Object.fromEntries(Object.entries(signals).map(([k, v]) => [k, { grid: v.grid, values: v.values.map(round) }]))
    };
  };

  // re-create module states for a version that was persisted without them (fits are deterministic)
  O.refit = async function (patient, version) {
    const snapshot = S.snapshot(patient, version.asOfDay);
    const mods = R.order().filter(m => version.modules[m.meta.id]);
    const ctxBase = { asOfDay: version.asOfDay, horizon: O.defaults.horizon, nParticles: O.defaults.nParticles, nForecast: O.defaults.nForecast, snapshot };
    const states = {}; let signals = {};
    const basePlan = version.plans[version.activePlanId] ? version.plans[version.activePlanId].plan : S.defaultPlan();
    for (const mod of mods) {
      const state = await mod.fit(snapshot, { ...ctxBase, signals: { ...signals } });
      const fc = await mod.forecast(state, basePlan, { ...ctxBase, signals: { ...signals } });
      Object.assign(signals, fc.signals || {}); states[mod.meta.id] = state;
    }
    return states;
  };
  O.stateCache = {};
  O.statesFor = async function (patient, version) {
    if (version.states) return version.states;
    if (!O.stateCache[version.id]) O.stateCache[version.id] = await O.refit(patient, version);
    return O.stateCache[version.id];
  };

  // what-if: evaluate a custom plan on an existing version (no persistence)
  O.whatIf = async function (patient, version, plan, weights) {
    const snapshot = S.snapshot(patient, version.asOfDay);
    const states = await O.statesFor(patient, version);
    const mods = R.order().filter(m => states[m.meta.id]);
    const ctxBase = { asOfDay: version.asOfDay, horizon: O.defaults.horizon, nForecast: O.defaults.nForecast, snapshot };
    const r = await chainForecast(mods, states, plan, ctxBase);
    const u = P.utility(r.outcomes, plan, weights || version.weights);
    const fcs = {}; for (const id in r.forecasts) fcs[id] = compactForecast(r.forecasts[id]);
    return { plan, forecasts: fcs, outcomes: r.outcomes, score: round(u.score), components: u.components };
  };

  // re-rank stored plans with new weights (no recomputation needed)
  O.rerank = function (version, weights) {
    const ids = Object.keys(version.plans);
    const scored = ids.map(id => { const pr = version.plans[id]; const u = P.utility(pr.outcomes, pr.plan, weights); return { id, score: u.score, components: u.components }; });
    scored.sort((a, b) => b.score - a.score);
    return scored;
  };

  // forecast skill of a previous version against observations that arrived afterwards
  O.evaluateForecast = function (prevVersion, patient, newAsOfDay, planId) {
    const pr = prevVersion.plans[planId || prevVersion.activePlanId];
    if (!pr) return null;
    const out = {};
    for (const modId in pr.forecasts) {
      const fc = pr.forecasts[modId];
      for (const v in fc.variables) {
        const obs = patient.observations.filter(o => o.variable === v && o.day > prevVersion.asOfDay && o.day <= newAsOfDay);
        if (!obs.length) continue;
        let se = 0, cov = 0, n = 0;
        for (const o of obs) {
          const i = fc.grid.indexOf(o.day); if (i < 0) continue;
          const pred = fc.variables[v].q50[i];
          const floor = v === 'plt' ? 2 : 0.05;
          const e = Math.log(o.value + floor) - Math.log(pred + floor);
          se += e * e; n++;
          if (o.value >= fc.variables[v].q05[i] && o.value <= fc.variables[v].q95[i]) cov++;
        }
        if (n) out[v] = { n, rmsle: Math.sqrt(se / n), coverage90: cov / n };
      }
    }
    return out;
  };
})(typeof window !== 'undefined' ? window : global);
