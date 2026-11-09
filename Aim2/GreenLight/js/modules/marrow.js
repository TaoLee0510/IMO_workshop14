/* GreenLight · modules/marrow.js
 * Aim 2 marrow module — browser surrogate of the Aim 2 mechanistic model
 * (H1 reserve / H2 inflammation / H3 clone). One ODE, three mechanisms:
 *   stem pools  S_W (normal), S_M (clone) sharing one niche K = 1
 *   neutrophil lineage  P -> T1 -> T2 -> T3 -> N (ANC)
 *   platelet lineage    MP -> M1 -> M2 -> M3 -> PLT
 *   inflammation        I (CRP-like) -> Im1 -> Im (marrow delay chain)
 * Learning: particle posterior over patient parameters
 *   R   marrow reserve (S(0) / healthy steady state)
 *   g   inflammatory gain
 *   f0  clone fraction of the stem pool at baseline
 *   tau maturation-delay scale
 *   N0, PLT0  set points (nuisance)
 * The real Aim 2 model (R/C++, Code/) can replace this module through the
 * same interface, locally or via modules/remote.js.
 */
(function (root) {
  const GL = root.GL = root.GL || {};
  const U = GL.util, S = GL.schema;
  GL.modules = GL.modules || {};

  // ---- population constants -------------------------------------------------
  const POP = {
    a: 0.036, d: 0.0072,         // stem self-renewal / differentiation  -> S_ss = 1 - d/a = 0.8
    cloneAdv: 1.3,               // clone self-renewal advantage
    delta: 0.1,                  // clone output per stem cell relative to normal
    kP: 0.55, kT: 0.48, kN: 2.2, // neutrophil lineage rates (1/day)
    gamma: 0.55, fbMax: 6,       // feedback (N0/N)^gamma capped
    kMP: 0.35, kM: 0.33, kPLT: 0.13, gammaP: 0.5,
    ld: { days: [-5, -3], kProg: 2.2, kProgT2: 0.6, kStem: 0.11, cloneResist: 0.2 },
    inflam: { Ibase: 5, kI: 0.32, kernelPeak: 6.5, tauM: 4.0, hill: 2.5, eProd: 0.95, qKill: 0.12, cloneResist: 0.2,
              tociBoost: 2.5, tociDays: 7, sterBoost: 1.2, sterGain: 0.6, anakGain: 0.7 },
    gcsf: { prodBoost: 2.2, transitBoost: 1.6 },
    tpo: { prodBoost: 2.0 },
    cd34: { add: 0.3 },
    pltTx: 25,
    crsThr: { 1: 40, 2: 90, 3: 160 }
  };
  const SS = 1 - POP.d / POP.a; // healthy stem steady state (fraction of niche)

  // ---- helpers -----------------------------------------------------------------
  function kernelD(t) { // default CAR T driver when no upstream signal: peak 1 at kernelPeak
    if (t <= 0) return 0;
    const p = POP.inflam.kernelPeak, x = t / p;
    return x * x * Math.exp(2 * (1 - x));
  }
  function ldFractions(snapshot, plan) {
    // locked LD stage -> recorded doses; else plan
    if (snapshot.asOfDay > -5 && snapshot.ld) {
      const f = snapshot.ld.flu_total !== undefined ? snapshot.ld.flu_total / 90 : (snapshot.ld.dose_reduced ? 0.75 : 1);
      const c = snapshot.ld.cy_total !== undefined ? snapshot.ld.cy_total / 1500 : (snapshot.ld.dose_reduced ? 0.75 : 1);
      return { flu: f, cy: c };
    }
    return { flu: plan.stages.ld.flu_frac, cy: plan.stages.ld.cy_frac };
  }

  // Build daily "input schedule" from recorded events (t <= asOfDay).
  // Policies for t > asOfDay are applied online during simulation.
  function recordedInputs(snapshot) {
    const inp = { gcsf: [], toci: [], steroid: [], anakinra: [], tpo: [], pltTx: [], cd34: [] };
    for (const e of snapshot.events) {
      const s = e.day_start, en = (e.day_end === null || e.day_end === undefined) ? s : e.day_end;
      switch (e.type) {
        case 'gcsf': inp.gcsf.push([s, en + 1]); break;
        case 'tocilizumab': inp.toci.push([s, s + POP.inflam.tociDays]); break;
        case 'steroid': inp.steroid.push([s, en + 1]); break;
        case 'anakinra': inp.anakinra.push([s, en + 1]); break;
        case 'tpo_ra': inp.tpo.push([s, Math.max(en, s + 60)]); break;
        case 'plt_tx': inp.pltTx.push(s); break;
        case 'cd34_boost': inp.cd34.push(s); break;
        default: break;
      }
    }
    return inp;
  }
  const inWin = (wins, t) => { for (const w of wins) if (t >= w[0] && t < w[1]) return true; return false; };

  // ---- the simulator ----------------------------------------------------------------
  // returns daily records from tStart..tEnd
  function simulate(p, snapshot, plan, ctx, tEnd) {
    const t0 = -6, dt = 0.25;
    const asOf = snapshot.asOfDay;
    const ld = ldFractions(snapshot, plan);
    const ldDose = 0.5 * ld.flu + 0.5 * ld.cy;
    const inp = recordedInputs(snapshot);
    const cz = plan.stages.crs, cy = plan.stages.cytopenia;
    const sig = ctx.signals && ctx.signals.car_expansion;
    const D = sig ? (t => Math.max(0, U.interp(sig.grid, sig.values, t))) : kernelD;

    // state indices
    const SW = 0, SM = 1, P = 2, T1 = 3, T2 = 4, T3 = 5, N = 6, I = 7, IM1 = 8, IM = 9, MP = 10, M1 = 11, M2 = 12, M3 = 13, PLT = 14;
    const y0 = new Float64Array(15);
    const S0 = U.clamp(p.R, 0.02, 1.25) * SS;
    y0[SW] = S0 * (1 - p.f0); y0[SM] = S0 * p.f0;
    // lineage initialised at the quasi steady state implied by the current output
    const out0 = (y0[SW] + POP.delta * y0[SM]) / SS;
    const n0 = p.N0 * Math.pow(out0, 1 / (1 + POP.gamma));
    const kinN = POP.kN * p.N0 / SS, kinP = POP.kPLT * p.PLT0 / SS;
    const fluxN = POP.kN * n0;
    y0[N] = n0; y0[T3] = fluxN / POP.kT; y0[T2] = y0[T3]; y0[T1] = y0[T3]; y0[P] = fluxN / POP.kP;
    const plt0 = p.PLT0 * Math.pow(out0, 1 / (1 + POP.gammaP));
    const fluxP = POP.kPLT * plt0;
    y0[PLT] = plt0; y0[M3] = fluxP / POP.kM; y0[M2] = y0[M3]; y0[M1] = y0[M3]; y0[MP] = fluxP / POP.kMP;
    y0[I] = POP.inflam.Ibase; y0[IM1] = POP.inflam.Ibase; y0[IM] = POP.inflam.Ibase;

    const kT = POP.kT * p.tau, kPp = POP.kP * p.tau;
    // online policy state
    const pol = { gcsfUntil: -1e9, gcsfCourses: 0, tociUntil: -1e9, tociDoses: 0, sterUntil: -1e9, anakUntil: -1e9, tpoFrom: 1e9, lastPltTx: -1e9, pltTx: 0, boosted: false, boosts: 0 };
    if (cz.anakinra_prophylaxis && asOf < 0) pol.anakUntil = 14; // prophylaxis window day 0-14 (future only)
    const rec = { days: [], anc: [], plt: [], crp: [], I: [], sw: [], sm: [] };
    let cur = { gcsf: 0, toci: 0, ster: 0, anak: 0, tpo: 0 };

    const rhs = function (t, y, dy) {
      // LD kill kernel: full effect on dosing days, 40% tail for 3 more days (drug persistence)
      const ldK = (t >= POP.ld.days[0] && t < POP.ld.days[1] + 1) ? 1 : (t >= POP.ld.days[1] + 1 && t < POP.ld.days[1] + 4) ? 0.4 : 0;
      const kLDp = POP.ld.kProg * ldDose * ldK, kLDs = POP.ld.kStem * ldDose * (ldK === 1 ? 1 : 0);
      const Stot = y[SW] + y[SM];
      const x = Math.pow(p.q * y[IM], POP.inflam.hill);
      const H = x / (1 + x);
      const kill = POP.inflam.qKill * H;
      dy[SW] = POP.a * y[SW] * (1 - Stot) - kLDs * y[SW] - kill * y[SW] - POP.d * y[SW];
      dy[SM] = POP.a * POP.cloneAdv * y[SM] * (1 - Stot) - POP.ld.cloneResist * kLDs * y[SM] - POP.inflam.cloneResist * kill * y[SM] - POP.d * y[SM];
      let fb = Math.pow(p.N0 / Math.max(y[N], 0.02), POP.gamma); if (fb > POP.fbMax) fb = POP.fbMax;
      let fbp = Math.pow(p.PLT0 / Math.max(y[PLT], 2), POP.gammaP); if (fbp > POP.fbMax) fbp = POP.fbMax;
      const prodSupp = 1 - POP.inflam.eProd * H;
      const gcsfK = cur.gcsf ? POP.gcsf.prodBoost : 1, gcsfT = cur.gcsf ? POP.gcsf.transitBoost : 1;
      const tpoK = cur.tpo ? POP.tpo.prodBoost : 1;
      const stemOut = y[SW] + POP.delta * y[SM];
      dy[P] = kinN * stemOut * fb * prodSupp * gcsfK - kPp * y[P] - kLDp * y[P];
      dy[T1] = kPp * y[P] - kT * gcsfT * y[T1] - kLDp * y[T1];
      dy[T2] = kT * gcsfT * (y[T1] - y[T2]) - POP.ld.kProgT2 * kLDp * y[T2];
      dy[T3] = kT * gcsfT * (y[T2] - y[T3]);
      dy[N] = kT * gcsfT * y[T3] - POP.kN * y[N];
      // inflammation
      const gEff = p.g * (cur.ster ? POP.inflam.sterGain : 1) * (cur.anak ? POP.inflam.anakGain : 1);
      const elim = POP.inflam.kI * (1 + (cur.toci ? POP.inflam.tociBoost : 0) + (cur.ster ? POP.inflam.sterBoost : 0));
      dy[I] = gEff * D(t) - elim * (y[I] - POP.inflam.Ibase);
      dy[IM1] = (y[I] - y[IM1]) / POP.inflam.tauM;
      dy[IM] = (y[IM1] - y[IM]) / POP.inflam.tauM;
      // platelets
      dy[MP] = kinP * stemOut * fbp * prodSupp * tpoK - POP.kMP * y[MP] - kLDp * y[MP];
      dy[M1] = POP.kMP * y[MP] - POP.kM * y[M1] - kLDp * y[M1];
      dy[M2] = POP.kM * (y[M1] - y[M2]);
      dy[M3] = POP.kM * (y[M2] - y[M3]);
      dy[PLT] = POP.kM * y[M3] - POP.kPLT * y[PLT];
    };

    const onStep = function (t, y, s) {
      const isDay = Math.abs(t - Math.round(t)) < 1e-9;
      if (isDay) {
        const d = Math.round(t);
        if (d <= asOf) {
          // recorded history
          cur.gcsf = inWin(inp.gcsf, d) ? 1 : 0; cur.toci = inWin(inp.toci, d) ? 1 : 0;
          cur.ster = inWin(inp.steroid, d) ? 1 : 0; cur.anak = inWin(inp.anakinra, d) ? 1 : 0; cur.tpo = inWin(inp.tpo, d) ? 1 : 0;
          if (inp.pltTx.includes(d)) y[PLT] += POP.pltTx;
          if (inp.cd34.includes(d)) { y[SW] += POP.cd34.add; pol.boosted = true; }
          // carry ongoing windows into policy state so they continue after asOf
          if (cur.toci) pol.tociUntil = Math.max(pol.tociUntil, d + 1);
          if (cur.gcsf) pol.gcsfUntil = Math.max(pol.gcsfUntil, d + 1);
          if (cur.tpo) pol.tpoFrom = Math.min(pol.tpoFrom, d);
        } else {
          // plan policies (future)
          const thrT = POP.crsThr[cz.toci_grade], thrS = POP.crsThr[cz.steroid_grade];
          if (d <= 28 && thrT && y[I] >= thrT && d >= pol.tociUntil && pol.tociDoses < 2) { pol.tociUntil = d + POP.inflam.tociDays; pol.tociDoses++; }
          if (d <= 28 && thrS && y[I] >= thrS && d >= pol.sterUntil) { pol.sterUntil = d + 3; }
          if (cy.gcsf === 'prophylactic' && d >= cy.gcsf_start_day && d <= 60 && d >= pol.gcsfUntil && y[N] < 1.0 && pol.gcsfCourses < 6) { pol.gcsfUntil = d + cy.gcsf_days; pol.gcsfCourses++; }
          if (cy.gcsf === 'reactive' && d >= cy.gcsf_start_day && d <= 90 && d >= pol.gcsfUntil && y[N] < cy.gcsf_anc_thr && pol.gcsfCourses < 5) { pol.gcsfUntil = d + cy.gcsf_days; pol.gcsfCourses++; }
          if (cy.tpo_ra && d >= 21 && y[PLT] < cy.tpo_plt_thr) pol.tpoFrom = Math.min(pol.tpoFrom, d);
          if (d >= 3 && y[PLT] < cy.plt_tx_thr && d - pol.lastPltTx >= 3) { y[PLT] += POP.pltTx; pol.lastPltTx = d; pol.pltTx++; }
          if (cy.cd34_backup && cy.cd34_boost === 'if_persistent' && d === cy.cd34_boost_day && y[N] < 0.5 && !pol.boosted) { y[SW] += POP.cd34.add; pol.boosted = true; pol.boosts++; }
          cur.gcsf = d < pol.gcsfUntil ? 1 : 0; cur.toci = d < pol.tociUntil ? 1 : 0; cur.ster = d < pol.sterUntil ? 1 : 0;
          cur.anak = d < pol.anakUntil ? 1 : 0; cur.tpo = d >= pol.tpoFrom ? 1 : 0;
        }
        rec.days.push(d); rec.anc.push(y[N]); rec.plt.push(y[PLT]); rec.crp.push(y[I]); rec.I.push(y[I]); rec.sw.push(y[SW]); rec.sm.push(y[SM]);
      }
    };
    GL.ode.integrate({ t0, t1: tEnd, dt, y0, rhs, onStep });
    rec.pol = pol;
    return rec;
  }

  // ---- prior ---------------------------------------------------------------------
  function priorParticles(snapshot, n, seed) {
    const r = U.rng(seed);
    const h = snapshot.history || {};
    // covariate-informed prior mean for reserve
    let mR = 0.85 - 0.04 * (h.prior_lines || 0) - 0.10 * (h.prior_asct ? 1 : 0) - 0.10 * (h.prior_bendamustine ? 1 : 0) - 0.10 * (h.marrow_involvement ? 1 : 0);
    if (h.months_since_last_chemo !== undefined && h.months_since_last_chemo < 3) mR -= 0.05;
    mR = U.clamp(mR, 0.3, 0.9);
    const mut = snapshot.mutations || { tested: false, variants: [] };
    const maxVaf = mut.variants && mut.variants.length ? Math.max(...mut.variants.map(v => v.vaf)) : 0;
    const gMed = snapshot.costim === 'CD28' ? 45 : 35;
    const parts = [];
    for (let i = 0; i < n; i++) {
      const R = U.clamp(mR + 0.2 * r.normal(), 0.05, 1.2);
      const g = r.lognormal(gMed, 0.85);
      let f0;
      if (maxVaf > 0) f0 = U.clamp(2 * maxVaf * Math.exp(0.35 * r.normal()), 0.005, 0.9);
      else if (mut.tested) f0 = r.next() < 0.9 ? r.uniform(0, 0.02) : r.uniform(0.02, 0.3);
      else f0 = r.next() < 0.75 ? r.uniform(0, 0.03) : r.uniform(0.03, 0.45);
      const tau = r.lognormal(1.0, 0.15);
      const q = r.lognormal(1 / 90, 0.7);
      const N0 = r.lognormal(3.6, 0.3), PLT0 = r.lognormal(230, 0.3);
      parts.push({ R, g, q, f0, tau, N0, PLT0 });
    }
    return parts;
  }

  // ---- likelihood -------------------------------------------------------------------
  const SIG = { anc: 0.35, plt: 0.32, crp: 0.55 };
  const JITTER = { R: { sd: 0.06, lo: 0.05, hi: 1.2 }, g: { log: true, sd: 0.25, lo: 1, hi: 1000 }, q: { log: true, sd: 0.2, lo: 1 / 1000, hi: 0.2 },
                   f0: { log: true, sd: 0.3, lo: 0.0005, hi: 0.9 }, tau: { log: true, sd: 0.05, lo: 0.5, hi: 2 }, N0: { log: true, sd: 0.1, lo: 1, hi: 12 }, PLT0: { log: true, sd: 0.1, lo: 40, hi: 600 } };
  function loglik(rec, snapshot) {
    let ll = 0, n = 0;
    const idx = {}; rec.days.forEach((d, i) => { idx[d] = i; });
    for (const o of snapshot.observations) {
      if (!(o.variable in SIG)) continue;
      const i = idx[o.day]; if (i === undefined) continue;
      const pred = o.variable === 'anc' ? rec.anc[i] : o.variable === 'plt' ? rec.plt[i] : rec.crp[i];
      const floor = o.variable === 'plt' ? 2 : 0.05;
      const z = (Math.log(o.value + floor) - Math.log(pred + floor)) / SIG[o.variable];
      ll += -0.5 * z * z; n++;
    }
    return { ll, n };
  }
  // ---- module ------------------------------------------------------------------------
  const M = {
    meta: { id: 'marrow', nameKey: 'mod.marrow.name', descKey: 'mod.marrow.desc', owner: 'Aim 2 · Team Green', version: '0.3-surrogate', kind: 'local', placeholder: false, engine: 'JS RK4 surrogate · particle posterior' },
    consumes: { variables: ['anc', 'plt', 'crp'], planKeys: ['ld.flu_frac', 'ld.cy_frac', 'crs.toci_grade', 'crs.steroid_grade', 'crs.anakinra_prophylaxis', 'cytopenia.*'], signals: ['car_expansion'] },
    provides: {
      signals: ['inflammation'],
      params: [
        { id: 'R', labelKey: 'param.R', unit: '× healthy', lo: 0, hi: 1.2, direction: 'higher', descKey: 'param.R.desc' },
        { id: 'g', labelKey: 'param.g', unit: 'mg/L·d⁻¹', log: true, lo: 2, hi: 600, direction: 'lower', descKey: 'param.g.desc' },
        { id: 'q', labelKey: 'param.q', unit: 'L/mg', log: true, lo: 0.002, hi: 0.1, direction: 'lower', descKey: 'param.q.desc' },
        { id: 'f0', labelKey: 'param.f0', unit: 'fraction', lo: 0, hi: 0.9, direction: 'lower', descKey: 'param.f0.desc' },
        { id: 'tau', labelKey: 'param.tau', unit: '× population', lo: 0.5, hi: 1.8, direction: 'neutral', descKey: 'param.tau.desc' }
      ],
      outcomes: [
        { id: 'p_quick', labelKey: 'out.p_quick', kind: 'prob', direction: 'higher' },
        { id: 'p_second_wave', labelKey: 'out.p_second_wave', kind: 'prob', direction: 'lower' },
        { id: 'p_persistent', labelKey: 'out.p_persistent', kind: 'prob', direction: 'lower' },
        { id: 'sev_neut_days', labelKey: 'out.sev_neut_days', kind: 'value', direction: 'lower', unit: 'd' },
        { id: 'anc_nadir_late', labelKey: 'out.anc_nadir_late', kind: 'value', direction: 'higher', unit: '×10⁹/L' },
        { id: 'p_plt20', labelKey: 'out.p_plt20', kind: 'prob', direction: 'lower' },
        { id: 'p_crs3', labelKey: 'out.p_crs3', kind: 'prob', direction: 'lower' },
        { id: 'gcsf_courses', labelKey: 'out.gcsf_courses', kind: 'value', direction: 'lower', unit: 'n' },
        { id: 'plt_tx', labelKey: 'out.plt_tx', kind: 'value', direction: 'lower', unit: 'n' }
      ],
      variables: ['anc', 'plt', 'crp']
    },
    POP, simulate, priorParticles,

    async fit(snapshot, ctx) {
      const n = (ctx && ctx.nParticles) || 500;
      const seed = (ctx && ctx.seed) || U.hashString(snapshot.patient_id);
      const prior = priorParticles(snapshot, n, seed);
      const plan = S.defaultPlan();
      const tEnd = Math.max(snapshot.asOfDay, -5);
      const simLL = (p) => loglik(simulate(p, snapshot, plan, ctx || {}, tEnd), snapshot);
      const fit = GL.smc.fit(prior, simLL, JITTER, seed);
      const params = {};
      for (const k of ['R', 'g', 'q', 'f0', 'tau', 'N0', 'PLT0']) params[k] = GL.smc.summarise(fit.particles, fit.weights, k);
      return { moduleId: 'marrow', asOfDay: snapshot.asOfDay, seed, n, nObs: fit.nObs, beta: fit.beta, ess: fit.ess, particles: fit.particles, weights: fit.weights, params, logEvidence: fit.logEvidence };
    },

    async forecast(state, plan, ctx) {
      const m = (ctx && ctx.nForecast) || 80, horizon = (ctx && ctx.horizon) || 120;
      const snapshot = ctx.snapshot;
      const parts = GL.smc.resample(state.particles, state.weights, m, state.seed + 7);
      const recs = parts.map(p => simulate(p, snapshot, plan, ctx, horizon));
      const grid = recs[0].days;
      const q = (arrs, sd) => GL.smc.predictiveBands(grid, arrs, sd);
      const variables = { anc: q(recs.map(r => r.anc), 0.22), plt: q(recs.map(r => r.plt), 0.22), crp: q(recs.map(r => r.crp), 0.45) };
      // outcomes per particle
      const ph = { quick: 0, second_wave: 0, persistent: 0 };
      const sev = [], nadir = [], plt20 = [], crs3 = [], gcsf = [], ptx = [], boosts = [];
      for (const r of recs) {
        ph[S.classifyPhenotype(r.days, r.anc)]++;
        let s = 0, nd = Infinity, p20 = 0, peakI = 0;
        for (let i = 0; i < r.days.length; i++) {
          const d = r.days[i];
          if (d >= 0 && d <= 90 && r.anc[i] < 0.5) s++;
          if (d >= 14 && d <= 60 && r.anc[i] < nd) nd = r.anc[i];
          if (d >= 10 && d <= 90 && r.plt[i] < 20) p20 = 1;
          if (r.I[i] > peakI) peakI = r.I[i];
        }
        sev.push(s); nadir.push(nd); plt20.push(p20); crs3.push(peakI >= POP.crsThr[3] ? 1 : 0);
        gcsf.push(r.pol.gcsfCourses); ptx.push(r.pol.pltTx); boosts.push(r.pol.boosts);
      }
      const val = GL.smc.val;
      const outcomes = {
        p_quick: { value: ph.quick / m }, p_second_wave: { value: ph.second_wave / m }, p_persistent: { value: ph.persistent / m },
        sev_neut_days: val(sev), anc_nadir_late: val(nadir), p_plt20: { value: U.mean(plt20) }, p_crs3: { value: U.mean(crs3) },
        gcsf_courses: val(gcsf), plt_tx: val(ptx), cd34_boosts: val(boosts)
      };
      const signals = { inflammation: { grid, values: variables.crp.q50 } };
      return { grid, variables, outcomes, signals, extra: { stem: { sw: recs.map(r => r.sw[r.sw.length - 1]), sm: recs.map(r => r.sm[r.sm.length - 1]) } } };
    }
  };
  GL.modules.marrow = M;
  if (GL.registry) GL.registry.register(M);
})(typeof window !== 'undefined' ? window : global);
