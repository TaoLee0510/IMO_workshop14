/* GreenLight · engine/model3_core.js
 * JavaScript port of the Model3 engine (MarrowTwin; Model3/Code/src/model3.cpp, 32 states, same
 * Dormand-Prince RK45 step control), of make_segments() / build_par() / eta_to_indiv()
 * (Model3/Code/R/model.R) and of the phenotype rule and phenotype_prob()
 * (Model3/Code/R/phenotype.R). Model2 = Model1 + chemotherapy kill of the proliferating chain
 * precursors (T1, T2 at kappa_c * th_aP * C) and of circulating monocytes (th_am * C * B_Mo);
 * Model3 = Model2 + (N) a neutrophil chain of n_tr_n compartments (up to 8: TN1..TN8), the first nk_n
 * proliferating (killed at kappa_c * th_aP * C), the rest the post-mitotic storage pool whose outflow
 * G-CSF multiplies by rel_N (ktr_n = n_tr_n / MTT_N); (P2) the platelet MPP component uses its own
 * suppression e_z_pl = plogis(qlogis(e_z) + d_zpl); (P1) -th_cp * D(t) * B_P in the platelet blood pool
 * (th_cp fixed at 0 in Model3, kept so the engine is identical).
 * Used in the browser for what-if forecasts and new data; the displayed as-of results are
 * precomputed by the R engine (step G).
 * Tested against R/C++ by tools/test_model3_ode.js (daily grid, < 1e-8 relative).
 */
(function (root) {
  const GL = root.GL = root.GL || {};
  const C = GL.model3 = {};

  // ---- layouts (m3_param_names / m3_state_names) -----------------------------------------
  C.PAR = ['r_s', 'f0', 'th_as', 'th_ap', 'th_bs', 'rho', 'a0',
    'r_a', 'r_b', 'delta', 'eta_syn', 'th_zs',
    'k_p', 'sr_p', 'gmax', 'gamma',
    'kappa_c', 'th_am', 'th_cp',
    'ktr_n', 'ktr_mo', 'ktr_pl',
    'c_n', 'c_mo', 'c_pl', 'c_l',
    'k_i', 'emax', 'ic50', 'hill',
    'th_al',
    'k_e', 'gcsf_fold', 'tpo_fold', 'k_g', 'gcsf_rel',
    'k_t', 'k_ag', 'kappa', 'k_z', 'z50', 'e_z',
    'e_z_pl', 'n_tr_n', 'nk_n'];
  const P = {}; C.PAR.forEach((k, i) => { P[k] = i; });
  C.STATES = ['SW', 'SX', 'PN', 'PMO', 'PPL', 'PL', 'C', 'I1', 'I',
    'TN1', 'TN2', 'TN3', 'TN4', 'TN5', 'TN6', 'TN7', 'TN8', 'BN',
    'TM1', 'TM2', 'TM3', 'BM', 'TP1', 'TP2', 'TP3', 'BP', 'BL', 'CT', 'AG', 'Z1', 'Z', 'GC'];
  const ix = (k) => C.STATES.indexOf(k);
  const [SW, SX, PN, PMO, PPL, PL, CX, I1, I2, TN1, BN, TM1, BM, TP1, BP, BL, CT, AG, Z1, Z2, GC] =
    ['SW', 'SX', 'PN', 'PMO', 'PPL', 'PL', 'C', 'I1', 'I', 'TN1', 'BN', 'TM1', 'BM', 'TP1', 'BP', 'BL', 'CT', 'AG', 'Z1', 'Z', 'GC'].map(ix);
  C.N_TR_MAX = 8;   // TN1..TN8
  const NS = C.STATES.length;

  // ---- model ids (model_spec) -------------------------------------------------------------------
  C.spec = function (id) {
    if (id === 'MF') return { id, family: 'mf', q2: null, clone: 'none', route: 'none', re: ['slope', 'mtt_scale'] };
    const m = /^(M0|Ma|Mb|Mab|Mabx)(C_(ddr|dta))?(_Z(h|m))?(_nc)?(_no(N|P1|P2))?$/.exec(id);
    if (!m) throw new Error('unknown model id ' + id);
    const q2 = m[1], clone = m[3] || 'none', route = m[5] === 'h' ? 'hsc' : m[5] === 'm' ? 'mpp' : 'both';
    const re = ['r_s', 'th_ap'];
    if (['Ma', 'Mab', 'Mabx'].includes(q2)) re.push('th_as');
    if (['Mb', 'Mab', 'Mabx'].includes(q2)) re.push('th_bs');
    if (clone !== 'none') re.push('f0');
    re.push('rho');
    return { id, family: 'm3', q2, clone, route, re, z_hsc: ['Mb', 'Mab', 'Mabx'].includes(q2) && route !== 'mpp', z_mpp: route !== 'hsc',
      f0_nocov: !!m[6], drop: m[8] || 'none' };
  };

  // eta (named random effects + b_* baseline deviations) -> natural-scale individual values
  C.etaToIndiv = function (eta, spec) {
    const out = { r_s: Math.exp(eta.r_s), f0: 0, th_as: 0, th_ap: Math.exp(eta.th_ap), th_bs: 0, rho: 0 };
    if ('th_as' in eta) out.th_as = Math.exp(eta.th_as);
    if ('th_bs' in eta) out.th_bs = Math.exp(eta.th_bs);
    if ('f0' in eta) out.f0 = 1 / (1 + Math.exp(-eta.f0));
    if ('rho' in eta) out.rho = Math.exp(eta.rho);
    return out;
  };

  // Model2 values of the three Model3 changes (config$model$drop_null; model ids _noN / _noP1 / _noP2)
  C.dropNull = { N: { n_tr_n: 3, nk_n: 2 }, P1: { th_cp: 0 }, P2: { d_zpl: 0 } };
  const plogis = (x) => 1 / (1 + Math.exp(-x)), qlogis = (p) => Math.log(p / (1 - p));

  // full parameter vector (build_par); struct on the natural scale incl. mtt_*, th_cp, d_zpl, n_tr_n, nk_n
  C.buildPar = function (ind, s0, spec) {
    const s = (spec.drop && spec.drop !== 'none') ? Object.assign({}, s0, C.dropNull[spec.drop]) : s0;
    if (!(s.sr_p < s.k_p)) throw new Error('self-renewal rate must be below the exit rate');
    const nN = s.n_tr_n, nk = s.nk_n;
    if (!(Number.isInteger(nN) && nN >= 3 && nN <= C.N_TR_MAX && Number.isInteger(nk) && nk >= 1 && nk <= nN - 1)) throw new Error('neutrophil chain: n_tr_n in 3..8, nk_n in 1..n_tr_n - 1');
    const on = (k) => spec.re.includes(k) ? ind[k] : 0;
    const ez = spec.z_mpp ? s.e_z : 0;
    const v = {
      r_s: ind.r_s, f0: on('f0'), th_as: on('th_as'), th_ap: ind.th_ap, th_bs: on('th_bs'),
      rho: ind.rho === undefined ? 0 : ind.rho, a0: ind.a0 === undefined ? 1 : ind.a0,
      r_a: s.r_a, r_b: s.r_b, delta: s.delta, eta_syn: spec.q2 === 'Mabx' ? s.eta_syn : 0,
      th_zs: spec.z_hsc ? s.th_zs : 0,
      k_p: s.k_p, sr_p: s.sr_p, gmax: s.gmax, gamma: s.gamma,
      kappa_c: s.kappa_c, th_am: s.th_am, th_cp: s.th_cp,
      ktr_n: nN / s.mtt_n, ktr_mo: 3 / s.mtt_mo, ktr_pl: 3 / s.mtt_pl,
      c_n: s.c_n, c_mo: s.c_mo, c_pl: s.c_pl, c_l: s.c_l,
      k_i: s.k_i, emax: s.emax, ic50: s.ic50, hill: s.hill, th_al: s.th_al,
      k_e: s.k_e, gcsf_fold: s.gcsf_fold, tpo_fold: s.tpo_fold, k_g: s.k_g, gcsf_rel: s.gcsf_rel,
      k_t: s.k_t, k_ag: s.k_ag, kappa: s.kappa, k_z: s.k_z, z50: s.z50, e_z: ez,
      // P2: platelet MPP component, logit offset d_zpl from e_z (0 without the MPP route)
      e_z_pl: ez > 0 ? plogis(qlogis(ez) + s.d_zpl) : 0,
      n_tr_n: nN, nk_n: nk
    };
    for (const k of C.PAR) if (typeof v[k] !== 'number' || Number.isNaN(v[k])) throw new Error('parameter ' + k + ' missing');
    return Float64Array.from(C.PAR.map(k => v[k]));
  };
  C.a0 = (s, ldh_c) => Math.exp(s.beta_a * (ldh_c === undefined || ldh_c === null ? 0 : ldh_c));

  // ---- inputs (make_segments; all breakpoints are whole days) -------------------------------------
  // design: { t0, egfr, bridging: 'none'|'rgo'|'cy', apheresis_day, ld_scale, no_car,
  //           gcsf: [[start, end)], tpo: [[start, end)], cd34: [{day, add}] }
  C.bridgingWindow = [2, -7];
  C.segments = function (design, s, tEnd) {
    const br = design.bridging === 'rgo' ? s.u_br_rgo : design.bridging === 'cy' ? s.u_br_cy : 0;
    const bw = br > 0 ? [design.apheresis_day + C.bridgingWindow[0], C.bridgingWindow[1]] : null;
    const ld = Math.pow(90 / design.egfr, s.beta_renal) * (design.ld_scale === undefined ? 1 : design.ld_scale);
    const gc = design.gcsf || [], tp = design.tpo || [], cd = design.cd34 || [];
    let bp = [design.t0, tEnd, -5, -2, 0];
    if (bw) bp.push(bw[0], bw[1]);
    gc.forEach(w => bp.push(w[0], w[1])); tp.forEach(w => bp.push(w[0], w[1])); cd.forEach(c => bp.push(c.day));
    bp = Array.from(new Set(bp.filter(x => x >= design.t0 && x <= tEnd))).sort((a, b) => a - b);
    const inWin = (ws, x) => ws.some(w => x >= w[0] && x < w[1]);
    const out = [];
    for (let i = 0; i < bp.length - 1; i++) {
      const st = bp[i], en = bp[i + 1], mid = (st + en) / 2;
      let dose = (mid >= -5 && mid < -2) ? ld : 0;
      if (bw && mid >= bw[0] && mid < bw[1]) dose += br;
      let sAdd = 0; cd.forEach(c => { if (c.day === st) sAdd += c.add; });
      out.push({ start: st, end: en, dose, gcsf: inWin(gc, mid) ? 1 : 0, tpo: inWin(tp, mid) ? 1 : 0, s_add: sAdd,
        car: (!design.no_car && st === 0) ? s.car_t0 : 0 });
    }
    return out;
  };

  C.drvAt = function (drv, t) {
    const x = (t - drv.t0) / drv.dt, n = drv.v.length;
    if (x <= 0) return drv.v[0];
    if (x >= n - 1) return drv.v[n - 1];
    const i = Math.floor(x), w = x - i;
    return (1 - w) * drv.v[i] + w * drv.v[i + 1];
  };
  C.kernel = (t, A, tp, h) => t > 0 ? A * Math.pow(t / tp, h) * Math.exp(h * (1 - t / tp)) : 0;

  // ---- right-hand side (model3.cpp Parent::rhs) ----------------------------------------------------
  // kc: chemotherapy kill of the proliferating precursors (T1, T2)
  function chain(y, dy, t1, k, c, inflow, kc, rel) {
    dy[t1] = k * (inflow - y[t1]) - kc * y[t1];
    dy[t1 + 1] = k * (y[t1] - y[t1 + 1]) - kc * y[t1 + 1];
    dy[t1 + 2] = k * (y[t1 + 1] - rel * y[t1 + 2]);
    dy[t1 + 3] = c * (rel * y[t1 + 2] - y[t1 + 3]);
  }
  // neutrophil chain (model3.cpp Parent::chain_n): n compartments (3..8), the first nk proliferating
  // (killed at kc), the rest the post-mitotic storage pool, whose outflow G-CSF multiplies by rel;
  // with n < 8 the unused compartments stay constant; blood gets the outflow of the last used one
  function chainN(y, dy, n, nk, k, c, inflow, kc, rel) {
    let prev = inflow;
    for (let j = 0; j < C.N_TR_MAX; j++) {
      const x = TN1 + j;
      if (j >= n) { dy[x] = 0; continue; }
      const out = (j >= nk ? rel : 1) * y[x];
      dy[x] = k * (prev - out) - (j < nk ? kc * y[x] : 0);
      prev = out;
    }
    dy[BN] = c * (prev - y[BN]);
  }
  C.rhs = function (t, y, dy, p, drv, u, dScale) {
    const mx = Math.max;
    const sw = mx(y[SW], 0), sx = mx(y[SX], 0), s = sw + sx;
    const cx = mx(y[CX], 0), i1 = mx(y[I1], 0), im = mx(y[I2], 0);
    let d = mx(C.drvAt(drv, t), 0);
    if (dScale && t > dScale.after) d *= dScale.factor;
    const bn = mx(y[BN], 1e-6), bm = mx(y[BM], 1e-6), bp = mx(y[BP], 1e-6), bl = mx(y[BL], 1e-6);
    const ct = mx(y[CT], 0), ag = mx(y[AG], 0), zz = mx(y[Z2], 0), gc = mx(y[GC], 0);
    const occ = ag / (ag + p[P.k_ag]), act = ct * occ;
    dy[CT] = (p[P.rho] * occ - p[P.k_t]) * ct;
    dy[AG] = -p[P.kappa] * act;
    dy[Z1] = p[P.k_z] * (act - y[Z1]);
    dy[Z2] = p[P.k_z] * (y[Z1] - y[Z2]);
    dy[GC] = p[P.k_g] * (u.gcsf - y[GC]);
    const regrow = p[P.r_s] * (1 - s);
    const killA = p[P.th_as] * cx;
    const killB = (p[P.th_bs] * im + p[P.th_zs] * zz) * (1 + p[P.eta_syn] * mx(1 - s, 0));
    dy[SW] = regrow * sw - (killA + killB) * sw;
    dy[SX] = regrow * sx - (p[P.r_a] * killA + p[P.r_b] * killB) * sx;
    const F = (sw + p[P.delta] * sx) / (1 - p[P.f0] + p[P.delta] * p[P.f0]);
    const ih = Math.pow(im, p[P.hill]);
    const e = p[P.emax] * ih / (ih + Math.pow(p[P.ic50], p[P.hill]));
    const sup = (1 - e) * (1 - p[P.e_z] * zz / (zz + p[P.z50]));
    const supPl = (1 - e) * (1 - p[P.e_z_pl] * zz / (zz + p[P.z50]));
    const kin = p[P.k_p] - p[P.sr_p], kp = p[P.th_ap] * cx;
    const mpp = (x, b, fb, su) => {
      const g = (fb ? Math.min(p[P.gmax], Math.pow(b, -p[P.gamma])) : 1) * su;
      const px = mx(y[x], 0);
      dy[x] = kin * F + (p[P.sr_p] * g - p[P.k_p]) * px - kp * px;
    };
    mpp(PN, bn, true, sup); mpp(PMO, bm, true, sup); mpp(PPL, bp, true, supPl); mpp(PL, bl, false, sup);
    const kl = p[P.th_al] * cx;
    dy[BL] = p[P.c_l] * (mx(y[PL], 0) - y[BL]) - kl * y[BL];
    dy[CX] = u.dose - p[P.k_e] * cx;
    dy[I1] = p[P.k_i] * (d - i1);
    dy[I2] = p[P.k_i] * (i1 - im);
    const gf = 1 + (p[P.gcsf_fold] - 1) * gc, relN = 1 + p[P.gcsf_rel] * gc, tf = u.tpo > 0 ? p[P.tpo_fold] : 1;
    const kc = p[P.kappa_c] * kp;
    // static_cast<int>(x + 0.5) in the C++ engine
    chainN(y, dy, Math.floor(p[P.n_tr_n] + 0.5), Math.floor(p[P.nk_n] + 0.5), p[P.ktr_n], p[P.c_n], mx(y[PN], 0) * gf, kc, relN);
    chain(y, dy, TM1, p[P.ktr_mo], p[P.c_mo], mx(y[PMO], 0), kc, 1);
    chain(y, dy, TP1, p[P.ktr_pl], p[P.c_pl], mx(y[PPL], 0) * tf, kc, 1);
    dy[BM] -= p[P.th_am] * cx * y[BM];
    // P1: inflammation-driven consumption of circulating platelets (th_cp = 0 in Model3)
    dy[BP] -= p[P.th_cp] * d * y[BP];
  };

  C.init = function (p) {
    const y = new Float64Array(NS).fill(1);
    y[SW] = 1 - p[P.f0]; y[SX] = p[P.f0];
    y[CX] = 0; y[I1] = 0; y[I2] = 0; y[CT] = 0; y[AG] = p[P.a0]; y[Z1] = 0; y[Z2] = 0; y[GC] = 0;
    return y;
  };

  // ---- Dormand-Prince RK45 (same coefficients and step control as the C++ engine) ------------------
  const c2 = 1 / 5, c3 = 3 / 10, c4 = 4 / 5, c5 = 8 / 9;
  const a21 = 1 / 5, a31 = 3 / 40, a32 = 9 / 40, a41 = 44 / 45, a42 = -56 / 15, a43 = 32 / 9;
  const a51 = 19372 / 6561, a52 = -25360 / 2187, a53 = 64448 / 6561, a54 = -212 / 729;
  const a61 = 9017 / 3168, a62 = -355 / 33, a63 = 46732 / 5247, a64 = 49 / 176, a65 = -5103 / 18656;
  const b1 = 35 / 384, b3 = 500 / 1113, b4 = 125 / 192, b5 = -2187 / 6784, b6 = 11 / 84;
  const e1 = 71 / 57600, e3 = -71 / 16695, e4 = 71 / 1920, e5 = -17253 / 339200, e6 = 22 / 525, e7 = -1 / 40;
  const K = Array.from({ length: 7 }, () => new Float64Array(NS)), YT = new Float64Array(NS), YN = new Float64Array(NS);

  function integrateTo(ctl, tEnd, y, p, drv, u, dScale) {
    const n = NS, [k1, k2, k3, k4, k5, k6, k7] = K;
    while (ctl.t < tEnd - 1e-12) {
      if (ctl.steps++ > ctl.maxSteps) return false;
      const hh = Math.min(Math.min(ctl.h, ctl.hmax), tEnd - ctl.t), t = ctl.t;
      C.rhs(t, y, k1, p, drv, u, dScale);
      for (let i = 0; i < n; i++) YT[i] = y[i] + hh * a21 * k1[i];
      C.rhs(t + c2 * hh, YT, k2, p, drv, u, dScale);
      for (let i = 0; i < n; i++) YT[i] = y[i] + hh * (a31 * k1[i] + a32 * k2[i]);
      C.rhs(t + c3 * hh, YT, k3, p, drv, u, dScale);
      for (let i = 0; i < n; i++) YT[i] = y[i] + hh * (a41 * k1[i] + a42 * k2[i] + a43 * k3[i]);
      C.rhs(t + c4 * hh, YT, k4, p, drv, u, dScale);
      for (let i = 0; i < n; i++) YT[i] = y[i] + hh * (a51 * k1[i] + a52 * k2[i] + a53 * k3[i] + a54 * k4[i]);
      C.rhs(t + c5 * hh, YT, k5, p, drv, u, dScale);
      for (let i = 0; i < n; i++) YT[i] = y[i] + hh * (a61 * k1[i] + a62 * k2[i] + a63 * k3[i] + a64 * k4[i] + a65 * k5[i]);
      C.rhs(t + hh, YT, k6, p, drv, u, dScale);
      for (let i = 0; i < n; i++) YN[i] = y[i] + hh * (b1 * k1[i] + b3 * k3[i] + b4 * k4[i] + b5 * k5[i] + b6 * k6[i]);
      C.rhs(t + hh, YN, k7, p, drv, u, dScale);
      let err = 0;
      for (let i = 0; i < n; i++) {
        const ei = hh * (e1 * k1[i] + e3 * k3[i] + e4 * k4[i] + e5 * k5[i] + e6 * k6[i] + e7 * k7[i]);
        const sc = ctl.atol + ctl.rtol * Math.max(Math.abs(y[i]), Math.abs(YN[i]));
        err += (ei / sc) * (ei / sc);
      }
      err = Math.sqrt(err / n);
      if (!isFinite(err)) return false;
      if (err <= 1) { ctl.t += hh; for (let i = 0; i < n; i++) y[i] = YN[i]; }
      const fac = err === 0 ? 5 : Math.min(5, Math.max(0.2, 0.9 * Math.pow(err, -0.2)));
      ctl.h = Math.max(hh * fac, 1e-8);
    }
    return true;
  }

  // Simulate on whole-day steps from design.t0 to tEnd (records every day, as the C++ engine does
  // with a daily observation grid). opt: { struct, design, drv, tEnd, dScale: {after, factor} }
  // returns { days, anc_rel, plt_rel, amc_rel, alc_rel, S, f_clone, M, C, I, z, cart } or null
  C.simulate = function (par, opt) {
    const segs = C.segments(opt.design, opt.struct, opt.tEnd);
    const y = C.init(par);
    const ctl = { t: opt.design.t0, h: 0.05, hmax: opt.hmax || 0.5, rtol: opt.rtol || 1e-6, atol: opt.atol || 1e-9, steps: 0, maxSteps: 200000 };
    const keys = ['days', 'anc_rel', 'plt_rel', 'amc_rel', 'alc_rel', 'S', 'f_clone', 'M', 'C', 'I', 'z', 'cart'];
    const out = {}; keys.forEach(k => { out[k] = []; });
    const rec = (t) => {
      const s = y[SW] + y[SX];
      out.days.push(t); out.anc_rel.push(y[BN]); out.plt_rel.push(y[BP]); out.amc_rel.push(y[BM]); out.alc_rel.push(y[BL]);
      out.S.push(s); out.f_clone.push(s > 1e-12 ? y[SX] / s : 0); out.M.push((y[PN] + y[PMO] + y[PPL]) / 3);
      out.C.push(y[CX]); out.I.push(y[I2]); out.z.push(y[Z2]); out.cart.push(y[CT]);
    };
    const u = { dose: 0, gcsf: 0, tpo: 0 };
    let si = 0;
    const t0 = opt.design.t0;
    for (let day = t0; day < opt.tEnd; day++) {
      while (si < segs.length - 1 && segs[si].end <= day) si++;
      const sg = segs[si];
      u.dose = sg.dose; u.gcsf = sg.gcsf; u.tpo = sg.tpo;
      rec(day);
      if (sg.start === day) { if (sg.s_add) y[SW] += sg.s_add; if (sg.car) y[CT] += sg.car; }
      if (!integrateTo(ctl, day + 1, y, par, opt.drv, u, opt.dScale)) return null;
    }
    rec(opt.tEnd);
    return out;
  };

  // ---- recovery phenotype (phenotype.R) -------------------------------------------------------------
  C.rule = { window_end: 30, nadir_window: [0, 14], rec_thr: 1.5, rec_n: 1, dip_thr: 1.0, dip_n: 1, dip_after: 0, drop_spikes: true, cens_value: 0.25 };
  function firstRun(cond, k, from) {
    for (let i = from; i < cond.length; i++) {
      if (i + k - 1 >= cond.length) return -1;
      let ok = true; for (let j = i; j < i + k; j++) if (!cond[j]) { ok = false; break; }
      if (ok) return i;
    }
    return -1;
  }
  // series: [{day, value, cens}]; returns 'quick'|'second_wave'|'delayed'|null
  C.phenotype = function (series, rule) {
    const r = rule || C.rule;
    const s = series.filter(o => o.value !== null && o.value !== undefined && !Number.isNaN(o.value)).slice().sort((a, b) => a.day - b.day);
    let day = [], x = [];
    s.forEach(o => { if (o.day >= 0 && o.day <= r.window_end) { day.push(o.day); x.push(o.cens ? r.cens_value : o.value); } });
    if (!x.length) return null;
    if (r.drop_spikes && x.length >= 2) {
      const drop = new Set();
      for (let i = 1; i < x.length; i++) {
        const nxt = i + 1 < x.length ? x[i + 1] : -Infinity;
        if (x[i - 1] < 1.0 && x[i] >= 1.5 && x[i] >= 2 * x[i - 1] && nxt < 1.0) drop.add(i);
      }
      if (drop.size) { day = day.filter((_, i) => !drop.has(i)); x = x.filter((_, i) => !drop.has(i)); }
    }
    let nad = -1;
    for (let i = 0; i < x.length; i++) if (day[i] >= r.nadir_window[0] && day[i] <= r.nadir_window[1] && (nad < 0 || x[i] < x[nad])) nad = i;
    if (nad < 0) return null;
    const rec = firstRun(x.map(v => v >= r.rec_thr), r.rec_n, nad + 1);
    if (rec < 0) return 'delayed';
    const dip = firstRun(x.map((v, i) => v < r.dip_thr && day[i] >= r.dip_after), r.dip_n, rec + r.rec_n);
    return dip >= 0 ? 'second_wave' : 'quick';
  };
  const Phi = (x) => { // standard normal cdf (Abramowitz-Stegun 7.1.26 on erf, |error| < 1.5e-7)
    const z = Math.abs(x) / Math.SQRT2, t = 1 / (1 + 0.3275911 * z);
    const erf = 1 - (((((1.061405429 * t - 1.453152027) * t) + 1.421413741) * t - 0.284496736) * t + 0.254829592) * t * Math.exp(-z * z);
    return x >= 0 ? 0.5 * (1 + erf) : 0.5 * (1 - erf);
  };
  // phenotype_prob(): margin-based soft classification of a model ANC series
  C.phenotypeProb = function (day, mu, rule, s, eps) {
    const r = rule || C.rule; eps = eps || 0;
    const idx = day.map((_, i) => i).sort((a, b) => day[a] - day[b]);
    const dd = idx.map(i => day[i]), lx = idx.map(i => Math.log(Math.max(mu[i], 1e-8))), n = lx.length;
    const nw = []; dd.forEach((d, i) => { if (d >= r.nadir_window[0] && d <= r.nadir_window[1]) nw.push(i); });
    if (!nw.length) return null;
    let nad = nw[0]; nw.forEach(i => { if (lx[i] < lx[nad]) nad = i; });
    const runStat = (from, k, f) => { const o = []; for (let i = from; i <= n - k; i++) o.push(f(lx.slice(i, i + k))); return o; };
    const rr = runStat(nad + 1, r.rec_n, a => Math.min(...a)).map(v => v - Math.log(r.rec_thr));
    let p;
    if (!rr.length) p = { quick: 0, second_wave: 0, delayed: 1 };
    else {
      const mR = Math.max(...rr);
      const fi = rr.findIndex(v => v >= 0);
      const rStar = nad + 1 + (fi >= 0 ? fi : rr.indexOf(mR));
      const from = rStar + r.rec_n;
      const dm = runStat(from, r.dip_n, a => Math.max(...a)).map((v, j) => ({ v: Math.log(r.dip_thr) - v, d: dd[from + j] })).filter(o => o.d >= r.dip_after);
      const pR = Phi(mR / s), pD = dm.length ? Phi(Math.max(...dm.map(o => o.v)) / s) : 0;
      p = { quick: pR * (1 - pD), second_wave: pR * pD, delayed: 1 - pR };
    }
    for (const k in p) p[k] = (1 - 3 * eps) * p[k] + eps;
    return p;
  };
  C.rank = (ph) => ({ quick: 1, second_wave: 2, delayed: 3 }[ph] || 0);
})(typeof window !== 'undefined' ? window : global);
