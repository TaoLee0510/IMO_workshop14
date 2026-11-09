/* GreenLight · modules/aim4_core.js
 * Aim 4 parent marrow model v1.2 in JavaScript — a line-by-line port of
 * Aim4/Code/src/aim4_model.cpp (same 24 states, same parameter layout as
 * a4_param_names(), same adaptive Dormand–Prince RK45 integrator), plus the
 * inputs (Aim4/Code/R/model.R make_segments), the CAR T-phase driver and the
 * recovery-phenotype rule (Aim4/Code/R/phenotype.R, primary rule + late class).
 * Pure functions: used by modules/marrow4.js in the browser and by the node
 * consistency test (tools/test_aim4_ode.js) against R/C++ reference output.
 * With opt.track the age / generation moments of the C++ model (track = TRUE) are integrated as
 * well: mean time and mean number of divisions since MPP entry of the circulating cells.
 *
 * Model v1.2 (lineage tree of Aim4/image (2).png; every compartment normalised to 1)
 *   C   dC/dt = dose(t) - k_e C                                 (bridging + flu-cy)
 *   I1, I  dI1/dt = k_I (D - I1),  dI/dt = k_I (I1 - I);  E(I) = Emax I^h / (I^h + IC50^h)
 *   S_W, S_X  dS_j/dt = r_S S_j (1-S)(1 + alpha (1-MPP)+) - [ra_j th_aS C + rb_j th_bS I (1 + eta (1-S)+)] S_j
 *   X in {MPP, CMP, GMP, MEP, CLP}:  dX/dt = (k_X - sr_X) U_X + (sr_X G_X - k_X) X - th_aP C X
 *     U: MPP <- F, CMP, CLP <- MPP, GMP, MEP <- CMP; G_X = min(gmax, s_X^-gamma_M) (1 - e_sr E) for
 *     MPP, CMP, GMP, MEP (s_X = geometric mean of the progeny's blood counts), G_CLP from ALC
 *   chains l: P_l -> T1 -> T2 -> T3 -> B_l (P_N = P_Mo = GMP, P_Pl = MEP), ktr = 3 / MTT
 *     dT1/dt = ktr [P_l B_l^-gam_l (1 - e_l E) u_l - T1]
 *   lymphocytes: dL/dt = (k_L - sr_L) CLP (1 - e_L E) + (sr_L G_L - k_L) L - th_aL C L,
 *     dB_L/dt = c_L (L - B_L) - th_aL C B_L
 */
(function (root) {
  const GL = root.GL = root.GL || {};
  const C = GL.aim4 = {};

  C.PAR = ['r_s', 'f0', 'th_as', 'th_ap', 'th_bs', 'alpha',
    'k_mpp', 'sr_mpp', 'k_cmp', 'sr_cmp', 'k_gmp', 'sr_gmp', 'k_mep', 'sr_mep', 'k_clp', 'sr_clp', 'k_ly', 'sr_ly',
    'gmax', 'gamma_m', 'fb_mpp', 'fb_cmp', 'fb_gmp', 'fb_mep', 'fb_clp', 'fb_ly',
    'r_a', 'r_b', 'delta', 'ktr_n', 'ktr_mo', 'ktr_pl', 'c_n', 'c_mo', 'c_pl', 'c_l',
    'gam_n', 'gam_mo', 'gam_pl', 'k_i', 'emax', 'ic50', 'hill', 'e_mo', 'e_pl', 'e_l', 'eta_syn', 'e_sr',
    'th_al', 'k_e', 'gcsf_fold', 'tpo_fold', 'd_n', 'd_mo', 'd_pl'];
  const P = {}; C.PAR.forEach((k, i) => { P[k] = i; });
  C.IX = P;
  const SW = 0, SX = 1, MPP = 2, CMP = 3, GMP = 4, MEP = 5, CLP = 6, LY = 7, CX = 8, I1 = 9, I2 = 10,
    TN1 = 11, BN = 14, TM1 = 15, BM = 18, TP1 = 19, BP = 22, BL = 23, NS = 24;
  C.NS = NS;
  C.IDX = { SW, SX, MPP, CMP, GMP, MEP, CLP, LY, CX, I1, I: I2, BN, BM, BP, BL };
  // compartments with age / generation moments (track): age mass at NS + q, generation mass at NS + NTRK + q
  const Q_MPP = 0, Q_CMP = 1, Q_GMP = 2, Q_MEP = 3, Q_CLP = 4, Q_LY = 5, Q_TN1 = 6, Q_BN = 9, Q_TM1 = 10, Q_BM = 13,
    Q_TP1 = 14, Q_BP = 17, Q_BL = 18, NTRK = 19;
  C.NTRK = NTRK;
  C.RE = ['r_s', 'f0', 'th_as', 'th_ap', 'th_bs'];

  // ---- model specification (mirrors model_spec() in model.R) ---------------------------
  C.spec = function (id) {
    const m = /^(M0|Ma|Mb|Mab|Mabx)(C_(ddr|dta))?$/.exec(id || 'MabC_ddr');
    if (!m) throw new Error('unknown Aim 4 model id ' + id);
    const q2 = m[1], clone = m[3] || 'none';
    const re = ['r_s', 'th_ap'];
    if (['Ma', 'Mab', 'Mabx'].includes(q2)) re.push('th_as');
    if (['Mb', 'Mab', 'Mabx'].includes(q2)) re.push('th_bs');
    if (clone !== 'none') re.push('f0');
    return { id, q2, clone, re };
  };

  // full parameter vector; terms the model does not contain are zero (nesting switches)
  C.buildPar = function (indiv, struct, spec) {
    const s = struct, on = (k) => spec.re.includes(k) ? indiv[k] : 0;
    ['mpp', 'cmp', 'gmp', 'mep', 'clp', 'ly'].forEach(x => {
      if (!(s['sr_' + x] < s['k_' + x])) throw new Error('self-renewal rate must be below the exit rate: ' + x);
    });
    const v = Object.assign({}, s, {
      r_s: indiv.r_s, f0: on('f0'), th_as: on('th_as'), th_ap: indiv.th_ap, th_bs: on('th_bs'),
      ktr_n: 3 / s.mtt_n, ktr_mo: 3 / s.mtt_mo, ktr_pl: 3 / s.mtt_pl,
      eta_syn: spec.q2 === 'Mabx' ? s.eta_syn : 0
    });
    const p = new Float64Array(C.PAR.length);
    C.PAR.forEach((k, i) => {
      if (typeof v[k] !== 'number' || !isFinite(v[k])) throw new Error('missing Aim 4 parameter ' + k);
      p[i] = v[k];
    });
    return p;
  };

  // ---- inputs (mirrors make_segments() in model.R; all breakpoints are whole days) ------
  // design: { t0, egfr, bridging: 'none'|'rgo'|'cy', apheresis_day, ld_scale,
  //           gcsf: [[start, end)], tpo: [[start, end)], cd34: [{day, add}] }
  C.bridgingWindow = [2, -7];
  C.segments = function (design, struct, tEnd) {
    const br = design.bridging === 'rgo' ? struct.u_br_rgo : design.bridging === 'cy' ? struct.u_br_cy : 0;
    const bw = br > 0 ? [design.apheresis_day + C.bridgingWindow[0], C.bridgingWindow[1]] : null;
    const ld = Math.pow(90 / design.egfr, struct.beta_renal) * (design.ld_scale === undefined ? 1 : design.ld_scale);
    const gc = design.gcsf || [], tp = design.tpo || [], cd = design.cd34 || [];
    let bp = [design.t0, tEnd, -5, -2];
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
      out.push({ start: st, end: en, dose, gcsf: inWin(gc, mid) ? 1 : 0, tpo: inWin(tp, mid) ? 1 : 0, s_add: sAdd });
    }
    return out;
  };

  // ---- driver ------------------------------------------------------------------------
  C.drvAt = function (drv, t) {
    const x = (t - drv.t0) / drv.dt, n = drv.v.length;
    if (x <= 0) return drv.v[0];
    if (x >= n - 1) return drv.v[n - 1];
    const i = Math.floor(x), w = x - i;
    return (1 - w) * drv.v[i] + w * drv.v[i + 1];
  };
  C.kernel = (t, A, tp, h) => t > 0 ? A * Math.pow(t / tp, h) * Math.exp(h * (1 - t / tp)) : 0;

  // ---- right-hand side (aim4_model.cpp Parent::rhs, track = false) ------------------------
  function chain(y, dy, t1, k, c, inflow) {
    dy[t1] = k * (inflow - y[t1]);
    dy[t1 + 1] = k * (y[t1] - y[t1 + 1]);
    dy[t1 + 2] = k * (y[t1 + 1] - y[t1 + 2]);
    dy[t1 + 3] = c * (y[t1 + 2] - y[t1 + 3]);
  }
  // self-renewal feedback from the blood counts of the compartment's progeny
  function feedback(p, signal, on) {
    return on > 0 ? Math.min(p[P.gmax], Math.pow(Math.max(signal, 1e-6), -p[P.gamma_m])) : 1;
  }
  C.rhs = function (t, y, dy, p, drv, u, dScale, track) {
    const sw = Math.max(y[SW], 0), sx = Math.max(y[SX], 0), s = sw + sx;
    const mpp = Math.max(y[MPP], 0), cmp = Math.max(y[CMP], 0), gmp = Math.max(y[GMP], 0);
    const mep = Math.max(y[MEP], 0), clp = Math.max(y[CLP], 0), ly = Math.max(y[LY], 0);
    const cx = Math.max(y[CX], 0), i1 = Math.max(y[I1], 0), im = Math.max(y[I2], 0);
    let d = Math.max(C.drvAt(drv, t), 0);
    if (dScale && t > dScale.after) d *= dScale.factor;
    const bn = Math.max(y[BN], 1e-6), bm = Math.max(y[BM], 1e-6), bp = Math.max(y[BP], 1e-6), bl = Math.max(y[BL], 1e-6);
    // HSC
    const regrow = p[P.r_s] * (1 - s) * (1 + p[P.alpha] * Math.max(1 - mpp, 0));
    const killA = p[P.th_as] * cx;
    const killB = p[P.th_bs] * im * (1 + p[P.eta_syn] * Math.max(1 - s, 0));
    dy[SW] = regrow * sw - (killA + killB) * sw;
    dy[SX] = regrow * sx - (p[P.r_a] * killA + p[P.r_b] * killB) * sx;
    const F = (sw + p[P.delta] * sx) / (1 - p[P.f0] + p[P.delta] * p[P.f0]);
    // marrow inflammation
    const ih = Math.pow(im, p[P.hill]);
    const e = p[P.emax] * ih / (ih + Math.pow(p[P.ic50], p[P.hill]));
    // self-renewing progenitors
    const my = Math.cbrt(bn * bm * bp), srInf = 1 - p[P.e_sr] * e;
    const gMpp = feedback(p, my, p[P.fb_mpp]) * srInf, gCmp = feedback(p, my, p[P.fb_cmp]) * srInf;
    const gGmp = feedback(p, Math.sqrt(bn * bm), p[P.fb_gmp]) * srInf, gMep = feedback(p, bp, p[P.fb_mep]) * srInf;
    const gClp = feedback(p, bl, p[P.fb_clp]), gLy = feedback(p, bl, p[P.fb_ly]);
    const kp = p[P.th_ap] * cx;
    dy[MPP] = (p[P.k_mpp] - p[P.sr_mpp]) * F + (p[P.sr_mpp] * gMpp - p[P.k_mpp]) * mpp - kp * mpp;
    dy[CMP] = (p[P.k_cmp] - p[P.sr_cmp]) * mpp + (p[P.sr_cmp] * gCmp - p[P.k_cmp]) * cmp - kp * cmp;
    dy[GMP] = (p[P.k_gmp] - p[P.sr_gmp]) * cmp + (p[P.sr_gmp] * gGmp - p[P.k_gmp]) * gmp - kp * gmp;
    dy[MEP] = (p[P.k_mep] - p[P.sr_mep]) * cmp + (p[P.sr_mep] * gMep - p[P.k_mep]) * mep - kp * mep;
    dy[CLP] = (p[P.k_clp] - p[P.sr_clp]) * mpp + (p[P.sr_clp] * gClp - p[P.k_clp]) * clp - kp * clp;
    // lymphocytes
    const kl = p[P.th_al] * cx, fLy = 1 - p[P.e_l] * e;
    dy[LY] = (p[P.k_ly] - p[P.sr_ly]) * clp * fLy + (p[P.sr_ly] * gLy - p[P.k_ly]) * ly - kl * ly;
    dy[BL] = p[P.c_l] * (ly - y[BL]) - kl * y[BL];
    // exposure and marrow lag
    dy[CX] = u.dose - p[P.k_e] * cx;
    dy[I1] = p[P.k_i] * (d - i1);
    dy[I2] = p[P.k_i] * (i1 - im);
    // maturation chains
    const gf = u.gcsf > 0 ? p[P.gcsf_fold] : 1, tf = u.tpo > 0 ? p[P.tpo_fold] : 1;
    const fN = Math.pow(bn, -p[P.gam_n]) * (1 - e) * gf, fMo = Math.pow(bm, -p[P.gam_mo]) * (1 - p[P.e_mo] * e);
    const fPl = Math.pow(bp, -p[P.gam_pl]) * (1 - p[P.e_pl] * e) * tf;
    chain(y, dy, TN1, p[P.ktr_n], p[P.c_n], gmp * fN);
    chain(y, dy, TM1, p[P.ktr_mo], p[P.c_mo], gmp * fMo);
    chain(y, dy, TP1, p[P.ktr_pl], p[P.c_pl], mep * fPl);
    if (!track) return;

    // ---- age (time since MPP entry) and generation moments (aim4_model.cpp, track = true) ----
    const A0 = NS, G0 = NS + NTRK;
    // self-renewing compartment q (size x): inflow coefficient `inc` times the upstream masses and size;
    // one division on entry (ent = 1, 0 for the MPP)
    const renewing = (q, x, k, sr, gx, kill, inc, aU, gU, Uz, ent) => {
      const r = sr * gx, A = y[A0 + q], G = y[G0 + q];
      dy[A0 + q] = x + inc * aU + (r - k - kill) * A;
      dy[G0 + q] = inc * (gU + ent * Uz) + r * (G + 2 * x) - (k + kill) * G;
    };
    const a = (q) => y[A0 + q], g = (q) => y[G0 + q];
    renewing(Q_MPP, mpp, p[P.k_mpp], p[P.sr_mpp], gMpp, kp, p[P.k_mpp] - p[P.sr_mpp], 0, 0, F, 0);
    renewing(Q_CMP, cmp, p[P.k_cmp], p[P.sr_cmp], gCmp, kp, p[P.k_cmp] - p[P.sr_cmp], a(Q_MPP), g(Q_MPP), mpp, 1);
    renewing(Q_GMP, gmp, p[P.k_gmp], p[P.sr_gmp], gGmp, kp, p[P.k_gmp] - p[P.sr_gmp], a(Q_CMP), g(Q_CMP), cmp, 1);
    renewing(Q_MEP, mep, p[P.k_mep], p[P.sr_mep], gMep, kp, p[P.k_mep] - p[P.sr_mep], a(Q_CMP), g(Q_CMP), cmp, 1);
    renewing(Q_CLP, clp, p[P.k_clp], p[P.sr_clp], gClp, kp, p[P.k_clp] - p[P.sr_clp], a(Q_MPP), g(Q_MPP), mpp, 1);
    renewing(Q_LY, ly, p[P.k_ly], p[P.sr_ly], gLy, kl, (p[P.k_ly] - p[P.sr_ly]) * fLy, a(Q_CLP), g(Q_CLP), clp, 1);
    dy[A0 + Q_BL] = y[BL] + p[P.c_l] * a(Q_LY) - (p[P.c_l] + kl) * a(Q_BL);
    dy[G0 + Q_BL] = p[P.c_l] * g(Q_LY) - (p[P.c_l] + kl) * g(Q_BL);
    // maturation chain: entry division (+1), d/2 divisions on each of T1->T2, T2->T3
    const mchain = (q1, t1, k, c, f, qp, Pz, dl) => {
      dy[A0 + q1] = y[t1] + k * f * a(qp) - k * a(q1);
      dy[G0 + q1] = k * f * (g(qp) + Pz) - k * g(q1);
      dy[A0 + q1 + 1] = y[t1 + 1] + k * a(q1) - k * a(q1 + 1);
      dy[G0 + q1 + 1] = k * (g(q1) + 0.5 * dl * y[t1]) - k * g(q1 + 1);
      dy[A0 + q1 + 2] = y[t1 + 2] + k * a(q1 + 1) - k * a(q1 + 2);
      dy[G0 + q1 + 2] = k * (g(q1 + 1) + 0.5 * dl * y[t1 + 1]) - k * g(q1 + 2);
      dy[A0 + q1 + 3] = y[t1 + 3] + c * a(q1 + 2) - c * a(q1 + 3);
      dy[G0 + q1 + 3] = c * g(q1 + 2) - c * g(q1 + 3);
    };
    mchain(Q_TN1, TN1, p[P.ktr_n], p[P.c_n], fN, Q_GMP, gmp, p[P.d_n]);
    mchain(Q_TM1, TM1, p[P.ktr_mo], p[P.c_mo], fMo, Q_GMP, gmp, p[P.d_mo]);
    mchain(Q_TP1, TP1, p[P.ktr_pl], p[P.c_pl], fPl, Q_MEP, mep, p[P.d_pl]);
  };

  // pre-treatment steady state; with track the moments start at their closed forms
  //   renewing compartment: age = a_in + 1/(k - sr), generation = g_in + entry + 2 sr/(k - sr)
  C.init = function (p, track) {
    const y = new Float64Array(track ? NS + 2 * NTRK : NS).fill(1);
    y[SW] = 1 - p[P.f0]; y[SX] = p[P.f0]; y[CX] = 0; y[I1] = 0; y[I2] = 0;
    if (!track) return y;
    const A0 = NS, G0 = NS + NTRK;
    const renewing = (q, k, sr, aU, gU, ent) => { y[A0 + q] = aU + 1 / (k - sr); y[G0 + q] = gU + ent + 2 * sr / (k - sr); };
    renewing(Q_MPP, p[P.k_mpp], p[P.sr_mpp], 0, 0, 0);
    renewing(Q_CMP, p[P.k_cmp], p[P.sr_cmp], y[A0 + Q_MPP], y[G0 + Q_MPP], 1);
    renewing(Q_GMP, p[P.k_gmp], p[P.sr_gmp], y[A0 + Q_CMP], y[G0 + Q_CMP], 1);
    renewing(Q_MEP, p[P.k_mep], p[P.sr_mep], y[A0 + Q_CMP], y[G0 + Q_CMP], 1);
    renewing(Q_CLP, p[P.k_clp], p[P.sr_clp], y[A0 + Q_MPP], y[G0 + Q_MPP], 1);
    renewing(Q_LY, p[P.k_ly], p[P.sr_ly], y[A0 + Q_CLP], y[G0 + Q_CLP], 1);
    y[A0 + Q_BL] = y[A0 + Q_LY] + 1 / p[P.c_l]; y[G0 + Q_BL] = y[G0 + Q_LY];
    const mchain = (q1, k, c, qp, dl) => {
      y[A0 + q1] = y[A0 + qp] + 1 / k; y[G0 + q1] = y[G0 + qp] + 1;
      y[A0 + q1 + 1] = y[A0 + q1] + 1 / k; y[G0 + q1 + 1] = y[G0 + q1] + 0.5 * dl;
      y[A0 + q1 + 2] = y[A0 + q1 + 1] + 1 / k; y[G0 + q1 + 2] = y[G0 + q1 + 1] + 0.5 * dl;
      y[A0 + q1 + 3] = y[A0 + q1 + 2] + 1 / c; y[G0 + q1 + 3] = y[G0 + q1 + 2];
    };
    mchain(Q_TN1, p[P.ktr_n], p[P.c_n], Q_GMP, p[P.d_n]);
    mchain(Q_TM1, p[P.ktr_mo], p[P.c_mo], Q_GMP, p[P.d_mo]);
    mchain(Q_TP1, p[P.ktr_pl], p[P.c_pl], Q_MEP, p[P.d_pl]);
    return y;
  };

  // ---- Dormand–Prince RK45 (same tableau and step control as the C++ integrator) -----------
  const c2 = 1 / 5, c3 = 3 / 10, c4 = 4 / 5, c5 = 8 / 9;
  const a21 = 1 / 5, a31 = 3 / 40, a32 = 9 / 40, a41 = 44 / 45, a42 = -56 / 15, a43 = 32 / 9;
  const a51 = 19372 / 6561, a52 = -25360 / 2187, a53 = 64448 / 6561, a54 = -212 / 729;
  const a61 = 9017 / 3168, a62 = -355 / 33, a63 = 46732 / 5247, a64 = 49 / 176, a65 = -5103 / 18656;
  const b1 = 35 / 384, b3 = 500 / 1113, b4 = 125 / 192, b5 = -2187 / 6784, b6 = 11 / 84;
  const e1 = 71 / 57600, e3 = -71 / 16695, e4 = 71 / 1920, e5 = -17253 / 339200, e6 = 22 / 525, e7 = -1 / 40;
  const NMAX = NS + 2 * NTRK;
  const K = Array.from({ length: 7 }, () => new Float64Array(NMAX)), YT = new Float64Array(NMAX), YN = new Float64Array(NMAX);

  // integrate y in place from ctl.t to tEnd; returns false on failure
  function integrateTo(ctl, tEnd, y, p, drv, u, dScale) {
    const n = y.length, tr = n > NS, [k1, k2, k3, k4, k5, k6, k7] = K;
    while (ctl.t < tEnd - 1e-12) {
      if (ctl.steps++ > ctl.maxSteps) return false;
      const hh = Math.min(ctl.h, ctl.hmax, tEnd - ctl.t), t = ctl.t;
      C.rhs(t, y, k1, p, drv, u, dScale, tr);
      for (let i = 0; i < n; i++) YT[i] = y[i] + hh * a21 * k1[i];
      C.rhs(t + c2 * hh, YT, k2, p, drv, u, dScale, tr);
      for (let i = 0; i < n; i++) YT[i] = y[i] + hh * (a31 * k1[i] + a32 * k2[i]);
      C.rhs(t + c3 * hh, YT, k3, p, drv, u, dScale, tr);
      for (let i = 0; i < n; i++) YT[i] = y[i] + hh * (a41 * k1[i] + a42 * k2[i] + a43 * k3[i]);
      C.rhs(t + c4 * hh, YT, k4, p, drv, u, dScale, tr);
      for (let i = 0; i < n; i++) YT[i] = y[i] + hh * (a51 * k1[i] + a52 * k2[i] + a53 * k3[i] + a54 * k4[i]);
      C.rhs(t + c5 * hh, YT, k5, p, drv, u, dScale, tr);
      for (let i = 0; i < n; i++) YT[i] = y[i] + hh * (a61 * k1[i] + a62 * k2[i] + a63 * k3[i] + a64 * k4[i] + a65 * k5[i]);
      C.rhs(t + hh, YT, k6, p, drv, u, dScale, tr);
      for (let i = 0; i < n; i++) YN[i] = y[i] + hh * (b1 * k1[i] + b3 * k3[i] + b4 * k4[i] + b5 * k5[i] + b6 * k6[i]);
      C.rhs(t + hh, YN, k7, p, drv, u, dScale, tr);
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

  // Simulate on whole-day steps from design.t0 to tEnd.
  // opt: { struct, design, drv, tEnd, rtol, atol, hmax, dScale:{after,factor},
  //        hook(day, y, u) -> called at the start of every day > hookFrom (may switch u.gcsf/u.tpo, add to y) }
  // opt.track: also { age_anc, age_plt, age_amc, age_alc, gen_anc, gen_plt, gen_amc, gen_alc } (mean days
  //   and divisions since MPP entry of the circulating cells)
  // returns { days, anc_rel, plt_rel, amc_rel, alc_rel, S, f_clone, M (= MPP), C, D, I }  (null on failure)
  C.simulate = function (par, opt) {
    const segs = C.segments(opt.design, opt.struct, opt.tEnd);
    const y = C.init(par, !!opt.track);
    const ctl = { t: opt.design.t0, h: 0.05, hmax: opt.hmax || 0.5, rtol: opt.rtol || 1e-6, atol: opt.atol || 1e-9, steps: 0, maxSteps: 200000 };
    const out = { days: [], anc_rel: [], plt_rel: [], amc_rel: [], alc_rel: [], S: [], f_clone: [], M: [], C: [], D: [], I: [] };
    const TRK = [['anc', Q_BN, BN], ['plt', Q_BP, BP], ['amc', Q_BM, BM], ['alc', Q_BL, BL]];
    if (opt.track) TRK.forEach(([v]) => { out['age_' + v] = []; out['gen_' + v] = []; });
    const rec = (t) => {
      const s = y[SW] + y[SX];
      out.days.push(t); out.anc_rel.push(y[BN]); out.plt_rel.push(y[BP]); out.amc_rel.push(y[BM]); out.alc_rel.push(y[BL]);
      out.S.push(s); out.f_clone.push(s > 1e-12 ? y[SX] / s : 0); out.M.push(y[MPP]); out.C.push(y[CX]);
      let d = C.drvAt(opt.drv, t); if (opt.dScale && t > opt.dScale.after) d *= opt.dScale.factor;
      out.D.push(d); out.I.push(y[I2]);
      if (opt.track) TRK.forEach(([v, q, xb]) => { const x = Math.max(y[xb], 1e-12); out['age_' + v].push(y[NS + q] / x); out['gen_' + v].push(y[NS + NTRK + q] / x); });
    };
    const u = { dose: 0, gcsf: 0, tpo: 0 };
    let si = 0;
    const t0 = opt.design.t0;
    for (let day = t0; day < opt.tEnd; day++) {
      while (si < segs.length - 1 && segs[si].end <= day) si++;
      const sg = segs[si];
      u.dose = sg.dose; u.gcsf = sg.gcsf; u.tpo = sg.tpo;
      // record first (as the C++ integrator does at a segment boundary), then apply jumps
      if (day === t0 || !opt.recordFrom || day >= opt.recordFrom) rec(day);
      if (sg.start === day && sg.s_add) y[SW] += sg.s_add;
      if (opt.hook) opt.hook(day, y, u);
      if (!integrateTo(ctl, day + 1, y, par, opt.drv, u, opt.dScale)) return null;
    }
    rec(opt.tEnd);
    return out;
  };

  // ---- CAR T-phase driver from ferritin ------------------------------------------------------
  // Whittaker smoother (second differences) of log(ferritin / baseline) on a daily grid from
  // day -3 to the last observation; positive part, zero before day 0; after the last
  // observation either exponential decay (decay) or the Gamma kernel continued smoothly.
  function whittaker(days, vals, lambda) {
    const lo = Math.min(...days), hi = Math.max(...days), n = hi - lo + 1;
    const w = new Float64Array(n), z = new Float64Array(n);
    for (let i = 0; i < days.length; i++) { w[days[i] - lo] += 1; z[days[i] - lo] += vals[i]; }
    for (let i = 0; i < n; i++) if (w[i] > 0) z[i] /= w[i];
    if (n < 3) return { lo, v: Array.from(z) };
    // A = diag(w) + lambda D'D (pentadiagonal) -> dense solve (n <= ~40)
    const A = Array.from({ length: n }, () => new Float64Array(n)), b = new Float64Array(n);
    for (let i = 0; i < n; i++) { A[i][i] += w[i]; b[i] = w[i] * z[i]; }
    for (let i = 0; i + 2 < n; i++) {
      const d = [1, -2, 1];
      for (let a = 0; a < 3; a++) for (let c = 0; c < 3; c++) A[i + a][i + c] += lambda * d[a] * d[c];
    }
    for (let k = 0; k < n; k++) {
      let piv = k; for (let i = k + 1; i < n; i++) if (Math.abs(A[i][k]) > Math.abs(A[piv][k])) piv = i;
      [A[k], A[piv]] = [A[piv], A[k]]; [b[k], b[piv]] = [b[piv], b[k]];
      const akk = A[k][k] || 1e-12;
      for (let i = k + 1; i < n; i++) { const f = A[i][k] / akk; if (!f) continue; for (let j = k; j < n; j++) A[i][j] -= f * A[k][j]; b[i] -= f * b[k]; }
    }
    const x = new Float64Array(n);
    for (let i = n - 1; i >= 0; i--) { let s = b[i]; for (let j = i + 1; j < n; j++) s -= A[i][j] * x[j]; x[i] = s / (A[i][i] || 1e-12); }
    return { lo, v: Array.from(x) };
  }
  C.smoothFerritin = function (obs, base, upTo) {
    const pts = obs.filter(o => o.day >= -3 && o.day <= Math.min(35, upTo) && o.value > 0);
    if (pts.length < 2) return null;
    const days = pts.map(o => Math.round(o.day)), vals = pts.map(o => Math.log(o.value / base));
    const sm = whittaker(days, vals, 4);
    return { lo: sm.lo, hi: sm.lo + sm.v.length - 1, v: sm.v };
  };
  // least-squares Gamma kernel through the observed log-ratios (days > 0); null if < 3 points
  C.fitKernel = function (obs, base, upTo, h) {
    const pts = obs.filter(o => o.day > 0 && o.day <= upTo && o.value > 0).map(o => ({ t: o.day, r: Math.log(o.value / base) }));
    if (pts.length < 3) return null;
    let best = null;
    for (let tp = 2; tp <= 25; tp += 0.25) {
      let sk = 0, skk = 0;
      for (const q of pts) { const k = C.kernel(q.t, 1, tp, h); sk += q.r * k; skk += k * k; }
      const A = Math.max(0.02, skk > 0 ? sk / skk : 0.02);
      let sse = 0; for (const q of pts) { const e = q.r - C.kernel(q.t, A, tp, h); sse += e * e; }
      if (!best || sse < best.sse) best = { A, tp, sse };
    }
    return best;
  };
  // driver grid from t0 to tEnd (dt 0.25): smoothed observations up to upTo, kernel (A, tp) after
  C.buildDriver = function (opt) {
    const dt = 0.25, n = Math.round((opt.tEnd - opt.t0) / dt) + 1, v = new Float64Array(n);
    const sm = opt.smooth, h = opt.h || 2;
    const smAt = (t) => { if (!sm) return 0; const x = t - sm.lo; if (x <= 0) return sm.v[0]; if (x >= sm.v.length - 1) return sm.v[sm.v.length - 1]; const i = Math.floor(x), w = x - i; return (1 - w) * sm.v[i] + w * sm.v[i + 1]; };
    const last = sm ? Math.min(sm.hi, opt.upTo) : 0;
    const dLast = sm && last > 0 ? Math.max(0, smAt(last)) : 0;
    const kLast = opt.kernel ? C.kernel(last, opt.kernel.A, opt.kernel.tp, h) : 0;
    for (let i = 0; i < n; i++) {
      const t = opt.t0 + i * dt;
      if (t <= 0) { v[i] = 0; continue; }
      if (sm && t <= last) { v[i] = Math.max(0, smAt(t)); continue; }
      if (opt.kernel) v[i] = Math.max(0, C.kernel(t, opt.kernel.A, opt.kernel.tp, h) + (dLast - kLast) * Math.exp(-0.3 * (t - last)));
      else v[i] = dLast * Math.exp(-(opt.decay || 0.07) * (t - last));
    }
    return { t0: opt.t0, dt, v };
  };

  // ---- recovery phenotype (phenotype.R primary rule; late class for simulations) ----------
  C.rule = { window_end: 30, nadir_window: [0, 14], rec_thr: 1.5, rec_n: 1, dip_thr: 1.0, dip_n: 1, dip_after: 0, drop_spikes: true, cens_value: 0.25 };
  C.persistent = { window: [45, 75], thr: 0.5 };
  function firstRun(cond, k, from) {
    for (let i = from; i < cond.length; i++) {
      if (i + k - 1 >= cond.length) return -1;
      let ok = true; for (let j = i; j < i + k; j++) if (!cond[j]) { ok = false; break; }
      if (ok) return i;
    }
    return -1;
  }
  const median = (a) => { const s = a.slice().sort((x, y) => x - y), n = s.length; return n % 2 ? s[(n - 1) / 2] : (s[n / 2 - 1] + s[n / 2]) / 2; };
  // series: [{day, value, cens}] (any order); returns 'quick'|'second_wave'|'delayed'|'persistent'|null
  C.phenotype = function (series, withPersistent) {
    const r = C.rule;
    const s = series.filter(o => o.value !== null && o.value !== undefined && !Number.isNaN(o.value)).slice().sort((a, b) => a.day - b.day);
    if (withPersistent) {
      const w = s.filter(o => o.day >= C.persistent.window[0] && o.day <= C.persistent.window[1]).map(o => o.value);
      if (w.length >= 3 && median(w) < C.persistent.thr) return 'persistent';
    }
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
  C.rank = (ph) => ({ quick: 1, second_wave: 2, delayed: 3, persistent: 4 }[ph] || 0);
  C.whenClass = function (act, noA, noB, noAB) {
    if (act === 'quick') return 'none';
    const ra = C.rank(act), ia = C.rank(noA) < ra, ib = C.rank(noB) < ra, iab = C.rank(noAB) < ra;
    if (ia && !ib) return 'a'; if (ib && !ia) return 'b'; if (ia && ib) return 'ab_joint'; if (iab) return 'ab_redundant';
    return 'other';
  };
  C.whoClass = function (act, noClone, reserve) {
    if (act === 'quick') return null;
    const ra = C.rank(act), ic = C.rank(noClone) < ra, ir = C.rank(reserve) < ra;
    return ic && ir ? 'both' : ic ? 'clone' : ir ? 'reserve' : 'other';
  };
})(typeof window !== 'undefined' ? window : global);
