/* GreenLight · modules/model3.js
 * MarrowTwin (Model3) as a registered module: the model card (equations per nested model, symbols,
 * outputs) from which the researcher portal is generated, the precomputed as-of prediction
 * (step G, R engine) and the in-browser what-if forecast (JS port, js/engine/model3_core.js).
 * Teal marks what Model3 adds to Model2: (N) the neutrophil chain of n_N compartments, the first k_N
 * dividing (killed by chemotherapy), the rest a post-mitotic storage pool whose release G-CSF speeds;
 * (P2) the platelet MPP component's own CAR T suppression e_Z,Pl. (P1) platelet consumption th_cP x D
 * was screened and not included (th_cP = 0): shown grey (switched off).
 */
(function (root) {
  const GL = root.GL = root.GL || {};
  const C = GL.model3;
  const OFF = (s) => `\\textcolor{#a3aab0}{${s}}`;
  // new in Model3 (neutrophil storage pool N, platelet e_Z P2), in every nested model that contains the term
  const NEW = (s) => `\\textcolor{#0e9fb5}{${s}}`;
  const on = (cond, s) => cond ? s : OFF(s);
  const onNew = (cond, s) => cond ? NEW(s) : OFF(s);
  // design constants of the neutrophil chain and th_cP of the exported population (Model3: 7, 4, 0)
  const pstruct = () => { const i = GL.loader && GL.loader.index && GL.loader.index(); return (i && i.population && i.population.struct) || {}; };
  C.chainN = () => { const st = pstruct(); return { n: st.n_tr_n === undefined ? 7 : +[].concat(st.n_tr_n)[0], k: st.nk_n === undefined ? 4 : +[].concat(st.nk_n)[0], thcp: st.th_cp === undefined ? 0 : +[].concat(st.th_cp)[0] }; };

  // ---- symbols (LaTeX, unit); meaning texts are i18n keys sym.<id> ----------------------------
  const SYM = {
    r_s: ['r_S', '1/d'], f0: ['f_0', '–'], th_as: ['\\theta_{aS}', '1/(d·C)'], th_ap: ['\\theta_{aP}', '1/(d·C)'], th_bs: ['\\theta_{bS}', '1/(d·I)'], rho: ['\\rho_T', '1/d'],
    r_a: ['\\rho^a_X', '–'], r_b: ['\\rho^b_X', '–'], delta: ['\\delta', '–'], eta_syn: ['\\eta', '–'], th_zs: ['\\theta_{zS}', '1/(d·Z)'],
    k_p: ['k_P', '1/d'], sr_p: ['\\sigma_P', '1/d'], gmax: ['G_{\\max}', '–'], gamma: ['\\gamma', '–'],
    kappa_c: ['\\kappa_c', '–'], th_am: ['\\theta_{aM}', '1/(d·C)'], th_cp: ['\\theta_{cP}', '1/(d·D)'],
    n_tr_n: ['n_N', '–'], nk_n: ['k_N', '–'], d_zpl: ['d_{Z,Pl}', '–'],
    mtt_n: ['\\mathrm{MTT}_N', 'd'], mtt_mo: ['\\mathrm{MTT}_{Mo}', 'd'], mtt_pl: ['\\mathrm{MTT}_{Pl}', 'd'],
    c_n: ['c_N', '1/d'], c_mo: ['c_{Mo}', '1/d'], c_pl: ['c_{Pl}', '1/d'], c_l: ['c_L', '1/d'], th_al: ['\\theta_{aL}', '1/(d·C)'],
    k_i: ['k_I', '1/d'], emax: ['E_{\\max}', '–'], ic50: ['\\mathrm{IC}_{50}', 'I'], hill: ['h', '–'],
    k_e: ['k_e', '1/d'], beta_renal: ['\\beta_{\\mathrm{renal}}', '–'], u_br_rgo: ['u_{\\mathrm{br,RGO}}', 'C/d'], u_br_cy: ['u_{\\mathrm{br,Cy}}', 'C/d'],
    gcsf_fold: ['u_G', '–'], gcsf_rel: ['u_R', '–'], k_g: ['k_G', '1/d'], tpo_fold: ['u_T', '–'],
    car_t0: ['C_{T,0}', '10⁹/L'], k_ag: ['K_A', '–'], kappa: ['\\kappa', 'L/(10⁹·d)'], k_t: ['k_T', '1/d'], beta_a: ['\\beta_A', '–'],
    k_z: ['k_Z', '1/d'], z50: ['Z_{50}', 'Z'], e_z: ['e_Z', '–']
  };
  const SOURCE = {
    th_ap: 'calib', th_al: 'calib', th_am: 'calib', kappa_c: 'mpp_equal', c_l: 'calib', emax: 'calib', kappa: 'calib', rho: 'calib', z50: 'calib', th_zs: 'calib',
    k_p: 'aim4', sr_p: 'aim4', gmax: 'aim4', gamma: 'aim4', mtt_n: 'lit', mtt_mo: 'lit', mtt_pl: 'lit', c_n: 'lit', c_mo: 'lit', c_pl: 'lit',
    r_a: 'scenario', r_b: 'scenario', delta: 'scenario',
    n_tr_n: 'screen', nk_n: 'screen', th_cp: 'screen', d_zpl: 'ez_centre'
  };

  // ---- equations of a nested model (spec from C.spec) -------------------------------------------------
  function equations(spec) {
    if (spec.family === 'mf') return [
      { id: 'mf_prol', titleKey: 'eq.mf_prol', descKey: 'eqd.mf_prol', latex: '\\frac{dP_\\ell}{dt}=k\\,P_\\ell\\,(1-\\mathrm{Slope}_\\ell\\,C)\\left(\\frac{B_0}{B_\\ell}\\right)^{\\gamma}-k\\,P_\\ell,\\quad k=\\frac{4}{\\mathrm{MTT}_\\ell\\,s_{\\mathrm{MTT}}}' },
      { id: 'mf_chain', titleKey: 'eq.mf_chain', descKey: 'eqd.mf_chain', latex: '\\frac{dT_{\\ell j}}{dt}=k\\,(T_{\\ell,j-1}-T_{\\ell j}),\\quad \\frac{dB_\\ell}{dt}=k\\,(T_{\\ell 3}-B_\\ell)' },
      { id: 'obs', titleKey: 'eq.obs', descKey: 'eqd.obs_mf', latex: '\\log y = \\log\\big(b_v\\,B_v(t)\\big)+\\varepsilon,\\ \\ \\varepsilon\\sim N(0,\\sigma_v^2),\\ v\\in\\{\\mathrm{ANC},\\mathrm{PLT},\\mathrm{AMC}\\}' },
      { id: 'pheno', titleKey: 'eq.pheno', descKey: 'eqd.pheno', latex: phenoTex() }
    ];
    const a = ['Ma', 'Mab', 'Mabx'].includes(spec.q2), b = ['Mb', 'Mab', 'Mabx'].includes(spec.q2), x = spec.q2 === 'Mabx', cl = spec.clone !== 'none';
    const zh = spec.z_hsc, zm = spec.z_mpp;
    return [
      { id: 'hsc', titleKey: 'eq.hsc', descKey: 'eqd.hsc', latex: `\\frac{dS_j}{dt}=r_S\\,S_j\\,(1-S)-\\Big[${on(cl, '\\rho^a_j')}${on(a, '\\theta_{aS}C')}+${on(cl, '\\rho^b_j')}\\big(${on(b, '\\theta_{bS}I')}+${on(zh, '\\theta_{zS}Z')}\\big)${on(x, '\\big(1+\\eta(1-S)_+\\big)')}\\Big]S_j,\\ \\ S=S_W+${on(cl, 'S_X')}` },
      { id: 'output', titleKey: 'eq.output', descKey: 'eqd.output', latex: `F=\\frac{S_W+${on(cl, '\\delta S_X')}}{1-${on(cl, 'f_0')}+${on(cl, '\\delta f_0')}}` },
      { id: 'mpp', titleKey: 'eq.mpp', descKey: 'eqd.mpp', latex: '\\frac{dP_\\ell}{dt}=(k_P-\\sigma_P)\\,F+(\\sigma_P\\,G_\\ell-k_P)\\,P_\\ell-\\theta_{aP}\\,C\\,P_\\ell,\\quad \\ell\\in\\{N,Mo,Pl,L\\}' },
      { id: 'feedback', titleKey: 'eq.feedback', descKey: 'eqd.feedback', latex: feedbackTex(zm) },
      { id: 'inflam', titleKey: 'eq.inflam', descKey: 'eqd.inflam', latex: '\\frac{dI_1}{dt}=k_I(D-I_1),\\ \\ \\frac{dI}{dt}=k_I(I_1-I),\\ \\ E(I)=\\frac{E_{\\max}I^h}{I^h+\\mathrm{IC}_{50}^h}' },
      { id: 'exposure', titleKey: 'eq.exposure', descKey: 'eqd.exposure', latex: '\\frac{dC}{dt}=u(t)-k_e\\,C,\\quad u_{\\mathrm{LD}}=\\Big(\\frac{90}{\\mathrm{eGFR}}\\Big)^{\\beta_{\\mathrm{renal}}}' },
      { id: 'chain_n', titleKey: 'eq.chain_n', descKey: 'eqd.chain_n', latex: chainNTex() },
      { id: 'chain', titleKey: 'eq.chain', descKey: 'eqd.chain', latex: chainTex() },
      { id: 'lymph', titleKey: 'eq.lymph', descKey: 'eqd.lymph', latex: '\\frac{dB_L}{dt}=c_L\\,(P_L-B_L)-\\theta_{aL}\\,C\\,B_L' },
      { id: 'cart', titleKey: 'eq.cart', descKey: 'eqd.cart', latex: '\\frac{dC_T}{dt}=\\Big(\\rho_T\\frac{A}{A+K_A}-k_T\\Big)C_T,\\ \\frac{dA}{dt}=-\\kappa\\,C_T\\frac{A}{A+K_A},\\ \\frac{dZ_1}{dt}=k_Z\\Big(C_T\\frac{A}{A+K_A}-Z_1\\Big),\\ \\frac{dZ}{dt}=k_Z(Z_1-Z)' },
      { id: 'gcsf', titleKey: 'eq.gcsf', descKey: 'eqd.gcsf', latex: '\\frac{dE_G}{dt}=k_G\\big(g(t)-E_G\\big),\\ \\ u_N=1+(u_G-1)E_G,\\ \\ r_N=1+u_R\\,E_G' },
      { id: 'obs', titleKey: 'eq.obs', descKey: 'eqd.obs', latex: '\\log y=\\log Y_v(t)+\\varepsilon,\\ \\varepsilon\\sim N(0,\\sigma_v^2),\\ \\ Y_v=b_v B_v,\\ \\ Y_{\\mathrm{ALC}}=b_{\\mathrm{ALC}}B_L+C_T' },
      { id: 'pheno', titleKey: 'eq.pheno', descKey: 'eqd.pheno', latex: phenoTex() }
    ];
  }
  // self-renewal feedback; the platelet component has its own CAR T suppression e_Z,Pl (Model3 P2)
  function feedbackTex(zm) {
    const zf = (e) => `\\Big(1-${e}\\frac{Z}{Z+Z_{50}}\\Big)`;
    return '\\begin{gathered}' +
      `G_\\ell=\\min\\!\\big(G_{\\max},B_\\ell^{-\\gamma}\\big)\\,\\big(1-E(I)\\big)\\,${on(zm, zf('e_{Z,\\ell}'))},\\quad G_L=\\big(1-E(I)\\big)\\,${on(zm, zf('e_Z'))}\\\\` +
      `${on(zm, 'e_{Z,N}=e_{Z,Mo}=e_Z')},\\quad ${onNew(zm, 'e_{Z,Pl}=\\mathrm{logit}^{-1}\\big(\\mathrm{logit}\\,e_Z+d_{Z,Pl}\\big)')}` +
      '\\end{gathered}';
  }
  // neutrophil chain (model3.cpp Parent::chain_n): n_N compartments, the first k_N dividing (killed at
  // kappa_c th_aP C), the rest the post-mitotic storage pool whose outflow G-CSF multiplies by r_N(t)
  function chainNTex() {
    const c = C.chainN();
    return '\\begin{gathered}' +
      `\\frac{dT_{Nj}}{dt}=\\kappa_N\\,(O_{j-1}-O_j)-${NEW('\\mathbf{1}[j\\le k_N]')}\\,\\kappa_c\\,\\theta_{aP}\\,C\\,T_{Nj},\\ \\ j=1,\\dots,n_N,\\quad \\frac{dB_N}{dt}=c_N\\,\\big(O_{n_N}-B_N\\big)\\\\` +
      `O_0=P_N\\,u_N(t),\\quad O_j=T_{Nj}\\ (j\\le k_N),\\quad ${NEW('O_j=r_N(t)\\,T_{Nj}\\ (j>k_N)')},\\quad \\kappa_N=\\frac{n_N}{\\mathrm{MTT}_N},\\ \\ ${NEW(`n_N=${c.n},\\ k_N=${c.k}`)}` +
      '\\end{gathered}';
  }
  // monocyte and megakaryocyte chains: T1, T2 proliferating (killed at kappa_c th_aP C), T3 non-dividing;
  // circulating monocytes killed at th_aM C; platelet consumption th_cP D (Model3 P1: screened, th_cP = 0)
  function chainTex() {
    const kc = '-\\kappa_c\\,\\theta_{aP}\\,C';
    return '\\begin{gathered}' +
      `\\frac{dT_{\\ell1}}{dt}=\\kappa_\\ell\\big(P_\\ell\\,u_\\ell(t)-T_{\\ell1}\\big)${kc}\\,T_{\\ell1},\\quad \\frac{dT_{\\ell2}}{dt}=\\kappa_\\ell(T_{\\ell1}-T_{\\ell2})${kc}\\,T_{\\ell2}\\\\` +
      `\\frac{dT_{\\ell3}}{dt}=\\kappa_\\ell\\big(T_{\\ell2}-T_{\\ell3}\\big),\\quad \\kappa_\\ell=\\frac{3}{\\mathrm{MTT}_\\ell},\\ \\ \\ell\\in\\{Mo,Pl\\}\\\\` +
      `\\frac{dB_\\ell}{dt}=c_\\ell\\big(T_{\\ell3}-B_\\ell\\big)-\\mathbf{1}[\\ell=\\mathrm{Mo}]\\,\\theta_{aM}\\,C\\,B_{\\mathrm{Mo}}${on(C.chainN().thcp > 0, '-\\mathbf{1}[\\ell=\\mathrm{Pl}]\\,\\theta_{cP}\\,D(t)\\,B_{\\mathrm{Pl}}')}` +
      '\\end{gathered}';
  }
  function phenoTex() {
    return '-\\,w\\log P\'(c\\mid\\psi),\\ P\'=(1-3\\epsilon)P+\\epsilon,\\ P(\\mathrm{delayed})=1-\\Phi\\Big(\\frac{m_R}{s}\\Big),\\ P(\\mathrm{2nd})=\\Phi\\Big(\\frac{m_R}{s}\\Big)\\Phi\\Big(\\frac{m_D}{s}\\Big)';
  }

  // ---- what-if forecast in the browser ------------------------------------------------------------------
  // a: as-of entry of a patient file; plan: { gcsf: {start, days} | null, reactive: bool,
  //   anti: {from, factor} | null, tpo: {start, end} | null, cd34: {day, add} | null }
  function drawPsi(a, n, seed) {
    const rng = GL.util.rng(seed), k = a.names.length, Lm = [];
    let c = 0;
    for (let j = 0; j < k; j++) for (let i = j; i < k; i++) { (Lm[i] = Lm[i] || [])[j] = a.chol[c++]; }
    const map = a.names.map(nm => [].concat(a.map[nm])[0]);
    const out = [];
    for (let d = 0; d < n; d++) {
      const z = a.names.map(() => d === 0 ? 0 : rng.normal());
      const psi = {};
      a.names.forEach((nm, i) => {
        let v = map[i]; for (let j = 0; j <= i; j++) v += (Lm[i][j] || 0) * z[j];
        psi[nm] = Math.min(a.hi[i], Math.max(a.lo[i], v));
      });
      out.push(psi);
    }
    return out;
  }
  function designAsOf(p, L, plan) {
    const d = p.design;
    const gc = (d.gcsf || []).filter(w => w[0] < L).map(w => [w[0], Math.min(w[1], L)]);
    const des = { t0: d.t0, egfr: d.egfr, bridging: d.bridging, apheresis_day: d.apheresis_day, gcsf: gc, tpo: [], cd34: [] };
    const dec = Math.max(Math.floor(L) + 1, 1);
    if (plan && plan.gcsf) des.gcsf = des.gcsf.concat([[Math.max(plan.gcsf.start, dec), Math.max(plan.gcsf.start, dec) + plan.gcsf.days]]);
    if (plan && plan.tpo) des.tpo = [[Math.max(plan.tpo.start, dec), plan.tpo.end]];
    if (plan && plan.cd34) des.cd34 = [{ day: Math.max(plan.cd34.day, dec), add: plan.cd34.add }];
    return des;
  }
  function observedUpTo(p, L) {
    return (p.obs.anc || []).filter(o => o[3] === 1 && o[0] >= 0 && o[0] <= L).map(o => ({ day: o[0], value: o[1], cens: o[2] === 1 }));
  }
  function simDraw(p, psi, spec, st, des, drv, dScale, tEnd) {
    const ind = C.etaToIndiv(psi, spec);
    ind.a0 = C.a0(st, p.design.ldh_c);
    const par = C.buildPar(ind, st, spec);
    return C.simulate(par, { design: des, struct: st, drv, tEnd, dScale });
  }
  function whatIf(p, a, plan, opt) {
    const idx = GL.loader.index(), pop = idx.population;
    const spec = C.spec(pop.model), st = Object.assign({}, pop.struct);
    const n = (opt && opt.n) || 24, L = a.L, tEnd = 90;
    const drv = { t0: a.drv.t0, dt: a.drv.dt, v: Float64Array.from(a.drv.v) };
    const dScale = plan && plan.anti ? { after: Math.max(L, 0, plan.anti.from - 1e-9), factor: plan.anti.factor } : null;
    const obs = observedUpTo(p, L), draws = drawPsi(a, n, 1009 + Math.round(L * 7));
    const cls = [], ancs = [], plts = [];
    let plt20 = 0;
    for (const psi of draws) {
      const des = designAsOf(p, L, plan);
      const b = { anc: p.base.anc * Math.exp(psi.b_anc || 0), plt: p.base.plt * Math.exp(psi.b_plt || 0) };
      let m = simDraw(p, psi, spec, st, des, drv, dScale, tEnd);
      if (m && plan && plan.reactive) {
        for (let course = 0; course < 2 && m; course++) {
          const from = Math.max(15, Math.floor(L) + 1, ...des.gcsf.map(w => w[1]));
          const i = m.days.findIndex((d, j) => d >= from && d <= 60 && m.anc_rel[j] * b.anc < 1.0);
          if (i < 0) break;
          des.gcsf.push([m.days[i], m.days[i] + 5]);
          m = simDraw(p, psi, spec, st, des, drv, dScale, tEnd);
        }
      }
      if (!m) continue;
      const anc = m.anc_rel.map(x => x * b.anc), plt = m.plt_rel.map(x => x * b.plt);
      const series = obs.concat(m.days.map((d, j) => ({ day: d, value: anc[j], cens: false })).filter(o => o.day > L && o.day <= 30));
      cls.push(C.phenotype(series));
      ancs.push({ days: m.days, v: anc }); plts.push({ days: m.days, v: plt });
      if (m.days.some((d, j) => d >= Math.max(14, L + 1) && d <= 30 && plt[j] < 20)) plt20++;
    }
    const nn = cls.length || 1;
    const med = (arr) => { const days = arr[0] ? arr[0].days : []; return { days, v: days.map((_, j) => GL.util.median(arr.map(x => x.v[j]))) }; };
    return { n: cls.length, p: { quick: cls.filter(c => c === 'quick').length / nn, second_wave: cls.filter(c => c === 'second_wave').length / nn, delayed: cls.filter(c => c === 'delayed').length / nn },
      plt20: plt20 / nn, anc: ancs.length ? med(ancs) : null, plt: plts.length ? med(plts) : null };
  }

  const card = {
    equations, symbols: SYM, source: SOURCE,
    outputs: [
      { id: 'p_quick', labelKey: 'out.p_quick', unit: '%', direction: 'higher' },
      { id: 'p_second_wave', labelKey: 'out.p_second_wave', unit: '%', direction: 'lower' },
      { id: 'p_delayed', labelKey: 'out.p_delayed', unit: '%', direction: 'lower' },
      { id: 'rec_day', labelKey: 'out.rec_day', unit: 'd', direction: 'lower' },
      { id: 'when', labelKey: 'out.when', unit: '%', direction: 'neutral' },
      { id: 'who', labelKey: 'out.who', unit: '%', direction: 'neutral' }],
    models: ['M0', 'Ma', 'Mb', 'Mab', 'Mabx', 'MabC_ddr', 'MabC_dta', 'Mab_Zh', 'Mab_Zm', 'MabxC_ddr', 'MF'],
    // what this version adds (N, P2) and what was screened and not included (P1); search-box limits from stability
    changes: { added: ['N: neutrophil storage pool (n_tr_n, nk_n)', 'P2: platelet e_Z (d_zpl)'], screened_not_included: ['P1: platelet consumption (th_cp = 0)'],
      stability_limits: { gamma: 1.0, mtt_n: 15 } },
    validation: () => (GL.loader.index() || {}).results
  };
  GL.registry.register({
    id: 'marrowtwin', version: 'MarrowTwin (Model3 v1.0)', type: 'mechanistic', nameKey: 'mod.model3', descKey: 'mod.model3.d',
    source: 'MarrowTwin = Model3/Code (R/C++), GreenLight port js/engine/model3_core.js', card,
    // precomputed as-of prediction (R engine)
    predict(snap, p) {
      const a = p && p.asof ? GL.model3.asofEntry(p, snap.asOf) : null;
      return a ? { p: { quick: a.p[0], second_wave: a.p[1], delayed: a.p[2] } } : null;
    },
    whatIf
  });
  GL.model3.asofEntry = function (p, day) {
    let best = null;
    for (const a of p.asof) if (a.L <= day + 1e-9) best = a;
    return best || p.asof[0];
  };
  GL.model3.whatIf = whatIf;
  GL.model3.equations = equations;
})(typeof window !== 'undefined' ? window : global);
