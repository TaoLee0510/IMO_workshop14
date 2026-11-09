/* GreenLight · ui/researcher.js
 * Researcher portal, generated from the registered model cards (GL.registry) and the results
 * bundle (data/model1_index.js): Model (clickable structure diagram), Equations (per nested
 * model; switched-off terms greyed; KaTeX, vendored), Symbols & parameters, Data, Fit &
 * accuracy, Digital twins vs truth (digital patients), Versions (as-of comparison), Interfaces.
 */
(function (root) {
  const GL = root.GL = root.GL || {};
  const t = (k, p) => GL.i18n.t(k, p), esc = (s) => GL.util.escapeHtml(s), U = GL.util;
  const TABS = ['model', 'equations', 'symbols', 'data', 'accuracy', 'twins', 'versions', 'interfaces'];
  const V = GL.views.research = { ui: { q2: 'Mab', clone: 'ddr', route: 'both', mf: false, node: 'hsc', q: '', twinPar: 'th_bs' } };
  const tex = (s, display) => { try { return katex.renderToString(s, { throwOnError: false, displayMode: !!display }); } catch (e) { return esc(s); } };
  const pct = (x) => (x === null || x === undefined || isNaN(x)) ? '–' : Math.round(100 * x) + '%';
  const dl = (d) => d === 0 ? 'D0' : (d > 0 ? 'D+' + d : 'D' + d);
  const card = () => GL.registry.get('model1').card;

  // ---- structure diagram ------------------------------------------------------------------------------
  const NODES = [
    { id: 'in_c', x: 20, y: 14, w: 200, h: 54, cls: 'dg-in', eq: ['exposure'], par: ['k_e', 'beta_renal', 'u_br_rgo', 'u_br_cy'] },
    { id: 'in_i', x: 260, y: 14, w: 210, h: 54, cls: 'dg-in', eq: ['inflam'], par: ['k_i', 'emax', 'ic50', 'hill'] },
    { id: 'in_z', x: 510, y: 14, w: 210, h: 54, cls: 'dg-in', eq: ['cart'], par: ['rho', 'k_t', 'k_ag', 'kappa', 'beta_a', 'car_t0', 'k_z', 'z50'] },
    { id: 'in_g', x: 760, y: 14, w: 180, h: 54, cls: 'dg-in', eq: ['gcsf'], par: ['k_g', 'gcsf_fold', 'gcsf_rel'] },
    { id: 'hsc', x: 20, y: 200, w: 180, h: 92, cls: 'dg-hsc', eq: ['hsc', 'output'], par: ['r_s', 'f0', 'th_as', 'th_bs', 'th_zs', 'eta_syn', 'r_a', 'r_b', 'delta'] },
    { id: 'mpp_n', x: 270, y: 130, w: 150, h: 50, cls: 'dg-mpp', eq: ['mpp', 'feedback'], par: ['k_p', 'sr_p', 'gmax', 'gamma', 'th_ap', 'emax', 'e_z'] },
    { id: 'mpp_mo', x: 270, y: 200, w: 150, h: 50, cls: 'dg-mpp', eq: ['mpp', 'feedback'], par: ['k_p', 'sr_p', 'gmax', 'gamma', 'th_ap'] },
    { id: 'mpp_pl', x: 270, y: 270, w: 150, h: 50, cls: 'dg-mpp', eq: ['mpp', 'feedback'], par: ['k_p', 'sr_p', 'gmax', 'gamma', 'th_ap'] },
    { id: 'mpp_l', x: 270, y: 340, w: 150, h: 50, cls: 'dg-mpp', eq: ['mpp', 'feedback'], par: ['k_p', 'sr_p', 'th_ap'] },
    { id: 'ch_n', x: 470, y: 130, w: 190, h: 50, cls: '', eq: ['chain', 'gcsf'], par: ['mtt_n', 'c_n', 'gcsf_fold', 'gcsf_rel'] },
    { id: 'ch_mo', x: 470, y: 200, w: 190, h: 50, cls: '', eq: ['chain'], par: ['mtt_mo', 'c_mo'] },
    { id: 'ch_pl', x: 470, y: 270, w: 190, h: 50, cls: '', eq: ['chain'], par: ['mtt_pl', 'c_pl', 'tpo_fold'] },
    { id: 'b_n', x: 700, y: 130, w: 110, h: 50, cls: 'dg-blood', eq: ['obs'], par: [] },
    { id: 'b_mo', x: 700, y: 200, w: 110, h: 50, cls: 'dg-blood', eq: ['obs'], par: [] },
    { id: 'b_pl', x: 700, y: 270, w: 110, h: 50, cls: 'dg-blood', eq: ['obs'], par: [] },
    { id: 'b_l', x: 700, y: 340, w: 110, h: 50, cls: 'dg-blood', eq: ['lymph', 'obs'], par: ['c_l', 'th_al'] },
    { id: 'lik', x: 840, y: 200, w: 110, h: 92, cls: 'dg-in', eq: ['obs', 'pheno'], par: [] }
  ];
  function diagram(sel) {
    const N = {}; NODES.forEach(n => { N[n.id] = n; });
    const c = (n, side) => side === 'r' ? [n.x + n.w, n.y + n.h / 2] : side === 'l' ? [n.x, n.y + n.h / 2] : side === 'b' ? [n.x + n.w / 2, n.y + n.h] : [n.x + n.w / 2, n.y];
    const edge = (a, b, cls) => { const [x1, y1] = c(N[a], 'r'), [x2, y2] = c(N[b], 'l'); const mx = (x1 + x2) / 2;
      return `<path class="${cls}" d="M${x1},${y1} C${mx},${y1} ${mx},${y2} ${x2},${y2}"/>`; };
    let s = `<defs><marker id="arr" viewBox="0 0 10 10" refX="9" refY="5" markerWidth="6" markerHeight="6" orient="auto-start-reverse"><path d="M0,0L10,5L0,10z" fill="var(--faint)"/></marker>
      <marker id="arrRed" viewBox="0 0 10 10" refX="9" refY="5" markerWidth="6" markerHeight="6" orient="auto-start-reverse"><path d="M0,0L10,5L0,10z" fill="var(--red)"/></marker>
      <marker id="arrTwin" viewBox="0 0 10 10" refX="9" refY="5" markerWidth="6" markerHeight="6" orient="auto-start-reverse"><path d="M0,0L10,5L0,10z" fill="var(--twin)"/></marker></defs>`;
    ['mpp_n', 'mpp_mo', 'mpp_pl', 'mpp_l'].forEach(m => { s += edge('hsc', m, 'dg-edge'); });
    s += `<text class="dg-lbl-twin" x="236" y="236" text-anchor="middle" style="fill:var(--muted)">F</text>`;
    [['mpp_n', 'ch_n'], ['mpp_mo', 'ch_mo'], ['mpp_pl', 'ch_pl'], ['ch_n', 'b_n'], ['ch_mo', 'b_mo'], ['ch_pl', 'b_pl'], ['mpp_l', 'b_l']].forEach(([a, b]) => { s += edge(a, b, 'dg-edge'); });
    ['b_n', 'b_mo', 'b_pl', 'b_l'].forEach(b => { s += edge(b, 'lik', 'dg-edge'); });
    // count feedback: blood -> own MPP component, in the gap below each row
    [['b_n', 'mpp_n'], ['b_mo', 'mpp_mo'], ['b_pl', 'mpp_pl']].forEach(([b, m]) => { const y = N[b].y + N[b].h + 9;
      s += `<path class="dg-fb" d="M${N[b].x + 30},${N[b].y + N[b].h} L${N[b].x + 30},${y} L${N[m].x + N[m].w - 30},${y} L${N[m].x + N[m].w - 30},${N[m].y + N[m].h + 1}"/>`; });
    s += `<text class="dg-lbl-twin" x="565" y="${N.b_n.y + N.b_n.h + 7}" text-anchor="middle" style="font-size:10.5px">${esc(t('dg.feedback'))}</text>`;
    // chemotherapy (a): stem cells and MPP
    s += `<path class="dg-hit" d="M80,68 L80,198"/><text class="dg-lbl" x="86" y="150">θaS·C</text>`;
    s += `<path class="dg-hit" d="M190,68 C190,110 250,120 268,145"/><text class="dg-lbl" x="200" y="112">θaP·C</text>`;
    // CAR T phase (b): inflammation I and CAR T signal Z join, hit the stem cells and suppress MPP self-renewal
    s += `<path class="dg-hit" style="marker-end:none" d="M365,68 L365,96"/><path class="dg-hit" style="marker-end:none" d="M615,68 L615,96 L365,96"/>`;
    s += `<path class="dg-hit" d="M365,96 L160,96 L160,198"/><text class="dg-lbl" x="236" y="90">b: θbS·I + θzS·Z</text>`;
    s += `<path class="dg-hit" d="M365,96 L365,128"/><text class="dg-lbl" x="372" y="118">E(I), e_Z</text>`;
    // G-CSF on the neutrophil chain
    s += `<path class="dg-hit" style="stroke:var(--blue)" d="M850,68 L850,100 L600,100 L600,128"/><text class="dg-lbl" style="fill:var(--blue)" x="720" y="94">u_G, u_R</text>`;
    s += `<text class="dg-lbl" x="755" y="${N.b_l.y + N.b_l.h + 16}" text-anchor="middle">θaL·C</text>`;
    NODES.forEach(n => {
      s += `<g class="dg-node ${n.cls} ${sel === n.id ? 'sel' : ''}" data-node="${n.id}" role="button" tabindex="0"><rect x="${n.x}" y="${n.y}" width="${n.w}" height="${n.h}" rx="9"/>
        <text x="${n.x + n.w / 2}" y="${n.y + n.h / 2 - 3}" text-anchor="middle">${esc(t('node.' + n.id))}</text>
        <text class="sub" x="${n.x + n.w / 2}" y="${n.y + n.h / 2 + 13}" text-anchor="middle">${esc(t('node.' + n.id + '.s'))}</text></g>`;
    });
    return `<svg viewBox="0 0 960 412" role="img" aria-label="${esc(t('dg.alt'))}">${s}</svg>`;
  }
  function paramValue(id) {
    const idx = GL.loader.index(), st = idx.population.struct;
    if (id in st) return st[id];
    const re = idx.params.re.find(r => r.id === id); return re ? (re.value || re.centre) : null;
  }
  function nodeDetail(id) {
    const n = NODES.find(x => x.id === id) || NODES[4];
    const spec = GL.model1.spec('MabxC_ddr'), eqs = card().equations(spec).filter(e => n.eq.includes(e.id));
    const sym = card().symbols;
    return `<div class="rgrid2"><div><h2>${esc(t('node.' + n.id))}</h2><p class="muted">${esc(t('node.' + n.id + '.d'))}</p>
      <div class="eqs">${eqs.map(e => `<div class="eq"><div class="eqt">${esc(t(e.titleKey))}</div><div class="eqx">${tex(e.latex, true)}</div></div>`).join('')}</div></div>
      <div>${n.par.length ? `<table class="tbl"><thead><tr><th>${esc(t('sym.col.sym'))}</th><th>${esc(t('sym.col.meaning'))}</th><th>${esc(t('sym.col.value'))}</th></tr></thead><tbody>
        ${n.par.map(k => `<tr><td>${sym[k] ? tex(sym[k][0]) : esc(k)}</td><td>${esc(t('sym.' + k))}</td><td>${U.fmt(paramValue(k), 3)} <span class="unit">${esc(sym[k] ? sym[k][1] : '')}</span></td></tr>`).join('')}</tbody></table>` : `<div class="muted small">${esc(t('rm.nopar'))}</div>`}</div></div>`;
  }

  // ---- tabs ------------------------------------------------------------------------------------------------
  const tabsHtml = {
    model() {
      return `<div class="rhead"><div><h1>${esc(t('rt.model'))}</h1><div class="muted">${esc(t('rm.lead'))}</div></div><span class="chip chip-twin">${esc(GL.loader.index().model_version)}</span></div>
        <div class="rcard diagram">${diagram(V.ui.node)}<div class="tiny muted">${esc(t('rm.click'))}</div></div><div class="rcard" id="nodeDetail">${nodeDetail(V.ui.node)}</div>
        <div class="rgrid3">${['rm.k1', 'rm.k2', 'rm.k3'].map(k => `<div class="rcard"><h2>${esc(t(k + '.t'))}</h2><div class="muted">${esc(t(k))}</div></div>`).join('')}</div>`;
    },
    equations() {
      const id = V.ui.mf ? 'MF' : V.ui.q2 + (V.ui.clone !== 'none' ? 'C_' + V.ui.clone : '') + (V.ui.route === 'hsc' ? '_Zh' : V.ui.route === 'mpp' ? '_Zm' : '');
      const spec = GL.model1.spec(id), eqs = card().equations(spec);
      const opt = (name, vals, cur) => `<select data-sel="${name}" ${V.ui.mf ? 'disabled' : ''}>${vals.map(v => `<option value="${v}" ${v === cur ? 'selected' : ''}>${esc(t('eqsel.' + name + '.' + v))}</option>`).join('')}</select>`;
      return `<div class="rhead"><div><h1>${esc(t('rt.equations'))}</h1><div class="muted">${esc(t('re.lead'))}</div></div><span class="chip chip-twin">${esc(t('re.model'))}: <b>${esc(id)}</b></span></div>
        <div class="rcard"><div class="selrow">${opt('q2', ['M0', 'Ma', 'Mb', 'Mab', 'Mabx'], V.ui.q2)}${opt('clone', ['none', 'ddr', 'dta'], V.ui.clone)}${opt('route', ['both', 'hsc', 'mpp'], V.ui.route)}
          <label class="sw-toggle"><input type="checkbox" id="eqMF" ${V.ui.mf ? 'checked' : ''}> ${esc(t('re.mf'))}</label><span class="tiny muted">${esc(t('re.grey'))}</span></div></div>
        <div class="eqs">${eqs.map(e => `<div class="eq"><div class="eqt">${esc(t(e.titleKey))}</div><div class="eqx">${tex(e.latex, true)}</div><div class="eqd">${esc(t(e.descKey))}</div></div>`).join('')}</div>`;
    },
    symbols() {
      const idx = GL.loader.index(), sym = card().symbols, src = card().source, q = V.ui.q.toLowerCase();
      const rows = [];
      idx.params.re.forEach(r => rows.push({ id: r.id, kind: 'individual', value: r.value || r.centre, extra: `ω ${r.omega} · ${esc(t('sym.box'))} [${r.box.join(', ')}]${r.covariates && r.covariates.length ? ' · ' + esc(t('sym.cov')) + ': ' + r.covariates.join(', ') : ''}`, src: 'prior' }));
      idx.params.struct.forEach(r => rows.push({ id: r.id, kind: r.kind, value: r.value, extra: r.kind === 'estimated' ? `${esc(t('sym.prior_sd'))} ${r.prior_sd} (${r.scale}) · ${esc(t('sym.box'))} [${r.box.join(', ')}]` : '', src: src[r.id] || 'fixed' }));
      const extra = [['sigma', '\\sigma_v', idx.population.sigma], ['tau', '\\tau_v', idx.population.tau]];
      const ph = idx.phenotype.term;
      const filt = rows.filter(r => !q || r.id.includes(q) || t('sym.' + r.id).toLowerCase().includes(q));
      return `<div class="rhead"><div><h1>${esc(t('rt.symbols'))}</h1><div class="muted">${esc(t('rs.lead'))}</div></div><input class="search" id="symQ" placeholder="${esc(t('rs.search'))}" value="${esc(V.ui.q)}"></div>
        <div class="rcard"><table class="tbl"><thead><tr><th>${esc(t('sym.col.sym'))}</th><th>${esc(t('sym.col.meaning'))}</th><th>${esc(t('sym.col.unit'))}</th><th>${esc(t('sym.col.value'))}</th><th>${esc(t('sym.col.kind'))}</th><th>${esc(t('sym.col.prior'))}</th><th>${esc(t('sym.col.source'))}</th></tr></thead><tbody>
          ${filt.map(r => `<tr><td>${sym[r.id] ? tex(sym[r.id][0]) : esc(r.id)}</td><td>${esc(t('sym.' + r.id))}<div class="tiny muted">${esc(r.id)}</div></td><td>${esc(sym[r.id] ? sym[r.id][1] : '')}</td><td><b>${U.fmt(r.value, 3)}</b></td><td><span class="kind kind-${r.kind}">${esc(t('kind.' + r.kind))}</span></td><td class="small">${r.extra}</td><td class="small">${esc(t('src.' + r.src))}</td></tr>`).join('')}
          ${!q ? extra.map(([k, s, o]) => `<tr><td>${tex(s)}</td><td>${esc(t('sym.' + k))}</td><td>log</td><td class="small">${Object.entries(o || {}).map(([v, x]) => `${v.toUpperCase()} ${U.fmt([].concat(x)[0], 3)}`).join(' · ')}</td><td><span class="kind kind-estimated">${esc(t('kind.estimated'))}</span></td><td></td><td class="small">${esc(t('src.mstep'))}</td></tr>`).join('') +
            `<tr><td>${tex('w,\\ s,\\ \\epsilon')}</td><td>${esc(t('sym.pheno'))}</td><td>–</td><td><b>${ph.w} · ${ph.s} · ${ph.eps}</b></td><td><span class="kind kind-fixed">${esc(t('kind.fixed'))}</span></td><td class="small">${esc(t('sym.pheno_sens'))}</td><td class="small">${esc(t('src.plan'))}</td></tr>` : ''}
        </tbody></table></div>`;
    },
    data() {
      // real cohort: research payloads only (pseudonyms, binned fields; small cells "<5"), opened at login
      const D = GL.session.cohort;
      if (!D) return `<div class="rhead"><div><h1>${esc(t('rt.data'))}</h1></div></div><div class="pending">${esc(t('rd.sealed'))}</div>`;
      const ph = D.phenotype, col = { quick: 'var(--ph-quick)', second_wave: 'var(--ph-sw)', delayed: 'var(--ph-del)' };
      const stat = (v, k) => `<div class="stat"><b>${v}</b><span class="muted small">${esc(t(k))}</span></div>`;
      const nObs = D.per_var.reduce((s, x) => s + x.n, 0);
      const xt = (k) => { const tb = D.crosstab[k], cols = Object.keys(tb);
        return `<table class="tbl"><thead><tr><th>${esc(t('rd.xt.' + k))}</th>${cols.map(c => `<th>${esc(t('rd.lvl.' + c))}</th>`).join('')}</tr></thead><tbody>
          ${['quick', 'second_wave', 'delayed'].map(r => `<tr><td>${esc(t('ph.' + r))}</td>${cols.map(c => `<td>${esc(tb[c][r])}</td>`).join('')}</tr>`).join('')}</tbody></table>`; };
      const recs = Object.values(GL.session.research).sort((x, y) => x.id.localeCompare(y.id));
      const nat = (a, k) => a && a.map[k] !== undefined ? Math.exp([].concat(a.map[k])[0]) : null;
      const prow = (r) => { const a = GL.model1.asofEntry(r, 30), s = r.summary;
        return `<tr><td><a href="#/research/versions" data-vp="${esc(r.id)}"><b>${esc(r.id)}</b></a></td><td>${esc(s.sex)}</td><td>${esc(s.age_band)}</td><td>${esc(s.lines_band)}</td><td>${s.asct ? esc(t('ui.yes')) : esc(t('ui.no'))}</td><td>${esc(t('rd.lvl.' + s.bridging_class))}</td><td>${s.crs === null || s.crs === undefined ? '–' : s.crs}</td>
          <td class="small">${['egfr', 'ferritin', 'ldh', 'crp'].map(k => `${k.toUpperCase()} ${esc(s.labs[k] || '–')}`).join(' · ')}</td><td>${esc(t('ph.' + (r.label && r.label.primary)))}</td><td>${pct(1 - a.p[0])}</td>
          <td class="small">${['r_s', 'th_as', 'th_bs', 'th_ap'].map(k => nat(a, k) === null ? '' : `${k} ${U.fmt(nat(a, k), 3)}`).filter(Boolean).join(' · ')}</td></tr>`; };
      return `<div class="rhead"><div><h1>${esc(t('rt.data'))}</h1><div class="muted">${esc(t('rd.lead'))}</div></div><span class="chip chip-warn">${esc(t('rd.local'))}</span></div>
        <div class="rgrid3" style="grid-template-columns:repeat(6,minmax(0,1fr))">${stat(D.n_patients, 'rd.patients')}${stat(nObs, 'rd.counts')}${stat(`${ph.quick}/${ph.second_wave}/${ph.delayed}`, 'rd.phenotypes')}${stat(D.gcsf_windows, 'rd.gcsf')}${stat(D.excluded.map(x => x.N).join(' + '), 'rd.excluded')}${stat(`${D.bridging_class.chemo}/${D.bridging_class.none}`, 'rd.bridging')}</div>
        <div class="rgrid2"><div class="rcard"><h2>${esc(t('rd.coverage'))} ${GL.help.btn('coverage')}</h2>${GL.charts.coverage({ rows: D.coverage.map(r => ({ id: r.id, color: col[r.label] || 'var(--faint)', days: r.days })) })}
          <div class="phlegend"><span><i class="sw ph-q"></i>${esc(t('ph.quick'))}</span><span><i class="sw ph-s"></i>${esc(t('ph.second_wave'))}</span><span><i class="sw ph-d"></i>${esc(t('ph.delayed'))}</span><span>${esc(t('rd.cov_note'))}</span></div></div>
          <div class="rcard"><h2>${esc(t('rd.crosstab'))} ${GL.help.btn('small_cells')}</h2>${xt('bridging_class')}<div style="height:8px"></div>${xt('lines_band')}
            <h2 style="margin-top:14px">${esc(t('rd.vars'))}</h2><table class="tbl"><thead><tr><th>${esc(t('rd.var'))}</th><th>${esc(t('rd.n'))}</th><th>${esc(t('rd.cens'))}</th><th>${esc(t('rd.pts'))}</th></tr></thead><tbody>
            ${D.per_var.map(x => `<tr><td><b>${x.var.toUpperCase()}</b></td><td>${x.n}</td><td>${x.censored}</td><td>${x.patients}</td></tr>`).join('')}</tbody></table></div></div>
        <div class="rcard"><h2>${esc(t('rd.patients_t'))} ${GL.help.btn('research_view')}</h2><div class="tiny muted">${esc(t('rd.patients_h'))}</div>
          <table class="tbl" id="researchTable"><thead><tr><th>${esc(t('rd.pseudo'))}</th><th>${esc(t('rd.sex'))}</th><th>${esc(t('sum.age'))}</th><th>${esc(t('sum.lines'))}</th><th>${esc(t('sum.asct'))}</th><th>${esc(t('sum.bridging'))}</th><th>${esc(t('sum.crs'))}</th><th>${esc(t('rd.labs'))}</th><th>${esc(t('rtw.recorded'))}</th><th>${esc(t('out.p_nonquick'))}</th><th>ψ̂ (D+30)</th></tr></thead><tbody>${recs.map(prow).join('')}</tbody></table></div>
        <div class="rcard"><h2>${esc(t('rd.rules'))}</h2><ol class="small" style="padding-left:18px;margin:0">${['blank', 'manual', 'gcsf', 'plt', 'follow', 'setpoint', 'phenotype', 'nottarget', 'research'].map(k => `<li>${esc(t('rd.rule.' + k))}</li>`).join('')}</ol></div>`;
    },
    accuracy() {
      const idx = GL.loader.index(), K = idx.k1 || {}, R = idx.results || {};
      const ut = K.unit_tests || [], pass = ut.filter(x => x.passed === true || x.passed === 'TRUE').length;
      const cal = (K.calibration || [])[0] || {};
      const pend = (k) => `<div class="pending"><span class="badge-pending">${esc(t('acc.pending'))}</span> ${esc(t(k))}</div>`;
      const fits = (R.fits || []).slice().sort((a, b) => a.BIC - b.BIC);
      return `<div class="rhead"><div><h1>${esc(t('rt.accuracy'))}</h1><div class="muted">${esc(t('ra.lead'))}</div></div></div>
        <div class="rgrid2">
          <div class="rcard"><h2>${esc(t('acc.k1'))} <span class="${pass === ut.length && ut.length ? 'ok-badge' : 'bad-badge'}">${pass}/${ut.length}</span></h2>
            <table class="tbl"><tbody>${ut.map(x => `<tr><td class="small">${esc(x.test)}</td><td>${x.passed === true || x.passed === 'TRUE' ? '<span class="ok-badge">✓</span>' : '<span class="bad-badge">✕</span>'}</td></tr>`).join('')}</tbody></table></div>
          <div class="rcard"><h2>${esc(t('acc.calib'))}</h2>
            <div class="kv"><span>${esc(t('acc.cal.thap'))}</span><b>${U.fmt(cal.th_ap_centre, 3)}</b></div><div class="kv"><span>${esc(t('acc.cal.thal'))}</span><b>${U.fmt(cal.th_al, 3)}</b></div>
            <div class="kv"><span>${esc(t('acc.cal.cl'))}</span><b>${U.fmt(cal.c_l, 3)}</b></div>
            <h2 style="margin-top:12px">${esc(t('acc.aim4'))}</h2><table class="tbl"><thead><tr><th>${esc(t('acc.engine'))}</th><th>${esc(t('rd.var'))}</th><th>${esc(t('acc.nadir'))}</th><th>${esc(t('acc.recovery'))}</th></tr></thead><tbody>
            ${(K.aim4_comparison || []).map(x => `<tr><td class="small">${esc(x.engine)}</td><td>${esc(x.lineage)}</td><td>${U.fmt(x.nadir_rel, 3)} · ${dl(x.nadir_day)}</td><td>${x.recovery_day === null ? '–' : dl(x.recovery_day)}</td></tr>`).join('')}</tbody></table></div>
        </div>
        <div class="rcard"><h2>${esc(t('acc.bic'))} ${GL.help.btn('bic')}</h2>${fits.length ? `<table class="tbl"><thead><tr><th>${esc(t('re.model'))}</th><th>BIC</th><th>ΔBIC</th><th>${esc(t('acc.npar'))}</th><th>${esc(t('acc.rounds'))}</th><th>${esc(t('acc.conv'))}</th></tr></thead><tbody>${fits.map(f => `<tr><td><b>${esc(f.id)}</b></td><td>${U.fmt(f.BIC, 1)}</td><td>${U.fmt(f.BIC - fits[0].BIC, 1)}</td><td>${f.n_par}</td><td>${f.rounds}</td><td>${f.converged ? '✓' : '✕'}</td></tr>`).join('')}</tbody></table>` : pend('acc.bic.p')}</div>
        <div class="rgrid2"><div class="rcard"><h2>${esc(t('acc.v1'))}</h2>${pend('acc.v1.p')}</div><div class="rcard"><h2>${esc(t('acc.vpc'))}</h2>${pend('acc.vpc.p')}</div></div>
        <div class="rcard"><h2>${esc(t('acc.syn'))}</h2><div class="muted small">${esc(t('acc.syn.lead'))}</div>
          <table class="tbl"><thead><tr><th>${esc(t('acc.cohort'))}</th><th>${esc(t('acc.gen'))}</th><th>${esc(t('ph.quick'))}</th><th>${esc(t('ph.second_wave'))}</th><th>${esc(t('ph.delayed'))}</th><th>${esc(t('acc.fit'))}</th></tr></thead><tbody>
          ${(R.cohorts || []).map(c => `<tr><td><b>${esc(c.cohort)}</b></td><td>${esc(t('acc.gen.' + c.generator))}</td><td>${c.quick}</td><td>${c.second_wave}</td><td>${c.delayed}</td><td>${fits.find(f => f.id === c.cohort) ? '✓' : `<span class="badge-pending">${esc(t('acc.running'))}</span>`}</td></tr>`).join('')}</tbody></table>
          ${pend('acc.syn.p')}</div>
        <div class="rcard"><h2>${esc(t('acc.resid'))}</h2>${pend('acc.resid.p')}</div>`;
    },
    twins() {
      return `<div class="rhead"><div><h1>${esc(t('rt.twins'))}</h1><div class="muted">${esc(t('rtw.lead'))}</div></div>
        <select id="twinPar">${['r_s', 'th_as', 'th_bs', 'th_ap', 'rho'].map(k => `<option value="${k}" ${k === V.ui.twinPar ? 'selected' : ''}>${esc(t('sym.' + k))}</option>`).join('')}</select></div>
        <div id="twinBody" class="muted">${esc(t('ui.loading'))}</div>`;
    },
    versions() {
      const idx = GL.loader.index();
      const ros = GL.session.roster('research');
      if (!V.ui.vp || !ros.some(r => r.id === V.ui.vp)) V.ui.vp = ros[0].id;
      const days = idx.asof_days;
      V.ui.v1 = V.ui.v1 === undefined ? 0 : V.ui.v1; V.ui.v2 = V.ui.v2 === undefined ? 14 : V.ui.v2;
      const sel = (id, cur) => `<select id="${id}">${days.map(d => `<option value="${d}" ${d === cur ? 'selected' : ''}>${dl(d)}</option>`).join('')}</select>`;
      return `<div class="rhead"><div><h1>${esc(t('rt.versions'))}</h1><div class="muted">${esc(t('rv.lead'))}</div></div>
        <div class="selrow"><select id="verP">${ros.map(r => `<option value="${r.id}" ${r.id === V.ui.vp ? 'selected' : ''}>${esc(r.id)}${r.type === 'digital' ? ' (DP)' : ''}</option>`).join('')}</select>${sel('verA', V.ui.v1)} → ${sel('verB', V.ui.v2)}</div></div>
        <div id="verBody" class="muted">${esc(t('ui.loading'))}</div>`;
    },
    interfaces() {
      const idx = GL.loader.index(), mods = GL.registry.modules;
      V.ui.cardId = V.ui.cardId || 'model1';
      const m = GL.registry.get(V.ui.cardId) || mods[0];
      const cardJson = JSON.stringify({ id: m.id, version: m.version, type: m.type, source: m.source,
        card: Object.assign({}, m.card, { equations: m.card && m.card.equations ? m.card.equations(GL.model1.spec('MabxC_ddr')).map(e => ({ id: e.id, latex: e.latex })) : undefined, validation: undefined }) }, null, 2);
      return `<div class="rhead"><div><h1>${esc(t('rt.interfaces'))}</h1><div class="muted">${esc(t('ri.lead'))}</div></div></div>
        <div class="rgrid2"><div class="rcard"><h2>${esc(t('ri.modules'))}</h2><table class="tbl"><thead><tr><th>ID</th><th>${esc(t('ri.type'))}</th><th>${esc(t('ri.version'))}</th><th>${esc(t('ri.outputs'))}</th><th></th></tr></thead><tbody>
          ${mods.map(x => `<tr><td><b>${esc(x.id)}</b><div class="tiny muted">${esc(t(x.nameKey))}</div></td><td>${esc(t('mtype.' + x.type))}</td><td class="small">${esc(x.version)}</td><td class="small">${((x.card && x.card.outputs) || []).map(o => esc(t(o.labelKey))).join(', ')}</td><td><button class="btn btn-xs" data-card="${esc(x.id)}">${esc(t('ri.card'))}</button></td></tr>`).join('')}</tbody></table>
          <h2 style="margin-top:14px">${esc(t('ri.remote'))}</h2><div class="small muted">${esc(t('ri.remote.d'))}</div>
          <div class="selrow" style="margin-top:6px"><input class="input" id="remUrl" placeholder="http://localhost:8000" style="flex:1"><button class="btn" id="remAdd">${esc(t('ri.add'))}</button></div>
          <h2 style="margin-top:14px">${esc(t('ri.bundle'))}</h2>
          <div class="kv"><span>schema</span><b>${esc(idx.schema)}</b></div><div class="kv"><span>${esc(t('ri.generated'))}</span><b>${esc(idx.generated)}</b></div>
          <div class="kv"><span>${esc(t('pop.label'))}</span><b>${esc(t('pop.' + idx.population.source))} · ${esc(idx.population.model)}</b></div><div class="kv"><span>${esc(t('ri.files'))}</span><b>${idx.roster.length}</b></div><div class="kv"><span>${esc(t('ri.sealed'))}</span><b>${idx.real_cohort ? idx.real_cohort.n : 0}</b></div>
          <h2 style="margin-top:14px">${esc(t('ri.docs'))}</h2><ul class="small"><li>docs/MODULE_API.md — ${esc(t('ri.doc.module'))}</li><li>docs/ML_API.md — ${esc(t('ri.doc.ml'))}</li><li>docs/DATA_API.md — ${esc(t('ri.doc.data'))}</li><li>api/plumber_stub.R, api/fastapi_stub.py — ${esc(t('ri.doc.stub'))}</li></ul></div>
          <div class="rcard"><h2>${esc(t('ri.card'))}: ${esc(m.id)}</h2><pre class="json">${esc(cardJson)}</pre></div></div>`;
    }
  };

  async function twinsFill() {
    const idx = GL.loader.index(), dps = idx.roster.filter(r => r.type === 'digital');
    const ps = await Promise.all(dps.map(r => GL.loader.patient(r.id)));
    const k = V.ui.twinPar;
    const nat = (k, x) => k === 'f0' ? 1 / (1 + Math.exp(-x)) : Math.exp(x);
    const rows = ps.map(p => { const a = p.asof[p.asof.length - 1], a30 = GL.model1.asofEntry(p, 30);
      return { p, est: a.map[k] !== undefined ? nat(k, [].concat(a.map[k])[0]) : null, truth: p.truth[k], p30: a30.p }; });
    const vals = rows.flatMap(r => [r.est, r.truth]).filter(x => x > 0);
    const dom = vals.length ? [Math.min(...vals) * 0.6, Math.max(...vals) * 1.6] : [0.01, 1];
    document.getElementById('twinBody').innerHTML = `<div class="rgrid2"><div class="rcard"><table class="tbl"><thead><tr><th>DP</th><th>${esc(t('rtw.mech'))}</th><th>${esc(t('rtw.latent'))}</th><th>${esc(t('rtw.recorded'))}</th><th>${esc(t('rtw.twin30'))}</th><th>${esc(t('rtw.truth'))}</th><th>${esc(t('rtw.est'))}</th></tr></thead><tbody>
      ${rows.map(r => `<tr><td><a href="#/clinic/${r.p.id}"><b>${r.p.id}</b></a> ${GL.app.typeMark(r.p)}</td><td class="small">${esc(t('mech.' + r.p.truth.mechanism))}</td><td>${esc(t('ph.' + r.p.truth.latent_phenotype))}</td><td>${esc(t('ph.' + r.p.label.primary))}</td><td style="min-width:90px"><div class="stack stack-sm"><i class="ph-q" style="width:${100 * r.p30[0]}%"></i><i class="ph-s" style="width:${100 * r.p30[1]}%"></i><i class="ph-d" style="width:${100 * r.p30[2]}%"></i></div></td><td>${U.fmt(r.truth, 3)}</td><td><b>${U.fmt(r.est, 3)}</b></td></tr>`).join('')}</tbody></table>
      <div class="tiny muted">${esc(t('rtw.note'))}</div></div>
      <div class="rcard"><h2>${esc(t('sym.' + k))}: ${esc(t('rtw.scatter'))}</h2>${GL.charts.scatter({ pts: rows.filter(r => r.est > 0 && r.truth > 0).map(r => ({ x: r.truth, y: r.est, label: r.p.id })), dom, log: true, xLabel: t('rtw.truth'), yLabel: t('rtw.est') })}
      <div class="tiny muted">${esc(t('rtw.scatter_note'))}</div></div></div>`;
  }
  async function versionsFill() {
    const p = await GL.loader.patient(V.ui.vp);
    if (p.kind === 'research') GL.audit.log('open', p.id, 'versions');
    const A = GL.model1.asofEntry(p, V.ui.v1), B = GL.model1.asofEntry(p, V.ui.v2);
    const sd = (a, i) => { let s = 0, c = 0; for (let j = 0; j < a.names.length; j++) for (let k = j; k < a.names.length; k++) { if (k === i) s += a.chol[c] * a.chol[c]; c++; } return Math.sqrt(s); };
    const nat = (k, x) => k.startsWith('b_') ? Math.exp(x) : k === 'f0' ? 1 / (1 + Math.exp(-x)) : Math.exp(x);
    const rowP = (k) => { const i = A.names.indexOf(k), xa = [].concat(A.map[k])[0], xb = [].concat(B.map[k])[0];
      return `<tr><td><b>${esc(k)}</b> <span class="tiny muted">${esc(k.startsWith('b_') ? t('rv.base') : t('sym.' + k))}</span></td><td>${U.fmt(nat(k, xa), 3)} <span class="tiny muted">±${U.fmt(sd(A, i), 2)}</span></td><td>${U.fmt(nat(k, xb), 3)} <span class="tiny muted">±${U.fmt(sd(B, B.names.indexOf(k)), 2)}</span></td><td>${U.fmt(nat(k, xb) / nat(k, xa), 2)}×</td></tr>`; };
    const pp = (a) => a.p.map(x => pct(x)).join(' / ');
    document.getElementById('verBody').innerHTML = `<div class="rgrid2"><div class="rcard"><h2>${esc(t('rv.params'))}</h2><table class="tbl"><thead><tr><th></th><th>${dl(A.L)}</th><th>${dl(B.L)}</th><th>${esc(t('rv.ratio'))}</th></tr></thead><tbody>${A.names.map(rowP).join('')}</tbody></table><div class="tiny muted">${esc(t('rv.note'))}</div></div>
      <div class="rcard"><h2>${esc(t('rv.outcomes'))}</h2><table class="tbl"><thead><tr><th></th><th>${dl(A.L)}</th><th>${dl(B.L)}</th></tr></thead><tbody>
        <tr><td>${esc(t('rv.p'))}</td><td>${pp(A)}</td><td>${pp(B)}</td></tr>
        <tr><td>${esc(t('est.recovery'))}</td><td>${A.rec[0] !== null && isFinite(A.rec[0]) ? dl(A.rec[0]) : '–'}</td><td>${B.rec[0] !== null && isFinite(B.rec[0]) ? dl(B.rec[0]) : '–'}</td></tr>
        ${A.top3 && B.top3 ? `<tr><td>${esc(t('rv.top'))}</td><td>${esc(t('plan.' + A.top3[0]))}</td><td>${esc(t('plan.' + B.top3[0]))}</td></tr>` : ''}
        <tr><td>${esc(t('rv.n'))}</td><td>${p.obs.anc.filter(o => o[0] <= A.L && o[4] === 0).length}</td><td>${p.obs.anc.filter(o => o[0] <= B.L && o[4] === 0).length}</td></tr></tbody></table></div></div>`;
  }

  V.render = async function (el, args) {
    const allowed = TABS.filter(k => GL.session.tabs().includes(k));
    const tab = allowed.includes(args[0]) ? args[0] : allowed[0];
    if (args[0] !== tab) { location.replace('#/research/' + tab); return; }
    V.ui.tab = tab;
    if (tab === 'data' && GL.session.cohort) GL.audit.log('open', 'cohort', 'data');
    el.innerHTML = `<div class="research"><nav class="rnav">${allowed.map(k => `<a href="#/research/${k}" class="${k === tab ? 'active' : ''}" data-tab="${k}">${esc(t('rt.' + k))}</a>`).join('')}<div class="foot">${esc(t('rnav.foot'))}</div></nav><section class="rmain" id="rmain">${tabsHtml[tab]()}</section></div>`;
    bind(el);
    if (tab === 'twins') await twinsFill();
    if (tab === 'versions') await versionsFill();
  };
  function bind(el) {
    el.querySelectorAll('[data-node]').forEach(n => { const go = () => { V.ui.node = n.dataset.node; el.querySelectorAll('.dg-node').forEach(x => x.classList.toggle('sel', x === n)); document.getElementById('nodeDetail').innerHTML = nodeDetail(V.ui.node); }; n.addEventListener('click', go); n.addEventListener('keydown', (e) => { if (e.key === 'Enter') go(); }); });
    el.querySelectorAll('[data-sel]').forEach(s => s.addEventListener('change', () => { V.ui[s.dataset.sel] = s.value; GL.app.render(); }));
    const mf = el.querySelector('#eqMF'); if (mf) mf.addEventListener('change', () => { V.ui.mf = mf.checked; GL.app.render(); });
    const q = el.querySelector('#symQ'); if (q) q.addEventListener('input', () => { V.ui.q = q.value; const pos = q.selectionStart; document.getElementById('rmain').innerHTML = tabsHtml.symbols(); bind(el); const q2 = el.querySelector('#symQ'); q2.focus(); q2.setSelectionRange(pos, pos); });
    const tp = el.querySelector('#twinPar'); if (tp) tp.addEventListener('change', () => { V.ui.twinPar = tp.value; twinsFill(); });
    ['verP', 'verA', 'verB'].forEach(id => { const s = el.querySelector('#' + id); if (s) s.addEventListener('change', () => { V.ui.vp = el.querySelector('#verP').value; V.ui.v1 = +el.querySelector('#verA').value; V.ui.v2 = +el.querySelector('#verB').value; versionsFill(); }); });
    el.querySelectorAll('[data-vp]').forEach(a => a.addEventListener('click', () => { V.ui.vp = a.dataset.vp; }));
    el.querySelectorAll('[data-card]').forEach(b => b.addEventListener('click', () => { V.ui.cardId = b.dataset.card; GL.app.render(); }));
    const ra = el.querySelector('#remAdd'); if (ra) ra.addEventListener('click', () => { const u = el.querySelector('#remUrl').value.trim(); if (u) { GL.remote.add(u); GL.app.render(); } });
  }
})(typeof window !== 'undefined' ? window : global);
