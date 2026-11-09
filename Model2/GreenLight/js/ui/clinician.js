/* GreenLight · ui/clinician.js
 * Clinician portal: one screen per patient (fits 1440x900 without scrolling). Left: patient list
 * with type marks (status colour + origin symbol) and filters. Main: patient summary and data
 * available; ANC/PLT/AMC/ALC with the digital twin's band (fixed axes); ferritin driver and
 * CAR T; how the twin's forecast evolved with the data; right: the twin's forecast as of the
 * selected day, why, top-3 strategies, other predictors. Bottom: pinned day scrubber.
 * What-if drawer: browser simulation of a custom plan (JS port of Model2).
 */
(function (root) {
  const GL = root.GL = root.GL || {};
  const t = (k, p) => GL.i18n.t(k, p), esc = (s) => GL.util.escapeHtml(s), U = GL.util;
  const LIN = [['anc', 'var(--l-anc)', [{ v: 0.5, label: '0.5' }, { v: 1, label: '1.0' }, { v: 1.5, label: '1.5' }], 0.05],
    ['plt', 'var(--l-plt)', [{ v: 20, label: '20' }, { v: 50, label: '50' }], 3],
    ['amc', 'var(--l-amc)', [{ v: 0.1, label: '0.1' }], 0.005],
    ['alc', 'var(--l-alc)', [{ v: 0.2, label: '0.2' }], 0.005]];
  const BATCHES = [-3, 0, 3, 5, 7, 10, 12, 14, 17, 21, 24, 28, 30, 37, 44, 51, 58, 65, 72, 79, 86];
  const pct = (x) => (x === null || x === undefined || isNaN(x)) ? '–' : Math.round(100 * x) + '%';
  const dl = (d) => { d = Math.round(d); return d === 0 ? 'D0' : (d > 0 ? 'D+' + d : 'D' + d); };
  const V = GL.views.clinic = { ui: { filter: GL.store.get('rosterFilter', 'all'), q: '' } };

  function rosterHtml(idx, pid) {
    const f = V.ui.filter, q = V.ui.q.toLowerCase(), all = GL.session.roster('clinic');
    const list = all.filter(p => (f === 'all' || (f === 'active' && p.status === 'active') || (f === 'completed' && p.status === 'completed') || (f === 'digital' && p.type === 'digital')) && (!q || p.id.toLowerCase().includes(q)))
      .sort((a, b) => (a.status === b.status ? 0 : a.status === 'active' ? -1 : 1) || a.id.localeCompare(b.id));
    const items = list.map(p => {
      const last = p.evolution.filter(e => p.status !== 'active' || e[0] <= (GL.store.patient(p).horizon))
        .slice(-1)[0] || [0, 0, 0, 0];
      const sub = p.status === 'active' ? `${t('roster.data_to')} ${dl(GL.store.patient(p).horizon)}` : (p.label ? t('ph.' + p.label) : '');
      return `<div class="ritem ${p.id === pid ? 'active' : ''}" data-pid="${esc(p.id)}" role="button" tabindex="0">
        <span>${GL.app.typeMark(p)}</span>
        <div><div class="rid">${esc(p.id)}</div><div class="rsub">${esc(sub)}</div></div>
        <div class="rph" title="${esc(t('roster.twin_now'))}"><i class="ph-q" style="width:${100 * last[1]}%"></i><i class="ph-s" style="width:${100 * last[2]}%"></i><i class="ph-d" style="width:${100 * last[3]}%"></i></div>
      </div>`;
    }).join('');
    const tab = (k) => `<button data-f="${k}" class="${f === k ? 'active' : ''}">${esc(t('filter.' + k))} <span class="muted">${all.filter(p => k === 'all' || (k === 'active' && p.status === 'active') || (k === 'completed' && p.status === 'completed') || (k === 'digital' && p.type === 'digital')).length}</span></button>`;
    return `<div class="roster-head">
        <div class="legend"><span>${GL.app.typeMark({ type: 'real', status: 'active' })} ${esc(t('type.active'))}</span><span>${GL.app.typeMark({ type: 'real', status: 'completed' })} ${esc(t('type.completed'))}</span>
          <span><svg class="pmark" viewBox="0 0 12 12"><circle cx="6" cy="6" r="4.5" fill="var(--slate)"/></svg> ${esc(t('type.real'))}</span><span><svg class="pmark" viewBox="0 0 12 12"><path d="M6 1l5 5-5 5-5-5z" fill="none" stroke="var(--slate)" stroke-width="2"/></svg><span class="dp-tag">DP</span> ${esc(t('type.digital'))}</span></div>
        <div class="ftabs">${['all', 'active', 'completed', 'digital'].map(tab).join('')}</div>
        <input class="search" id="rosterQ" placeholder="${esc(t('roster.search'))}" value="${esc(V.ui.q)}">
      </div>
      <div class="roster-list" id="rosterList">${items || `<div class="muted small" style="padding:10px">${esc(t('roster.none'))}</div>`}</div>
      <div style="padding:6px 10px;border-top:1px solid var(--border);display:flex;justify-content:space-between;align-items:center"><span class="tiny muted">${esc(t('roster.real_note'))}</span><button class="btn btn-xs" id="btnReset" title="${esc(t('roster.reset_h'))}">${esc(t('roster.reset'))}</button></div>`;
  }

  function ydom(p) {
    if (p._ydom) return p._ydom;
    const out = {};
    LIN.forEach(([v, , , floor]) => {
      let lo = Infinity, hi = 0;
      p.asof.forEach(a => { a.band[v][0].forEach(x => { if (x > 0) lo = Math.min(lo, x); }); a.band[v][2].forEach(x => { hi = Math.max(hi, x); }); });
      (p.obs[v] || []).forEach(o => { if (o[0] >= -10 && o[0] <= 90) { lo = Math.min(lo, o[1]); hi = Math.max(hi, o[1]); } });
      out[v] = [Math.max(floor, lo * 0.8), Math.max(hi * 1.25, floor * 10)];
    });
    return (p._ydom = out);
  }
  const facts = (p) => {
    const s = p.summary, c = (k, v) => `<span class="chip">${esc(t(k))}: <b>${esc(v)}</b></span>`;
    return [c('sum.age', `${s.age} · ${s.female ? t('sum.f') : t('sum.m')}`), c('sum.lines', s.prior_lines), c('sum.asct', s.prior_asct ? t('ui.yes') : t('ui.no')),
      c('sum.bridging', t('br.' + s.bridging)), c('sum.ferritin', s.ferritin_ld), c('sum.ldh', s.ldh_ld), c('sum.egfr', Math.round(s.egfr)),
      c('sum.crs', s.max_crs === null || s.max_crs === undefined ? '–' : s.max_crs), s.hematotox !== undefined && s.hematotox !== null ? c('sum.hx', s.hematotox) : ''].join('');
  };
  function dataFacts(p, asOf) {
    const n = {}; let tot = 0;
    for (const v in p.obs) { n[v] = p.obs[v].filter(o => o[0] <= asOf && o[4] === 0 && o[3] === 1).length; tot += n[v]; }
    const anc = p.obs.anc.filter(o => o[0] <= asOf && o[3] === 1).slice(-1)[0];
    const gw = (p.design.gcsf || []).filter(w => w[0] < asOf).length;
    return `<div><b>${tot}</b> ${esc(t('data.counts_to', { day: dl(asOf) }))} ${GL.help.btn('data_points')}</div>
      <div>ANC ${n.anc} · PLT ${n.plt} · AMC ${n.amc} · ALC ${n.alc}</div>
      <div>${anc ? esc(t('data.last_anc', { v: anc[1], day: dl(anc[0]) })) : ''}${gw ? ' · ' + esc(t('data.gcsf', { n: gw })) : ''}</div>`;
  }

  function stackBar(p3, cls) { return `<div class="stack ${cls || ''}"><i class="ph-q" style="width:${100 * p3[0]}%"></i><i class="ph-s" style="width:${100 * p3[1]}%"></i><i class="ph-d" style="width:${100 * p3[2]}%"></i></div>`; }
  const light = (p3) => (p3[2] + 0.5 * p3[1]) >= 0.4 ? 'red' : (1 - p3[0]) >= 0.3 ? 'amber' : 'green';

  function panelHtml(p, a, prev, asOf) {
    const p3 = a.p, nq = 1 - p3[0], dnq = prev ? nq - (1 - prev.p[0]) : null;
    const rec = a.rec, recTxt = isFinite(rec[0]) && rec[0] !== null && rec[0] <= 90
      ? t('est.rec_by', { day: dl(rec[0]), lo: dl(rec[1]), hi: isFinite(rec[2]) && rec[2] !== null ? dl(rec[2]) : '>' + dl(90) })
      : t('est.rec_none');
    const recorded = p.status === 'completed' && asOf >= 30 && p.label && p.label.primary ? `<div class="kv"><span>${esc(t('est.recorded'))}</span><b>${esc(t('ph.' + p.label.primary))}</b></div>` : '';
    const when = a.when || {}, who = a.who || {}, fx = a.f || {};
    const wrow = (k, v) => `<div class="bar-row"><span>${esc(t('when.' + k))}</span><span class="track"><span class="fill" style="width:${100 * (v || 0)}%;display:block"></span></span><span class="v">${pct(v)}</span></div>`;
    const hrow = (k, v) => `<div class="bar-row"><span>${esc(t('who.' + k))}</span><span class="track"><span class="fill" style="width:${100 * (v || 0)}%;display:block;background:var(--dp)"></span></span><span class="v">${pct(v)}</span></div>`;
    const plans = a.plans.slice(0, 3).map((pl, i) => `<div class="plan-row ${i === 0 ? 'best' : ''}" data-plan="${pl.id}">
        <span class="rk">${i + 1}</span><div><div class="pname">${esc(t('plan.' + pl.id))} ${GL.help.btn('plan.' + pl.id)}</div>${stackBar(pl.p, 'stack-sm')}</div>
        <span class="small">${esc(t('plan.nonquick'))} <b>${pct(1 - pl.p[0])}</b><br>PLT&lt;20 <b>${pct(pl.plt20)}</b></span><span class="sc">${pl.score === null ? '–' : U.fmt(pl.score, 2)}</span></div>`).join('');
    const mods = GL.registry.predictAll(p, asOf);
    const hx = mods.car_hematotox, lg = mods.counts_logistic;
    return `<div class="card">
        <div class="card-head"><h3><span class="light ${light(p3)}"></span>${esc(t('est.title', { day: dl(asOf) }))} ${GL.help.btn('twin_forecast')}</h3><span class="chip chip-twin">${esc(t('est.twin'))}</span></div>
        ${stackBar(p3)}
        <div class="phlegend"><span><i class="sw ph-q"></i>${esc(t('ph.quick'))} <b>${pct(p3[0])}</b></span><span><i class="sw ph-s"></i>${esc(t('ph.second_wave'))} <b>${pct(p3[1])}</b></span><span><i class="sw ph-d"></i>${esc(t('ph.delayed'))} <b>${pct(p3[2])}</b></span></div>
        <div class="kv" style="margin-top:4px"><span>${esc(t('est.nonquick'))} ${GL.help.btn('nonquick')}</span><span><b>${pct(nq)}</b> ${dnq !== null && Math.abs(dnq) >= 0.005 ? `<span class="delta ${dnq > 0 ? 'up' : 'down'}">${dnq > 0 ? '▲' : '▼'} ${Math.round(Math.abs(100 * dnq))} ${esc(t('est.pts'))}</span>` : ''}</span></div>
        <div class="kv"><span>${esc(t('est.recovery'))} ${GL.help.btn('rec_day')}</span><b>${esc(recTxt)}</b></div>
        <div class="kv"><span>${esc(t('est.nadir'))}</span><b>${a.nadir && a.nadir[1] !== null ? `${U.fmt(a.nadir[1], 2)} · ${dl(a.nadir[0])}` : '–'}</b></div>
        ${recorded}
      </div>
      <div class="card">
        <div class="card-head"><h3>${esc(t('why.title'))} ${GL.help.btn('why')}</h3><span class="sub">${esc(t('why.sub'))}</span></div>
        <div class="why2"><div><h4>${esc(t('why.when'))}</h4><div class="bars">${['a', 'b', 'ab_joint', 'ab_redundant', 'other'].map(k => wrow(k, when[k])).join('')}</div></div>
          <div><h4>${esc(t('why.who'))}</h4>${nq >= 0.1 ? `<div class="bars">${['reserve', 'clone', 'both'].map(k => hrow(k, who[k])).join('')}</div>` : `<div class="tiny muted">${esc(t('why.who_quick'))}</div>`}</div></div>
        <div class="tiny muted" style="margin-top:4px">${esc(t('why.loss', { a: pct(fx.a), bi: pct(fx.bi), bz: pct(fx.bz) }))}</div>
      </div>
      <div class="card">
        <div class="card-head"><h3>${esc(t('plans.title'))} ${GL.help.btn('plans')}</h3><span class="sub">${esc(t('plans.illustrative'))}</span><button class="btn btn-xs" id="btnWhatIf">${GL.app.icon('sliders')} ${esc(t('plans.whatif'))}</button></div>
        <div class="plans">${plans}</div>
      </div>
      <div class="card">
        <div class="card-head"><h3>${esc(t('mods.title'))} ${GL.help.btn('modules')}</h3></div>
        <div class="mods">
          <div class="mod"><span class="tp">${esc(t('mtype.mechanistic'))}</span>${esc(t('mod.model2.short'))}<b>${pct(nq)}</b><span class="tiny muted">${esc(t('out.p_nonquick'))}</span></div>
          <div class="mod"><span class="tp">${esc(t('mtype.score'))}</span>CAR-HEMATOTOX<b>${hx ? hx.score : '–'}</b><span class="tiny muted">${hx ? esc(t('hx.' + hx.group)) : esc(t('ui.na'))}</span></div>
          <div class="mod"><span class="tp">${esc(t('mtype.ml'))}</span>${esc(t('mod.logistic.short'))}<b>${lg ? pct(lg.p_nonquick) : '–'}</b><span class="tiny muted">${lg ? esc(t('mods.landmark', { day: dl(lg.landmark) })) : esc(t('mods.from_d0'))}</span></div>
        </div>
      </div>`;
  }

  function drawerHtml(p, a) {
    const L = a.L, dec = Math.max(Math.floor(L) + 1, 1);
    const w = V.ui.wi || { gcsf: false, gstart: Math.max(dec, 15), gdays: 5, reactive: false, anti: false, afrom: Math.max(dec, 3), afac: 0.5, tpo: false, tstart: Math.max(dec, 7), cd34: false, cday: Math.max(dec, 7), cadd: 0.3 };
    V.ui.wi = w;
    const rng = (id, v, lo, hi, st, lbl, unit) => `<div class="field"><label for="${id}">${esc(t(lbl))}</label><span class="fv" id="${id}_v">${v}${unit || ''}</span><input type="range" id="${id}" min="${lo}" max="${hi}" step="${st}" value="${v}"></div>`;
    const tog = (id, v, lbl) => `<label class="sw-toggle"><input type="checkbox" id="${id}" ${v ? 'checked' : ''}> ${esc(t(lbl))} ${GL.help.btn('wi.' + id)}</label>`;
    return `<div class="drawer-head"><h3>${GL.app.icon('sliders')} ${esc(t('wi.title', { day: dl(L) }))}</h3><button class="btn btn-xs" id="wiClose">✕</button></div>
      <div class="drawer-body">
        <div class="small muted">${esc(t('wi.intro'))}</div>
        <div class="selrow">${a.plans.slice(0, 7).map(pl => `<button class="btn btn-xs" data-preset="${pl.id}">${esc(t('plan.' + pl.id))}</button>`).join('')}</div>
        <div class="card">${tog('wGcsf', w.gcsf, 'wi.gcsf')}${rng('wGstart', w.gstart, dec, 40, 1, 'wi.start', '')}${rng('wGdays', w.gdays, 1, 10, 1, 'wi.days', ' d')}</div>
        <div class="card">${tog('wReactive', w.reactive, 'wi.reactive')}</div>
        <div class="card">${tog('wAnti', w.anti, 'wi.anti')}${rng('wAfrom', w.afrom, Math.max(0, dec), 30, 1, 'wi.from', '')}${rng('wAfac', w.afac, 0.1, 1, 0.05, 'wi.factor', '×')}</div>
        <div class="card">${tog('wTpo', w.tpo, 'wi.tpo')}${rng('wTstart', w.tstart, dec, 40, 1, 'wi.start', '')}</div>
        <div class="card">${tog('wCd34', w.cd34, 'wi.cd34')}${rng('wCday', w.cday, dec, 40, 1, 'wi.day', '')}${rng('wCadd', w.cadd, 0.05, 0.6, 0.05, 'wi.size', '')}</div>
        <div class="selrow"><button class="btn btn-primary" id="wiRun">${esc(t('wi.run'))}</button><button class="btn" id="wiSave">${esc(t('wi.save'))}</button><button class="btn btn-ghost" id="wiClear">${esc(t('wi.clear'))}</button></div>
        <div id="wiOut"></div>
        <div id="wiSaved"></div>
      </div>`;
  }
  function planFromUi() {
    const w = V.ui.wi;
    return { gcsf: w.gcsf ? { start: w.gstart, days: w.gdays } : null, reactive: w.reactive, anti: w.anti ? { from: w.afrom, factor: w.afac } : null,
      tpo: w.tpo ? { start: w.tstart, end: 60 } : null, cd34: w.cd34 ? { day: w.cday, add: w.cadd } : null };
  }
  function presetToUi(id, dec) {
    const w = V.ui.wi;
    Object.assign(w, { gcsf: false, reactive: false, anti: false, tpo: false, cd34: false });
    if (id === 'gcsf_early') Object.assign(w, { gcsf: true, gstart: dec, gdays: 5 });
    if (id === 'gcsf_reactive') w.reactive = true;
    if (id === 'anti_inflam') Object.assign(w, { anti: true, afrom: Math.max(0, dec - 1), afac: 0.5 });
    if (id === 'tpo_ra') Object.assign(w, { tpo: true, tstart: Math.max(dec, 7) });
    if (id === 'cd34_boost') Object.assign(w, { cd34: true, cday: Math.max(dec, 7), cadd: 0.3 });
    if (id === 'combo') Object.assign(w, { anti: true, afrom: Math.max(0, dec - 1), afac: 0.5, reactive: true });
  }

  V.drawCharts = function () {
    const p = V.ui.p, a = V.ui.a, asOf = V.ui.asOf, idx = GL.loader.index(), dom = ydom(p), hz = V.ui.horizon;
    LIN.forEach(([v, col, thr]) => {
      const box = document.querySelector(`#ch_${v} .plot`); if (!box) return;
      const wr = V.ui.whatif && V.ui.whatif.L === a.L ? V.ui.whatif.res[v] : null;
      const wi = wr ? { days: wr.days.filter(d => d >= -10), v: wr.v.filter((_, j) => wr.days[j] >= -10) } : null;
      box.innerHTML = GL.charts.lineage({ w: box.clientWidth || 300, h: box.clientHeight || 160, yDomain: dom[v], log: true, grid: idx.grid_days,
        band: { lo: a.band[v][0], med: a.band[v][1], hi: a.band[v][2] }, obs: p.obs[v], asOf, horizon: hz, color: col, thr, whatif: wi,
        gcsf: v === 'anc' ? (p.design.gcsf || []).filter(w => w[0] <= asOf) : null, future: p.status === 'completed' ? 'hollow' : 'hide' });
    });
    const fb = document.querySelector('#ch_fer .plot');
    if (fb) {
      const ymax = Math.max(1, ...p.ferritin.map(f => Math.log(f[1] / p.summary.ferritin_ld)), ...p.asof.map(x => Math.max(...x.drv.v))) * 1.1;
      fb.innerHTML = GL.charts.ferritin({ w: fb.clientWidth || 300, h: fb.clientHeight || 50, ferritin: p.ferritin, base: p.summary.ferritin_ld, drv: a.drv, asOf, horizon: hz, yMax: ymax });
    }
    const eb = document.querySelector('#ch_evo .plot');
    if (eb) eb.innerHTML = GL.charts.evolution({ w: eb.clientWidth || 300, h: eb.clientHeight || 80, evolution: p.asof.map(x => [x.L].concat(x.p)), asOf, horizon: hz,
      label: p.status === 'completed' && p.label ? t('evo.recorded', { ph: t('ph.' + p.label.primary) }) : '' });
  };

  let scrubLog = 0;
  V.update = function (asOf) {
    const p = V.ui.p;
    if (!p) return;
    if (V.ui.real && V.ui.asOf !== undefined && asOf !== V.ui.asOf && Date.now() - scrubLog > 5000) { scrubLog = Date.now(); GL.audit.log('scrub', p.id, asOf); }
    V.ui.asOf = asOf;
    const i = p.asof.findIndex(x => x.L === GL.model2.asofEntry(p, asOf).L);
    V.ui.a = p.asof[i]; V.ui.prev = i > 0 ? p.asof[i - 1] : null;
    const st = GL.store.patient(p); st.asOf = asOf; GL.store.savePatient(p.id, st);
    V.drawCharts();
    document.getElementById('cright').innerHTML = panelHtml(p, V.ui.a, V.ui.prev, asOf);
    document.getElementById('dataFacts').innerHTML = dataFacts(p, asOf);
    bindPanel();
  };

  function bindPanel() {
    const b = document.getElementById('btnWhatIf'); if (b) b.addEventListener('click', openDrawer);
  }
  function openDrawer() {
    const dr = document.getElementById('drawer'), a = V.ui.a, p = V.ui.p;
    dr.innerHTML = drawerHtml(p, a); dr.classList.add('open');
    const dec = Math.max(Math.floor(a.L) + 1, 1);
    const sync = () => {
      const w = V.ui.wi, g = (id) => document.getElementById(id);
      w.gcsf = g('wGcsf').checked; w.reactive = g('wReactive').checked; w.anti = g('wAnti').checked; w.tpo = g('wTpo').checked; w.cd34 = g('wCd34').checked;
      w.gstart = +g('wGstart').value; w.gdays = +g('wGdays').value; w.afrom = +g('wAfrom').value; w.afac = +g('wAfac').value; w.tstart = +g('wTstart').value; w.cday = +g('wCday').value; w.cadd = +g('wCadd').value;
      [['wGstart', ''], ['wGdays', ' d'], ['wAfrom', ''], ['wAfac', '×'], ['wTstart', ''], ['wCday', ''], ['wCadd', '']].forEach(([id, u]) => { g(id + '_v').textContent = g(id).value + u; });
    };
    dr.querySelectorAll('input').forEach(i => i.addEventListener('input', sync));
    dr.querySelector('#wiClose').addEventListener('click', () => dr.classList.remove('open'));
    dr.querySelectorAll('[data-preset]').forEach(b => b.addEventListener('click', () => { presetToUi(b.dataset.preset, dec); openDrawer(); }));
    dr.querySelector('#wiRun').addEventListener('click', runWhatIf);
    dr.querySelector('#wiClear').addEventListener('click', () => { V.ui.whatif = null; document.getElementById('wiOut').innerHTML = ''; V.drawCharts(); });
    dr.querySelector('#wiSave').addEventListener('click', async () => {
      if (!V.ui.whatif) await runWhatIf();
      await GL.store.saveScenario({ patient: p.id, L: a.L, plan: planFromUi(), p: V.ui.whatif.res.p, base: V.ui.whatif.base.p, saved: new Date().toISOString() });
      listSaved();
    });
    listSaved();
  }
  async function listSaved() {
    const box = document.getElementById('wiSaved'); if (!box) return;
    const sc = await GL.store.scenarios(V.ui.p.id);
    box.innerHTML = sc.length ? `<h3>${esc(t('wi.saved'))}</h3>` + sc.map(s => `<div class="kv small"><span>${dl(s.L)} · ${esc(Object.keys(s.plan).filter(k => s.plan[k]).map(k => t('wi.k.' + k)).join(' + ') || t('plan.standard'))}</span><b>${pct(1 - s.p.quick)}</b></div>`).join('') : '';
  }
  async function runWhatIf() {
    const out = document.getElementById('wiOut'), p = V.ui.p, a = V.ui.a;
    if (V.ui.real && !GL.session.can('whatif.real')) return;
    if (V.ui.real) GL.audit.log('what-if', p.id, a.L);
    out.innerHTML = `<div class="muted small">${esc(t('wi.running'))}</div>`;
    await new Promise(r => setTimeout(r, 20));
    const t0 = performance.now();
    const base = GL.model2.whatIf(p, a, null, { n: 24 }), res = GL.model2.whatIf(p, a, planFromUi(), { n: 24 });
    const ms = Math.round(performance.now() - t0);
    res.amc = null; res.alc = null;
    V.ui.whatif = { L: a.L, res, base };
    const row = (lbl, r) => `<div class="kv"><span>${esc(lbl)}</span><span><b>${pct(1 - r.p.quick)}</b> ${esc(t('plan.nonquick'))} · PLT&lt;20 ${pct(r.plt20)}</span></div>${stackBar([r.p.quick, r.p.second_wave, r.p.delayed], 'stack-sm')}`;
    out.innerHTML = `<div class="card">${row(t('wi.base'), base)}<div style="height:6px"></div>${row(t('wi.plan'), res)}<div class="tiny muted" style="margin-top:6px">${esc(t('wi.note', { n: res.n, ms }))}</div></div>`;
    V.drawCharts();
  }

  V.render = async function (el, args) {
    GL.ui.scrubber.stopAll();
    const idx = GL.loader.index();
    const roster = GL.session.roster('clinic');
    const last = GL.store.lastPid();
    const pid = args[0] && roster.some(r => r.id === args[0]) ? args[0] : (roster.some(r => r.id === last) ? last : roster[0].id);
    if (args[0] !== pid) { location.replace('#/clinic/' + encodeURIComponent(pid)); return; }
    GL.store.lastPid(pid);
    el.innerHTML = `<div class="clinic"><aside class="roster" id="roster">${rosterHtml(idx, pid)}</aside><section class="cmain"><div class="boot">${esc(t('ui.loading'))}</div></section></div><div class="drawer" id="drawer"></div>`;
    bindRoster(el, idx, pid);
    let p;
    try { p = await GL.loader.patient(pid); } catch (e) { el.querySelector('.cmain').innerHTML = `<div class="boot">${esc(t('ui.load_error', { id: pid }))}</div>`; return; }
    const meta = roster.find(r => r.id === pid) || p;
    V.ui.real = !GL.loader.isDigital(pid);
    if (V.ui.real) GL.audit.log('open', pid);
    const st = GL.store.patient(p);
    const horizon = p.status === 'active' ? st.horizon : 90;
    const days = idx.asof_days;
    const selectable = days.filter(d => d <= horizon);
    let asOf = st.asOf !== undefined && st.asOf !== null && st.asOf <= horizon ? st.asOf : (p.status === 'active' ? selectable[selectable.length - 1] : 30);
    if (!selectable.includes(asOf)) asOf = selectable.reduce((b, x) => Math.abs(x - asOf) < Math.abs(b - asOf) ? x : b, selectable[0]);
    V.ui.p = p; V.ui.horizon = horizon; V.ui.whatif = null; V.ui.wi = null; V.ui.asOf = undefined;
    const chart = (v) => `<div class="card chart-card" id="ch_${v}"><div class="card-head"><h3>${esc(t('lin.' + v))} <span class="unit">${esc(t('lin.' + v + '.u'))}</span> ${GL.help.btn('lin.' + v)}</h3><span class="sub">${esc(t('chart.band'))}</span></div><div class="plot"></div></div>`;
    el.querySelector('.cmain').innerHTML = `
      <div class="psum"><div class="pid">${esc(p.id)} ${GL.app.typeMark(meta, true)}</div><div class="facts">${facts(p)}</div><div class="datafacts" id="dataFacts"></div></div>
      <div class="cgrid">
        <div class="cleft">
          <div class="charts4">${LIN.map(l => chart(l[0])).join('')}</div>
          <div class="card chart-card" id="ch_fer"><div class="card-head"><h3>${esc(t('fer.title'))} ${GL.help.btn('ferritin')}</h3><span class="sub">${esc(t('fer.sub'))}</span></div><div class="plot"></div></div>
          <div class="card chart-card" id="ch_evo"><div class="card-head"><h3>${esc(t('evo.title'))} ${GL.help.btn('evolution')}</h3><span class="sub">${esc(t('evo.sub'))}</span></div><div class="plot"></div></div>
        </div>
        <div class="cright" id="cright"></div>
      </div>
      <div class="scrub" id="scrub"></div>`;
    const next = BATCHES.find(b => b > horizon);
    V.ui.scrub = GL.ui.scrubber(document.getElementById('scrub'), {
      days, horizon, asOf, active: p.status === 'active', gcsf: p.design.gcsf || [], canIngest: next !== undefined,
      onChange: (d) => V.update(d),
      onPlay: () => { if (V.ui.real) GL.audit.log('play', p.id); },
      onIngest: () => { const s2 = GL.store.patient(p); s2.horizon = next; s2.asOf = next; GL.store.savePatient(p.id, s2); V.render(el, args); }
    });
    V.update(asOf);
  };
  function bindRoster(el, idx, pid) {
    el.querySelectorAll('.ritem').forEach(r => { const go = () => { location.hash = '#/clinic/' + encodeURIComponent(r.dataset.pid); }; r.addEventListener('click', go); r.addEventListener('keydown', (e) => { if (e.key === 'Enter') go(); }); });
    el.querySelectorAll('.ftabs [data-f]').forEach(b => b.addEventListener('click', () => { V.ui.filter = b.dataset.f; GL.store.set('rosterFilter', V.ui.filter); document.getElementById('roster').innerHTML = rosterHtml(idx, pid); bindRoster(el, idx, pid); }));
    const q = el.querySelector('#rosterQ'); if (q) q.addEventListener('input', () => { V.ui.q = q.value; const pos = q.selectionStart; document.getElementById('roster').innerHTML = rosterHtml(idx, pid); bindRoster(el, idx, pid); const q2 = el.querySelector('#rosterQ'); q2.focus(); q2.setSelectionRange(pos, pos); });
    const rs = el.querySelector('#btnReset'); if (rs) rs.addEventListener('click', () => { GL.store.resetDemo(); GL.app.render(); });
  }
  V.onResize = function () { if (V.ui.p && document.getElementById('scrub')) { V.drawCharts(); if (V.ui.scrub) V.ui.scrub.redraw(); } };
})(typeof window !== 'undefined' ? window : global);
