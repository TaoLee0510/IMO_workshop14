/* GreenLight · app.js — shell, router, patient/version lifecycle, toasts. */
(function (root) {
  const GL = root.GL = root.GL || {};
  const U = GL.util, S = GL.schema, ST = GL.store, O = GL.orchestrator;
  const t = (k, p) => GL.i18n.t(k, p);
  const A = GL.app = { state: { adapter: null, patients: [], patient: null, pid: null, view: 'roster', customResults: {}, changes: null, busy: false, theme: 'light' }, views: {} };
  GL.views = A.views;

  // ---- lifecycle -------------------------------------------------------------
  A.init = async function () {
    GL.i18n.init();
    GL.help.init();
    ST.load();
    try { A.state.theme = localStorage.getItem('gl.theme') || 'light'; } catch (e) { /* ignore */ }
    document.documentElement.setAttribute('data-theme', A.state.theme);
    A.state.adapter = GL.adapters.DemoDataAdapter;
    try { A.state.adapterInfo = await A.state.adapter.connect(); A.state.patients = await A.state.adapter.listPatients(); }
    catch (e) { A.state.adapterInfo = { ok: false, error: String(e) }; A.state.patients = []; }
    GL.i18n.onChange(() => A.render());
    window.addEventListener('hashchange', () => A.route());
    await A.route();
  };
  A.setAdapter = async function (adapter, args) {
    A.progress(true, t('ui.connecting'));
    try {
      const info = await adapter.connect(...(args || []));
      A.state.adapter = adapter; A.state.adapterInfo = info; A.state.patients = await adapter.listPatients();
      A.state.patient = null; A.state.pid = null; A.toast(t('ui.connected', { n: A.state.patients.length }));
    } catch (e) { A.toast(t('ui.connect_failed') + ': ' + e.message, 'error'); }
    A.progress(false); A.render();
  };
  A.route = async function () {
    const h = (location.hash || '#/roster').replace(/^#\/?/, '').split('/');
    const view = h[0] || 'roster', pid = h[1] ? decodeURIComponent(h[1]) : null;
    A.state.view = A.views[view] ? view : 'roster';
    if (pid && pid !== A.state.pid) await A.loadPatient(pid);
    if (!pid && ['twin', 'planner', 'history'].includes(A.state.view) && A.state.pid) { location.hash = '#/' + A.state.view + '/' + encodeURIComponent(A.state.pid); return; }
    A.render();
  };
  A.go = (view, pid) => { location.hash = '#/' + view + (pid ? '/' + encodeURIComponent(pid) : ''); };

  A.loadPatient = async function (pid) {
    const p = await A.state.adapter.getPatient(pid);
    if (!p) { A.toast(t('ui.patient_not_found', { id: pid }), 'error'); return; }
    A.state.patient = p; A.state.pid = pid; A.state.changes = null;
    ST.state(pid, p);
    await A.ensureVersions(pid);
  };
  A.checkpointDay = (i) => A.state.patient.demo.checkpoints[i];
  A.isActive = (p) => !!p && (p.demo.status || 'active') !== 'completed';
  A.currentDay = () => A.checkpointDay(ST.state(A.state.pid).viewIndex);      // day being displayed
  A.horizonDay = () => A.checkpointDay(ST.state(A.state.pid).horizonIndex);   // last day with data
  // build the versions for the given checkpoint indices if missing (sequential, with progress)
  A.ensureIndices = async function (pid, indices) {
    const p = A.state.patient, st = ST.state(pid, p);
    const need = indices.filter((i, k, arr) => i >= 0 && i < p.demo.checkpoints.length && arr.indexOf(i) === k && !ST.versionAt(pid, p.demo.checkpoints[i])).map(i => p.demo.checkpoints[i]).sort((a, b) => a - b);
    if (!need.length) return;
    A.progress(true, t('ui.building_twin', { n: need.length }));
    for (let k = 0; k < need.length; k++) {
      A.progress(true, t('ui.building_twin_step', { day: U.dayLabel(need[k]), k: k + 1, n: need.length }));
      await new Promise(r => setTimeout(r, 20));
      const v = await O.buildVersion(p, need[k], { weights: st.weights });
      ST.addVersion(v);
    }
    A.progress(false);
  };
  A.ensureVersions = async function (pid) { const st = ST.state(pid, A.state.patient); await A.ensureIndices(pid, [st.viewIndex, st.viewIndex - 1]); };
  A.ensureAllVersions = async function (pid) { const st = ST.state(pid, A.state.patient); await A.ensureIndices(pid, Array.from({ length: st.horizonIndex + 1 }, (_, i) => i)); };
  A.currentVersion = () => A.state.pid ? ST.versionAt(A.state.pid, A.currentDay()) : null;
  A.previousVersion = () => { const st = ST.state(A.state.pid); return st.viewIndex > 0 ? ST.versionAt(A.state.pid, A.checkpointDay(st.viewIndex - 1)) : null; };
  A.ingestNext = async function () {
    const pid = A.state.pid, p = A.state.patient, st = ST.state(pid);
    if (!A.isActive(p)) { A.toast(t('ui.record_complete')); return; }
    if (st.horizonIndex >= p.demo.checkpoints.length - 1) { A.toast(t('ui.no_more_data')); return; }
    const prevIdx = st.horizonIndex;
    ST.setState(pid, { horizonIndex: st.horizonIndex + 1, viewIndex: st.horizonIndex + 1 });
    await A.ensureIndices(pid, [prevIdx, prevIdx + 1]);
    const cur = A.currentVersion(), prev = ST.versionAt(pid, A.checkpointDay(prevIdx));
    A.state.changes = A.computeChanges(prev, cur);
    A.toast(t('ui.ingested', { day: U.dayLabel(cur.asOfDay), v: cur.version }));
    A.go('twin', pid); A.render();
  };
  // move the viewing point along the timeline (never beyond the data horizon)
  A.setView = async function (index) {
    const pid = A.state.pid, st = ST.state(pid);
    index = Math.max(0, Math.min(st.horizonIndex, index));
    if (index === st.viewIndex) return;
    ST.setState(pid, { viewIndex: index }); A.state.changes = null;
    await A.ensureIndices(pid, [index, index - 1]);
    A.render();
  };
  A.rewind = (index) => A.setView(index);
  A.resetDemo = async function () {
    ST.reset(); A.state.customResults = {}; A.state.changes = null;
    if (A.state.pid) { await A.loadPatient(A.state.pid); }
    A.toast(t('ui.reset_done')); A.render();
  };

  // ---- plans & ranking -------------------------------------------------------
  A.planName = (plan) => plan.custom ? plan.name : t(plan.nameKey);
  // display order: real modules first, placeholders last
  A.moduleOrder = () => GL.registry.order().slice().sort((a, b) => (a.meta.placeholder ? 1 : 0) - (b.meta.placeholder ? 1 : 0));
  A.rankedPlans = function (version, weights) {
    weights = weights || ST.state(A.state.pid).weights;
    const list = [];
    for (const id in version.plans) { const pr = version.plans[id]; const u = GL.plans.utility(pr.outcomes, pr.plan, weights); list.push({ id, plan: pr.plan, outcomes: pr.outcomes, forecasts: pr.forecasts, score: u.score, components: u.components, custom: false }); }
    const cr = A.state.customResults[version.id] || {};
    for (const id in cr) { const pr = cr[id]; const u = GL.plans.utility(pr.outcomes, pr.plan, weights); list.push({ id, plan: pr.plan, outcomes: pr.outcomes, forecasts: pr.forecasts, score: u.score, components: u.components, custom: true }); }
    list.sort((a, b) => b.score - a.score);
    list.forEach((x, i) => { x.rank = i + 1; });
    return list;
  };
  A.planResult = function (version, planId) { return A.rankedPlans(version).find(x => x.id === planId) || A.rankedPlans(version)[0]; };
  A.evaluateCustom = async function (version, plan) {
    const r = await O.whatIf(A.state.patient, version, plan, ST.state(A.state.pid).weights);
    (A.state.customResults[version.id] = A.state.customResults[version.id] || {})[plan.id] = r;
    return r;
  };
  A.ensureCustomResults = async function (version) {
    const plans = ST.customPlansFor(A.state.pid); const have = A.state.customResults[version.id] || {};
    for (const p of plans) if (!have[p.id]) await A.evaluateCustom(version, p);
  };

  // ---- change summary between two versions -----------------------------------
  A.computeChanges = function (prev, cur) {
    if (!prev || !cur) return null;
    const pid = A.state.pid, patient = A.state.patient;
    const rp = A.rankedPlans(prev), rc = A.rankedPlans(cur);
    const soc = (v) => (v.plans.soc && v.plans.soc.outcomes.marrow) || {};
    const out = { prev, cur, phen: {}, topPrev: rp[0], topCur: rc[0], rankMoves: [], params: {}, skill: O.evaluateForecast(prev, patient, cur.asOfDay, prev.activePlanId) };
    for (const k of ['p_quick', 'p_second_wave', 'p_persistent', 'p_crs3']) out.phen[k] = { from: soc(prev)[k] ? soc(prev)[k].value : null, to: soc(cur)[k] ? soc(cur)[k].value : null };
    const ti = (v) => (v.plans.soc && v.plans.soc.outcomes.tumor_immune) || {};
    out.phen.p_response = { from: ti(prev).p_response ? ti(prev).p_response.value : null, to: ti(cur).p_response ? ti(cur).p_response.value : null };
    rc.forEach(x => { const o = rp.find(y => y.id === x.id); out.rankMoves.push({ id: x.id, plan: x.plan, from: o ? o.rank : null, to: x.rank, score: x.score, scorePrev: o ? o.score : null }); });
    for (const mid in cur.modules) { out.params[mid] = {}; const pc = cur.modules[mid].params, pp = prev.modules[mid] ? prev.modules[mid].params : {}; for (const k in pc) out.params[mid][k] = { cur: pc[k], prev: pp[k] || null }; }
    out.newObs = patient.observations.filter(o => o.day > prev.asOfDay && o.day <= cur.asOfDay).length;
    return out;
  };

  // ---- shell -------------------------------------------------------------------
  A.render = function () {
    const app = document.getElementById('app');
    if (!app) return;
    const pid = A.state.pid, p = A.state.patient;
    const nav = [['roster', 'nav.roster', icon('users')], ['twin', 'nav.twin', icon('twin')], ['planner', 'nav.planner', icon('sliders')], ['history', 'nav.history', icon('clock')], ['integration', 'nav.integration', icon('plug')]];
    const navHtml = nav.map(([v, k, ic]) => {
      const needsPid = ['twin', 'planner', 'history'].includes(v);
      const dis = needsPid && !pid;
      return `<a class="nav-item ${A.state.view === v ? 'active' : ''} ${dis ? 'disabled' : ''}" href="${dis ? '#' : '#/' + v + (needsPid ? '/' + encodeURIComponent(pid) : '')}">${ic}<span>${t(k)}</span></a>`;
    }).join('');
    const langOpts = GL.i18n.languages.map(l => `<option value="${l.code}" ${GL.i18n.lang === l.code ? 'selected' : ''}>${l.label}</option>`).join('');
    const patientOpts = A.state.patients.map(x => `<option value="${x.patient_id}" ${x.patient_id === pid ? 'selected' : ''}>${x.label || x.patient_id}</option>`).join('');
    const st = pid ? ST.state(pid) : null;
    const lastIdx = p ? p.demo.checkpoints.length - 1 : 0;
    app.innerHTML = `
      <aside class="sidebar">
        <div class="brand">${A.logo()}<div class="brand-text"><div class="brand-name"><span class="b1">Green</span><span class="b2">Light</span></div><div class="brand-sub">${t('brand.tagline')}</div></div></div>
        <nav class="nav">${navHtml}</nav>
        ${A.sidebarArt()}
        <div class="sidebar-foot">
          <div class="foot-row"><span class="dot ${A.state.adapterInfo && A.state.adapterInfo.ok ? 'ok' : 'bad'}"></span>${t('adapter.' + A.state.adapter.id)} · ${A.state.patients.length} ${t('ui.patients')}</div>
          <div class="foot-row muted">${t('brand.org')}</div>
          <div class="foot-row muted">v0.5 · ${t('ui.demo_build')}</div>
        </div>
      </aside>
      <div class="main">
        <header class="topbar">
          <div class="topbar-left">
            <label class="field-inline"><span>${t('ui.patient')}</span>
              <select id="patientSelect"><option value="">${t('ui.select_patient')}</option>${patientOpts}</select></label>
            ${p ? `<span class="pill ${A.isActive(p) ? 'pill-active' : 'pill-done'}">${t(A.isActive(p) ? 'status.active' : 'status.completed')}</span><span class="chip chip-day">${t('ui.as_of')} <b>${U.dayLabel(A.currentDay())}</b> ${GL.help.btn('asof', 'help-xs')}</span><span class="chip">${t('ui.twin_version')} <b>v${(A.currentVersion() || {}).version || '–'}</b> ${GL.help.btn('version', 'help-xs')}</span>` : ''}
          </div>
          <div class="topbar-right">
            ${p && A.isActive(p) ? `<button class="btn btn-primary" id="btnIngest" ${st.horizonIndex >= lastIdx ? 'disabled' : ''}>${icon('download')} ${t('ui.ingest_next')}${st.horizonIndex < lastIdx ? ' (' + U.dayLabel(A.checkpointDay(st.horizonIndex + 1)) + ')' : ''}</button>` : ''}
            <select id="langSelect" class="lang" title="${t('ui.language')}">${langOpts}</select>
            <button class="btn btn-ghost" id="btnTheme" title="${t('ui.theme')}">${A.state.theme === 'dark' ? icon('sun') : icon('moon')}</button>
          </div>
        </header>
        <main class="content" id="view"></main>
        <footer class="disclaimer">${t('ui.disclaimer')}</footer>
      </div>
      <div id="toasts" class="toasts"></div>
      <div id="progress" class="progress-overlay" style="display:none"><div class="progress-card"><div class="spinner"></div><div id="progressText"></div></div></div>`;
    document.getElementById('patientSelect').addEventListener('change', (e) => { if (e.target.value) A.go(['twin', 'planner', 'history'].includes(A.state.view) ? A.state.view : 'twin', e.target.value); });
    document.getElementById('langSelect').addEventListener('change', (e) => GL.i18n.set(e.target.value));
    document.getElementById('btnTheme').addEventListener('click', () => { A.state.theme = A.state.theme === 'dark' ? 'light' : 'dark'; try { localStorage.setItem('gl.theme', A.state.theme); } catch (e) { /* ignore */ } document.documentElement.setAttribute('data-theme', A.state.theme); A.render(); });
    const bi = document.getElementById('btnIngest'); if (bi) bi.addEventListener('click', () => A.ingestNext());
    const view = A.views[A.state.view];
    const vc = document.getElementById('view');
    if (['twin', 'planner', 'history'].includes(A.state.view) && !p) { vc.innerHTML = `<div class="empty">${t('ui.pick_patient_first')}</div>`; return; }
    Promise.resolve(view.render(vc)).catch(e => { vc.innerHTML = `<div class="empty error">${U.escapeHtml(String(e && e.stack || e))}</div>`; console.error(e); });
  };

  A.toast = function (msg, kind) {
    const box = document.getElementById('toasts'); if (!box) return;
    const d = document.createElement('div'); d.className = 'toast ' + (kind || ''); d.textContent = msg; box.appendChild(d);
    setTimeout(() => { d.classList.add('out'); setTimeout(() => d.remove(), 400); }, 3200);
  };
  A.progress = function (show, text) {
    const o = document.getElementById('progress'); if (!o) return;
    o.style.display = show ? 'flex' : 'none'; const tx = document.getElementById('progressText'); if (tx && text) tx.textContent = text;
  };
  A.logo = () => `<svg class="logo" viewBox="0 0 48 48" aria-hidden="true"><defs><linearGradient id="gl-logo-g" x1="0" y1="0" x2="1" y2="1"><stop offset="0" stop-color="#4ade80"/><stop offset=".55" stop-color="#16a34a"/><stop offset="1" stop-color="#0e7490"/></linearGradient></defs><rect x="2" y="2" width="44" height="44" rx="12" fill="url(#gl-logo-g)"/><circle cx="24" cy="24" r="16" fill="none" stroke="#fff" stroke-width="1.6" opacity=".35"/><path d="M24 8a16 16 0 0 1 16 16" fill="none" stroke="#fff" stroke-width="2.6" stroke-linecap="round"/><circle cx="24" cy="24" r="9.5" fill="none" stroke="#fff" stroke-width="2.4" opacity=".9"/><circle cx="24" cy="24" r="4.2" fill="#fff"/><circle cx="24" cy="24" r="4.2" fill="#fff" opacity=".5"><animate attributeName="r" values="4.2;8;4.2" dur="2.8s" repeatCount="indefinite"/><animate attributeName="opacity" values=".5;0;.5" dur="2.8s" repeatCount="indefinite"/></circle></svg>`;
  // sidebar illustration depends on the current view
  A.sidebarArt = function () {
    const v = A.state.view, p = A.state.patient;
    let variant = 'pair', caption = t('ui.art_caption'), opts = {};
    if (v === 'twin' && p) {
      const ver = A.currentVersion(); const m = ver && ver.plans.soc ? ver.plans.soc.outcomes.marrow : null;
      const risk = m ? (m.p_persistent.value >= 0.4 ? 'red' : (m.p_second_wave.value + m.p_persistent.value) >= 0.4 ? 'amber' : 'green') : 'green';
      variant = 'model'; opts = { risk, completed: !A.isActive(p) };
      caption = t('art.caption.model', { id: p.label || p.patient_id, v: ver ? ver.version : '–', day: U.dayLabel(A.currentDay()) });
    } else if (v === 'planner') { variant = 'planner'; caption = t('art.caption.planner'); }
    else if (v === 'history') { variant = 'history'; caption = t('art.caption.history'); }
    else if (v === 'integration') { variant = 'integration'; caption = t('art.caption.integration'); }
    return `<div class="sidebar-art art-${variant}" aria-hidden="true">${GL.figure.svg({ id: 'glf-side', variant, width: 204, height: 186, ...opts })}<div class="art-caption">${caption}</div></div>`;
  };
  function icon(name) {
    const paths = {
      users: '<path d="M16 21v-2a4 4 0 0 0-4-4H6a4 4 0 0 0-4 4v2"/><circle cx="9" cy="7" r="4"/><path d="M22 21v-2a4 4 0 0 0-3-3.87M16 3.13a4 4 0 0 1 0 7.75"/>',
      twin: '<circle cx="12" cy="12" r="3"/><path d="M12 2v3M12 19v3M2 12h3M19 12h3M5 5l2 2M17 17l2 2M5 19l2-2M17 7l2-2"/>',
      sliders: '<path d="M4 21v-7M4 10V3M12 21v-9M12 8V3M20 21v-5M20 12V3M1 14h6M9 8h6M17 16h6"/>',
      clock: '<circle cx="12" cy="12" r="10"/><path d="M12 6v6l4 2"/>',
      plug: '<path d="M12 22v-5M9 8V2M15 8V2M6 8h12v4a6 6 0 0 1-12 0z"/>',
      download: '<path d="M21 15v4a2 2 0 0 1-2 2H5a2 2 0 0 1-2-2v-4M7 10l5 5 5-5M12 15V3"/>',
      sun: '<circle cx="12" cy="12" r="5"/><path d="M12 1v2M12 21v2M4.2 4.2l1.4 1.4M18.4 18.4l1.4 1.4M1 12h2M21 12h2M4.2 19.8l1.4-1.4M18.4 5.6l1.4-1.4"/>',
      moon: '<path d="M21 12.8A9 9 0 1 1 11.2 3a7 7 0 0 0 9.8 9.8z"/>'
    };
    return `<svg class="ic" viewBox="0 0 24 24" fill="none" stroke="currentColor" stroke-width="2" stroke-linecap="round" stroke-linejoin="round" aria-hidden="true">${paths[name] || ''}</svg>`;
  }
  A.icon = icon;

  // ---- shared view helpers -----------------------------------------------------
  A.h = {
    // dir: 'higher' (increase is good), 'lower' (increase is bad) or undefined (neutral)
    delta(a, b, isPct, dir) {
      if (a === null || a === undefined || b === null || b === undefined) return '';
      const d = b - a, mag = Math.max(Math.abs(a), Math.abs(b));
      const dec = isPct ? 0 : mag >= 20 ? 0 : mag >= 2 ? 1 : 2;
      if (Math.abs(d) < (isPct ? 0.005 : Math.pow(10, -dec) / 2)) return `<span class="delta flat">=</span>`;
      const cls = !dir ? 'neutral' : (d > 0) === (dir === 'higher') ? 'good' : 'bad';
      return `<span class="delta ${cls}">${d > 0 ? '▲' : '▼'} ${isPct ? U.pct(Math.abs(d)) : U.fmt(Math.abs(d), dec)}</span>`;
    },
    patientHeader(p) {
      const h = p.history || {}; const mut = p.mutations || {};
      const chips = [];
      chips.push(`${p.demographics.age} ${t('ui.years')} · ${t('sex.' + p.demographics.sex)}`);
      chips.push(`${t('ui.product')} ${p.product} (${p.costim})`);
      chips.push(`${t('ui.prior_lines')}: ${h.prior_lines}`);
      if (h.prior_asct) chips.push(t('ui.prior_asct'));
      if (h.prior_bendamustine) chips.push(t('ui.prior_bendamustine'));
      if (h.marrow_involvement) chips.push(t('ui.marrow_involvement'));
      chips.push(`${t('ui.bridging')}: ${t('bridging.' + (p.bridging || 'none'))}`);
      if (p.ld && p.ld.dose_reduced) chips.push(t('ui.ld_reduced'));
      if (p.crs && A.state.pid === p.patient_id && A.currentVersion() && A.currentVersion().dataSummary.crsKnown) chips.push(`${t('ui.crs_grade')} ${p.crs.grade} (${U.dayLabel(p.crs.onset_day)}${p.crs.end_day ? '–' + U.dayLabel(p.crs.end_day) : ''})`);
      if (mut.tested) chips.push(mut.variants.length ? mut.variants.map(v => `${v.gene} ${v.protein_change} (VAF ${U.pct(v.vaf, 1)})`).join(', ') : t('ui.ch_negative')); else chips.push(t('ui.ch_untested'));
      return chips.map(c => `<span class="chip chip-sm">${c}</span>`).join('');
    },
    outcomeLabel(id) { return t('out.' + id); },
    // label + unit + help button for an outcome definition
    outcomeHead(def) { const unit = def.kind === 'prob' ? '%' : (def.unit || ''); return `${t(def.labelKey)}${unit ? ` <span class="unit">(${unit})</span>` : ''} ${GL.help.btn('out.' + def.id)}`; },
    fmtOutcome(def, o) { if (!o) return '–'; if (def.kind === 'prob') return U.pct(o.value); const v = U.fmt(o.value, def.unit === 'd' ? 0 : undefined); return o.lo !== undefined ? `${v} <span class="muted">[${U.fmt(o.lo, def.unit === 'd' ? 0 : undefined)}–${U.fmt(o.hi, def.unit === 'd' ? 0 : undefined)}]</span>` : v; }
  };
})(typeof window !== 'undefined' ? window : global);
