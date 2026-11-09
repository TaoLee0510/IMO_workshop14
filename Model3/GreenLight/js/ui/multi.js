/* GreenLight · ui/multi.js
 * Multi-patient view of the clinician portal (#/clinic/multi). Left: the same patient list with tick boxes;
 * "Show selected" puts the ticked patients into the drawing area. Main: tabs ANC / platelets / monocytes /
 * lymphocytes + CAR T; under the tab one chart row per patient (one column; the same y axis for all rows of a
 * count and the time axis of the single-patient view), and the day scrubber under them. Right: every shown
 * patient's predictions as of the selected day. The day is shared; a patient in treatment is shown at
 * min(day, its data horizon). Click or drag on a chart row = move the scrubber (GL.ui.timePick).
 * The selection lives in memory only (real patient IDs never go to browser storage).
 */
(function (root) {
  const GL = root.GL = root.GL || {};
  const t = (k, p) => GL.i18n.t(k, p), esc = (s) => GL.util.escapeHtml(s), U = GL.util;
  const V = GL.views.clinic, M = V.multi = { pts: [] };
  const L = () => V.lib;

  M.render = async function (el) {
    const idx = GL.loader.index(), ids = new Set(GL.session.roster('clinic').map(r => r.id));
    V.ui.sel = new Set(Array.from(V.ui.sel || []).filter(id => ids.has(id)));
    V.ui.shown = (V.ui.shown || []).filter(id => ids.has(id));
    el.innerHTML = `<div class="clinic"><aside class="roster" id="roster">${L().rosterHtml(idx, null, true)}</aside><section class="cmain" id="mmain"></section></div>`;
    L().bindRoster(el, idx, null, true); L().restoreRoster(el);
    await M.show();
  };

  const rowHtml = (x) => `<div class="card chart-card mrow" data-pid="${esc(x.p.id)}">
      <div class="mrow-head"><a class="mrow-id" href="#/clinic/${encodeURIComponent(x.p.id)}" title="${esc(t('multi.open'))}">${esc(x.p.id)}</a>${GL.app.typeMark(x.meta)}
        <span class="tiny muted mrow-sub"></span><span class="mrow-ph"></span></div>
      <div class="plot"></div></div>`;

  M.show = async function () {
    const main = document.getElementById('mmain'); if (!main) return;
    const idx = GL.loader.index(), roster = GL.session.roster('clinic'), ids = V.ui.shown || [];
    M.pts = []; M.started = false;
    if (!ids.length) {
      main.innerHTML = `<div class="mempty">${GL.app.icon('list')}<b>${esc(t('multi.title'))}</b><span>${esc(t('multi.empty'))}</span></div>`;
      return;
    }
    main.innerHTML = `<div class="boot">${esc(t('ui.loading'))}</div>`;
    let pts;
    try { pts = await Promise.all(ids.map(id => GL.loader.patient(id))); } catch (e) { main.innerHTML = `<div class="boot">${esc(t('ui.load_error', { id: ids.join(', ') }))}</div>`; return; }
    if (V.ui.mode !== 'multi' || document.getElementById('mmain') !== main) return;
    M.pts = pts.map(p => ({ p, meta: roster.find(r => r.id === p.id) || p, real: !GL.loader.isDigital(p.id), hz: p.status === 'active' ? GL.store.patient(p).horizon : 90 }));
    M.pts.filter(x => x.real).forEach(x => GL.audit.log('open', x.p.id));
    const days = idx.asof_days, hzMax = Math.max(...M.pts.map(x => x.hz)), sel = days.filter(d => d <= hzMax);
    const want = V.ui.masOf !== undefined ? V.ui.masOf : Math.min(30, hzMax);
    const asOf = sel.reduce((b, d) => Math.abs(d - want) < Math.abs(b - want) ? d : b, sel[0]);
    M.tab = L().LIN.some(l => l[0] === GL.store.get('multiTab', 'anc')) ? GL.store.get('multiTab', 'anc') : 'anc';
    main.innerHTML = `
      <div class="psum mhead">
        <div class="pid">${esc(t('multi.title'))} <span class="chip">${esc(t('multi.n', { n: M.pts.length }))}</span> ${GL.help.btn('multi')}</div>
        <div class="mtabs" role="tablist">${L().LIN.map(([v]) => `<button role="tab" data-tab="${v}" class="${v === M.tab ? 'on' : ''}" aria-selected="${v === M.tab}">${esc(t('lin.' + v))}</button>`).join('')}</div>
        <div class="msub tiny muted">${esc(t('multi.sub'))}</div>
      </div>
      <div class="cgrid">
        <div class="cleft mleft">
          <div class="mrows" id="mrows">${M.pts.map(rowHtml).join('')}</div>
          <div class="mscrub"><div class="card scrub" id="scrub"></div></div>
        </div>
        <div class="cright" id="cright"></div>
      </div>`;
    main.querySelectorAll('.mtabs [data-tab]').forEach(b => b.addEventListener('click', () => {
      M.tab = b.dataset.tab; GL.store.set('multiTab', M.tab);
      main.querySelectorAll('.mtabs [data-tab]').forEach(x => { x.classList.toggle('on', x === b); x.setAttribute('aria-selected', x === b); });
      drawRows();
    }));
    M.scrub = GL.ui.scrubber(document.getElementById('scrub'), { days, horizon: hzMax, asOf, active: hzMax < 90, gcsf: [], onChange: (d) => M.update(d), onIngest: null, note: t('multi.scrub') });
    main.querySelectorAll('.mrow .plot').forEach(box => GL.ui.timePick(box, { days: () => sel, set: (d) => M.scrub.set(d) }));
    // a row and its table line light up together
    const hl = (id, on) => main.querySelectorAll(`[data-pid="${CSS.escape(id)}"]`).forEach(n => n.classList.toggle('hl', on));
    main.addEventListener('pointerover', (e) => { const n = e.target.closest && e.target.closest('[data-pid]'); if (n) hl(n.dataset.pid, true); });
    main.addEventListener('pointerout', (e) => { const n = e.target.closest && e.target.closest('[data-pid]'); if (n) hl(n.dataset.pid, false); });
    M.update(asOf);
  };

  let scrubLog = 0;
  M.update = function (asOf) {
    if (!M.pts.length) return;
    if (M.started && asOf !== V.ui.masOf && Date.now() - scrubLog > 5000 && M.pts.some(x => x.real)) { scrubLog = Date.now(); M.pts.filter(x => x.real).forEach(x => GL.audit.log('scrub', x.p.id, asOf)); }
    M.started = true; V.ui.masOf = asOf;
    M.pts.forEach(x => { x.eff = Math.min(asOf, x.hz); x.a = GL.model3.asofEntry(x.p, x.eff); });
    drawRows();
    const { stackBar, pct, dl } = L();
    M.pts.forEach(x => {
      const row = document.querySelector(`.mrow[data-pid="${CSS.escape(x.p.id)}"]`); if (!row) return;
      row.querySelector('.mrow-sub').textContent = x.p.status === 'active' ? `${t('roster.data_to')} ${dl(x.hz)}`
        : (x.p.label && x.p.label.primary ? `${t('est.recorded')}: ${t('ph.' + x.p.label.primary)}` : '');
      row.querySelector('.mrow-ph').innerHTML = `${stackBar(x.a.p, 'stack-sm')}<span class="tiny">${esc(t('plan.nonquick'))} <b>${pct(1 - x.a.p[0])}</b></span>`;
    });
    document.getElementById('cright').innerHTML = tableHtml(asOf);
  };

  function drawRows() {
    const v = M.tab, lin = L().LIN.find(l => l[0] === v), idx = GL.loader.index();
    const doms = M.pts.map(x => L().ydom(x.p)[v]), dom = [Math.min(...doms.map(d => d[0])), Math.max(...doms.map(d => d[1]))];
    M.pts.forEach(x => {
      const box = document.querySelector(`.mrow[data-pid="${CSS.escape(x.p.id)}"] .plot`); if (!box || !x.a) return;
      box.innerHTML = GL.charts.lineage({ w: box.clientWidth || 300, h: box.clientHeight || 100, yDomain: dom, log: true, grid: idx.grid_days,
        band: { lo: x.a.band[v][0], med: x.a.band[v][1], hi: x.a.band[v][2] }, obs: x.p.obs[v], asOf: x.eff, horizon: x.hz, color: lin[1], thr: lin[2],
        gcsf: v === 'anc' ? (x.p.design.gcsf || []).filter(w => w[0] <= x.hz).map(w => [w[0], Math.min(w[1], x.hz)]) : null, future: 'hollow' });
    });
  }

  function tableHtml(asOf) {
    const { stackBar, pct, dl } = L();
    const mean = M.pts.reduce((s, x) => s + 1 - x.a.p[0], 0) / M.pts.length;
    const rows = M.pts.map(x => {
      const rec = x.a.rec || [], recTxt = isFinite(rec[0]) && rec[0] !== null && rec[0] <= 90 ? dl(rec[0]) : '>' + dl(90);
      return `<tr data-pid="${esc(x.p.id)}"><td><a href="#/clinic/${encodeURIComponent(x.p.id)}" title="${esc(t('multi.open'))}">${esc(x.p.id)}</a>${x.eff < asOf ? ` <span class="tiny muted" title="${esc(t('multi.scrub'))}">${dl(x.eff)}</span>` : ''}</td>
        <td>${stackBar(x.a.p, 'stack-sm')}</td><td class="num"><b>${pct(1 - x.a.p[0])}</b></td><td class="num">${recTxt}</td>
        <td class="num">${x.a.nadir && x.a.nadir[1] !== null ? U.fmt(x.a.nadir[1], 2) : '–'}</td></tr>`;
    }).join('');
    return `<div class="card">
        <div class="card-head"><h3>${esc(t('est.title', { day: dl(asOf) }))}</h3><span class="chip chip-twin">${esc(t('est.twin'))}</span></div>
        <div class="kv small"><span>${esc(t('multi.mean'))}</span><b>${pct(mean)}</b></div>
        <table class="tbl mtbl"><thead><tr><th>${esc(t('multi.col.patient'))}</th><th>${esc(t('multi.col.ph'))}</th><th class="num">${esc(t('plan.nonquick'))}</th><th class="num">${esc(t('multi.col.rec'))}</th><th class="num">${esc(t('multi.col.nadir'))}</th></tr></thead><tbody>${rows}</tbody></table>
        <div class="phlegend"><span><i class="sw ph-q"></i>${esc(t('ph.quick'))}</span><span><i class="sw ph-s"></i>${esc(t('ph.second_wave'))}</span><span><i class="sw ph-d"></i>${esc(t('ph.delayed'))}</span></div>
      </div>`;
  }

  M.redraw = function () { if (M.pts.length && document.getElementById('mrows')) { drawRows(); if (M.scrub) M.scrub.redraw(); } };
})(typeof window !== 'undefined' ? window : global);
