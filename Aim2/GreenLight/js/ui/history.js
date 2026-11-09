/* GreenLight · ui/history.js — version timeline and diff. */
(function (root) {
  const GL = root.GL = root.GL || {};
  const U = GL.util, S = GL.schema, ST = GL.store, A = GL.app, R = GL.registry, O = GL.orchestrator;
  const t = (k, p) => GL.i18n.t(k, p);
  const V = A.views.history = {
    ui: { a: null, b: null, variable: 'anc' },
    async render(el) {
      const p = A.state.patient, pid = A.state.pid, st = ST.state(pid);
      await A.ensureAllVersions(pid);
      const versions = ST.listVersions(pid);
      if (!versions.length) { el.innerHTML = `<div class="empty">${t('ui.no_version')}</div>`; return; }
      const cur = A.currentVersion();
      let b = versions.find(x => x.version === V.ui.b) || cur, a = versions.find(x => x.version === V.ui.a);
      if (!a || a.version >= b.version) a = versions.filter(x => x.version < b.version).pop() || null;
      V.ui.a = a ? a.version : null; V.ui.b = b.version;
      const timeline = versions.map(v => {
        const m = v.plans.soc.outcomes.marrow; const top = A.rankedPlans(v)[0];
        const nObs = Object.values(v.dataSummary.nObs).reduce((x, y) => x + y, 0);
        return `<div class="tl-item ${v.version === b.version ? 'sel-b' : ''} ${a && v.version === a.version ? 'sel-a' : ''} ${v.asOfDay === cur.asOfDay ? 'is-current' : ''}" data-v="${v.version}">
          <div class="tl-dot"></div><div class="tl-ver">v${v.version}</div><div class="tl-day">${U.dayLabel(v.asOfDay)}</div>
          <div class="ph-mini" data-p='${JSON.stringify({ quick: m.p_quick.value, second_wave: m.p_second_wave.value, persistent: m.p_persistent.value })}'></div>
          <div class="tl-meta">${nObs} ${t('roster.obs')} · ${A.planName(top.plan)}</div><div class="tl-meta muted">${GL.i18n.fmtDate(v.createdAt)}</div>
          <div class="tl-actions"><button class="btn btn-xs" data-a="${v.version}">A</button><button class="btn btn-xs" data-b="${v.version}">B</button><button class="btn btn-xs btn-ghost" data-view="${p.demo.checkpoints.indexOf(v.asOfDay)}">${t('hist.view')}</button></div></div>`;
      }).join('');
      el.innerHTML = `
        <div class="page-head"><div><h1>${t('hist.title')} <span class="muted">· ${p.label || pid}</span></h1><p class="muted">${t('hist.subtitle')}</p></div></div>
        <div class="card"><div class="timeline">${timeline}</div></div>
        ${a ? V.diff(a, b, p) : `<div class="card"><div class="empty">${t('hist.need_two')}</div></div>`}`;
      el.querySelectorAll('.ph-mini').forEach(c => GL.charts.phenotypeBar(c, JSON.parse(c.dataset.p), { height: 12 }));
      el.querySelectorAll('[data-a]').forEach(x => x.addEventListener('click', () => { V.ui.a = Number(x.dataset.a); if (V.ui.a >= V.ui.b) V.ui.b = Math.min(versions.length, V.ui.a + 1); A.render(); }));
      el.querySelectorAll('[data-b]').forEach(x => x.addEventListener('click', () => { V.ui.b = Number(x.dataset.b); if (V.ui.a !== null && V.ui.a >= V.ui.b) V.ui.a = V.ui.b - 1 || null; A.render(); }));
      el.querySelectorAll('[data-view]').forEach(x => x.addEventListener('click', async () => { await A.setView(Number(x.dataset.view)); A.go('twin', pid); }));
      if (a) { V.drawSkill(el, a, b, p); el.querySelector('#skTabs').addEventListener('click', (e) => { const bt = e.target.closest('[data-var]'); if (!bt) return; V.ui.variable = bt.dataset.var; el.querySelectorAll('#skTabs .tab').forEach(z => z.classList.toggle('active', z.dataset.var === V.ui.variable)); V.drawSkill(el, a, b, p); }); }
    },
    diff(a, b, p) {
      const ra = A.rankedPlans(a), rb = A.rankedPlans(b);
      const newObs = p.observations.filter(o => o.day > a.asOfDay && o.day <= b.asOfDay);
      const byVar = {}; newObs.forEach(o => { byVar[o.variable] = (byVar[o.variable] || 0) + 1; });
      // parameters
      const prow = [];
      A.moduleOrder().forEach(mod => { if (!b.modules[mod.meta.id]) return; mod.provides.params.forEach(pd => { const pa = a.modules[mod.meta.id] && a.modules[mod.meta.id].params[pd.id], pb = b.modules[mod.meta.id].params[pd.id]; const f = (x) => pd.log && x < 0.01 ? x.toExponential(1) : U.fmt(x); prow.push(`<tr><td><span class="muted small">${t(mod.meta.nameKey)}</span><br>${t(pd.labelKey)}${pd.unit ? ` <span class="unit">(${pd.unit})</span>` : ''} ${GL.help.btn('param.' + pd.id, 'help-xs')}</td><td>${pa ? `${f(pa.median)} <span class="muted">[${f(pa.lo)}–${f(pa.hi)}]</span>` : '–'}</td><td><b>${f(pb.median)}</b> <span class="muted">[${f(pb.lo)}–${f(pb.hi)}]</span></td><td>${pa ? A.h.delta(pa.median, pb.median) : ''}<div class="muted small">${pa ? (() => { const w = 100 * (1 - (pb.hi - pb.lo) / ((pa.hi - pa.lo) || 1)); return w >= 0 ? t('hist.width', { w: U.fmt(w, 0) }) : t('hist.wider', { w: U.fmt(-w, 0) }); })() : ''}</div></td></tr>`); }); });
      // outcomes under SOC and under each version's top plan
      const orow = [];
      A.moduleOrder().forEach(mod => mod.provides.outcomes.forEach(def => { const oa = a.plans.soc.outcomes[mod.meta.id], ob = b.plans.soc.outcomes[mod.meta.id]; if (!ob) return; const va = oa && oa[def.id] ? oa[def.id].value : null, vb = ob[def.id] ? ob[def.id].value : null; orow.push(`<tr><td>${A.h.outcomeHead(def)}</td><td>${oa ? A.h.fmtOutcome(def, oa[def.id]) : '–'}</td><td><b>${A.h.fmtOutcome(def, ob[def.id])}</b></td><td>${A.h.delta(va, vb, def.kind === 'prob', def.direction)}</td></tr>`); }));
      // ranking moves
      const phCell = (x) => { const m = x && x.outcomes.marrow; return m && m.p_quick ? `<div class="ph-mini" data-p='${JSON.stringify({ quick: m.p_quick.value, second_wave: m.p_second_wave.value, persistent: m.p_persistent.value })}'></div>` : '–'; };
      const rrow = rb.map(x => { const o = ra.find(y => y.id === x.id); const mv = o ? o.rank - x.rank : null; return `<tr class="${x.rank <= 3 ? 'row-top' : ''}"><td>#${x.rank}</td><td>${A.planName(x.plan)} ${x.custom ? `<span class="badge">${t('ui.custom')}</span>` : ''}</td><td>${phCell(o)}</td><td>${phCell(x)}</td><td>${o ? '#' + o.rank : '–'}</td><td>${mv === null ? '' : mv > 0 ? `<span class="delta good">▲ ${mv}</span>` : mv < 0 ? `<span class="delta bad">▼ ${-mv}</span>` : '<span class="delta flat">=</span>'}</td><td>${o ? U.fmt(o.score, 2) : '–'}</td><td><b>${U.fmt(x.score, 2)}</b></td></tr>`; }).join('');
      // forecast skill
      const sk = O.evaluateForecast(a, p, b.asOfDay, a.activePlanId) || {};
      const skrow = Object.entries(sk).map(([k, s]) => `<tr><td>${t(S.variables[k].labelKey)}</td><td>${s.n}</td><td>${U.fmt(s.rmsle, 2)}</td><td>${U.pct(s.coverage90)}</td></tr>`).join('');
      const vars = Object.keys(sk); if (vars.length && !vars.includes(V.ui.variable)) V.ui.variable = vars[0];
      return `
        <div class="diff-head"><span class="chip chip-a">A · v${a.version} · ${U.dayLabel(a.asOfDay)}</span><span class="arrow">→</span><span class="chip chip-b">B · v${b.version} · ${U.dayLabel(b.asOfDay)}</span><span class="muted">${t('hist.new_data', { n: newObs.length })}: ${Object.entries(byVar).map(([k, n]) => t(S.variables[k].labelKey) + ' ' + n).join(', ') || '–'}</span></div>
        <div class="grid-1-1">
          <div class="card"><div class="card-head"><h3>${t('hist.params')} ${GL.help.btn('param_bars')}</h3></div><table class="table compact"><thead><tr><th>${t('hist.param')}</th><th>A</th><th>B</th><th>Δ</th></tr></thead><tbody>${prow.join('')}</tbody></table></div>
          <div class="card"><div class="card-head"><h3>${t('hist.outcomes')} <span class="muted small">· ${t('plan.soc')}</span></h3></div><table class="table compact"><thead><tr><th>${t('planner.outcome')}</th><th>A</th><th>B</th><th>Δ</th></tr></thead><tbody>${orow.join('')}</tbody></table></div>
        </div>
        <div class="grid-1-1">
          <div class="card"><div class="card-head"><h3>${t('hist.ranking')} ${GL.help.btn('score')}</h3></div><table class="table compact"><thead><tr><th>B</th><th>${t('hist.plan')}</th><th>${t('twin.phenotype_risk')} A</th><th>${t('twin.phenotype_risk')} B</th><th>A</th><th>${t('hist.move')}</th><th>${t('twin.score')} A</th><th>${t('twin.score')} B</th></tr></thead><tbody>${rrow}</tbody></table></div>
          <div class="card"><div class="card-head"><h3>${t('hist.skill')} ${GL.help.btn('skill')}</h3><div class="tabs" id="skTabs">${vars.map(x => `<button class="tab ${x === V.ui.variable ? 'active' : ''}" data-var="${x}">${t(S.variables[x].labelKey)}</button>`).join('')}</div></div>
            <p class="muted small">${t('hist.skill_note', { a: 'v' + a.version, plan: A.planName(a.plans[a.activePlanId].plan) })}</p>
            <div id="skChart" class="chart-box"></div>
            <table class="table compact"><thead><tr><th>${t('hist.variable')}</th><th>n</th><th>RMSLE</th><th>${t('hist.coverage')}</th></tr></thead><tbody>${skrow || `<tr><td colspan="4" class="muted">${t('hist.no_new')}</td></tr>`}</tbody></table></div>
        </div>`;
    },
    drawSkill(el, a, b, p) {
      const vid = V.ui.variable, vd = S.variables[vid]; if (!vd) return;
      const pr = a.plans[a.activePlanId]; let bands = null, grid = null, bb = null;
      for (const mid in pr.forecasts) if (pr.forecasts[mid].variables[vid]) { bands = pr.forecasts[mid].variables[vid]; grid = pr.forecasts[mid].grid; }
      const pb = b.plans[b.activePlanId] || b.plans.soc; for (const mid in pb.forecasts) if (pb.forecasts[mid].variables[vid]) bb = pb.forecasts[mid].variables[vid];
      GL.charts.trajectory(el.querySelector('#skChart'), { grid: grid || [], bands: bb, prevBands: bands, obs: p.observations.filter(o => o.variable === vid && o.day <= a.asOfDay), obsFuture: p.observations.filter(o => o.variable === vid && o.day > a.asOfDay && o.day <= b.asOfDay), asOfDay: a.asOfDay, events: p.events.filter(e => e.day_start <= b.asOfDay), variable: { id: vid, label: t(vd.labelKey), unit: vd.unit, log: vd.log, thresholds: vd.thresholds }, height: 260, xRange: [-6, Math.max(b.asOfDay + 30, 60)] });
    }
  };
})(typeof window !== 'undefined' ? window : global);
