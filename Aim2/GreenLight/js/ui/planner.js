/* GreenLight · ui/planner.js — physician plan builder (what-if). */
(function (root) {
  const GL = root.GL = root.GL || {};
  const U = GL.util, S = GL.schema, ST = GL.store, A = GL.app, P = GL.plans;
  const t = (k, p) => GL.i18n.t(k, p);
  const V = A.views.planner = {
    ui: { basePlanId: null, draft: null, result: null, variable: 'anc', timer: null },
    async render(el) {
      const p = A.state.patient, pid = A.state.pid, st = ST.state(pid), v = A.currentVersion();
      if (!v) { el.innerHTML = `<div class="empty">${t('ui.no_version')}</div>`; return; }
      await A.ensureCustomResults(v);
      const ranked = A.rankedPlans(v);
      const active = ranked.find(x => x.id === st.activePlanId) || ranked[0];
      if (!V.ui.draft || V.ui.draftVersion !== v.id) {
        const base = ranked.find(x => x.id === V.ui.basePlanId) || active;
        V.ui.draft = U.deepClone(base.plan); V.ui.draft.id = 'custom_' + Date.now().toString(36); V.ui.draft.custom = true; V.ui.draft.name = t('planner.custom_name', { base: A.planName(base.plan) }); V.ui.draft.basedOn = base.id;
        V.ui.draftVersion = v.id; V.ui.result = null;
      }
      const d = V.ui.draft;
      const stages = S.stageOrder.map(stg => {
        const locked = S.stageLocked(stg, v.asOfDay);
        const fields = P.fields.filter(f => f.stage === stg).map(f => V.field(f, d.stages[stg][f.key], locked)).join('');
        const win = S.stageWindows[stg];
        return `<fieldset class="stage ${locked ? 'locked' : ''}"><legend>${t('stage.' + stg)} ${GL.help.btn('stage.' + stg, 'help-xs')} <span class="muted small">${U.dayLabel(win[0])}${win[1] !== win[0] ? '–' + U.dayLabel(win[1]) : ''}</span> ${locked ? `<span class="badge">${t('planner.administered')}</span>` : stg === 'crs' || stg === 'cytopenia' ? `<span class="badge badge-ok">${t('planner.policy')}</span>` : ''}</legend><div class="fields">${fields}</div>${locked ? `<div class="muted small">${t('planner.locked_note')}</div>` : ''}</fieldset>`;
      }).join('');
      el.innerHTML = `
        <div class="page-head"><div><h1>${t('planner.title')} <span class="muted">· ${p.label || pid}</span></h1><p class="muted">${t('planner.subtitle')}</p></div>
          <div class="page-actions"><label class="field-inline"><span>${t('planner.start_from')}</span><select id="baseSel">${ranked.map(x => `<option value="${x.id}" ${x.id === d.basedOn ? 'selected' : ''}>#${x.rank} ${A.planName(x.plan)}</option>`).join('')}</select></label></div></div>
        <div class="grid-1-1">
          <div class="card"><div class="card-head"><h3>${t('planner.parameters')}</h3><input id="planName" class="input" value="${U.escapeHtml(d.name)}"></div><form id="planForm">${stages}</form>
            <div class="card-foot"><button class="btn btn-primary" id="btnEval">${t('planner.evaluate')}</button><button class="btn" id="btnSave" ${V.ui.result ? '' : 'disabled'}>${t('planner.save')}</button><span class="muted small" id="evalStatus"></span></div></div>
          <div class="card"><div class="card-head"><h3>${t('planner.predicted')}</h3><div class="tabs" id="pvTabs">${['anc', 'plt', 'crp', 'car_copies', 'tumor_burden'].map(x => `<button class="tab ${x === V.ui.variable ? 'active' : ''}" data-var="${x}">${t(S.variables[x].labelKey)}</button>`).join('')}</div></div>
            <div id="plChart" class="chart-box"></div>
            <div class="legend"><span><i class="lg lg-fc"></i>${t('planner.custom')}</span><span><i class="lg lg-prev"></i>${t('ui.active_plan')}: ${A.planName(active.plan)}</span></div>
            <div id="plTable"></div></div>
        </div>`;
      const form = el.querySelector('#planForm');
      form.addEventListener('input', (e) => { V.readForm(form); V.updateLabels(form); V.schedule(el, v, active); });
      form.addEventListener('change', (e) => { V.readForm(form); V.updateLabels(form); V.schedule(el, v, active, true); });
      el.querySelector('#planName').addEventListener('input', (e) => { d.name = e.target.value; });
      el.querySelector('#baseSel').addEventListener('change', (e) => { V.ui.basePlanId = e.target.value; V.ui.draft = null; A.render(); });
      el.querySelector('#btnEval').addEventListener('click', () => V.evaluate(el, v, active));
      el.querySelector('#btnSave').addEventListener('click', async () => { const plan = U.deepClone(d); plan.id = 'custom_' + Date.now().toString(36); ST.addCustomPlan(pid, plan); await A.evaluateCustom(v, plan); A.toast(t('planner.saved')); A.go('twin', pid); });
      el.querySelector('#pvTabs').addEventListener('click', (e) => { const b = e.target.closest('[data-var]'); if (!b) return; V.ui.variable = b.dataset.var; el.querySelectorAll('#pvTabs .tab').forEach(x => x.classList.toggle('active', x.dataset.var === V.ui.variable)); V.draw(el, v, active); });
      if (V.ui.result) V.draw(el, v, active); else V.evaluate(el, v, active);
    },
    field(f, val, locked) {
      const id = 'f_' + f.stage + '_' + f.key, dis = locked ? 'disabled' : '';
      const hb = GL.help.btn('field.' + f.stage + '.' + f.key, 'help-xs');
      if (f.type === 'bool') return `<div class="fld check"><input type="checkbox" id="${id}" data-stage="${f.stage}" data-key="${f.key}" ${val ? 'checked' : ''} ${dis}><label for="${id}">${t(f.labelKey)}</label> ${hb}</div>`;
      if (f.type === 'select') return `<label class="fld"><span>${t(f.labelKey)} ${hb}</span><select id="${id}" data-stage="${f.stage}" data-key="${f.key}" ${dis}>${f.options.map(o => `<option value="${o}" ${String(o) === String(val) ? 'selected' : ''}>${GL.i18n.has('opt.' + f.key + '.' + o) ? t('opt.' + f.key + '.' + o) : o}</option>`).join('')}</select></label>`;
      return `<label class="fld"><span>${t(f.labelKey)} ${hb} <b class="fval" data-for="${id}">${f.fmt(Number(val))}</b></span><input type="range" id="${id}" data-stage="${f.stage}" data-key="${f.key}" min="${f.min}" max="${f.max}" step="${f.step}" value="${val}" ${dis}></label>`;
    },
    readForm(form) {
      const d = V.ui.draft;
      form.querySelectorAll('[data-stage]').forEach(inp => {
        const f = P.fields.find(x => x.stage === inp.dataset.stage && x.key === inp.dataset.key);
        let val;
        if (f.type === 'bool') val = inp.checked; else if (f.type === 'range') val = Number(inp.value); else { val = inp.value; if (f.options.every(o => typeof o === 'number')) val = Number(val); }
        d.stages[f.stage][f.key] = val;
      });
    },
    updateLabels(form) { form.querySelectorAll('.fval').forEach(b => { const inp = form.querySelector('#' + b.dataset.for); const f = P.fields.find(x => x.stage === inp.dataset.stage && x.key === inp.dataset.key); b.textContent = f.fmt(Number(inp.value)); }); },
    schedule(el, v, active, now) { clearTimeout(V.ui.timer); V.ui.timer = setTimeout(() => V.evaluate(el, v, active), now ? 50 : 500); },
    async evaluate(el, v, active) {
      const stat = el.querySelector('#evalStatus'); if (stat) stat.textContent = t('planner.evaluating');
      const t0 = Date.now();
      const errs = S.validatePlan(V.ui.draft); if (errs.length) { if (stat) stat.textContent = errs.join('; '); return; }
      V.ui.result = await GL.orchestrator.whatIf(A.state.patient, v, V.ui.draft, ST.state(A.state.pid).weights);
      if (stat) stat.textContent = t('planner.evaluated_in', { ms: Date.now() - t0 });
      const bs = el.querySelector('#btnSave'); if (bs) bs.disabled = false;
      V.draw(el, v, active);
    },
    draw(el, v, active) {
      const r = V.ui.result; if (!r) return;
      const p = A.state.patient, vid = V.ui.variable, vd = S.variables[vid];
      let bands = null, grid = null, prevBands = null;
      for (const mid in r.forecasts) if (r.forecasts[mid].variables[vid]) { bands = r.forecasts[mid].variables[vid]; grid = r.forecasts[mid].grid; }
      for (const mid in active.forecasts) if (active.forecasts[mid].variables[vid]) prevBands = active.forecasts[mid].variables[vid];
      GL.charts.trajectory(el.querySelector('#plChart'), { grid: grid || [], bands, prevBands, obs: p.observations.filter(o => o.variable === vid && o.day <= v.asOfDay), obsFuture: [], asOfDay: v.asOfDay, events: p.events.filter(e => e.day_start <= v.asOfDay), variable: { id: vid, label: t(vd.labelKey), unit: vd.unit, log: vd.log, thresholds: vd.thresholds }, height: 300 });
      const ranked = A.rankedPlans(v); const top = ranked[0];
      const rows = [];
      const phCell = (x) => { const m = x.outcomes.marrow; return m && m.p_quick ? `<div class="ph-mini ph-row" data-p='${JSON.stringify({ quick: m.p_quick.value, second_wave: m.p_second_wave.value, persistent: m.p_persistent.value })}'></div>` : '–'; };
      rows.push(`<tr class="row-pheno"><td>${t('twin.phenotype_risk')} <span class="unit">(%)</span> ${GL.help.btn('phenotype', 'help-xs')}</td><td>${phCell(r)}</td><td>${phCell(active)}</td><td>${phCell(top)}</td></tr>`);
      A.moduleOrder().forEach(mod => mod.provides.outcomes.forEach(def => {
        const g = (x) => x.outcomes[mod.meta.id] ? A.h.fmtOutcome(def, x.outcomes[mod.meta.id][def.id]) : '–';
        rows.push(`<tr><td>${A.h.outcomeHead(def)}</td><td><b>${g(r)}</b></td><td>${g(active)}</td><td>${g(top)}</td></tr>`);
      }));
      const comp = (k, help) => `<tr class="row-comp"><td>${t('twin.c_' + k)} ${GL.help.btn(help, 'help-xs')}</td><td><b>${U.fmt(r.components[k], 2)}</b></td><td>${U.fmt(active.components[k], 2)}</td><td>${U.fmt(top.components[k], 2)}</td></tr>`;
      rows.push(comp('efficacy', 'efficacy'), comp('recovery', 'recovery'), comp('toxicity', 'toxicity'), comp('burden', 'burden'));
      rows.push(`<tr class="row-score"><td>${t('twin.score')} ${GL.help.btn('score', 'help-xs')}</td><td><b>${U.fmt(r.score, 2)}</b></td><td>${U.fmt(active.score, 2)}</td><td>${U.fmt(top.score, 2)}</td></tr>`);
      const box = el.querySelector('#plTable');
      box.innerHTML = `<table class="table compact"><thead><tr><th>${t('planner.outcome')}</th><th>${t('planner.custom')}</th><th>${t('ui.active_plan')}</th><th>#1 ${A.planName(top.plan)}</th></tr></thead><tbody>${rows.join('')}</tbody></table>`;
      box.querySelectorAll('.ph-mini').forEach(c => GL.charts.phenotypeBar(c, JSON.parse(c.dataset.p), { height: 16 }));
    }
  };
})(typeof window !== 'undefined' ? window : global);
