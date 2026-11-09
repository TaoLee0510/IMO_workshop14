/* GreenLight · ui/twin.js — per-patient twin dashboard. */
(function (root) {
  const GL = root.GL = root.GL || {};
  const U = GL.util, S = GL.schema, ST = GL.store, A = GL.app, R = GL.registry;
  const t = (k, p) => GL.i18n.t(k, p);
  const V = A.views.twin = {
    ui: { variable: 'anc', planId: null, showPrev: true },
    async render(el) {
      const p = A.state.patient, pid = A.state.pid, st = ST.state(pid), v = A.currentVersion(), prev = A.previousVersion();
      if (!v) { el.innerHTML = `<div class="empty">${t('ui.no_version')}</div>`; return; }
      await A.ensureCustomResults(v);
      const ranked = A.rankedPlans(v), top3 = ranked.slice(0, 3);
      const active = ranked.find(x => x.id === st.activePlanId) || ranked[0];
      const shown = ranked.find(x => x.id === V.ui.planId) || active;
      V.ui.planId = shown.id;
      const m = GL.plans.marrowOut(active.outcomes), ti = active.outcomes.tumor_immune || null;
      const pm = prev ? (A.planResult(prev, active.id) || {}).outcomes || {} : {};
      const pmm = GL.plans.marrowOut(pm), pti = pm.tumor_immune || {};
      const ch = A.state.changes;
      const g = (o, k) => o && o[k] ? o[k].value : null;
      const isActive = A.isActive(p); const nObs = Object.values(v.dataSummary.nObs).reduce((a, b) => a + b, 0);
      el.innerHTML = `
        <div class="card patient-banner">
          <div class="pb-top"><div class="pb-id"><div class="eyebrow">${t('eyebrow.patient')}</div><h1>${p.label || pid}</h1><span class="pill ${isActive ? 'pill-active' : 'pill-done'}">${t(A.statusKey(p))}</span></div>
            <div class="pb-chips chips">${A.h.patientHeader(p)}</div>
            <div class="pb-meta"><div class="pb-meta-item"><span class="muted">${t('ui.twin_version')}</span><b>v${v.version}</b></div><div class="pb-meta-item"><span class="muted">${t('ui.as_of')}</span><b>${U.dayLabel(v.asOfDay)}</b></div><div class="pb-meta-item"><span class="muted">${t('twin.obs_label')}</span><b>${nObs}</b></div><div class="pb-meta-item"><span class="muted">${t('ui.active_plan')}</span><b>${A.planName(active.plan)}</b></div></div></div>
          <div class="pb-timeline" id="tlBox"></div>
        </div>
        ${ch ? V.changesBanner(ch) : ''}
        <div class="tiles tiles-6">
          <div class="tile tile-pheno"><div class="tile-label">${t('twin.phenotype_risk')} <span class="unit">(%)</span> ${GL.help.btn('phenotype')} <span class="muted">· ${t('plan.short.' + (active.plan.custom ? 'custom' : active.id))}</span></div><div class="ph-big" id="phBig"></div>
            <div class="ph-legend">${Object.keys(GL.plans.phenoProbs(m)).map(k => `<span><i class="sw sw-${k}"></i>${t('pheno.' + k)} ${U.pct(g(m, 'p_' + k))} ${A.h.delta(g(pmm, 'p_' + k), g(m, 'p_' + k), true, k === 'quick' ? 'higher' : 'lower')}</span>`).join('')}</div></div>
          <div class="tile"><div class="tile-label">${t('out.p_response')} ${GL.help.btn('out.p_response')}</div><div class="tile-value">${ti ? U.pct(g(ti, 'p_response')) : '–'}</div><div class="tile-sub">${A.h.delta(g(pti, 'p_response'), g(ti, 'p_response'), true, 'higher')} ${t('mod.ti.short')}</div></div>
          <div class="tile"><div class="tile-label">${t('twin.age_tile')} ${GL.help.btn('age_since_mpp')}</div><div class="tile-value">${U.fmt(g(m, 'age_anc'), 1)} <span class="unit">d</span></div><div class="tile-sub">${t('twin.age_tile_sub', { gen: U.fmt(g(m, 'gen_anc'), 1) })} ${A.h.delta(g(pmm, 'age_anc'), g(m, 'age_anc'), false)}</div></div>
          <div class="tile"><div class="tile-label">${t('out.sev_neut_days')} ${GL.help.btn('out.sev_neut_days')}</div><div class="tile-value">${U.fmt(g(m, 'sev_neut_days'), 0)} <span class="unit">d</span></div><div class="tile-sub">${A.h.delta(g(pmm, 'sev_neut_days'), g(m, 'sev_neut_days'), false, 'lower')} · ${t('twin.d0_90')}</div></div>
          <div class="tile"><div class="tile-label">${t('out.peak_expansion')} ${GL.help.btn('out.peak_expansion')}</div><div class="tile-value">${ti ? U.fmt(g(ti, 'peak_expansion'), 1) : '–'}</div><div class="tile-sub">×10³ copies/µg</div></div>
          <div class="tile tile-accent"><div class="tile-label">${t('twin.recommended')} ${GL.help.btn('recommended')}</div><div class="tile-value small">${A.planName(top3[0].plan)}</div><div class="tile-sub">${t('twin.score')} ${U.fmt(top3[0].score, 2)} ${ch && ch.topPrev && ch.topPrev.id !== top3[0].id ? '<span class="delta neutral">' + t('twin.changed') + '</span>' : ''}</div></div>
        </div>
        ${V.whenWhy(m, p)}
        ${V.phenoOutlook(ranked, active, shown, v, st)}
        <div class="grid-2-1">
          <div class="card">
            <div class="eyebrow">${t('eyebrow.forecast')}</div>
            <div class="card-head"><div class="tabs" id="varTabs">${V.variableTabs(v)}</div>
              <div class="card-tools"><label class="field-inline"><span>${t('twin.plan_shown')}</span><select id="planShown">${ranked.map(x => `<option value="${x.id}" ${x.id === shown.id ? 'selected' : ''}>#${x.rank} ${A.planName(x.plan)}${x.id === active.id ? ' · ' + t('ui.active') : ''}</option>`).join('')}</select></label>
              <label class="check"><input type="checkbox" id="showPrev" ${V.ui.showPrev ? 'checked' : ''} ${prev ? '' : 'disabled'}><span>${t('twin.show_prev')}</span></label></div></div>
            <div id="trajChart" class="chart-box"></div>
            <div class="muted small var-note" id="varNote"></div>
            <div class="legend"><span><i class="lg lg-obs"></i>${t('legend.observed')}</span><span><i class="lg lg-future"></i>${t('legend.future_obs')}</span><span><i class="lg lg-fit"></i>${t('legend.fit')}</span><span><i class="lg lg-fc"></i>${t('legend.forecast')}</span><span><i class="lg lg-band"></i>${t('legend.band')}</span><span><i class="lg lg-prev"></i>${t('legend.prev')}</span><span><i class="lg lg-ld"></i>${t('legend.ld')}</span><span><i class="lg lg-ev"></i>${t('legend.events')}</span><span>${GL.help.btn('band')}</span></div>
          </div>
          <div class="card">
            <div class="eyebrow">${t('eyebrow.state')}</div>
            <div class="card-head"><h3>${t('twin.state')} ${GL.help.btn('param_bars')}</h3><span class="muted small">v${v.version} · ${U.dayLabel(v.asOfDay)}${prev ? ' · ' + t('twin.vs') + ' v' + prev.version : ''}</span></div>
            ${A.moduleOrder().filter(mod => v.modules[mod.meta.id]).map(mod => V.moduleState(mod, v, prev)).join('')}
          </div>
        </div>
        <div class="card">
          <div class="eyebrow">${t('eyebrow.plans')}</div>
          <div class="card-head"><h3>${t('twin.top3')} ${GL.help.btn('score')}</h3>
            <div class="card-tools weights"><label class="field-inline slider"><span>${t('twin.w_safety')}</span><input type="range" id="wSlider" min="0.2" max="1.8" step="0.1" value="${st.weights.eff}"><span>${t('twin.w_efficacy')}</span></label><label class="field-inline slider"><span>${t('twin.w_recovery')}</span><input type="range" id="wRec" min="0" max="2" step="0.1" value="${st.weights.rec}"><b id="wRecVal">${U.fmt(st.weights.rec, 1)}</b></label>${GL.help.btn('weights')}</div></div>
          <div class="plans-row">${top3.map(x => V.planCard(x, active, ranked[0])).join('')}</div>
          <div class="muted small" style="margin-top:8px">${t('twin.rank_note', { n: ranked.length })}</div>
        </div>`;
      GL.charts.phenotypeBar(el.querySelector('#phBig'), GL.plans.phenoProbs(m), { height: 22 });
      el.querySelectorAll('.ww-bar').forEach(c => V.drawShareBar(c, JSON.parse(c.dataset.p), c.dataset.kind));
      el.querySelectorAll('.ph-mini').forEach(c => GL.charts.phenotypeBar(c, JSON.parse(c.dataset.p), { height: c.classList.contains('ph-row') ? 18 : 14 }));
      el.querySelectorAll('[data-show]').forEach(b => b.addEventListener('click', () => { V.ui.planId = b.dataset.show; A.render(); }));
      el.querySelector('#wRec').addEventListener('input', (e) => { el.querySelector('#wRecVal').textContent = U.fmt(Number(e.target.value), 1); });
      el.querySelector('#wRec').addEventListener('change', (e) => { ST.setState(pid, { weights: Object.assign({}, st.weights, { rec: Number(e.target.value) }) }); A.render(); });
      GL.ui.timeline(el.querySelector('#tlBox'), { patient: p, checkpoints: p.demo.checkpoints, horizonIndex: st.horizonIndex, viewIndex: st.viewIndex, active: isActive, onSelect: (i) => A.setView(i) });
      V.drawTrajectory(el, v, prev, shown);
      el.querySelector('#varTabs').addEventListener('click', (e) => { const b = e.target.closest('[data-var]'); if (!b) return; V.ui.variable = b.dataset.var; el.querySelectorAll('#varTabs .tab').forEach(x => x.classList.toggle('active', x.dataset.var === V.ui.variable)); V.drawTrajectory(el, v, prev, shown); });
      el.querySelector('#planShown').addEventListener('change', (e) => { V.ui.planId = e.target.value; A.render(); });
      el.querySelector('#showPrev').addEventListener('change', (e) => { V.ui.showPrev = e.target.checked; V.drawTrajectory(el, v, prev, shown); });
      el.querySelector('#wSlider').addEventListener('change', (e) => { const w = Number(e.target.value); ST.setState(pid, { weights: Object.assign({}, st.weights, { eff: w, tox: 2 - w }) }); A.render(); });
      el.querySelectorAll('[data-adopt]').forEach(b => b.addEventListener('click', () => { ST.setState(pid, { activePlanId: b.dataset.adopt }); V.ui.planId = b.dataset.adopt; A.toast(t('twin.adopted')); A.render(); }));
      el.querySelectorAll('[data-builder]').forEach(b => b.addEventListener('click', () => { A.views.planner.ui.basePlanId = b.dataset.builder; A.go('planner', pid); }));
      const dismiss = el.querySelector('#dismissChanges'); if (dismiss) dismiss.addEventListener('click', () => { A.state.changes = null; A.render(); });
      const toHist = el.querySelector('#toHistory'); if (toHist) toHist.addEventListener('click', () => A.go('history', pid));
    },
    variableTabs(v) {
      const vars = []; R.order().forEach(mod => (mod.provides.variables || []).forEach(x => { if (!vars.includes(x)) vars.push(x); }));
      const order = Object.keys(S.variables); vars.sort((a, b) => order.indexOf(a) - order.indexOf(b));
      if (!vars.includes(V.ui.variable)) V.ui.variable = vars[0];
      return vars.map(x => `<button class="tab ${x === V.ui.variable ? 'active' : ''}" data-var="${x}">${t(S.variables[x].labelKey)}</button>`).join('');
    },
    drawTrajectory(el, v, prev, shown) {
      const p = A.state.patient, vid = V.ui.variable, vd = S.variables[vid];
      let bands = null, prevBands = null, grid = null;
      for (const mid in shown.forecasts) if (shown.forecasts[mid].variables[vid]) { bands = shown.forecasts[mid].variables[vid]; grid = shown.forecasts[mid].grid; }
      if (prev && V.ui.showPrev) { const pr = A.planResult(prev, shown.id); if (pr) for (const mid in pr.forecasts) if (pr.forecasts[mid].variables[vid]) prevBands = pr.forecasts[mid].variables[vid]; }
      const obs = p.observations.filter(o => o.variable === vid && o.day <= v.asOfDay), obsFuture = p.observations.filter(o => o.variable === vid && o.day > v.asOfDay);
      const events = p.events.filter(e => e.day_start <= v.asOfDay);
      GL.charts.trajectory(el.querySelector('#trajChart'), { grid: grid || [], bands, prevBands, obs, obsFuture, asOfDay: v.asOfDay, events, variable: { id: vid, label: t(vd.labelKey), unit: vd.unit, log: vd.log, thresholds: vd.thresholds }, height: 320 });
      const note = el.querySelector('#varNote'); if (note) note.textContent = vid === 'alc' ? t('twin.alc_note') : '';
    },
    moduleState(mod, v, prev) {
      const mi = v.modules[mod.meta.id], pi = prev ? prev.modules[mod.meta.id] : null;
      const params = mod.provides.params.map(pd => ({ id: pd.id, help: 'param.' + pd.id, label: t(pd.labelKey), median: mi.params[pd.id].median, lo: mi.params[pd.id].lo, hi: mi.params[pd.id].hi, prev: pi && pi.params[pd.id] ? pi.params[pd.id] : null, range: [pd.lo, pd.hi], log: pd.log, fmt: (x) => pd.log ? (x < 0.01 ? x.toExponential(1) : U.fmt(x)) : U.fmt(x) }));
      const cov = Object.entries(mi.coverage).map(([k, n]) => `<span class="chip chip-sm ${n ? '' : 'chip-off'}">${t(S.variables[k].labelKey)} ${n}</span>`).join('');
      const id = 'pb_' + mod.meta.id;
      setTimeout(() => { const c = document.getElementById(id); if (c) GL.charts.paramBars(c, params); }, 0);
      return `<div class="module-block"><div class="module-head"><b>${t(mod.meta.nameKey)}</b>${mod.meta.placeholder ? `<span class="badge badge-warn">${t('ui.placeholder')}</span>` : ''}<span class="badge">${mod.meta.kind}</span><span class="muted small">${t('twin.fit_ms', { ms: mi.fitMs })} · ESS ${U.fmt(mi.ess, 0)} ${GL.help.btn('ess', 'help-xs')}</span></div>
        <div class="chips">${cov}${GL.help.btn('coverage', 'help-xs')}</div><div id="${id}" class="param-box"></div>
        <div class="muted small">${t('twin.param_note')}</div></div>`;
    },
    planCard(x, active, top) {
      const soc = GL.schema.defaultPlan(); const d = GL.plans.diff(soc, x.plan);
      const diffs = d.length ? d.map(c => `<li><span class="muted">${t(c.labelKey)}</span> ${GL.plans.formatValue(c.field, c.to)}</li>`).join('') : `<li class="muted">${t('plan.same_as_soc')}</li>`;
      const m = GL.plans.marrowOut(x.outcomes), ti = x.outcomes.tumor_immune;
      const comp = x.components;
      const bar = (label, val, max, cls) => `<div class="cbar"><span>${label}</span><div class="cbar-track"><div class="cbar-fill ${cls}" style="width:${Math.min(100, 100 * (val || 0) / max)}%"></div></div><b>${U.fmt(val, 2)}</b></div>`;
      return `<div class="plan-card ${x.id === active.id ? 'is-active' : ''} ${x.rank === 1 ? 'is-top' : ''}">
        <div class="plan-rank">#${x.rank}</div>
        <div class="plan-name">${A.planName(x.plan)} ${x.custom ? `<span class="badge">${t('ui.custom')}</span>` : GL.help.btn('plan.' + x.id)} ${x.id === active.id ? `<span class="badge badge-ok">${t('ui.active')}</span>` : ''}</div>
        <div class="plan-score">${t('twin.score')} <b>${U.fmt(x.score, 2)}</b> ${GL.help.btn('score', 'help-xs')}</div>
        ${bar(t('twin.c_efficacy') + ' ' + GL.help.btn('efficacy', 'help-xs'), comp.efficacy, 1, 'c-eff')}${bar(t('twin.c_recovery') + ' ' + GL.help.btn('recovery', 'help-xs'), comp.recovery, 1, 'c-rec')}${bar(t('twin.c_toxicity') + ' ' + GL.help.btn('toxicity', 'help-xs'), comp.toxicity, 1, 'c-tox')}${bar(t('twin.c_burden') + ' ' + GL.help.btn('burden', 'help-xs'), comp.burden, 0.3, 'c-bur')}
        <div class="plan-ph-label muted small">${t('twin.phenotype_risk')} ${GL.help.btn('phenotype', 'help-xs')}</div>
        <div class="ph-mini plan-ph" data-p='${JSON.stringify(GL.plans.phenoProbs(m))}'></div>
        <div class="plan-outcomes">${Object.entries(GL.plans.phenoProbs(m)).map(([k, pv]) => `<span><i class="sw sw-${k}"></i>${t('pheno.' + k)} <b>${U.pct(pv)}</b></span>`).join('')}<span>${t('out.p_response')} <b>${ti ? U.pct(ti.p_response.value) : '–'}</b></span></div>
        <div class="plan-diff"><div class="muted small">${t('twin.vs_soc')}</div><ul>${diffs}</ul></div>
        <div class="plan-actions"><button class="btn btn-sm btn-primary" data-adopt="${x.id}" ${x.id === active.id ? 'disabled' : ''}>${t('twin.adopt')}</button><button class="btn btn-sm" data-builder="${x.id}">${t('twin.customise')}</button></div></div>`;
    },
    // Recovery phenotype under every candidate plan: which plan best avoids persistent cytopenia / second wave
    phenoOutlook(ranked, active, shown, v, st) {
      const g = (x, k) => { const mm = GL.plans.marrowOut(x.outcomes); return mm[k] ? mm[k].value : null; };
      const best = GL.plans.phenoBest(ranked);
      const keys = ranked.length ? Object.keys(GL.plans.phenoProbs(GL.plans.marrowOut(ranked[0].outcomes))) : ['quick', 'second_wave', 'persistent'];
      const strip = [['quick', 'twin.best_quick'], ['second_wave', 'twin.best_sw'], ['delayed', 'twin.best_del'], ['persistent', 'twin.best_pers']].filter(([k]) => keys.includes(k)).map(([k, lk]) => {
        const b = best[k]; if (!b) return '';
        return `<div class="best best-${k}"><span class="best-k"><i class="sw sw-${k}"></i>${t(lk)}</span><b>${A.planName(b.plan)}</b><span class="best-v">${t('pheno.' + k)} <b>${U.pct(b.value)}</b> · #${b.rank}</span></div>`;
      }).join('');
      const rows = ranked.map(x => {
        const ti = x.outcomes.tumor_immune;
        const ph = {}; keys.forEach(k => { ph[k] = g(x, 'p_' + k) || 0; });
        const badges = `${x.custom ? `<span class="badge">${t('ui.custom')}</span>` : GL.help.btn('plan.' + x.id, 'help-xs')}${x.rank === 1 ? ` <span class="badge badge-ok">${t('twin.recommended')}</span>` : ''}${x.id === active.id ? ` <span class="badge">${t('ui.active')}</span>` : ''}`;
        return `<tr class="${x.rank === 1 ? 'row-top' : ''} ${x.id === active.id ? 'row-active' : ''}"><td class="num"><b>#${x.rank}</b></td><td class="plan-cell">${A.planName(x.plan)} ${badges}</td>
          <td class="dist-cell"><div class="ph-mini ph-row" data-p='${JSON.stringify(ph)}'></div></td>
          ${keys.map(k => `<td class="num ph-${k}">${U.pct(ph[k])}</td>`).join('')}
          <td class="num">${ti && ti.p_response ? U.pct(ti.p_response.value) : '–'}</td><td class="num"><b>${U.fmt(x.score, 2)}</b></td>
          <td class="actions"><button class="btn btn-xs" data-show="${x.id}" ${x.id === shown.id ? 'disabled' : ''}>${t('twin.show')}</button><button class="btn btn-xs btn-primary" data-adopt="${x.id}" ${x.id === active.id ? 'disabled' : ''}>${t('twin.adopt')}</button></td></tr>`;
      }).join('');
      return `<div class="card outlook"><div class="eyebrow">${t('eyebrow.outlook')}</div>
        <div class="card-head"><h3>${t('twin.pheno_by_plan')} ${GL.help.btn('pheno_by_plan')}</h3><span class="muted small">${t('twin.pheno_as_of', { n: ranked.length, day: U.dayLabel(v.asOfDay) })}</span></div>
        <div class="best-strip">${strip}</div>
        <div class="table-wrap"><table class="table compact pheno-table"><thead><tr><th>#</th><th>${t('twin.col_plan')}</th><th>${t('twin.col_dist')} <span class="unit">(%)</span></th>${keys.map(k => `<th class="num"><i class="sw sw-${k}"></i>${t('pheno.' + k)}</th>`).join('')}<th class="num">${t('out.p_response')}</th><th class="num">${t('twin.score')}</th><th></th></tr></thead><tbody>${rows}</tbody></table></div>
        <div class="muted small" style="margin-top:8px">${t('twin.pheno_by_plan_note', { day: U.dayLabel(v.asOfDay) })}</div></div>`;
    },
    // "When was the stem cell pool hit · why" (Aim 4 module only): counterfactual readouts on the
    // current posterior; follows the timeline because it reads the displayed version
    whenWhy(m, p) {
      if (!m || !m.p_a_driven) return '';
      const v = A.currentVersion(); const pr = v && v.plans.soc ? v.plans.soc.forecasts : null;
      const at = (pr && pr.marrow4 && pr.marrow4.extra) ? pr.marrow4.extra.attribution : null;
      const st = v && v.states && v.states.marrow4 ? v.states.marrow4 : null;
      const attr = (st && st.attr) || at || null;
      const when = attr ? attr.when : { none: null, a: m.p_a_driven.value, b: m.p_b_driven.value, ab: m.p_ab_driven.value, other: null };
      const who = attr ? attr.who : { clone: m.p_clone.value, reserve: m.p_reserve.value, both: 0, other: null, none: null };
      const iv = (o) => o ? `${U.fmt(o.value, 2)} <span class="muted">[${U.fmt(o.lo, 2)}–${U.fmt(o.hi, 2)}]</span>` : '–';
      const iv1 = (o) => o ? `${U.fmt(o.value, 1)} <span class="muted">[${U.fmt(o.lo, 1)}–${U.fmt(o.hi, 1)}]</span>` : '–';
      // mean time / divisions since MPP entry of the circulating cells (model moments)
      const sm = attr && attr.sinceMpp ? attr.sinceMpp : null;
      const ageTab = sm ? `<div class="ww-age"><div class="ww-label">${t('twin.age_title', { day: U.dayLabel(sm.day) })} ${GL.help.btn('age_since_mpp', 'help-xs')}</div>
        <table class="table compact"><thead><tr><th>${t('twin.age_col_line')}</th><th class="num">${t('twin.age_col_now')}</th><th class="num">${t('twin.age_col_pre')}</th><th class="num">${t('twin.gen_col_now')}</th><th class="num">${t('twin.gen_col_pre')}</th></tr></thead>
        <tbody>${sm.lines.map(l => `<tr><td>${t(S.variables[l.id].labelKey)}</td><td class="num">${iv1(l.age)}</td><td class="num">${U.fmt(l.age0, 1)}</td><td class="num">${iv1(l.gen)}</td><td class="num">${U.fmt(l.gen0, 1)}</td></tr>`).join('')}</tbody></table>
        <div class="muted small">${t('twin.age_note')}</div></div>` : '';
      const leg = (obj, kind) => Object.keys(obj).filter(k => obj[k] !== null && obj[k] !== undefined).map(k => `<span><i class="sw ww-${kind}-${k}"></i>${t(kind + '.' + k)} <b>${U.pct(obj[k])}</b></span>`).join('');
      const pop = root.GL_AIM4_POP || {};
      const re = p.aim4 && p.aim4.r_engine ? p.aim4.r_engine : null;
      const reRow = re ? `<div class="ww-r muted small"><b>${t('twin.ww_r_engine')}</b> ${GL.help.btn('r_engine', 'help-xs')} · ${t('when.a')} ${U.pct(re.p_when_a)} · ${t('when.b')} ${U.pct(re.p_when_b)} · ${t('when.ab')} ${U.pct((re.p_when_ab_joint || 0) + (re.p_when_ab_redundant || 0))} · ${t('who.clone')} ${U.pct((re.p_who_clone || 0) + (re.p_who_both || 0))} · ${t('who.reserve')} ${U.pct((re.p_who_reserve || 0) + (re.p_who_both || 0))}</div>` : '';
      return `<div class="card ww-card"><div class="eyebrow">${t('eyebrow.ww')}</div>
        <div class="card-head"><h3>${t('twin.ww_title')} ${GL.help.btn('when_why')}</h3><span class="muted small">${t('twin.ww_meta', { model: pop.model || '–', src: t(pop.source === 'fitted' ? 'integ.pop_fitted_short' : 'integ.pop_prior_short'), day: U.dayLabel(A.currentDay()) })}</span></div>
        <div class="ww-grid">
          <div class="ww-col"><div class="ww-label">${t('twin.ww_when')} ${GL.help.btn('ww_when', 'help-xs')}</div><div class="ww-bar" data-kind="when" data-p='${JSON.stringify(when)}'></div><div class="ph-legend ww-legend">${leg(when, 'when')}</div></div>
          <div class="ww-col"><div class="ww-label">${t('twin.ww_who')} ${GL.help.btn('ww_who', 'help-xs')}</div><div class="ww-bar" data-kind="who" data-p='${JSON.stringify(who)}'></div><div class="ph-legend ww-legend">${leg(who, 'who')}</div></div>
          <div class="ww-col ww-nums"><div class="ww-num"><span>${t('out.f_a')} ${GL.help.btn('out.f_a', 'help-xs')}</span><b>${iv(m.f_a)}</b></div><div class="ww-num"><span>${t('out.f_b')} ${GL.help.btn('out.f_b', 'help-xs')}</span><b>${iv(m.f_b)}</b></div><div class="ww-num"><span>${t('out.f_clone30')} ${GL.help.btn('out.f_clone30', 'help-xs')}</span><b>${iv(m.f_clone30)}</b></div></div>
        </div>${reRow}${ageTab}
        <div class="muted small" style="margin-top:6px">${t('twin.ww_note')}</div></div>`;
    },
    drawShareBar(container, probs, kind) {
      container.innerHTML = '';
      const W = 300, H = 20, ns = 'http://www.w3.org/2000/svg';
      const svg = document.createElementNS(ns, 'svg'); svg.setAttribute('viewBox', `0 0 ${W} ${H}`); svg.setAttribute('preserveAspectRatio', 'none'); svg.setAttribute('class', 'gl-chart gl-pheno');
      let x = 0;
      for (const k of Object.keys(probs)) {
        const w = W * (probs[k] || 0); if (!(w > 0)) continue;
        const r = document.createElementNS(ns, 'rect'); r.setAttribute('x', x); r.setAttribute('y', 0); r.setAttribute('width', w); r.setAttribute('height', H); r.setAttribute('class', `ww-${kind}-${k}`);
        const ti = document.createElementNS(ns, 'title'); ti.textContent = t(kind + '.' + k) + ' ' + U.pct(probs[k]); r.appendChild(ti); svg.appendChild(r);
        if (w > 36) { const tx = document.createElementNS(ns, 'text'); tx.setAttribute('x', x + w / 2); tx.setAttribute('y', H / 2 + 4); tx.setAttribute('text-anchor', 'middle'); tx.setAttribute('class', 'gl-ph-text'); tx.textContent = U.pct(probs[k]); svg.appendChild(tx); }
        x += w;
      }
      container.appendChild(svg);
    },
    changesBanner(ch) {
      const ph = ch.phen; const dirs = { p_persistent: 'lower', p_delayed: 'lower', p_second_wave: 'lower', p_quick: 'higher', p_crs3: 'lower', p_response: 'higher' };
      const row = (k, lbl) => ph[k] && ph[k].to !== null ? `<span class="chg"><span class="muted">${lbl}</span> ${U.pct(ph[k].from)} → <b>${U.pct(ph[k].to)}</b> ${A.h.delta(ph[k].from, ph[k].to, true, dirs[k])}</span>` : '';
      const skill = ch.skill ? Object.entries(ch.skill).map(([k, s]) => `${t(S.variables[k].labelKey)}: ${t('hist.coverage')} ${U.pct(s.coverage90)}`).join(' · ') : '';
      const topChanged = ch.topPrev && ch.topCur && ch.topPrev.id !== ch.topCur.id;
      return `<div class="banner"><div class="banner-head"><b>${t('twin.changes_title', { a: 'v' + ch.prev.version, b: 'v' + ch.cur.version })}</b> ${GL.help.btn('changes', 'help-xs')} <span class="muted">· ${t('twin.new_obs', { n: ch.newObs })}</span><div class="spacer"></div><button class="btn btn-sm" id="toHistory">${t('twin.full_diff')}</button><button class="btn btn-ghost btn-sm" id="dismissChanges">✕</button></div>
        <div class="banner-row">${row('p_persistent', t('pheno.persistent'))}${row('p_delayed', t('pheno.delayed'))}${row('p_second_wave', t('pheno.second_wave'))}${row('p_quick', t('pheno.quick'))}${row('p_crs3', t('out.p_crs3'))}${row('p_response', t('out.p_response'))}</div>
        <div class="banner-row"><span class="chg"><span class="muted">${t('twin.recommended')}</span> ${A.planName(ch.topPrev.plan)} → <b>${A.planName(ch.topCur.plan)}</b> ${topChanged ? `<span class="delta neutral">${t('twin.changed')}</span>` : `<span class="delta flat">${t('twin.unchanged')}</span>`}</span>${skill ? `<span class="chg"><span class="muted">${t('twin.prev_skill')}</span> ${skill}</span>` : ''}</div></div>`;
    }
  };
})(typeof window !== 'undefined' ? window : global);
