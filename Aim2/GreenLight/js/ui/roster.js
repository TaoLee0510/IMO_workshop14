/* GreenLight · ui/roster.js — patient roster: patients in treatment vs completed records. */
(function (root) {
  const GL = root.GL = root.GL || {};
  const U = GL.util, ST = GL.store, A = GL.app;
  const t = (k, p) => GL.i18n.t(k, p);
  const V = A.views.roster = {
    ui: { tab: 'active' },
    async render(el) {
      const rows = A.state.patients;
      for (const r of rows) await ST.loadPatient(r.patient_id);
      const active = rows.filter(r => (r.demo.status || 'active') !== 'completed'), completed = rows.filter(r => (r.demo.status || 'active') === 'completed');
      const q = (A.state.rosterQuery || '').toLowerCase();
      const fitted = rows.filter(r => ST.latest(r.patient_id)).length;
      const highRisk = active.filter(r => { const v = ST.latest(r.patient_id); return v && v.plans.soc && v.plans.soc.outcomes.marrow && v.plans.soc.outcomes.marrow.p_persistent.value >= 0.4; }).length;
      const tab = V.ui.tab;
      el.innerHTML = `
        <div class="hero"><div class="hero-text"><div class="hero-kicker">${t('brand.org')}</div><h1>${t('roster.title')}</h1><p>${t('roster.subtitle')}</p><p class="hero-note">${t('roster.hero_note')}</p>
            <div class="hero-stats"><div><b>${active.length}</b><span>${t('roster.n_active')}</span></div><div><b>${completed.length}</b><span>${t('roster.n_completed')}</span></div><div><b>${fitted}</b><span>${t('roster.n_twins')}</span></div><div><b class="${highRisk ? 'warn' : ''}">${highRisk}</b><span>${t('roster.n_review')}</span></div></div></div>
          <div class="hero-art">${GL.figure.svg({ id: 'glf-hero', variant: 'pair', width: 300, height: 250 })}</div></div>
        <div class="page-head"><div class="tabs tabs-lg" id="rosterTabs"><button class="tab ${tab === 'active' ? 'active' : ''}" data-tab="active">${t('roster.tab_active', { n: active.length })}</button><button class="tab ${tab === 'completed' ? 'active' : ''}" data-tab="completed">${t('roster.tab_completed', { n: completed.length })}</button></div>
          <div class="page-actions"><input id="rosterSearch" class="search" placeholder="${t('roster.search')}" value="${U.escapeHtml(A.state.rosterQuery || '')}"><button class="btn btn-ghost" id="btnReset">${t('ui.reset_demo')}</button></div></div>
        ${tab === 'active' ? V.activeSection(active, q, highRisk) : V.completedSection(completed, q)}`;
      el.querySelector('#rosterSearch').addEventListener('input', (e) => { A.state.rosterQuery = e.target.value; A.render(); });
      el.querySelector('#btnReset').addEventListener('click', () => A.resetDemo());
      el.querySelector('#rosterTabs').addEventListener('click', (e) => { const b = e.target.closest('[data-tab]'); if (b) { V.ui.tab = b.dataset.tab; A.render(); } });
      el.querySelectorAll('[data-open]').forEach(b => b.addEventListener('click', () => A.go('twin', b.dataset.open)));
      el.querySelectorAll('.ph-mini').forEach(c => GL.charts.phenotypeBar(c, JSON.parse(c.dataset.p), { height: 14 }));
    },
    activeSection(list, q, highRisk) {
      return `
        <div class="section-head"><div class="eyebrow">${t('status.active')}</div><h2>${t('roster.active_title')}</h2><p class="muted">${t('roster.active_note')}</p></div>
        <div class="tiles">
          <div class="tile"><div class="tile-label">${t('roster.n_active')}</div><div class="tile-value">${list.length}</div></div>
          <div class="tile"><div class="tile-label">${t('roster.n_twins')} ${GL.help.btn('models_built')}</div><div class="tile-value">${list.filter(r => ST.latest(r.patient_id)).length}</div></div>
          <div class="tile"><div class="tile-label">${t('roster.n_review')} ${GL.help.btn('review')}</div><div class="tile-value ${highRisk ? 'warn' : ''}">${highRisk}</div><div class="tile-sub">${t('roster.n_review_sub')}</div></div>
          <div class="tile"><div class="tile-label">${t('roster.source')}</div><div class="tile-value small">${t('adapter.' + A.state.adapter.id)}</div><div class="tile-sub">${U.escapeHtml((A.state.adapterInfo && A.state.adapterInfo.source) || '')}</div></div>
        </div>
        <div class="card"><table class="table" id="rosterTable"><thead><tr>
          <th>${t('ui.patient')}</th><th>${t('roster.demo')}</th><th>${t('ui.product')}</th><th>${t('roster.phase')}</th><th>${t('tl.horizon')}</th><th>${t('roster.phenotype')} ${GL.help.btn('phenotype')}</th><th>${t('out.p_response')} <span class="unit">(%)</span> ${GL.help.btn('out.p_response')}</th><th>${t('roster.top_plan')}</th><th>${t('ui.twin_version')}</th><th></th>
        </tr></thead><tbody>${list.filter(r => !q || r.patient_id.toLowerCase().includes(q)).map(r => V.activeRow(r)).join('')}</tbody></table></div>`;
    },
    completedSection(list, q) {
      const ph = { quick: 0, second_wave: 0, persistent: 0 }; let crs3 = 0, relapse = 0, ch = 0;
      list.forEach(r => { const d = A.state.adapter.data ? A.state.adapter.data.patients.find(x => x.patient_id === r.patient_id) : null; if (!d) return; if (d.recorded && d.recorded.phenotype_clinical && ph[d.recorded.phenotype_clinical] !== undefined) ph[d.recorded.phenotype_clinical]++; if (d.crs && d.crs.grade >= 3) crs3++; if (d.recorded && d.recorded.relapse_day) relapse++; if (d.mutations && d.mutations.variants.length) ch++; });
      const n = list.length || 1;
      return `
        <div class="section-head"><div class="eyebrow">${t('status.completed')}</div><h2>${t('roster.cohort_title')}</h2><p class="muted">${t('roster.cohort_note')}</p></div>
        <div class="tiles">
          <div class="tile tile-pheno"><div class="tile-label">${t('roster.recorded_phenotype')} ${GL.help.btn('phenotype')}</div><div class="ph-big ph-mini" data-p='${JSON.stringify({ quick: ph.quick / n, second_wave: ph.second_wave / n, persistent: ph.persistent / n })}'></div><div class="ph-legend"><span><i class="sw sw-quick"></i>${t('pheno.quick')} ${ph.quick}</span><span><i class="sw sw-second_wave"></i>${t('pheno.second_wave')} ${ph.second_wave}</span><span><i class="sw sw-persistent"></i>${t('pheno.persistent')} ${ph.persistent}</span></div></div>
          <div class="tile"><div class="tile-label">${t('out.p_crs3')}</div><div class="tile-value">${U.pct(crs3 / n)}</div><div class="tile-sub">${crs3} / ${list.length}</div></div>
          <div class="tile"><div class="tile-label">${t('roster.relapse')}</div><div class="tile-value">${U.pct(relapse / n)}</div><div class="tile-sub">${relapse} / ${list.length}</div></div>
          <div class="tile"><div class="tile-label">${t('roster.ch_positive')}</div><div class="tile-value">${ch}</div><div class="tile-sub">${t('roster.ch_sub')}</div></div>
        </div>
        <div class="card"><table class="table" id="rosterTable"><thead><tr>
          <th>${t('ui.patient')}</th><th>${t('roster.demo')}</th><th>${t('ui.product')}</th><th>${t('roster.recorded_phenotype')}</th><th>${t('roster.crs')}</th><th>${t('roster.relapse')}</th><th>${t('roster.followup')}</th><th>${t('roster.model_fit')} ${GL.help.btn('model_fit')}</th><th>${t('roster.ch')}</th><th></th>
        </tr></thead><tbody>${list.filter(r => !q || r.patient_id.toLowerCase().includes(q)).map(r => V.completedRow(r)).join('')}</tbody></table></div>`;
    },
    activeRow(r) {
      const v = ST.latest(r.patient_id); const st = ST.patientState[r.patient_id];
      const day = st ? r.demo.checkpoints[st.horizonIndex !== undefined ? st.horizonIndex : st.checkpointIndex] : r.demo.checkpoints[r.demo.initialCheckpoint];
      const phase = day <= -5 ? t('phase.pre') : day <= -3 ? t('phase.ld') : day < 0 ? t('phase.post_ld') : day <= 28 ? t('phase.acute') : day <= 90 ? t('phase.early') : t('phase.late');
      const m = v && v.plans.soc ? v.plans.soc.outcomes.marrow : null; const ti = v && v.plans.soc ? v.plans.soc.outcomes.tumor_immune : null;
      const pheno = m ? `<div class="ph-mini" data-p='${JSON.stringify({ quick: m.p_quick.value, second_wave: m.p_second_wave.value, persistent: m.p_persistent.value })}'></div>` : `<span class="muted">${t('roster.not_fitted')}</span>`;
      const top = v ? A.planName(A.rankedPlans(v, (st && st.weights) || undefined)[0].plan) : '–';
      const light = m ? (m.p_persistent.value >= 0.4 ? 'red' : (m.p_second_wave.value + m.p_persistent.value) >= 0.4 ? 'amber' : 'green') : 'grey';
      return `<tr><td><span class="light ${light}"></span><b>${r.label || r.patient_id}</b></td><td>${r.demographics.age} · ${t('sex.' + r.demographics.sex)}</td><td>${r.product}</td><td>${phase}</td><td><b>${U.dayLabel(day)}</b></td><td class="w-180">${pheno}</td><td>${ti ? U.pct(ti.p_response.value) : '–'}</td><td>${top}</td><td>${v ? 'v' + v.version : '–'}</td><td><button class="btn btn-sm" data-open="${r.patient_id}">${t('roster.open')}</button></td></tr>`;
    },
    completedRow(r) {
      const d = A.state.adapter.data ? A.state.adapter.data.patients.find(x => x.patient_id === r.patient_id) : null;
      const rec = (d && d.recorded) || {}; const v = ST.latest(r.patient_id);
      const m = v && v.plans.soc ? v.plans.soc.outcomes.marrow : null;
      const fit = m ? `<div class="ph-mini" data-p='${JSON.stringify({ quick: m.p_quick.value, second_wave: m.p_second_wave.value, persistent: m.p_persistent.value })}'></div>` : `<span class="muted">${t('roster.not_fitted')}</span>`;
      const ph = rec.phenotype_clinical ? `<span class="pheno-tag ph-${rec.phenotype_clinical}">${t('pheno.' + rec.phenotype_clinical)}</span>` : '–';
      const chv = d && d.mutations && d.mutations.variants.length ? d.mutations.variants.map(x => x.gene).join(', ') : (d && d.mutations && d.mutations.tested ? t('ui.ch_negative') : '–');
      return `<tr><td><span class="light grey"></span><b>${r.label || r.patient_id}</b></td><td>${r.demographics.age} · ${t('sex.' + r.demographics.sex)}</td><td>${r.product}</td><td>${ph}</td><td>${d && d.crs ? d.crs.grade : '0'}</td><td>${rec.relapse_day ? U.dayLabel(rec.relapse_day) : '—'}</td><td>${rec.last_followup_day ? rec.last_followup_day + ' d' : '–'}</td><td class="w-180">${fit}</td><td class="small">${chv}</td><td><button class="btn btn-sm" data-open="${r.patient_id}">${t('roster.open')}</button></td></tr>`;
    }
  };
})(typeof window !== 'undefined' ? window : global);
