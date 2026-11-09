/* GreenLight · ui/integration.js — data sources, model modules, contracts. */
(function (root) {
  const GL = root.GL = root.GL || {};
  const U = GL.util, S = GL.schema, A = GL.app, R = GL.registry;
  const t = (k, p) => GL.i18n.t(k, p);
  const V = A.views.integration = {
    ui: { csvUrl: '../TestData/', restUrl: 'http://localhost:8787/api', restToken: '', modUrl: 'http://localhost:8788/module', modToken: '', status: {} },
    async render(el) {
      const cur = A.state.adapter.id, p = A.state.patient, v = A.currentVersion();
      const snap = p && v ? S.snapshot(p, v.asOfDay) : null;
      const mods = R.modules.map(mod => {
        const cov = snap ? R.coverage(mod, snap) : null;
        return `<tr class="${R.enabled[mod.meta.id] ? '' : 'row-off'}"><td><b>${t(mod.meta.nameKey)}</b>${mod.meta.placeholder ? ` <span class="badge badge-warn">${t('ui.placeholder')}</span>` : ''}<div class="muted small">${mod.meta.owner} · ${mod.meta.version}</div></td>
          <td><span class="badge">${mod.meta.kind}</span><div class="muted small">${U.escapeHtml(mod.meta.engine || '')}</div></td>
          <td>${(mod.consumes.variables || []).map(x => `<span class="chip chip-sm ${cov && !cov[x] ? 'chip-off' : ''}">${t(S.variables[x].labelKey)}${cov ? ' ' + cov[x] : ''}</span>`).join('')}<div class="muted small">${(mod.consumes.signals || []).map(s => '⇠ ' + s).join(' ')}</div></td>
          <td>${(mod.provides.signals || []).map(s => `<span class="chip chip-sm">⇢ ${s}</span>`).join('')}<div class="muted small">${mod.provides.outcomes.length} ${t('integ.outcomes')} · ${mod.provides.params.length} ${t('integ.params')}</div></td>
          <td>${v && v.modules[mod.meta.id] ? `${t('twin.fit_ms', { ms: v.modules[mod.meta.id].fitMs })}<div class="muted small">ESS ${U.fmt(v.modules[mod.meta.id].ess, 0)}</div>` : '–'}</td>
          <td><label class="switch"><input type="checkbox" data-mod="${mod.meta.id}" ${R.enabled[mod.meta.id] ? 'checked' : ''}><span></span></label></td></tr>`;
      }).join('');
      const order = R.order().map(m => t(m.meta.nameKey)).join('  →  ');
      const st = V.ui.status;
      el.innerHTML = `
        <div class="page-head"><div><h1>${t('integ.title')}</h1><p class="muted">${t('integ.subtitle')}</p></div></div>
        <div class="grid-3">
          <div class="card ${cur === 'demo' ? 'is-active' : ''}"><div class="card-head"><h3>${t('adapter.demo')}</h3>${cur === 'demo' ? `<span class="badge badge-ok">${t('integ.connected')}</span>` : ''}</div><p class="muted small">${t('integ.demo_desc')}</p><button class="btn btn-sm" id="useDemo" ${cur === 'demo' ? 'disabled' : ''}>${t('integ.use')}</button></div>
          <div class="card ${cur === 'csv' ? 'is-active' : ''}"><div class="card-head"><h3>${t('adapter.csv')}</h3>${cur === 'csv' ? `<span class="badge badge-ok">${t('integ.connected')}</span>` : ''}</div><p class="muted small">${t('integ.csv_desc')}</p><label class="fld"><span>${t('integ.base_url')}</span><input class="input" id="csvUrl" value="${U.escapeHtml(V.ui.csvUrl)}"></label><button class="btn btn-sm" id="useCsv">${t('integ.connect')}</button> <span class="muted small">${st.csv || ''}</span></div>
          <div class="card ${cur === 'rest' ? 'is-active' : ''}"><div class="card-head"><h3>${t('adapter.rest')}</h3>${cur === 'rest' ? `<span class="badge badge-ok">${t('integ.connected')}</span>` : ''}</div><p class="muted small">${t('integ.rest_desc')}</p><label class="fld"><span>${t('integ.base_url')}</span><input class="input" id="restUrl" value="${U.escapeHtml(V.ui.restUrl)}"></label><label class="fld"><span>${t('integ.token')}</span><input class="input" id="restToken" type="password" value="${U.escapeHtml(V.ui.restToken)}"></label><button class="btn btn-sm" id="useRest">${t('integ.connect')}</button> <span class="muted small">${st.rest || ''}</span></div>
        </div>
        <div class="card"><div class="card-head"><h3>${t('integ.modules')}</h3><span class="muted small">${t('integ.run_order')}: ${order}</span></div>
          <table class="table"><thead><tr><th>${t('integ.module')}</th><th>${t('integ.kind')}</th><th>${t('integ.consumes')}</th><th>${t('integ.provides')}</th><th>${t('integ.last_fit')}</th><th>${t('integ.enabled')}</th></tr></thead><tbody>${mods}</tbody></table>
          <div class="add-remote"><h4>${t('integ.add_remote')}</h4><p class="muted small">${t('integ.add_remote_desc')}</p>
            <div class="fields-row"><label class="fld"><span>${t('integ.base_url')}</span><input class="input" id="modUrl" value="${U.escapeHtml(V.ui.modUrl)}"></label><label class="fld"><span>${t('integ.token')}</span><input class="input" id="modToken" type="password" value="${U.escapeHtml(V.ui.modToken)}"></label><button class="btn btn-sm" id="addMod">${t('integ.register')}</button><span class="muted small">${st.mod || ''}</span></div></div></div>
        <div class="grid-1-1">
          <div class="card"><div class="card-head"><h3>${t('integ.module_contract')}</h3><span class="muted small">docs/MODULE_API.md</span></div><pre class="code">${U.escapeHtml(V.moduleContract())}</pre></div>
          <div class="card"><div class="card-head"><h3>${t('integ.data_contract')}</h3><span class="muted small">docs/DATA_API.md</span></div><pre class="code">${U.escapeHtml(V.dataContract())}</pre></div>
        </div>`;
      el.querySelector('#useDemo').addEventListener('click', () => A.setAdapter(GL.adapters.DemoDataAdapter));
      el.querySelector('#useCsv').addEventListener('click', async () => { V.ui.csvUrl = el.querySelector('#csvUrl').value; await A.setAdapter(GL.adapters.CsvDataAdapter, [V.ui.csvUrl]); });
      el.querySelector('#useRest').addEventListener('click', async () => { V.ui.restUrl = el.querySelector('#restUrl').value; V.ui.restToken = el.querySelector('#restToken').value; await A.setAdapter(GL.adapters.RestDataAdapter, [V.ui.restUrl, V.ui.restToken]); });
      el.querySelector('#addMod').addEventListener('click', async () => {
        V.ui.modUrl = el.querySelector('#modUrl').value; V.ui.modToken = el.querySelector('#modToken').value; V.ui.status.mod = t('ui.connecting'); A.render();
        try { const m = await GL.modules.makeRemote({ baseUrl: V.ui.modUrl, token: V.ui.modToken, id: 'remote_' + Date.now().toString(36) }).load(); R.register(m); V.ui.status.mod = t('integ.registered', { name: m.meta.nameKey }); A.toast(V.ui.status.mod); }
        catch (e) { V.ui.status.mod = t('ui.connect_failed') + ': ' + e.message; }
        A.render();
      });
      el.querySelectorAll('[data-mod]').forEach(sw => sw.addEventListener('change', () => { R.setEnabled(sw.dataset.mod, sw.checked); A.toast(t('integ.toggle_note')); A.render(); }));
    },
    moduleContract() {
      return `GET  /meta      -> { meta:{id,nameKey,owner,version,kind}, consumes:{variables,planKeys,signals}, provides:{signals,params,outcomes,variables} }
POST /fit       { snapshot, ctx:{asOfDay,horizon,nParticles,seed,signals} }
                -> state  { moduleId, asOfDay, nObs, ess, params:{id:{median,lo,hi}}, ...opaque }
POST /forecast  { state, plan, ctx:{asOfDay,horizon,signals,snapshot} }
                -> { grid:[day], variables:{v:{q05,q25,q50,q75,q95}},
                     outcomes:{id:{value,lo?,hi?}}, signals:{id:{grid,values}} }
plan = { id, stages:{ bridging:{type}, ld:{flu_frac,cy_frac}, infusion:{car_dose_frac},
         crs:{toci_grade,steroid_grade,anakinra_prophylaxis},
         cytopenia:{gcsf,gcsf_start_day,gcsf_anc_thr,cd34_backup,cd34_boost,cd34_boost_day,tpo_ra,plt_tx_thr},
         followup:{cbc_interval_days} } }`;
    },
    dataContract() {
      return `patient = { patient_id, demographics:{age,sex}, product, costim,
  history:{prior_lines,prior_asct,prior_bendamustine,months_since_last_chemo,marrow_involvement},
  bridging, ld:{flu_total,cy_total,dose_reduced}, car_dose_frac,
  crs:{grade,onset_day,end_day}|null,
  mutations:{tested,sample_day,variants:[{gene,protein_change,vaf,depth}]},
  observations:[{day,variable,value,unit?,source?}],   // long format, any variable
  events:[{type,day_start,day_end,dose}],
  demo:{checkpoints:[...],initialCheckpoint} }
variables: ${Object.keys(S.variables).join(', ')}
events:    ${S.eventTypes.join(', ')}
REST: GET /patients · GET /patients/{id} · GET /patients/{id}/observations?upto=D
      GET /patients/{id}/events?upto=D · GET /patients/{id}/mutations
CSV:  patients.csv, labs_long.csv, inflam_long.csv, events.csv, mutation_tested.csv, mutations_pre.csv`;
    }
  };
})(typeof window !== 'undefined' ? window : global);
