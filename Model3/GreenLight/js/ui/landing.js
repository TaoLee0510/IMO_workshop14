/* GreenLight · ui/landing.js
 * Home page: the GreenLight logo, the digital-twin idea in one picture, and the two entries
 * (clinician portal, researcher portal).
 */
(function (root) {
  const GL = root.GL = root.GL || {};
  const t = (k, p) => GL.i18n.t(k, p), esc = (s) => GL.util.escapeHtml(s);

  // patient <-> digital twin: data flow in, forecasts and decisions flow out
  function concept() {
    const person = (x, fill, stroke, dash) => `<g transform="translate(${x},0)"><circle cx="0" cy="58" r="20" fill="${fill}" stroke="${stroke}" stroke-width="2.4" ${dash ? 'stroke-dasharray="4 3"' : ''}/><path d="M-34 150c0-30 15-52 34-52s34 22 34 52z" fill="${fill}" stroke="${stroke}" stroke-width="2.4" ${dash ? 'stroke-dasharray="4 3"' : ''}/></g>`;
    let nodes = '';
    const pts = [[0, 58], [-14, 112], [14, 112], [0, 132], [-22, 140], [22, 140], [0, 92]];
    pts.forEach(([x, y]) => { nodes += `<circle cx="${640 + x}" cy="${y}" r="3.2" fill="var(--twin)"/>`; });
    for (let i = 0; i < pts.length; i++) for (let j = i + 1; j < pts.length; j++) if (Math.hypot(pts[i][0] - pts[j][0], pts[i][1] - pts[j][1]) < 50) nodes += `<line x1="${640 + pts[i][0]}" y1="${pts[i][1]}" x2="${640 + pts[j][0]}" y2="${pts[j][1]}" stroke="var(--twin)" stroke-width="1" opacity=".6"/>`;
    const dots = [0, 1, 2, 3, 4, 5].map(i => `<circle r="4" cx="300" cy="78" fill="var(--l-anc)"><animate attributeName="cx" from="300" to="560" dur="3s" begin="${i * 0.5}s" repeatCount="indefinite"/><animate attributeName="cy" values="78;70;78" dur="3s" begin="${i * 0.5}s" repeatCount="indefinite"/></circle>`).join('');
    const back = [0, 1, 2].map(i => `<circle r="4" cx="560" cy="126" fill="var(--brand-600)"><animate attributeName="cx" from="560" to="300" dur="3s" begin="${i * 1}s" repeatCount="indefinite"/><animate attributeName="cy" values="126;134;126" dur="3s" begin="${i * 1}s" repeatCount="indefinite"/></circle>`).join('');
    return `<svg viewBox="0 0 860 180" role="img" aria-label="${esc(t('land.concept_alt'))}">
      ${person(220, 'var(--surface)', 'var(--brand-600)', false)}
      ${person(640, 'var(--twin-100)', 'var(--twin)', true)}${nodes}
      <path d="M300 74 C 400 60, 460 60, 560 74" fill="none" stroke="var(--l-anc)" stroke-width="2" stroke-dasharray="2 5"/>${dots}
      <path d="M560 130 C 460 144, 400 144, 300 130" fill="none" stroke="var(--brand-600)" stroke-width="2" stroke-dasharray="2 5"/>${back}
      <text x="430" y="50" text-anchor="middle" style="font-size:12.5px;font-weight:700;fill:var(--l-anc)">${esc(t('land.flow_in'))}</text>
      <text x="430" y="168" text-anchor="middle" style="font-size:12.5px;font-weight:700;fill:var(--brand-600)">${esc(t('land.flow_out'))}</text>
      <text x="220" y="176" text-anchor="middle" style="font-size:12px;font-weight:700;fill:var(--text)">${esc(t('land.patient'))}</text>
      <text x="640" y="176" text-anchor="middle" style="font-size:12px;font-weight:700;fill:var(--twin)">${esc(t('land.twin'))}</text>
    </svg>`;
  }

  GL.views.landing = {
    async render(el) {
      const idx = GL.loader.index(), S = GL.session;
      const r = idx ? idx.roster : [];
      const n = (f) => r.filter(f).length;
      const nReal = idx && idx.real_cohort ? idx.real_cohort.n : 0;
      const scope = (k) => `<div class="tiny scope">${esc(t(k))}</div>`;
      const clinicScope = S.can('clinic.real') ? scope('land.scope.assigned') : scope('land.scope.digital');
      const researchScope = S.role === 'researcher' ? scope('land.scope.research') : S.role === 'clinician' ? scope('land.scope.model') : scope('land.scope.digital');
      el.innerHTML = `<section class="landing">
        <div class="hero">
          <h1 class="hero-brand">${GL.app.brand('bi-hero')}</h1>
          <div class="claim">${esc(t('land.claim'))}</div>
          <div class="lead">${esc(t('land.lead'))}</div>
        </div>
        <div class="entries" style="${S.can('clinic') && S.can('research') ? '' : 'grid-template-columns:minmax(280px,380px)'}">
          ${S.can('clinic') ? `<a class="entry entry-clinic" href="#/clinic" id="enterClinic">
            <div class="ehead"><span class="eic">${GL.app.icon('steth')}</span><div><h2>${esc(t('land.clinic.t'))}</h2><div class="small muted">${esc(t('land.clinic.s'))}</div></div></div>
            <ul><li>${esc(t('land.clinic.1'))}</li><li>${esc(t('land.clinic.2'))}</li><li>${esc(t('land.clinic.3'))}</li></ul>
            ${clinicScope}<span class="go">${esc(t('land.enter'))} ${GL.app.icon('arrow')}</span>
          </a>` : ''}
          ${S.can('research') ? `<a class="entry entry-research" href="#/research/model" id="enterResearch">
            <div class="ehead"><span class="eic">${GL.app.icon('flask')}</span><div><h2>${esc(t('land.research.t'))}</h2><div class="small muted">${esc(t('land.research.s'))}</div></div></div>
            <ul><li>${esc(t('land.research.1'))}</li><li>${esc(t('land.research.2'))}</li><li>${esc(t('land.research.3'))}</li></ul>
            ${researchScope}<span class="go">${esc(t('land.enter'))} ${GL.app.icon('arrow')}</span>
          </a>` : ''}
        </div>
        <div class="concept">${concept()}</div>
        <div class="land-foot">
          <span class="chip">${GL.app.typeMark({ type: 'real', status: 'completed' })} <b>${nReal}</b> ${esc(t('land.n_real'))} · ${GL.app.icon('lock')} ${esc(t('land.sealed'))}</span>
          <span class="chip">${GL.app.typeMark({ type: 'digital', status: 'active' })} <b>${n(p => p.status === 'active')}</b> ${esc(t('land.n_active'))}</span>
          <span class="chip">${GL.app.typeMark({ type: 'digital', status: 'completed' })} <b>${n(p => p.type === 'digital')}</b> ${esc(t('land.n_digital'))}</span>
          <span class="chip chip-twin">${esc(GL.loader.modelLabel())} · ${esc(t('pop.label'))}: ${esc(t('pop.' + (idx ? idx.population.source : 'prior')))}</span>
        </div>
        <div class="disclaimer">${esc(t('ui.disclaimer'))}</div>
      </section>`;
    }
  };
})(typeof window !== 'undefined' ? window : global);
