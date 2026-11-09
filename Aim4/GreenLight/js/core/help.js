/* GreenLight · core/help.js
 * "?" help buttons and a single popover. Any element with data-help="<id>"
 * opens an explanation resolved from the module/plan definitions plus the
 * help dictionary (js/i18n/help_dict.js): meaning, unit, preferred direction,
 * reference values.
 */
(function (root) {
  const GL = root.GL = root.GL || {};
  const H = GL.help = {};
  const t = (k, p) => GL.i18n.t(k, p);
  const esc = (s) => GL.util.escapeHtml(s);

  H.btn = (id, cls) => `<button type="button" class="help-btn ${cls || ''}" data-help="${id}" aria-label="${esc(t('help.open'))}" title="${esc(t('help.open'))}">?</button>`;
  // SVG variant (for charts): small circle with a question mark
  H.svgBtn = (id, x, y) => `<g class="help-svg" data-help="${id}" transform="translate(${x},${y})"><circle r="6.5"/><text y="3.2" text-anchor="middle">?</text></g>`;

  function findOutcome(id) { for (const m of GL.registry.modules) { const d = (m.provides.outcomes || []).find(o => o.id === id); if (d) return { def: d, mod: m }; } return null; }
  function findParam(id) { for (const m of GL.registry.modules) { const d = (m.provides.params || []).find(o => o.id === id); if (d) return { def: d, mod: m }; } return null; }

  H.resolve = function (id) {
    const has = (k) => GL.i18n.has(k);
    const body = has('help.' + id) ? t('help.' + id) : '';
    if (id.startsWith('out.')) {
      const f = findOutcome(id.slice(4)); if (!f) return null;
      const unit = f.def.kind === 'prob' ? '%' : (f.def.unit || '');
      return { title: t(f.def.labelKey), unit, direction: f.def.direction, body, source: t(f.mod.meta.nameKey) + (f.mod.meta.placeholder ? ' · ' + t('ui.placeholder') : '') };
    }
    if (id.startsWith('param.')) {
      const f = findParam(id.slice(6)); if (!f) return null;
      return { title: t(f.def.labelKey), unit: f.def.unit || '', direction: f.def.direction, body: body || (has(f.def.descKey) ? t(f.def.descKey) : ''), source: t(f.mod.meta.nameKey) };
    }
    if (id.startsWith('field.')) {
      const [, stage, key] = id.split('.'); const f = GL.plans.fields.find(x => x.stage === stage && x.key === key); if (!f) return null;
      return { title: t(f.labelKey), unit: f.unit || '', direction: null, body, source: t('stage.' + stage) };
    }
    if (id.startsWith('plan.')) return { title: t(id), unit: '', direction: null, body, source: t('twin.top3') };
    if (id.startsWith('stage.')) return { title: t(id), unit: '', direction: null, body, source: t('planner.title') };
    const title = has('help.' + id + '.t') ? t('help.' + id + '.t') : '';
    const dir = has('help.' + id + '.dir') ? t('help.' + id + '.dir') : null; // 'higher' | 'lower' | 'neutral'
    return { title, unit: has('help.' + id + '.unit') ? t('help.' + id + '.unit') : '', direction: dir, body, source: '' };
  };

  let pop = null, anchor = null;
  function ensurePop() {
    if (pop) return pop;
    pop = document.createElement('div'); pop.className = 'help-pop'; pop.setAttribute('role', 'dialog'); pop.style.display = 'none';
    document.body.appendChild(pop);
    return pop;
  }
  H.close = function () { if (pop) { pop.style.display = 'none'; } if (anchor) { anchor.classList.remove('is-open'); anchor = null; } };
  H.open = function (btn) {
    const id = btn.dataset.help; const r = H.resolve(id);
    const p = ensurePop();
    if (!r) { p.innerHTML = `<div class="help-body muted">${esc(id)}</div>`; }
    else {
      const dirCls = r.direction === 'higher' ? 'good' : r.direction === 'lower' ? 'bad' : 'neutral';
      const dirTxt = r.direction === 'higher' ? t('help.dir.higher') : r.direction === 'lower' ? t('help.dir.lower') : r.direction === 'neutral' ? t('help.dir.neutral') : '';
      p.innerHTML = `<div class="help-head"><b>${esc(r.title)}</b><button type="button" class="help-x" aria-label="${esc(t('help.close'))}">✕</button></div>
        <div class="help-meta">${r.unit ? `<span class="chip chip-sm">${esc(t('help.unit'))}: <b>${esc(r.unit)}</b></span>` : ''}${dirTxt ? `<span class="chip chip-sm dir-${dirCls}">${esc(dirTxt)}</span>` : ''}</div>
        <div class="help-body">${esc(r.body).replace(/\n/g, '<br>')}</div>${r.source ? `<div class="help-src muted small">${esc(r.source)}</div>` : ''}`;
    }
    p.style.display = 'block';
    if (anchor) anchor.classList.remove('is-open');
    anchor = btn; btn.classList.add('is-open');
    H.reposition();
  };
  H.reposition = function () {
    if (!pop || !anchor || pop.style.display === 'none') return;
    const b = anchor.getBoundingClientRect(), w = Math.min(340, window.innerWidth - 24);
    pop.style.width = w + 'px';
    let left = b.left + window.scrollX - 8, top = b.bottom + window.scrollY + 8;
    if (left + w > window.scrollX + window.innerWidth - 12) left = window.scrollX + window.innerWidth - 12 - w;
    if (left < 12) left = 12;
    pop.style.left = left + 'px'; pop.style.top = top + 'px';
    const ph = pop.getBoundingClientRect().height;
    if (b.bottom + 8 + ph > window.innerHeight && b.top - 8 - ph > 0) pop.style.top = (b.top + window.scrollY - 8 - ph) + 'px';
  };
  H.init = function () {
    document.addEventListener('click', (e) => {
      const x = e.target.closest && e.target.closest('.help-x'); if (x) { H.close(); return; }
      const btn = e.target.closest && e.target.closest('[data-help]');
      if (btn) { e.preventDefault(); e.stopPropagation(); if (anchor === btn) H.close(); else H.open(btn); return; }
      if (pop && pop.style.display !== 'none' && !pop.contains(e.target)) H.close();
    }, true);
    document.addEventListener('keydown', (e) => { if (e.key === 'Escape') H.close(); });
    window.addEventListener('resize', H.reposition);
    GL.i18n.onChange(H.close);
  };
})(typeof window !== 'undefined' ? window : global);
