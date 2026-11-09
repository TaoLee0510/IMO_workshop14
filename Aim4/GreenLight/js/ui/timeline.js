/* GreenLight · ui/timeline.js — treatment timeline scrubber.
 * Shows the course (LD window, infusion, interventions), the checkpoints, the data horizon
 * (for patients in treatment) and the viewing point. Click or drag to view the model as of
 * any checkpoint up to the horizon; keyboard ← → moves one checkpoint. */
(function (root) {
  const GL = root.GL = root.GL || {};
  const U = GL.util, A = GL.app;
  const t = (k, p) => GL.i18n.t(k, p);
  GL.ui = GL.ui || {};
  GL.ui.timeline = function (container, opt) {
    const cps = opt.checkpoints, horizon = opt.horizonIndex, view = opt.viewIndex, p = opt.patient, active = opt.active;
    const x0 = -8, x1 = 94, W = 1000, H = 86, m = { l: 18, r: 18 };
    const X = (d) => m.l + (d - x0) / (x1 - x0) * (W - m.l - m.r);
    const trackY = 44;
    const events = p.events.filter(e => active ? e.day_start <= cps[horizon] : true);
    const evHtml = events.map(e => `<rect class="tl-ev gl-ev-${e.type}" x="${X(e.day_start) - 1.5}" y="${trackY - 20}" width="${Math.max(3, X(Math.max(e.day_end || e.day_start, e.day_start) + 1) - X(e.day_start))}" height="7" rx="1.5"><title>${t('event.' + e.type)} ${U.dayLabel(e.day_start)}</title></rect>`).join('');
    const futureX = X(cps[horizon]);
    const cpHtml = cps.map((d, i) => {
      const dis = i > horizon, sel = i === view;
      const label = d === -5 ? U.dayLabel(d) + ' · ' + t('tl.baseline') : d === 0 ? 'D0 · ' + t('tl.infusion') : U.dayLabel(d);
      const above = d === -5; // avoid overlap between the baseline and infusion labels
      return `<g class="tl-cp ${dis ? 'is-disabled' : ''} ${sel ? 'is-selected' : ''} ${i === horizon && active ? 'is-horizon' : ''}" data-idx="${i}" transform="translate(${X(d)},${trackY})">
        <circle r="${sel ? 9 : 6}" class="tl-cp-dot"/>${sel ? '<circle r="3" class="tl-cp-core"/>' : ''}
        <text y="${above ? -14 : 26}" text-anchor="${above ? 'start' : 'middle'}" x="${above ? -6 : 0}" class="tl-cp-label">${label}</text></g>`;
    }).join('');
    container.innerHTML = `<div class="tl-head"><div class="tl-title">${t('tl.title')} ${GL.help.btn('timeline', 'help-xs')}</div>
      <div class="tl-status">${active ? `<span class="chip chip-sm chip-day">${t('tl.horizon')} <b>${U.dayLabel(cps[horizon])}</b></span>` : `<span class="chip chip-sm">${t(GL.app.statusKey(p))} · ${t('roster.followup')} ${p.recorded && p.recorded.last_followup_day ? p.recorded.last_followup_day + ' d' : '–'}</span>`}
      <span class="chip chip-sm"><span class="muted">${t('tl.viewing')}</span>&nbsp;<b>${U.dayLabel(cps[view])}</b></span>${view < horizon ? `<button class="btn btn-xs" data-tl-latest>${t('twin.view_latest')}</button>` : ''}</div></div>
    <svg class="tl-svg" viewBox="0 0 ${W} ${H}" tabindex="0" role="slider" aria-valuemin="0" aria-valuemax="${horizon}" aria-valuenow="${view}">
      <defs><pattern id="tl-hatch" width="8" height="8" patternUnits="userSpaceOnUse" patternTransform="rotate(45)"><line class="tl-hatch-line" x1="0" y1="0" x2="0" y2="8"/></pattern></defs>
      <rect class="tl-ld" x="${X(-5)}" y="${trackY - 26}" width="${X(-2) - X(-5)}" height="52" rx="3"/>
      ${active ? `<rect class="tl-future" x="${futureX}" y="${trackY - 24}" width="${X(x1) - futureX}" height="48" rx="4" fill="url(#tl-hatch)"/>` : ''}
      <line class="tl-track" x1="${X(x0)}" x2="${X(x1)}" y1="${trackY}" y2="${trackY}"/>
      <line class="tl-track-done" x1="${X(x0)}" x2="${X(cps[view])}" y1="${trackY}" y2="${trackY}"/>
      <line class="tl-infusion" x1="${X(0)}" x2="${X(0)}" y1="${trackY - 26}" y2="${trackY + 10}"/>
      ${evHtml}${cpHtml}
      ${active ? `<text x="${X(x1)}" y="${trackY - 30}" text-anchor="end" class="tl-future-label">${t('tl.future')} →</text>` : ''}
    </svg>
    <div class="tl-foot muted small">${view < horizon ? t('tl.retro') : (active ? t('tl.live_note') : t('tl.completed_note'))}</div>`;
    const svg = container.querySelector('.tl-svg');
    const pick = (clientX) => {
      const pt = svg.createSVGPoint(); pt.x = clientX; pt.y = 0; const sp = pt.matrixTransform(svg.getScreenCTM().inverse());
      const d = x0 + (sp.x - m.l) / (W - m.l - m.r) * (x1 - x0);
      let best = 0, bd = Infinity; cps.forEach((c, i) => { if (i <= horizon && Math.abs(c - d) < bd) { bd = Math.abs(c - d); best = i; } });
      return best;
    };
    let dragging = false;
    svg.addEventListener('pointerdown', (e) => { dragging = true; try { svg.setPointerCapture(e.pointerId); } catch (err) { /* synthetic pointer */ } const i = pick(e.clientX); if (i !== view) opt.onSelect(i); });
    svg.addEventListener('pointermove', (e) => { if (!dragging) return; const i = pick(e.clientX); if (i !== view) { dragging = false; opt.onSelect(i); } });
    svg.addEventListener('pointerup', () => { dragging = false; });
    svg.addEventListener('keydown', (e) => { if (e.key === 'ArrowRight') opt.onSelect(Math.min(horizon, view + 1)); if (e.key === 'ArrowLeft') opt.onSelect(Math.max(0, view - 1)); });
    const latest = container.querySelector('[data-tl-latest]'); if (latest) latest.addEventListener('click', () => opt.onSelect(horizon));
  };
})(typeof window !== 'undefined' ? window : global);
