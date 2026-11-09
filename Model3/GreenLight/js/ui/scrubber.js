/* GreenLight · ui/scrubber.js
 * Day scrubber: the digital twin "as of" any day (D-5..D+30 daily, weekly to D+90). A card under the
 * charts of the left column; its time axis has the margins of the full-width charts above (inflammation
 * driver, how the twin learned), so the days line up. Drag or click the track, step with the buttons or
 * the arrow keys (no playback). For a patient in treatment the track stops at the data horizon (hatched future).
 * G-CSF windows: blue up to the selected day (known to the twin then), grey after it.
 */
(function (root) {
  const GL = root.GL = root.GL || {};
  const t = (k, p) => GL.i18n.t(k, p), esc = (s) => GL.util.escapeHtml(s), I = (n) => GL.app.icon(n);
  GL.ui = GL.ui || {};
  // o: { days (all as-of days), horizon, asOf, active, gcsf, onChange(day), onIngest() | null, canIngest, note }
  GL.ui.scrubber = function (el, o) {
    const sel = o.days.filter(d => d <= o.horizon);
    const label = (d) => d === 0 ? 'D0' : (d > 0 ? 'D+' + d : 'D' + d);
    el.innerHTML = `<div class="card-head scrub-head">
        <h3>${esc(t('scrub.asof'))} ${GL.help.btn('asof')}</h3><b class="scrub-day" id="scrubDay">${label(o.asOf)}</b>
        <button class="btn btn-icon btn-xs" id="scPrev" title="${esc(t('scrub.prev'))}">${I('prev')}</button>
        <button class="btn btn-icon btn-xs" id="scNext" title="${esc(t('scrub.next'))}">${I('next')}</button>
        <span class="spacer"></span>
        ${o.active && o.onIngest ? `<button class="btn btn-xs" id="scIngest" ${o.canIngest ? '' : 'disabled'} title="${esc(t('scrub.ingest_h'))}">${I('ingest')} ${esc(t('scrub.ingest'))}</button>` : `<span class="tiny muted">${esc(o.note || t('scrub.completed'))}</span>`}
      </div>
      <svg id="scSvg" tabindex="0" role="slider" aria-label="${esc(t('scrub.asof'))}" aria-valuenow="${o.asOf}"></svg>`;
    const svg = el.querySelector('#scSvg');
    const draw = (asOf) => {
      const W = Math.max(300, svg.clientWidth || 600), H = 46, m = GL.charts.xMargin;
      const X = GL.charts.timeScale(m.l, W - m.r, -10, 90), y = 20;
      let s = `<defs><pattern id="hatch" width="7" height="7" patternUnits="userSpaceOnUse" patternTransform="rotate(45)"><line class="tl-hatch-line" x1="0" y1="0" x2="0" y2="7"/></pattern></defs>`;
      s += `<rect class="tl-ld" x="${X(-5)}" y="${y - 10}" width="${X(-2) - X(-5)}" height="20" rx="3"><title>${esc(t('scrub.ld'))}</title></rect>`;
      // G-CSF given up to the selected day in blue (the twin knows it), later G-CSF grey
      GL.charts.gcsfParts((o.gcsf || []).filter(w => w[0] <= o.horizon).map(w => [w[0], Math.min(w[1], o.horizon)]), asOf)
        .forEach(g => { s += `<rect class="tl-gcsf${g.future ? ' tl-gcsf-future' : ''}" x="${X(g.s)}" y="${y + 6}" width="${Math.max(2, X(g.e) - X(g.s))}" height="4"><title>G-CSF</title></rect>`; });
      if (o.active && o.horizon < 90) s += `<rect class="tl-future" x="${X(o.horizon)}" y="${y - 10}" width="${X(90) - X(o.horizon)}" height="20" rx="3"/>`;
      s += `<line class="tl-axis" x1="${X(-10)}" x2="${X(90)}" y1="${y}" y2="${y}"/><line class="tl-done" x1="${X(-10)}" x2="${X(asOf)}" y1="${y}" y2="${y}"/>`;
      GL.charts.xTicks.forEach(d => { s += `<g class="tl-tick"><line x1="${X(d)}" x2="${X(d)}" y1="${y + 8}" y2="${y + 11}"/><text x="${X(d)}" y="${y + 22}" text-anchor="${d === 90 ? 'end' : 'middle'}">${label(d)}</text></g>`; });
      sel.forEach(d => { s += `<circle cx="${X(d)}" cy="${y}" r="1.6" fill="var(--surface)" opacity=".9"/>`; });
      if (o.active) s += `<line class="tl-horizon" x1="${X(o.horizon)}" x2="${X(o.horizon)}" y1="${y - 13}" y2="${y + 11}"/><text class="tl-label" x="${X(o.horizon) + 4}" y="${y - 12}">${esc(t('scrub.horizon'))}</text>`;
      s += `<text class="tl-label" x="${(X(-5) + X(-2)) / 2}" y="${y - 12}" text-anchor="middle">LD</text>`;
      s += `<circle class="tl-handle" cx="${X(asOf)}" cy="${y}" r="8"/><circle cx="${X(asOf)}" cy="${y}" r="3" fill="var(--brand-600)"/>`;
      svg.setAttribute('viewBox', `0 0 ${W} ${H}`); svg.innerHTML = s;
      svg.setAttribute('aria-valuenow', asOf);
      el.querySelector('#scrubDay').textContent = label(asOf);
      return X;
    };
    let X = draw(o.asOf), cur = o.asOf;
    const nearest = (d) => sel.reduce((b, x) => Math.abs(x - d) < Math.abs(b - d) ? x : b, sel[0]);
    const pick = (clientX) => {
      const r = svg.getBoundingClientRect(), W = Math.max(300, svg.clientWidth || 600), px = (clientX - r.left) / r.width * W;
      let best = sel[0], bd = Infinity; sel.forEach(d => { const dd = Math.abs(X(d) - px); if (dd < bd) { bd = dd; best = d; } });
      return best;
    };
    const set = (d) => { d = nearest(d); if (d === cur) return; cur = d; X = draw(d); o.onChange(d); };
    let drag = false;
    svg.addEventListener('pointerdown', (e) => { drag = true; try { svg.setPointerCapture(e.pointerId); } catch (err) { /* synthetic */ } set(pick(e.clientX)); });
    svg.addEventListener('pointermove', (e) => { if (drag) set(pick(e.clientX)); });
    svg.addEventListener('pointerup', () => { drag = false; });
    const step = (k) => { const i = sel.indexOf(cur); const j = Math.min(sel.length - 1, Math.max(0, i + k)); set(sel[j]); };
    svg.addEventListener('keydown', (e) => { if (e.key === 'ArrowRight') step(1); if (e.key === 'ArrowLeft') step(-1); });
    el.querySelector('#scPrev').addEventListener('click', () => step(-1));
    el.querySelector('#scNext').addEventListener('click', () => step(1));
    const ing = el.querySelector('#scIngest'); if (ing) ing.addEventListener('click', () => o.onIngest());
    return { set, redraw: () => { X = draw(cur); }, get day() { return cur; } };
  };

  // click or drag on a time chart = move the scrubber to that day; hovering shows a guide at the nearest selectable day.
  // plot: a persistent chart container (its svg is redrawn); o: { days() -> selectable days, set(day) }
  GL.ui.timePick = function (plot, o) {
    plot.classList.add('tplot');
    const dayAt = (e) => {
      const r = plot.getBoundingClientRect(), m = GL.charts.xMargin;
      const d = GL.charts.timeInvert(m.l, r.width - m.r)(e.clientX - r.left), sel = o.days();
      return sel.length ? { day: sel.reduce((b, x) => Math.abs(x - d) < Math.abs(b - d) ? x : b, sel[0]), x: (r2) => GL.charts.timeScale(m.l, r2.width - m.r) } : null;
    };
    const guide = (e) => {
      const h = dayAt(e); if (!h) return;
      let g = plot.querySelector(':scope > .tguide');
      if (!g) { g = document.createElement('div'); g.className = 'tguide'; g.setAttribute('aria-hidden', 'true'); g.innerHTML = '<span></span>'; plot.appendChild(g); }
      g.style.left = h.x(plot.getBoundingClientRect())(h.day) + 'px';
      g.firstChild.textContent = h.day === 0 ? 'D0' : (h.day > 0 ? 'D+' + h.day : 'D' + h.day);
    };
    let drag = false;
    plot.addEventListener('pointerdown', (e) => { const h = dayAt(e); if (!h) return; drag = true; try { plot.setPointerCapture(e.pointerId); } catch (err) { /* synthetic */ } o.set(h.day); guide(e); });
    plot.addEventListener('pointermove', (e) => { if (drag) { const h = dayAt(e); if (h) o.set(h.day); } guide(e); });
    plot.addEventListener('pointerup', () => { drag = false; });
    plot.addEventListener('pointerleave', () => { if (!drag) { const g = plot.querySelector(':scope > .tguide'); if (g) g.remove(); } });
  };
})(typeof window !== 'undefined' ? window : global);
