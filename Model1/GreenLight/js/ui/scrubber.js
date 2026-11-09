/* GreenLight · ui/scrubber.js
 * Pinned day scrubber: the digital twin "as of" any day (D-5..D+30 daily, weekly to D+90).
 * Drag or click the track, step with the buttons or the arrow keys, play once / loop /
 * ping-pong (back and forth) at 0.5x-4x. For a patient in treatment the track stops at the data
 * horizon (hatched future). The time axis is the same piecewise axis as the charts.
 */
(function (root) {
  const GL = root.GL = root.GL || {};
  const t = (k, p) => GL.i18n.t(k, p), esc = (s) => GL.util.escapeHtml(s), I = (n) => GL.app.icon(n);
  GL.ui = GL.ui || {};
  // o: { days (all as-of days), horizon, asOf, active, gcsf, onChange(day), onPlay(), onIngest() | null, canIngest }
  GL.ui.scrubber = function (el, o) {
    const S = GL.ui.scrubber.state = GL.ui.scrubber.state || { playing: false, mode: GL.store.get('playMode', 'pingpong'), speed: GL.store.get('playSpeed', 1), dir: 1, timer: null };
    const sel = o.days.filter(d => d <= o.horizon);
    const label = (d) => d === 0 ? 'D0' : (d > 0 ? 'D+' + d : 'D' + d);
    el.innerHTML = `<div class="ctrls">
        <div class="asof"><div class="tiny muted">${esc(t('scrub.asof'))} ${GL.help.btn('asof')}</div><div class="big" id="scrubDay">${label(o.asOf)}</div></div>
        <button class="btn btn-icon" id="scPrev" title="${esc(t('scrub.prev'))}">${I('prev')}</button>
        <button class="btn btn-icon btn-primary" id="scPlay" title="${esc(t(S.playing ? 'scrub.pause' : 'scrub.play'))}">${I(S.playing ? 'pause' : 'play')}</button>
        <button class="btn btn-icon" id="scNext" title="${esc(t('scrub.next'))}">${I('next')}</button>
      </div>
      <svg id="scSvg" tabindex="0" role="slider" aria-label="${esc(t('scrub.asof'))}" aria-valuenow="${o.asOf}"></svg>
      <div class="rctrls">
        <div class="ctrls">
          <button class="btn btn-xs ${S.mode === 'once' ? 'on' : ''}" data-mode="once" title="${esc(t('scrub.once'))}">1×</button>
          <button class="btn btn-xs ${S.mode === 'loop' ? 'on' : ''}" data-mode="loop" title="${esc(t('scrub.loop'))}">${I('loop')}</button>
          <button class="btn btn-xs ${S.mode === 'pingpong' ? 'on' : ''}" data-mode="pingpong" id="scPing" title="${esc(t('scrub.pingpong'))}">${I('pingpong')}</button>
          <select id="scSpeed" title="${esc(t('scrub.speed'))}">${[0.5, 1, 2, 4].map(v => `<option value="${v}" ${+S.speed === v ? 'selected' : ''}>${v}×</option>`).join('')}</select>
        </div>
        ${o.active ? `<button class="btn btn-xs" id="scIngest" ${o.canIngest ? '' : 'disabled'} title="${esc(t('scrub.ingest_h'))}">${I('ingest')} ${esc(t('scrub.ingest'))}</button>` : `<span class="tiny muted">${esc(t('scrub.completed'))}</span>`}
      </div>`;
    const svg = el.querySelector('#scSvg');
    const draw = (asOf) => {
      const W = Math.max(300, svg.clientWidth || 600), H = 62, m = { l: 12, r: 12 };
      const X = GL.charts.timeScale(m.l, W - m.r, -10, 90), y = 30;
      let s = `<defs><pattern id="hatch" width="7" height="7" patternUnits="userSpaceOnUse" patternTransform="rotate(45)"><line class="tl-hatch-line" x1="0" y1="0" x2="0" y2="7"/></pattern></defs>`;
      s += `<rect class="tl-ld" x="${X(-5)}" y="${y - 14}" width="${X(-2) - X(-5)}" height="28" rx="3"><title>${esc(t('scrub.ld'))}</title></rect>`;
      (o.gcsf || []).forEach(w => { if (w[0] <= o.horizon) s += `<rect class="tl-gcsf" x="${X(w[0])}" y="${y + 8}" width="${Math.max(2, X(Math.min(w[1], o.horizon)) - X(w[0]))}" height="4"><title>G-CSF</title></rect>`; });
      if (o.active && o.horizon < 90) s += `<rect class="tl-future" x="${X(o.horizon)}" y="${y - 14}" width="${X(90) - X(o.horizon)}" height="28" rx="3"/>`;
      s += `<line class="tl-axis" x1="${X(-10)}" x2="${X(90)}" y1="${y}" y2="${y}"/><line class="tl-done" x1="${X(-10)}" x2="${X(asOf)}" y1="${y}" y2="${y}"/>`;
      [-5, 0, 7, 14, 21, 30, 44, 58, 72, 86].forEach(d => { s += `<g class="tl-tick"><line x1="${X(d)}" x2="${X(d)}" y1="${y + 6}" y2="${y + 10}"/><text x="${X(d)}" y="${y + 22}" text-anchor="middle">${label(d)}</text></g>`; });
      sel.forEach(d => { s += `<circle cx="${X(d)}" cy="${y}" r="1.6" fill="var(--surface)" opacity=".9"/>`; });
      if (o.active) s += `<line class="tl-horizon" x1="${X(o.horizon)}" x2="${X(o.horizon)}" y1="${y - 18}" y2="${y + 14}"/><text class="tl-label" x="${X(o.horizon)}" y="${y - 20}" text-anchor="middle">${esc(t('scrub.horizon'))}</text>`;
      s += `<text class="tl-label" x="${(X(-5) + X(-2)) / 2}" y="${y - 18}" text-anchor="middle">LD</text>`;
      s += `<circle class="tl-handle" cx="${X(asOf)}" cy="${y}" r="9"/><circle cx="${X(asOf)}" cy="${y}" r="3.2" fill="var(--brand-600)"/>`;
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
    svg.addEventListener('pointerdown', (e) => { drag = true; try { svg.setPointerCapture(e.pointerId); } catch (err) { /* synthetic */ } stop(); set(pick(e.clientX)); });
    svg.addEventListener('pointermove', (e) => { if (drag) set(pick(e.clientX)); });
    svg.addEventListener('pointerup', () => { drag = false; });
    const step = (k) => { const i = sel.indexOf(cur); const j = Math.min(sel.length - 1, Math.max(0, i + k)); set(sel[j]); };
    svg.addEventListener('keydown', (e) => { if (e.key === 'ArrowRight') { stop(); step(1); } if (e.key === 'ArrowLeft') { stop(); step(-1); } if (e.key === ' ') { e.preventDefault(); toggle(); } });
    el.querySelector('#scPrev').addEventListener('click', () => { stop(); step(-1); });
    el.querySelector('#scNext').addEventListener('click', () => { stop(); step(1); });
    function tick() {
      const i = sel.indexOf(cur);
      let j = i + S.dir;
      if (j >= sel.length || j < 0) {
        if (S.mode === 'once') { stop(); return; }
        if (S.mode === 'loop') j = S.dir > 0 ? 0 : sel.length - 1;
        else { S.dir = -S.dir; j = i + S.dir; }
      }
      set(sel[Math.min(sel.length - 1, Math.max(0, j))]);
    }
    function play() { stop(true); S.playing = true; if (o.onPlay) o.onPlay(); if (S.mode !== 'pingpong') S.dir = 1; if (S.mode === 'once' && cur === sel[sel.length - 1]) set(sel[0]); S.timer = setInterval(tick, 420 / S.speed); refreshBtn(); }
    function stop(silent) { if (S.timer) clearInterval(S.timer); S.timer = null; S.playing = false; if (!silent) refreshBtn(); }
    function toggle() { if (S.playing) stop(); else play(); }
    function refreshBtn() { const b = el.querySelector('#scPlay'); if (b) { b.innerHTML = I(S.playing ? 'pause' : 'play'); b.title = t(S.playing ? 'scrub.pause' : 'scrub.play'); } }
    el.querySelector('#scPlay').addEventListener('click', toggle);
    el.querySelectorAll('[data-mode]').forEach(b => b.addEventListener('click', () => { S.mode = b.dataset.mode; GL.store.set('playMode', S.mode); el.querySelectorAll('[data-mode]').forEach(x => x.classList.toggle('on', x === b)); }));
    el.querySelector('#scSpeed').addEventListener('change', (e) => { S.speed = +e.target.value; GL.store.set('playSpeed', S.speed); if (S.playing) play(); });
    const ing = el.querySelector('#scIngest'); if (ing) ing.addEventListener('click', () => { stop(); o.onIngest(); });
    if (S.playing) { S.playing = false; play(); }
    return { set, stop, redraw: () => { X = draw(cur); }, get day() { return cur; } };
  };
  GL.ui.scrubber.stopAll = function () { const S = GL.ui.scrubber.state; if (S && S.timer) { clearInterval(S.timer); S.timer = null; S.playing = false; } };
})(typeof window !== 'undefined' ? window : global);
