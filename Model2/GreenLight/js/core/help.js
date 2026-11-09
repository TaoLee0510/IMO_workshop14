/* GreenLight · core/help.js
 * "?" buttons and one popover. Texts come from the dictionaries:
 *   help.<id>.t (title), help.<id> (body), help.<id>.unit, help.<id>.dir ('higher'|'lower'|'neutral').
 */
(function (root) {
  const GL = root.GL = root.GL || {};
  const H = GL.help = {};
  const t = (k, p) => GL.i18n.t(k, p);
  const esc = (s) => GL.util.escapeHtml(s);
  H.btn = (id) => `<button type="button" class="help-btn" data-help="${id}" aria-label="${esc(t('help.open'))}" title="${esc(t('help.open'))}">?</button>`;
  let pop = null, anchor = null;
  H.close = function () { if (pop) pop.style.display = 'none'; if (anchor) anchor.classList.remove('is-open'); anchor = null; };
  H.open = function (btn) {
    const id = btn.dataset.help, has = (k) => GL.i18n.has(k);
    if (!pop) { pop = document.createElement('div'); pop.className = 'help-pop'; pop.setAttribute('role', 'dialog'); document.body.appendChild(pop); }
    const dir = has('help.' + id + '.dir') ? t('help.' + id + '.dir') : '';
    const dirTxt = dir === 'higher' ? t('help.dir.higher') : dir === 'lower' ? t('help.dir.lower') : dir === 'neutral' ? t('help.dir.neutral') : '';
    const unit = has('help.' + id + '.unit') ? t('help.' + id + '.unit') : '';
    pop.innerHTML = `<div class="help-head"><b>${esc(has('help.' + id + '.t') ? t('help.' + id + '.t') : id)}</b><button type="button" class="help-x" aria-label="${esc(t('help.close'))}">✕</button></div>
      <div class="help-meta">${unit ? `<span class="chip">${esc(t('help.unit'))}: <b>${esc(unit)}</b></span>` : ''}${dirTxt ? `<span class="chip ${dir === 'higher' ? 'dir-good' : dir === 'lower' ? 'dir-bad' : ''}">${esc(dirTxt)}</span>` : ''}</div>
      <div>${esc(has('help.' + id) ? t('help.' + id) : '').replace(/\n/g, '<br>')}</div>`;
    pop.style.display = 'block';
    if (anchor) anchor.classList.remove('is-open');
    anchor = btn; btn.classList.add('is-open');
    const b = btn.getBoundingClientRect(), w = Math.min(340, window.innerWidth - 24);
    pop.style.width = w + 'px';
    let left = Math.min(Math.max(12, b.left + window.scrollX - 8), window.scrollX + window.innerWidth - 12 - w);
    pop.style.left = left + 'px'; pop.style.top = (b.bottom + window.scrollY + 8) + 'px';
    const ph = pop.getBoundingClientRect().height;
    if (b.bottom + 8 + ph > window.innerHeight && b.top - 8 - ph > 0) pop.style.top = (b.top + window.scrollY - 8 - ph) + 'px';
  };
  H.init = function () {
    document.addEventListener('click', (e) => {
      if (e.target.closest && e.target.closest('.help-x')) { H.close(); return; }
      const b = e.target.closest && e.target.closest('[data-help]');
      if (b) { e.preventDefault(); e.stopPropagation(); if (anchor === b) H.close(); else H.open(b); return; }
      if (pop && !(e.target.closest && e.target.closest('.help-pop'))) H.close();
    });
    document.addEventListener('keydown', (e) => { if (e.key === 'Escape') H.close(); });
    window.addEventListener('resize', H.close);
  };
})(typeof window !== 'undefined' ? window : global);
