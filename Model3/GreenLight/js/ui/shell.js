/* GreenLight · ui/shell.js
 * App shell: router (#/ landing, #/clinic/<id>, #/research/<tab>[/<id>], #/admin), login and lock
 * screens, role guards (GL.session), top bar with logo, portal switch, user, population source,
 * language and theme, the access banner and the watermark over patient views. Views register in GL.views.
 */
(function (root) {
  const GL = root.GL = root.GL || {};
  const t = (k, p) => GL.i18n.t(k, p), esc = (s) => GL.util.escapeHtml(s);
  const A = GL.app = { state: { view: 'landing', args: [] } };
  GL.views = GL.views || {};

  // logo: the icon of assets/image.png (logo_mark.png); brand: the full logo with the wordmark, a light-text
  // variant for the dark theme (CSS shows one of the two)
  A.logo = (cls) => `<img class="${cls || 'logo'}" src="assets/logo_mark.png" alt="" aria-hidden="true">`;
  A.brand = (cls) => `<span class="brandimg ${cls || ''}"><img class="lt" src="assets/logo_full.png" alt="GreenLight"><img class="dk" src="assets/logo_full_dark.png" alt="GreenLight"></span>`;
  A.wordmark = () => `<span class="wordmark"><span class="b1">Green</span><span class="b2">Light</span></span>`;
  A.icon = function (name) {
    const p = {
      steth: '<path d="M6 3v6a4 4 0 0 0 8 0V3"/><path d="M10 13v3a5 5 0 0 0 10 0v-2"/><circle cx="20" cy="12" r="2"/>',
      flask: '<path d="M9 3h6"/><path d="M10 3v6L4.5 19a1.5 1.5 0 0 0 1.3 2h12.4a1.5 1.5 0 0 0 1.3-2L14 9V3"/><path d="M7 15h10"/>',
      play: '<path d="M7 5l12 7-12 7z"/>', pause: '<path d="M7 5h4v14H7zM13 5h4v14h-4z"/>', prev: '<path d="M15 6l-6 6 6 6"/>', next: '<path d="M9 6l6 6-6 6"/>',
      loop: '<path d="M17 2l4 4-4 4"/><path d="M3 11V9a3 3 0 0 1 3-3h15"/><path d="M7 22l-4-4 4-4"/><path d="M21 13v2a3 3 0 0 1-3 3H3"/>',
      pingpong: '<path d="M4 12h16"/><path d="M8 8l-4 4 4 4"/><path d="M16 8l4 4-4 4"/>',
      sun: '<circle cx="12" cy="12" r="4"/><path d="M12 2v2M12 20v2M4.9 4.9l1.4 1.4M17.7 17.7l1.4 1.4M2 12h2M20 12h2M4.9 19.1l1.4-1.4M17.7 6.3l1.4-1.4"/>',
      moon: '<path d="M21 12.8A9 9 0 1 1 11.2 3a7 7 0 0 0 9.8 9.8z"/>', home: '<path d="M3 11l9-8 9 8"/><path d="M5 10v10h14V10"/>',
      sliders: '<path d="M4 21v-7M4 10V3M12 21v-9M12 8V3M20 21v-5M20 12V3M1 14h6M9 8h6M17 16h6"/>', download: '<path d="M12 3v12M7 10l5 5 5-5M4 21h16"/>',
      list: '<path d="M8 6h13M8 12h13M8 18h13M3 6h.01M3 12h.01M3 18h.01"/>', lock: '<rect x="4" y="11" width="16" height="10" rx="2"/><path d="M8 11V7a4 4 0 0 1 8 0v4"/>', arrow: '<path d="M5 12h14M13 6l6 6-6 6"/>', ingest: '<path d="M12 3v12M7 10l5 5 5-5"/><rect x="3" y="17" width="18" height="4" rx="1"/>'
    }[name] || '';
    return `<svg class="ic" viewBox="0 0 24 24" fill="none" stroke="currentColor" stroke-width="2" stroke-linecap="round" stroke-linejoin="round" aria-hidden="true">${p}</svg>`;
  };

  // patient-type mark: status colour + origin symbol (real = filled dot, digital = hollow diamond + DP)
  A.typeMark = function (p, withText) {
    const st = p.status === 'active' ? 'active' : 'completed';
    const col = st === 'active' ? 'var(--brand-600)' : 'var(--slate)';
    const sym = p.type === 'digital'
      ? `<svg class="pmark" viewBox="0 0 12 12"><path d="M6 1l5 5-5 5-5-5z" fill="none" stroke="${col}" stroke-width="2"/></svg>`
      : `<svg class="pmark" viewBox="0 0 12 12"><circle cx="6" cy="6" r="4.5" fill="${col}"/></svg>`;
    const stat = st === 'active' ? '<span class="pulse" title="' + esc(t('type.active')) + '"></span>' : '<span class="check" title="' + esc(t('type.completed')) + '">✓</span>';
    return `<span class="ptype st-${st}">${sym}${p.type === 'digital' ? '<span class="dp-tag">DP</span>' : ''}${stat}${withText ? `<span class="small muted">${esc(t('type.' + st))} · ${esc(t(p.type === 'digital' ? 'type.digital' : 'type.real'))}</span>` : ''}</span>`;
  };

  A.theme = function (th) {
    if (th) { GL.store.set('theme', th); }
    document.documentElement.dataset.theme = GL.store.get('theme', 'light');
  };

  A.top = function () {
    const v = A.state.view, idx = GL.loader.index(), S = GL.session;
    const pop = idx ? idx.population : null;
    const langs = GL.i18n.languages.map(l => `<option value="${l.code}" ${l.code === GL.i18n.lang ? 'selected' : ''}>${l.label}</option>`).join('');
    const dark = GL.store.get('theme', 'light') === 'dark';
    const portals = S.active() && !['landing', 'login', 'admin'].includes(v) && !S.can('admin');
    const who = S.user ? `<span class="chip chip-user" id="userChip">${esc(S.user)} · <b>${esc(t('role.' + S.role))}</b></span>
        <button class="btn btn-xs" id="btnLock" title="${esc(t('auth.lock_h'))}">${esc(t('auth.lock'))}</button><button class="btn btn-xs" id="btnLogout">${esc(t('auth.logout'))}</button>`
      : S.guest ? `<span class="chip chip-user" id="userChip">${esc(t('role.guest'))}</span><button class="btn btn-xs" id="btnLogout">${esc(t('auth.signin'))}</button>` : '';
    return `<header class="top">
      <a class="brand" href="#/">${A.brand('bi-top')}</a>
      <span class="tag">${esc(t('brand.tag'))}</span>
      ${portals ? `<nav class="portal-switch">${S.can('clinic') ? `<a href="#/clinic" class="${v === 'clinic' ? 'active' : ''}" id="navClinic">${A.icon('steth')}${esc(t('portal.clinic'))}</a>` : ''}${S.can('research') ? `<a href="#/research/model" class="${v === 'research' ? 'active' : ''}" id="navResearch">${A.icon('flask')}${esc(t('portal.research'))}</a>` : ''}</nav>` : ''}
      <span class="spacer"></span>
      ${pop && v !== 'login' ? `<span class="chip chip-twin pop-src" title="${esc(t('pop.title'))}">${esc(t('pop.label'))}: <b>${esc(t('pop.' + pop.source))}</b> · ${esc(pop.model)}</span> ${GL.help.btn('population')}` : ''}
      ${who}
      <select id="langSel" class="lang" aria-label="${esc(t('ui.language'))}">${langs}</select>
      <button class="btn btn-icon" id="btnTheme" title="${esc(t('ui.theme'))}">${A.icon(dark ? 'sun' : 'moon')}</button>
    </header><div class="authbar" role="note">${A.icon('lock')}<span>${esc(t('auth.banner'))}</span></div>`;
  };

  // watermark over patient views: user name and time (deterrent, docs/PRIVACY.md)
  A.watermark = function (show) {
    let w = document.getElementById('wm');
    if (!show) { if (w) w.remove(); return; }
    if (!w) { w = document.createElement('div'); w.id = 'wm'; w.setAttribute('aria-hidden', 'true'); document.body.appendChild(w); }
    const txt = `${GL.session.user} · ${new Date().toISOString().slice(0, 16).replace('T', ' ')}`;
    const svg = `<svg xmlns="http://www.w3.org/2000/svg" width="420" height="220"><text x="20" y="130" transform="rotate(-22 210 110)" font-family="sans-serif" font-size="18" fill="#808080">${esc(txt)}</text></svg>`;
    w.style.backgroundImage = `url("data:image/svg+xml;utf8,${encodeURIComponent(svg)}")`;
    w.dataset.text = txt;
  };
  setInterval(() => { if (document.getElementById('wm') && GL.session.user) A.watermark(true); }, 60000);

  A.route = async function () {
    const S = GL.session;
    const h = (location.hash || '#/').replace(/^#\/?/, '');
    const parts = h.split('/').filter(Boolean).map(decodeURIComponent);
    let view = parts[0] && GL.views[parts[0]] ? parts[0] : 'landing';
    if (!S.active()) view = 'login';
    else if (S.can('admin')) { if (view !== 'admin') { location.replace('#/admin'); return; } }
    else if (view === 'admin' || view === 'login' || (view === 'clinic' && !S.can('clinic')) || (view === 'research' && !S.can('research'))) view = 'landing';
    A.state.view = view; A.state.args = parts.slice(1);
    await A.render();
  };
  A.render = async function () {
    const app = document.getElementById('app'), S = GL.session;
    document.body.classList.toggle('fit-screen', A.state.view === 'clinic');
    GL.help.close();
    app.innerHTML = A.top() + '<main id="view"></main>';
    document.getElementById('langSel').addEventListener('change', (e) => { GL.i18n.set(e.target.value); A.render(); });
    document.getElementById('btnTheme').addEventListener('click', () => { A.theme(GL.store.get('theme', 'light') === 'dark' ? 'light' : 'dark'); A.render(); });
    const lo = document.getElementById('btnLogout'); if (lo) lo.addEventListener('click', () => S.logout());
    const lk = document.getElementById('btnLock'); if (lk) lk.addEventListener('click', () => S.lock());
    A.watermark(!!S.user && !S.locked && ['clinic', 'research'].includes(A.state.view));
    if (S.locked) { GL.views.lock.render(document.getElementById('view')); return; }
    await GL.views[A.state.view].render(document.getElementById('view'), A.state.args);
  };

  // ---- login, lock and admin views ------------------------------------------------------------------
  GL.views.login = {
    async render(el) {
      const S = GL.session, none = !S.users().length;
      el.innerHTML = `<section class="login"><form class="login-card" id="loginForm" autocomplete="off">
        <div class="login-brand">${A.brand('bi-login')}</div>
        <div class="muted small">${esc(t('brand.tag'))}</div>
        <h2>${esc(t('auth.title'))}</h2>
        <label class="small" for="lgUser">${esc(t('auth.user'))}</label><input class="input" id="lgUser" autocomplete="username" ${none ? 'disabled' : ''}>
        <label class="small" for="lgPw">${esc(t('auth.password'))}</label><input class="input" id="lgPw" type="password" autocomplete="current-password" ${none ? 'disabled' : ''}>
        <div class="login-err" id="lgErr" role="alert">${none ? esc(t('auth.no_users')) : ''}${GL.crypto.available() ? '' : esc(t('auth.no_crypto'))}</div>
        <button class="btn btn-primary" id="lgGo" type="submit" ${none ? 'disabled' : ''}>${esc(t('auth.signin'))}</button>
        <div class="login-or"><span>${esc(t('auth.or'))}</span></div>
        <button class="btn" id="lgGuest" type="button">${GL.app.typeMark({ type: 'digital', status: 'completed' })} ${esc(t('auth.guest'))}</button>
        <div class="tiny muted">${esc(t('auth.guest_h'))}</div>
        <div class="tiny muted" style="margin-top:6px">${esc(t('auth.note'))}</div>
      </form></section>`;
      const f = el.querySelector('#loginForm');
      f.addEventListener('submit', async (e) => {
        e.preventDefault();
        const u = el.querySelector('#lgUser').value.trim(), pw = el.querySelector('#lgPw');
        const btn = el.querySelector('#lgGo'), err = el.querySelector('#lgErr');
        btn.disabled = true; err.textContent = t('auth.opening');
        try { await S.login(u, pw.value); pw.value = ''; location.hash = S.can('admin') ? '#/admin' : '#/'; await A.route(); }
        catch (x) { pw.value = ''; err.textContent = t('auth.failed'); btn.disabled = false; }
      });
      el.querySelector('#lgGuest').addEventListener('click', () => { S.guestMode(); location.hash = '#/'; A.route(); });
      const uf = el.querySelector('#lgUser'); if (uf && !uf.disabled) uf.focus();
    }
  };
  GL.views.lock = {
    render(el) {
      el.innerHTML = `<section class="login"><form class="login-card" id="lockForm" autocomplete="off">
        <div class="login-brand">${A.brand('bi-login')}</div>
        <h2>${A.icon('lock')} ${esc(t('auth.locked'))}</h2>
        <div class="muted small">${esc(t('auth.locked_h', { user: GL.session.user }))}</div>
        <label class="small" for="ulPw">${esc(t('auth.password'))}</label><input class="input" id="ulPw" type="password" autocomplete="current-password">
        <div class="login-err" id="ulErr" role="alert"></div>
        <button class="btn btn-primary" id="ulGo" type="submit">${esc(t('auth.unlock'))}</button>
        <button class="btn" id="ulOut" type="button">${esc(t('auth.logout'))}</button>
      </form></section>`;
      el.querySelector('#lockForm').addEventListener('submit', async (e) => {
        e.preventDefault();
        const pw = el.querySelector('#ulPw'), err = el.querySelector('#ulErr');
        err.textContent = t('auth.opening');
        try { await GL.session.unlock(pw.value); pw.value = ''; } catch (x) { pw.value = ''; err.textContent = t('auth.failed'); }
      });
      el.querySelector('#ulOut').addEventListener('click', () => GL.session.logout());
      el.querySelector('#ulPw').focus();
    }
  };
  GL.views.admin = {
    async render(el) {
      const users = GL.session.users(), log = await GL.audit.all(), ver = await GL.audit.verify(log);
      el.innerHTML = `<section class="rmain admin">
        <div class="rhead"><div><h1>${esc(t('admin.title'))}</h1><div class="muted">${esc(t('admin.lead'))}</div></div></div>
        <div class="rgrid2"><div class="rcard"><h2>${esc(t('admin.users'))}</h2><table class="tbl"><thead><tr><th>${esc(t('auth.user'))}</th><th>${esc(t('admin.role'))}</th><th>${esc(t('admin.keys'))}</th></tr></thead><tbody>
          ${users.map(u => `<tr><td><b>${esc(u.u)}</b></td><td>${esc(t('role.' + u.role))}</td><td>${u.files.length}</td></tr>`).join('') || `<tr><td colspan="3" class="muted">${esc(t('auth.no_users'))}</td></tr>`}</tbody></table>
          <div class="tiny muted" style="margin-top:6px">${esc(t('admin.cli'))}</div></div>
          <div class="rcard"><h2>${esc(t('admin.audit'))}</h2>
            <div class="kv"><span>${esc(t('admin.entries'))}</span><b id="auditN">${log.length}</b></div>
            <div class="kv"><span>${esc(t('admin.chain'))}</span><b id="auditChain" class="${ver.ok ? 'ok-badge' : 'bad-badge'}">${esc(ver.ok ? t('admin.chain_ok') : t('admin.chain_bad', { seq: ver.at }))}</b></div>
            <button class="btn btn-primary" id="auditExport" style="margin-top:8px">${A.icon('download')} ${esc(t('admin.export'))}</button>
            <div class="tiny muted" style="margin-top:6px">${esc(t('admin.local'))}</div></div></div>
        <div class="rcard"><h2>${esc(t('admin.recent'))}</h2><table class="tbl"><thead><tr><th>#</th><th>${esc(t('admin.time'))}</th><th>${esc(t('auth.user'))}</th><th>${esc(t('admin.role'))}</th><th>${esc(t('admin.subject'))}</th><th>${esc(t('admin.action'))}</th></tr></thead><tbody>
          ${log.slice(-100).reverse().map(e => `<tr><td>${e.seq}</td><td class="small">${esc(e.ts)}</td><td>${esc(e.user)}</td><td>${esc(e.role)}</td><td>${esc(e.subject)}</td><td>${esc(e.action)}</td></tr>`).join('')}</tbody></table></div>
      </section>`;
      el.querySelector('#auditExport').addEventListener('click', async () => {
        await GL.audit.log('export', 'audit-log');
        const all = await GL.audit.all(), v = await GL.audit.verify(all);
        const blob = new Blob([JSON.stringify({ exported: new Date().toISOString(), by: GL.session.user, chain_ok: v.ok, entries: all }, null, 1)], { type: 'application/json' });
        const a = document.createElement('a'); a.href = URL.createObjectURL(blob); a.download = 'greenlight_audit_' + Date.now() + '.json'; a.id = 'auditLink';
        document.body.appendChild(a); a.click(); setTimeout(() => { URL.revokeObjectURL(a.href); a.remove(); }, 1000);
        GL.views.admin.lastExport = { n: all.length, ok: v.ok };
      });
    }
  };

  A.init = async function () {
    GL.i18n.init();
    A.theme();
    GL.help.init();
    GL.remote.registerAll();
    await GL.store.purgeLegacy();
    await GL.audit.init();
    window.addEventListener('hashchange', () => A.route());
    let rt = null;
    window.addEventListener('resize', () => { clearTimeout(rt); rt = setTimeout(() => { const v = GL.views[A.state.view]; if (v && v.onResize) v.onResize(); }, 120); });
    await A.route();
  };
})(typeof window !== 'undefined' ? window : global);
