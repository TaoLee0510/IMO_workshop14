/* GreenLight · core/session.js
 * Login, roles and session security (Model1 plan 13.4–13.6).
 *   roles      clinician  -> clinician portal: assigned real patients (clinical payloads) + digital patients;
 *                            research tabs without patient data (Model, Equations, Symbols, Interfaces)
 *              researcher -> research portal incl. Data, Fit & accuracy, Twins vs truth, Versions with the
 *                            research payloads (pseudonyms, binned fields, no plans / what-if); clinician
 *                            portal with digital patients only
 *              steward    -> user list and audit-log export only (no patient data)
 *              guest      -> no login; both portals, digital patients only
 *   data       decrypted payloads live only in this object (memory); real-patient UI state and what-if
 *              scenarios too. wipe() drops them (logout, idle lock).
 *   idle lock  after idleMs (15 min) of no input: wipe + lock screen; the password reopens the bundle.
 */
(function (root) {
  const GL = root.GL = root.GL || {};
  const S = GL.session = { user: null, role: null, guest: false, locked: false, real: {}, research: {}, cohort: null,
    scenarios: [], state: {}, idleMs: 15 * 60 * 1000, timer: null };
  const ALL_TABS = ['model', 'equations', 'symbols', 'data', 'accuracy', 'twins', 'versions', 'interfaces'];
  S.users = () => (root.GLM1_USERS && root.GLM1_USERS.users) || [];
  S.active = () => !!S.user || S.guest;
  S.can = function (f) {
    const r = S.role;
    return { clinic: ['clinician', 'researcher', 'guest'].includes(r), 'clinic.real': r === 'clinician',
      research: ['clinician', 'researcher', 'guest'].includes(r), admin: r === 'steward', 'whatif.real': r === 'clinician' }[f] || false;
  };
  S.tabs = () => ({ clinician: ['model', 'equations', 'symbols', 'interfaces'], researcher: ALL_TABS,
    guest: ['model', 'equations', 'symbols', 'accuracy', 'twins', 'versions', 'interfaces'] }[S.role] || []);
  S.isReal = (id) => !!(S.real[id] || S.research[id]);

  async function load(user, pk) {
    const recs = await Promise.all(user.files.map(async (e) => {
      const sealed = await GL.loader.sealed(e.f);
      const p = await GL.crypto.openPayload(pk, e.k, sealed);
      GL.loader.dropSealed(e.f);
      return p;
    }));
    S.real = {}; S.research = {}; S.cohort = null;
    for (const p of recs) {
      if (p.kind === 'clinical' && user.role === 'clinician') S.real[p.id] = p;
      else if (p.kind === 'research' && user.role === 'researcher') S.research[p.id] = p;
      else if (p.kind === 'cohort' && user.role === 'researcher') S.cohort = p;
    }
  }
  // same error for an unknown user and a wrong password
  S.login = async function (username, password) {
    const user = S.users().find(u => u.u === username);
    let pk = null;
    try { if (user) pk = await GL.crypto.openPrivateKey(user, password); } catch (e) { pk = null; }
    if (!pk) throw new Error('credentials');
    await load(user, pk);
    S.user = user.u; S.role = user.role; S.guest = false; S.locked = false;
    startIdle();
    GL.audit.log('login', '-', Object.keys(S.real).length + Object.keys(S.research).length);
  };
  S.guestMode = function () { S.wipe(); S.user = null; S.role = 'guest'; S.guest = true; S.locked = false; };
  S.wipe = function () {
    S.real = {}; S.research = {}; S.cohort = null; S.scenarios = []; S.state = {};
    if (GL.ui && GL.ui.scrubber) GL.ui.scrubber.stopAll();
    const c = GL.views && GL.views.clinic; if (c) Object.assign(c.ui, { p: null, a: null, prev: null, whatif: null, wi: null, scrub: null });
    const r = GL.views && GL.views.research; if (r) Object.assign(r.ui, { vp: null });
    const d = typeof document !== 'undefined' && document.getElementById('drawer'); if (d) { d.innerHTML = ''; d.classList.remove('open'); }
  };
  S.lock = function () {
    if (!S.user || S.locked) return;
    GL.audit.log('lock', '-');
    S.wipe(); S.locked = true; stopIdle();
    GL.app.render();
  };
  S.unlock = async function (password) {
    const user = S.users().find(u => u.u === S.user);
    let pk = null;
    try { pk = await GL.crypto.openPrivateKey(user, password); } catch (e) { pk = null; }
    if (!pk) throw new Error('credentials');
    await load(user, pk);
    S.locked = false; startIdle();
    GL.audit.log('unlock', '-');
    await GL.app.render();
  };
  S.logout = async function () {
    if (S.user) await GL.audit.log('logout', '-');
    S.wipe(); S.user = null; S.role = null; S.guest = false; S.locked = false; stopIdle();
    if (location.hash !== '#/') location.hash = '#/'; else GL.app.route();
  };

  // ---- idle lock ---------------------------------------------------------------------------------------
  let lastReset = 0;
  function bump() { if (!S.user || S.locked) return; const now = Date.now(); if (now - lastReset < 1000) return; lastReset = now; startIdle(); }
  function startIdle() { stopIdle(); S.timer = setTimeout(() => S.lock(), S.idleMs); }
  function stopIdle() { if (S.timer) clearTimeout(S.timer); S.timer = null; }
  S.setIdle = function (ms) { S.idleMs = ms; if (S.user && !S.locked) startIdle(); };
  if (typeof window !== 'undefined') ['pointerdown', 'pointermove', 'keydown', 'wheel', 'touchstart'].forEach(ev => window.addEventListener(ev, bump, { passive: true }));

  // ---- rosters -------------------------------------------------------------------------------------------
  const toRoster = (rec) => ({ id: rec.id, type: 'real', status: 'completed', research: rec.kind === 'research', summary: rec.summary,
    label: rec.label && rec.label.primary, evolution: (rec.asof || []).map(a => [a.L].concat(a.p)) });
  // clinic: digital patients + (clinician) assigned real patients; research: digital + (researcher) pseudonyms
  S.roster = function (kind) {
    const idx = GL.loader.index(), dig = idx ? idx.roster : [];
    const real = kind === 'clinic' ? (S.can('clinic.real') ? Object.values(S.real) : []) : (S.role === 'researcher' ? Object.values(S.research) : []);
    return real.map(toRoster).sort((a, b) => a.id.localeCompare(b.id)).concat(dig);
  };
})(typeof window !== 'undefined' ? window : global);
