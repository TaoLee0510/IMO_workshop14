/* GreenLight · core/store.js
 * Browser storage under the access rules (Model1 plan 13.6, decision A6):
 *   localStorage (prefix glm1.*): UI preferences only (language, theme, filters, playback) and the demo
 *     state of digital patients; real-patient state lives in GL.session.state (memory).
 *   IndexedDB greenlight_model1: 'scenarios' = what-if scenarios of digital patients only (real-patient
 *     scenarios stay in memory and are cleared at logout / lock); 'audit' = hash-chained audit log.
 * Every access is guarded: private windows or blocked storage fall back to memory.
 */
(function (root) {
  const GL = root.GL = root.GL || {};
  const S = GL.store = {};
  const mem = {};
  const digital = (id) => /^DP\d+$/.test(id);
  S.get = function (k, dflt) {
    try { const v = localStorage.getItem('glm1.' + k); return v === null ? (k in mem ? mem[k] : dflt) : JSON.parse(v); } catch (e) { return k in mem ? mem[k] : dflt; }
  };
  S.set = function (k, v) { mem[k] = v; try { localStorage.setItem('glm1.' + k, JSON.stringify(v)); } catch (e) { /* memory only */ } };
  // per-patient state: data horizon of an in-treatment patient (demo ingest) and the viewed as-of day
  S.patient = function (p) {
    const st = (digital(p.id) ? S.get('p.' + p.id, null) : GL.session.state[p.id]) || {};
    if (p.status === 'active') { if (st.horizon === undefined || st.horizon === null || st.horizon < p.horizon) st.horizon = p.horizon; }
    else st.horizon = Infinity;
    return st;
  };
  S.savePatient = function (id, st) {
    if (digital(id)) S.set('p.' + id, Object.assign({}, st, { horizon: st.horizon === Infinity ? null : st.horizon }));
    else GL.session.state[id] = Object.assign({}, st);
  };
  S.lastPid = (id) => { if (id === undefined) return GL.session.state.__last || S.get('lastPid', null); if (digital(id)) S.set('lastPid', id); else GL.session.state.__last = id; };
  S.resetDemo = function () {
    try { Object.keys(localStorage).filter(k => k.startsWith('glm1.p.')).forEach(k => localStorage.removeItem(k)); } catch (e) { /* ignore */ }
    for (const k in mem) if (k.startsWith('p.')) delete mem[k];
  };

  // ---- IndexedDB -----------------------------------------------------------------------------------
  let dbp = null;
  S.db = function () {
    if (dbp) return dbp;
    dbp = new Promise((res) => {
      try {
        const r = indexedDB.open('greenlight_model1', 2);
        r.onupgradeneeded = () => { const d = r.result;
          if (!d.objectStoreNames.contains('scenarios')) d.createObjectStore('scenarios', { keyPath: 'key' });
          if (!d.objectStoreNames.contains('audit')) d.createObjectStore('audit', { keyPath: 'seq' }); };
        r.onsuccess = () => res(r.result); r.onerror = () => res(null);
      } catch (e) { res(null); }
    });
    return dbp;
  };
  const memScen = [];
  S.saveScenario = async function (sc) {
    sc.key = sc.patient + '|' + Date.now();
    if (!digital(sc.patient)) { GL.session.scenarios.push(sc); return sc; }
    const d = await S.db();
    if (!d) { memScen.push(sc); return sc; }
    await new Promise(res => { const tx = d.transaction('scenarios', 'readwrite'); tx.objectStore('scenarios').put(sc); tx.oncomplete = res; tx.onerror = res; });
    return sc;
  };
  S.scenarios = async function (patient) {
    if (!digital(patient)) return GL.session.scenarios.filter(s => s.patient === patient);
    const d = await S.db();
    if (!d) return memScen.filter(s => s.patient === patient);
    return new Promise(res => {
      const out = []; const r = d.transaction('scenarios', 'readonly').objectStore('scenarios').openCursor();
      r.onsuccess = () => { const c = r.result; if (c) { if (c.value.patient === patient) out.push(c.value); c.continue(); } else res(out); };
      r.onerror = () => res(out);
    });
  };
  // remove anything an earlier version may have stored about real patients
  S.purgeLegacy = async function () {
    try { Object.keys(localStorage).forEach(k => { if ((k.startsWith('glm1.p.') && !digital(k.slice(7))) || (k === 'glm1.lastPid' && !digital(S.get('lastPid', '')))) localStorage.removeItem(k); }); } catch (e) { /* ignore */ }
    const d = await S.db(); if (!d) return;
    await new Promise(res => { const tx = d.transaction('scenarios', 'readwrite'); const r = tx.objectStore('scenarios').openCursor();
      r.onsuccess = () => { const c = r.result; if (c) { if (!digital(c.value.patient)) c.delete(); c.continue(); } }; tx.oncomplete = res; tx.onerror = res; });
  };
})(typeof window !== 'undefined' ? window : global);
