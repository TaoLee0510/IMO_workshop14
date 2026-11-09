/* GreenLight · core/audit.js
 * Hash-chained audit log in IndexedDB (greenlight_model2 / audit): every entry
 *   { seq, ts, user, role, subject (patient id, research pseudonym or '-'), action, detail, prev, hash }
 * with hash = SHA-256(prev | seq | ts | user | role | subject | action | detail). Changing or deleting an
 * entry breaks the chain (verify). Actions: login, logout, lock, unlock, open, scrub, play, what-if, export.
 * The data steward exports the log as JSON. Limitation: the log is local; it can be deleted or exported
 * but not silently altered (docs/PRIVACY.md).
 */
(function (root) {
  const GL = root.GL = root.GL || {};
  const Au = GL.audit = { last: null, seq: 0, ready: null, memory: [] };
  const fields = (e) => [e.prev, e.seq, e.ts, e.user, e.role, e.subject, e.action, e.detail].join('|');
  async function tail() {
    const d = await GL.store.db();
    if (!d) return null;
    return new Promise(res => { const r = d.transaction('audit', 'readonly').objectStore('audit').openCursor(null, 'prev');
      r.onsuccess = () => res(r.result ? r.result.value : null); r.onerror = () => res(null); });
  }
  Au.init = async function () { Au.ready = tail().then(t => { Au.last = t ? t.hash : 'genesis'; Au.seq = t ? t.seq : 0; }); return Au.ready; };
  let q = Promise.resolve();
  Au.log = function (action, subject, detail) {
    const s = GL.session;
    const e = { ts: new Date().toISOString(), user: s.user || 'guest', role: s.role || 'guest', subject: subject || '-', action, detail: detail === undefined ? '' : String(detail) };
    q = q.then(async () => {
      await Au.ready;
      e.seq = ++Au.seq; e.prev = Au.last;
      e.hash = await GL.crypto.sha256hex(fields(e));
      Au.last = e.hash;
      const d = await GL.store.db();
      if (!d) { Au.memory.push(e); return; }
      await new Promise(res => { const tx = d.transaction('audit', 'readwrite'); tx.objectStore('audit').put(e); tx.oncomplete = res; tx.onerror = res; });
    });
    return q;
  };
  Au.all = async function () {
    await q;
    const d = await GL.store.db();
    if (!d) return Au.memory.slice();
    return new Promise(res => { const out = []; const r = d.transaction('audit', 'readonly').objectStore('audit').openCursor();
      r.onsuccess = () => { const c = r.result; if (c) { out.push(c.value); c.continue(); } else res(out); }; r.onerror = () => res(out); });
  };
  Au.verify = async function (entries) {
    let prev = 'genesis';
    for (const e of entries) {
      if (e.prev !== prev || await GL.crypto.sha256hex(fields(e)) !== e.hash) return { ok: false, at: e.seq };
      prev = e.hash;
    }
    return { ok: true, n: entries.length };
  };
})(typeof window !== 'undefined' ? window : global);
