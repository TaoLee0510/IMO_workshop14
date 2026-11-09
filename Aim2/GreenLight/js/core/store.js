/* GreenLight · core/store.js
 * Twin version store: immutable TwinVersions per patient, patient UI state
 * (current checkpoint, active plan, custom plans). Small state lives in
 * localStorage; versions (large) live in IndexedDB and are loaded per patient.
 * Every fit becomes a new numbered version; nothing is overwritten.
 */
(function (root) {
  const GL = root.GL = root.GL || {};
  const ST = GL.store = { versions: {}, patientState: {}, customPlans: {}, loaded: {}, _db: null };
  const KEY = 'gl.state.v1', DB = 'greenlight', STORE = 'versions';

  function idb() {
    if (ST._db) return ST._db;
    ST._db = new Promise((resolve) => {
      try {
        const req = indexedDB.open(DB, 1);
        req.onupgradeneeded = () => { const db = req.result; if (!db.objectStoreNames.contains(STORE)) db.createObjectStore(STORE, { keyPath: 'key' }).createIndex('pid', 'patient_id'); };
        req.onsuccess = () => resolve(req.result);
        req.onerror = () => resolve(null);
      } catch (e) { resolve(null); }
    });
    return ST._db;
  }
  async function idbPut(v) { const db = await idb(); if (!db) return; await new Promise((res) => { const tx = db.transaction(STORE, 'readwrite'); tx.objectStore(STORE).put({ key: v.patient_id + '#' + v.version, patient_id: v.patient_id, version: v }); tx.oncomplete = res; tx.onerror = res; }); }
  async function idbAll(pid) { const db = await idb(); if (!db) return []; return new Promise((res) => { const tx = db.transaction(STORE, 'readonly'); const req = tx.objectStore(STORE).index('pid').getAll(pid); req.onsuccess = () => res(req.result.map(r => r.version)); req.onerror = () => res([]); }); }
  async function idbClear() { const db = await idb(); if (!db) return; await new Promise((res) => { const tx = db.transaction(STORE, 'readwrite'); tx.objectStore(STORE).clear(); tx.oncomplete = res; tx.onerror = res; }); }

  ST.load = function () {
    try { const raw = localStorage.getItem(KEY); if (raw) { const o = JSON.parse(raw); ST.patientState = o.patientState || {}; ST.customPlans = o.customPlans || {}; } }
    catch (e) { ST.patientState = {}; ST.customPlans = {}; }
  };
  ST.save = function () { try { localStorage.setItem(KEY, JSON.stringify({ patientState: ST.patientState, customPlans: ST.customPlans })); } catch (e) { /* ignore */ } };
  ST.loadPatient = async function (pid) {
    if (ST.loaded[pid]) return;
    const vs = await idbAll(pid);
    ST.versions[pid] = vs.sort((a, b) => a.version - b.version); ST.loaded[pid] = true;
  };
  ST.reset = async function () { ST.versions = {}; ST.patientState = {}; ST.customPlans = {}; ST.loaded = {}; try { localStorage.removeItem(KEY); } catch (e) { /* ignore */ } await idbClear(); };

  ST.listVersions = (pid) => (ST.versions[pid] || []).slice().sort((a, b) => a.version - b.version);
  ST.getVersion = (pid, n) => (ST.versions[pid] || []).find(v => v.version === n) || null;
  ST.latest = (pid) => { const vs = ST.listVersions(pid); return vs.length ? vs[vs.length - 1] : null; };
  ST.versionAt = (pid, asOfDay) => { const vs = ST.listVersions(pid).filter(v => v.asOfDay === asOfDay); return vs.length ? vs[vs.length - 1] : null; };
  ST.addVersion = function (version) {
    const list = ST.versions[version.patient_id] = ST.versions[version.patient_id] || [];
    version.version = list.length ? Math.max(...list.map(v => v.version)) + 1 : 1;
    Object.freeze(version);
    list.push(version);
    idbPut(version);
    return version;
  };
  // horizonIndex: last checkpoint whose data have arrived (data horizon); viewIndex: checkpoint being displayed (<= horizon)
  ST.state = function (pid, patient) {
    if (!ST.patientState[pid]) {
      const init = patient ? patient.demo.initialCheckpoint : 0;
      ST.patientState[pid] = { horizonIndex: init, viewIndex: init, activePlanId: 'soc', weights: Object.assign({}, GL.plans.defaultWeights) };
    }
    const st = ST.patientState[pid];
    st.weights = GL.plans.normWeights(st.weights); // older saved states lack the recovery weight
    if (st.horizonIndex === undefined) { st.horizonIndex = st.checkpointIndex !== undefined ? st.checkpointIndex : (patient ? patient.demo.initialCheckpoint : 0); }
    if (st.viewIndex === undefined || st.viewIndex > st.horizonIndex) st.viewIndex = st.horizonIndex;
    return st;
  };
  ST.setState = function (pid, patch) { Object.assign(ST.state(pid), patch); ST.save(); };
  ST.addCustomPlan = function (pid, plan) { const l = ST.customPlans[pid] = ST.customPlans[pid] || []; l.push(plan); ST.save(); return plan; };
  ST.customPlansFor = (pid) => ST.customPlans[pid] || [];
  ST.removeCustomPlan = function (pid, id) { ST.customPlans[pid] = (ST.customPlans[pid] || []).filter(p => p.id !== id); ST.save(); };
})(typeof window !== 'undefined' ? window : global);
