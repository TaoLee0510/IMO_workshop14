/* GreenLight · core/registry.js
 * Module registry. Every module (mechanistic model, ML predictor, clinical score, remote
 * service) registers a model card (docs/MODULE_API.md, docs/ML_API.md):
 *   { id, version, type: 'mechanistic'|'ml'|'score', nameKey, descKey, source,
 *     card: { equations: [{id, titleKey, latex(spec), descKey, when(spec)}], symbols: [{sym, id, key, unit}],
 *             params: [...], outputs: [{id, labelKey, unit, direction}], validation: {...} },
 *     predict(snapshot) -> { p: {quick, second_wave, delayed} | null, score?, forecast? } }
 * The researcher portal is generated from the cards; the clinician dashboard shows every
 * module's prediction beside the digital twin.
 * Snapshot (input of predict): { patient, asOf, obs: {anc: [[day, value, cens]], ...}, summary, design }.
 */
(function (root) {
  const GL = root.GL = root.GL || {};
  const R = GL.registry = { modules: [] };
  R.register = function (m) {
    if (!m || !m.id || !m.type) throw new Error('module needs id and type');
    R.modules = R.modules.filter(x => x.id !== m.id).concat([m]);
    return m;
  };
  R.get = (id) => R.modules.find(m => m.id === id);
  R.byType = (type) => R.modules.filter(m => m.type === type);
  // snapshot of a patient record as of a day (observations up to that day only)
  R.snapshot = function (p, asOf) {
    const obs = {};
    for (const v in p.obs) obs[v] = p.obs[v].filter(o => o[0] <= asOf && o[3] === 1).map(o => [o[0], o[1], o[2]]);
    return { patient: p.id, asOf, obs, summary: p.summary, design: p.design };
  };
  R.predictAll = function (p, asOf) {
    const snap = R.snapshot(p, asOf), out = {};
    R.modules.filter(m => typeof m.predict === 'function').forEach(m => { try { out[m.id] = m.predict(snap, p); } catch (e) { out[m.id] = { error: String(e) }; } });
    return out;
  };
})(typeof window !== 'undefined' ? window : global);
