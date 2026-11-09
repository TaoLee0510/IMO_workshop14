/* GreenLight · core/registry.js
 * Module registry: modules register themselves; the orchestrator resolves
 * run order from provides/consumes signals; the UI renders whatever is here.
 */
(function (root) {
  const GL = root.GL = root.GL || {};
  const R = GL.registry = { modules: [], enabled: {} };

  R.register = function (mod) {
    const errs = GL.schema.validateModule(mod);
    if (errs.length) throw new Error('module ' + (mod && mod.meta && mod.meta.id) + ' invalid: ' + errs.join('; '));
    const i = R.modules.findIndex(m => m.meta.id === mod.meta.id);
    if (i >= 0) R.modules[i] = mod; else R.modules.push(mod);
    if (!(mod.meta.id in R.enabled)) R.enabled[mod.meta.id] = true;
    return mod;
  };
  R.get = (id) => R.modules.find(m => m.meta.id === id);
  R.setEnabled = (id, on) => { R.enabled[id] = !!on; };
  R.active = () => R.modules.filter(m => R.enabled[m.meta.id]);

  // topological order by signal dependencies (Kahn); modules whose required
  // signals are missing still run (they fall back to defaults) but after providers.
  R.order = function () {
    const mods = R.active();
    const providerOf = {};
    mods.forEach(m => (m.provides.signals || []).forEach(s => { providerOf[s] = m.meta.id; }));
    const indeg = {}, adj = {};
    mods.forEach(m => { indeg[m.meta.id] = 0; adj[m.meta.id] = []; });
    mods.forEach(m => (m.consumes.signals || []).forEach(s => {
      const p = providerOf[s];
      if (p && p !== m.meta.id) { adj[p].push(m.meta.id); indeg[m.meta.id]++; }
    }));
    const q = mods.filter(m => indeg[m.meta.id] === 0).map(m => m.meta.id), out = [];
    while (q.length) {
      const id = q.shift(); out.push(id);
      adj[id].forEach(n => { if (--indeg[n] === 0) q.push(n); });
    }
    mods.forEach(m => { if (!out.includes(m.meta.id)) out.push(m.meta.id); }); // cycles: append
    return out.map(id => R.get(id));
  };

  // data coverage for a module given a snapshot: {variable: nObs}
  R.coverage = function (mod, snapshot) {
    const cov = {};
    (mod.consumes.variables || []).forEach(v => { cov[v] = snapshot.observations.filter(o => o.variable === v).length; });
    return cov;
  };
})(typeof window !== 'undefined' ? window : global);
