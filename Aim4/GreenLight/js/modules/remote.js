/* GreenLight · modules/remote.js
 * RemoteModule: wraps any REST endpoint that implements the module contract
 * (docs/MODULE_API.md) so other teams' ODE models (R plumber, Python FastAPI,
 * Julia …) plug in without touching the UI.
 *   GET  {base}/meta                 -> { meta, consumes, provides }
 *   POST {base}/fit                  { snapshot, ctx }        -> state
 *   POST {base}/forecast             { state, plan, ctx }     -> forecast
 * Create with GL.modules.makeRemote({ baseUrl, token }) then await .load() and register.
 */
(function (root) {
  const GL = root.GL = root.GL || {};
  GL.modules = GL.modules || {};
  GL.modules.makeRemote = function (cfg) {
    const base = cfg.baseUrl.replace(/\/$/, '');
    const headers = () => { const h = { 'Content-Type': 'application/json' }; if (cfg.token) h.Authorization = 'Bearer ' + cfg.token; return h; };
    const post = async (path, body) => {
      const r = await fetch(base + path, { method: 'POST', headers: headers(), body: JSON.stringify(body) });
      if (!r.ok) throw new Error('remote module ' + path + ' HTTP ' + r.status);
      return r.json();
    };
    const mod = {
      meta: { id: cfg.id || 'remote', nameKey: cfg.name || 'Remote module', owner: cfg.owner || '', version: '?', kind: 'remote', placeholder: false, engine: base },
      consumes: { variables: [], planKeys: [], signals: [] },
      provides: { signals: [], params: [], outcomes: [], variables: [] },
      async load() {
        const r = await fetch(base + '/meta', { headers: headers() });
        if (!r.ok) throw new Error('remote module meta HTTP ' + r.status);
        const m = await r.json();
        Object.assign(mod.meta, m.meta || {}, { kind: 'remote', engine: base });
        mod.consumes = m.consumes || mod.consumes; mod.provides = m.provides || mod.provides;
        return mod;
      },
      async fit(snapshot, ctx) { const { snapshot: _s, ...c } = ctx || {}; return post('/fit', { snapshot, ctx: c }); },
      async forecast(state, plan, ctx) { const { snapshot, ...c } = ctx || {}; return post('/forecast', { state, plan, ctx: { ...c, snapshot } }); }
    };
    return mod;
  };
})(typeof window !== 'undefined' ? window : global);
