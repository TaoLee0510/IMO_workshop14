/* GreenLight · modules/remote.js
 * Remote module adapter (docs/MODULE_API.md): a model served elsewhere (R plumber or Python
 * FastAPI, api/*_stub.*) is registered with a base URL; GET {url}/card returns its model card,
 * POST {url}/predict takes a snapshot and returns { p: {quick, second_wave, delayed} } and/or
 * { p_nonquick } / { score }. Remote predictions are asynchronous and shown on request.
 */
(function (root) {
  const GL = root.GL = root.GL || {};
  GL.remote = {
    list: () => GL.store.get('remotes', []),
    add(url, id) { const r = GL.remote.list().filter(x => x.url !== url).concat([{ url, id: id || 'remote_' + GL.remote.list().length }]); GL.store.set('remotes', r); GL.remote.registerAll(); },
    remove(url) { GL.store.set('remotes', GL.remote.list().filter(x => x.url !== url)); GL.registry.modules = GL.registry.modules.filter(m => m.url !== url); },
    registerAll() {
      GL.remote.list().forEach(r => GL.registry.register({
        id: r.id, url: r.url, version: 'remote', type: 'ml', nameKey: 'mod.remote', descKey: 'mod.remote.d', source: r.url, card: { outputs: [] },
        async predictAsync(snap) {
          const res = await fetch(r.url.replace(/\/$/, '') + '/predict', { method: 'POST', headers: { 'Content-Type': 'application/json' }, body: JSON.stringify(snap) });
          if (!res.ok) throw new Error('HTTP ' + res.status);
          return res.json();
        },
        async fetchCard() { const res = await fetch(r.url.replace(/\/$/, '') + '/card'); return res.json(); }
      }));
    }
  };
})(typeof window !== 'undefined' ? window : global);
