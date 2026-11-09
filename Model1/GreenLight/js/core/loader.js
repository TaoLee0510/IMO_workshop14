/* GreenLight · core/loader.js
 * Data access. Public: data/model1_index.js (digital roster, aggregate model information) and the digital
 * patients data/patients/DP*.js (GLM1.addPatient). Sealed: data/sealed/<fid>.js (GLM1.addSealed, AES-GCM
 * ciphertext; opened by GL.session with the user's keys). Real patients are served from the decrypted
 * session only; there is no plaintext real-patient file. Script tags work from file:// (no fetch / CORS).
 */
(function (root) {
  const GL = root.GL = root.GL || {};
  const cache = {}, waiting = {}, sealedCache = {}, sealedWaiting = {};
  root.GLM1 = root.GLM1 || {};
  root.GLM1.addPatient = function (p) { if (!/^DP\d+$/.test(p.id)) return; cache[p.id] = p; (waiting[p.id] || []).forEach(fn => fn(p)); delete waiting[p.id]; };
  root.GLM1.addSealed = function (fid, iv, ct) { sealedCache[fid] = { iv, ct }; (sealedWaiting[fid] || []).forEach(fn => fn(sealedCache[fid])); delete sealedWaiting[fid]; };
  function script(src, onerror) { const s = document.createElement('script'); s.src = src; s.dataset.gl = '1'; s.onerror = onerror; document.head.appendChild(s); return s; }
  GL.loader = {
    schema: 'glm1-1.0',
    index: () => root.GLM1_INDEX,
    cached: (id) => cache[id] || (GL.session && (GL.session.real[id] || GL.session.research[id])),
    isDigital: (id) => /^DP\d+$/.test(id),
    patient(id) {
      const S = GL.session;
      if (S && S.real[id]) return Promise.resolve(S.real[id]);
      if (S && S.research[id]) return Promise.resolve(S.research[id]);
      if (!GL.loader.isDigital(id)) return Promise.reject(new Error('not available'));
      if (cache[id]) return Promise.resolve(cache[id]);
      return new Promise((res, rej) => {
        (waiting[id] = waiting[id] || []).push(res);
        if (waiting[id].length > 1) return;
        script('data/patients/' + encodeURIComponent(id) + '.js', () => { delete waiting[id]; rej(new Error('cannot load ' + id)); });
      });
    },
    sealed(fid) {
      if (sealedCache[fid]) return Promise.resolve(sealedCache[fid]);
      return new Promise((res, rej) => {
        (sealedWaiting[fid] = sealedWaiting[fid] || []).push(res);
        if (sealedWaiting[fid].length > 1) return;
        script('data/sealed/' + fid + '.js', () => { delete sealedWaiting[fid]; rej(new Error('cannot load a sealed file')); });
      });
    },
    dropSealed(fid) { delete sealedCache[fid]; document.querySelectorAll(`script[src="data/sealed/${fid}.js"]`).forEach(s => s.remove()); }
  };
})(typeof window !== 'undefined' ? window : global);
