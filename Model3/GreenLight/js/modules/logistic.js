/* GreenLight · modules/logistic.js
 * Example ML module (interface demonstration, docs/ML_API.md): ridge logistic regression of
 * non-quick recovery on log ANC and log PLT at the landmark (last value of the 3 days before),
 * one model per landmark (D0, D7, D10, D14, D21), coefficients from step G (in sample on the
 * real cohort; not a validated predictor - V1 gives the out-of-sample comparison).
 */
(function (root) {
  const GL = root.GL = root.GL || {};
  function last(obs, L) {
    const w = obs.filter(o => o[0] <= L && o[0] >= L - 3);
    const x = (w.length ? w : obs.filter(o => o[0] <= L)).sort((a, b) => a[0] - b[0]);
    return x.length ? Math.max(x[x.length - 1][1], 0.01) : null;
  }
  GL.registry.register({
    id: 'counts_logistic', version: 'example 1', type: 'ml', nameKey: 'mod.logistic', descKey: 'mod.logistic.d',
    source: 'step G (an_G_export.R g_logistic_module)', card: { outputs: [{ id: 'p_nonquick', labelKey: 'out.p_nonquick', unit: '%', direction: 'lower' }],
      inputs: ['log_anc', 'log_plt'] },
    predict(snap) {
      const ms = ((GL.loader.index() || {}).logistic || []).filter(m => m.L <= snap.asOf);
      if (!ms.length) return null;
      const m = ms[ms.length - 1], a = last(snap.obs.anc || [], m.L), p = last(snap.obs.plt || [], m.L);
      if (a === null || p === null) return null;
      const z = [(Math.log(a) - m.centre[0]) / m.scale[0], (Math.log(p) - m.centre[1]) / m.scale[1]];
      const eta = m.intercept + m.coef[0] * z[0] + m.coef[1] * z[1];
      return { p_nonquick: 1 / (1 + Math.exp(-eta)), landmark: m.L };
    }
  });
})(typeof window !== 'undefined' ? window : global);
