/* GreenLight · core/ode.js
 * Fixed-step RK4 integrator on a daily grid with impulse hooks.
 * sim = GL.ode.integrate({ t0, t1, dt, y0, rhs(t, y, dydt), onStep(t, y) })
 *  - rhs writes derivatives into dydt (preallocated Float64Array)
 *  - onStep is called at every grid time BEFORE the step (lets callers apply
 *    impulses / policy rules and record outputs)
 */
(function (root) {
  const GL = root.GL = root.GL || {};
  GL.ode = {};
  GL.ode.integrate = function (opt) {
    const n = opt.y0.length, dt = opt.dt;
    const y = Float64Array.from(opt.y0);
    const k1 = new Float64Array(n), k2 = new Float64Array(n), k3 = new Float64Array(n), k4 = new Float64Array(n), tmp = new Float64Array(n);
    const steps = Math.round((opt.t1 - opt.t0) / dt);
    let t = opt.t0;
    for (let s = 0; s <= steps; s++) {
      if (opt.onStep) opt.onStep(t, y, s);
      if (s === steps) break;
      opt.rhs(t, y, k1);
      for (let i = 0; i < n; i++) tmp[i] = y[i] + 0.5 * dt * k1[i];
      opt.rhs(t + 0.5 * dt, tmp, k2);
      for (let i = 0; i < n; i++) tmp[i] = y[i] + 0.5 * dt * k2[i];
      opt.rhs(t + 0.5 * dt, tmp, k3);
      for (let i = 0; i < n; i++) tmp[i] = y[i] + dt * k3[i];
      opt.rhs(t + dt, tmp, k4);
      for (let i = 0; i < n; i++) {
        y[i] += dt / 6 * (k1[i] + 2 * k2[i] + 2 * k3[i] + k4[i]);
        if (y[i] < 0) y[i] = 0; // all states are non-negative quantities
      }
      t = opt.t0 + (s + 1) * dt;
    }
    return y;
  };
})(typeof window !== 'undefined' ? window : global);
