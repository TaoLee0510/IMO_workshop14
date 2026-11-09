// Loads the browser-namespace core into node (no bundler). Usage: require('./load_core')
// Data files written as `window.X = ...` (data/aim4_pop.js, data/real_cohort.js) are evaluated
// with `window` bound to the node global.
const path = require('path'), fs = require('fs');
const files = ['core/util.js', 'core/ode.js', 'core/schema.js', 'core/registry.js', 'core/smc.js', 'modules/tumor_immune.js',
  'modules/aim4_core.js', 'modules/marrow4.js', 'core/plans.js', 'core/orchestrator.js'];
for (const f of ['data/aim4_pop.js', 'data/real_cohort.js']) {
  const p = path.join(__dirname, '..', f);
  if (fs.existsSync(p)) new Function('window', fs.readFileSync(p, 'utf8'))(global);
}
for (const f of files) { try { require(path.join(__dirname, '..', 'js', f)); } catch (e) { if (e.code !== 'MODULE_NOT_FOUND') throw e; } }
module.exports = global.GL;
