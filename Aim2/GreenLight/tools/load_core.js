// Loads the browser-namespace core into node (no bundler). Usage: require('./load_core')
const path = require('path');
const files = ['core/util.js', 'core/ode.js', 'core/schema.js', 'core/registry.js', 'core/smc.js', 'modules/marrow.js', 'modules/tumor_immune.js', 'core/plans.js', 'core/orchestrator.js'];
for (const f of files) { try { require(path.join(__dirname, '..', 'js', f)); } catch (e) { if (e.code !== 'MODULE_NOT_FOUND') throw e; } }
module.exports = global.GL;
