/* Drive the demo in headless Chrome over the DevTools protocol: load pages, ingest data,
 * switch language/theme, edit a plan, and capture screenshots + console errors.
 * Usage: node tools/browser_test.js <outDir> [chromePath]
 */
const { spawn } = require('child_process'); const fs = require('fs'); const path = require('path'); const http = require('http');
const out = process.argv[2] || path.join(__dirname, '..', '..', 'shots'); fs.mkdirSync(out, { recursive: true });
const chrome = process.argv[3] || '/Applications/Google Chrome.app/Contents/MacOS/Google Chrome';
const port = 9333, profile = '/tmp/gl-chrome-cdp-' + Date.now();
const url = 'file://' + path.resolve(__dirname, '..', 'index.html');
const proc = spawn(chrome, ['--headless=new', '--disable-gpu', '--no-sandbox', '--hide-scrollbars', '--window-size=1440,1500', '--remote-debugging-port=' + port, '--user-data-dir=' + profile, 'about:blank'], { stdio: 'ignore' });
const sleep = (ms) => new Promise(r => setTimeout(r, ms));
const getJson = (u) => new Promise((res, rej) => http.get(u, r => { let d = ''; r.on('data', c => d += c); r.on('end', () => res(JSON.parse(d))); }).on('error', rej));
(async () => {
  let targets = null;
  for (let i = 0; i < 150 && !(targets && targets.length); i++) { try { targets = await getJson(`http://localhost:${port}/json/list`); } catch (e) { targets = null; } if (!(targets && targets.length)) await sleep(300); }
  const page = targets.find(t => t.type === 'page');
  const ws = new WebSocket(page.webSocketDebuggerUrl);
  await new Promise(r => { ws.onopen = r; });
  let id = 0; const pending = {}; const errors = [];
  ws.onmessage = (m) => { const msg = JSON.parse(m.data); if (msg.id && pending[msg.id]) { pending[msg.id](msg); delete pending[msg.id]; }
    if (msg.method === 'Runtime.exceptionThrown') errors.push('EXC: ' + (msg.params.exceptionDetails.exception && msg.params.exceptionDetails.exception.description || msg.params.exceptionDetails.text));
    if (msg.method === 'Runtime.consoleAPICalled' && (msg.params.type === 'error' || msg.params.type === 'warning')) errors.push(msg.params.type.toUpperCase() + ': ' + msg.params.args.map(a => a.value || a.description).join(' ')); };
  const send = (method, params) => new Promise(res => { const i = ++id; pending[i] = res; ws.send(JSON.stringify({ id: i, method, params: params || {} })); });
  const evalJs = async (expr) => { const r = await send('Runtime.evaluate', { expression: expr, awaitPromise: true, returnByValue: true }); if (r.result && r.result.exceptionDetails) errors.push('EVAL: ' + JSON.stringify(r.result.exceptionDetails.exception && r.result.exceptionDetails.exception.description)); return r.result && r.result.result ? r.result.result.value : undefined; };
  const shot = async (name) => { await send('Emulation.setDeviceMetricsOverride', { width: 1440, height: 1500, deviceScaleFactor: 1, mobile: false }); const r = await send('Page.captureScreenshot', { format: 'png', captureBeyondViewport: true }); fs.writeFileSync(path.join(out, name + '.png'), Buffer.from(r.result.data, 'base64')); console.log('shot', name); };
  const waitFor = async (expr, ms) => { const t0 = Date.now(); while (Date.now() - t0 < (ms || 30000)) { if (await evalJs(expr)) return true; await sleep(250); } console.log('TIMEOUT waiting for', expr); return false; };
  await send('Runtime.enable'); await send('Page.enable');
  const nav = async (hash) => { await send('Page.navigate', { url: url + hash }); await sleep(800); await waitFor('!!(window.GL && GL.app && GL.app.state.patients.length)'); };

  // 0. help popovers (tile, SVG param button, plan card) in light and dark theme
  await nav("#/twin/P148"); await waitFor("!!GL.app.currentVersion()"); await sleep(300);
  await evalJs(`document.querySelector('[data-help="phenotype"]').click()`); await sleep(400); await shot("00a_help_phenotype");
  await evalJs(`GL.help.close(); document.querySelector('[data-help="param.R"]').dispatchEvent(new MouseEvent('click', { bubbles: true }))`); await sleep(400); await shot("00b_help_param");
  await evalJs(`GL.help.close(); document.getElementById('btnTheme').click()`); await sleep(300);
  await evalJs(`document.querySelector('[data-help="score"]').click()`); await sleep(400); await shot("00c_help_score_dark");
  await evalJs(`GL.help.close(); document.getElementById('btnTheme').click()`); await sleep(200);
  // 1. P148 pre-LD, then ingest twice, then scrub back on the timeline
  await shot('01_twin_P148_preLD');
  // 1a. recovery-phenotype outlook: one row per ranked plan; "Show" switches the forecast plan; recovery weight re-ranks
  const nRows = await evalJs("document.querySelectorAll('.pheno-table tbody tr').length"), nRanked = await evalJs('GL.app.rankedPlans(GL.app.currentVersion()).length');
  console.log('phenotype-outlook rows:', nRows, '/ ranked plans:', nRanked, nRows === nRanked ? 'OK' : 'MISMATCH');
  await evalJs(`document.querySelector('[data-help="pheno_by_plan"]').click()`); await sleep(400); await shot('00d_help_pheno_by_plan'); await evalJs('GL.help.close()');
  await evalJs(`document.querySelector('.pheno-table [data-show]:not([disabled])').click()`); await sleep(400);
  console.log('forecast plan after Show:', await evalJs('GL.app.views.twin.ui.planId'));
  const setRec = async (w) => { await evalJs(`(() => { const s = document.getElementById('wRec'); s.value = ${w}; s.dispatchEvent(new Event('change', { bubbles: true })); return true; })()`); await sleep(400); return evalJs("document.querySelector('.pheno-table tbody tr .plan-cell').firstChild.textContent.trim()"); };
  const top0 = await setRec(0); await shot('01b_twin_P148_recovery_weight_0'); const top2 = await setRec(2); console.log('top plan at recovery weight 0:', top0, '| at 2:', top2); await setRec(1.2);
  console.log('stored weights:', JSON.stringify(await evalJs('GL.store.state("P148").weights')));
  await evalJs('GL.app.ingestNext()'); await waitFor('GL.store.state("P148").horizonIndex === 1 && !!GL.app.currentVersion()'); await sleep(300); await shot('02_twin_P148_after_ingest');
  await evalJs('GL.app.ingestNext()'); await waitFor('GL.store.state("P148").horizonIndex === 2 && !!GL.app.currentVersion()'); await sleep(300);
  await evalJs(`document.querySelector('.tl-cp[data-idx="0"]').dispatchEvent(new PointerEvent('pointerdown', { bubbles: true, clientX: document.querySelector('.tl-cp[data-idx="0"]').getBoundingClientRect().left + 4 }))`);
  await waitFor('GL.store.state("P148").viewIndex === 0'); await sleep(400); await shot('02b_twin_P148_scrubbed_back');
  await evalJs(`GL.app.setView(2)`); await waitFor('GL.store.state("P148").viewIndex === 2'); await sleep(200);
  // 2. history diff
  await nav('#/history/P148'); await sleep(500); await shot('03_history_P148');
  // 3. planner: change plan and evaluate
  await nav('#/planner/P148'); await waitFor('!!(GL.app.views.planner.ui.result)'); await sleep(200);
  await evalJs(`(async () => { const f = document.querySelector('#f_ld_flu_frac'); f.value = 0.75; f.dispatchEvent(new Event('input', { bubbles: true })); const c = document.querySelector('#f_cytopenia_cd34_backup'); c.checked = true; c.dispatchEvent(new Event('change', { bubbles: true })); const b = document.querySelector('#f_cytopenia_cd34_boost'); b.value = 'if_persistent'; b.dispatchEvent(new Event('change', { bubbles: true })); await new Promise(r => setTimeout(r, 1500)); return true; })()`);
  await sleep(500); await shot('04_planner_P148_custom');
  await evalJs(`(async () => { document.querySelector('#btnSave').click(); await new Promise(r => setTimeout(r, 2500)); return true; })()`); await waitFor('location.hash.startsWith("#/twin")'); await sleep(500); await shot('05_twin_P148_with_custom');
  // 4. languages
  for (const l of ['zh-CN', 'ja', 'de']) { await evalJs(`GL.i18n.set('${l}')`); await sleep(400); await shot('06_twin_P148_' + l); }
  await evalJs(`GL.i18n.set('en')`);
  // 5. completed record: open P039 (full record), scrub to D+21, planner, dark theme
  await nav('#/twin/P039'); await waitFor('!!GL.app.currentVersion()'); await sleep(300); await shot('07_twin_P039_completed_D90');
  await evalJs(`GL.app.setView(4)`); await waitFor('GL.store.state("P039").viewIndex === 4 && !!GL.app.currentVersion()'); await sleep(400); await shot('07b_twin_P039_as_of_D21');
  await evalJs(`document.getElementById('btnTheme').click()`); await nav('#/twin/P010'); await waitFor('!!GL.app.currentVersion()'); await sleep(300); await shot('07c_twin_P010_dark');
  await evalJs(`document.getElementById('btnTheme').click()`);
  // 6. roster tabs; integration
  await nav('#/roster'); await sleep(500); await shot('08_roster_active');
  await evalJs(`document.querySelector('[data-tab="completed"]').click()`); await sleep(500); await shot('08b_roster_completed');
  await nav('#/integration'); await sleep(300); await shot('09_integration');
  // 7. reload persistence check
  await nav('#/history/P148'); await sleep(600); const nv = await evalJs('GL.store.listVersions("P148").length'); console.log('versions after reload for P148:', nv);
  console.log('console errors/warnings:', errors.length); errors.slice(0, 20).forEach(e => console.log('  ', e.slice(0, 300)));
  ws.close(); proc.kill(); await sleep(1500); try { fs.rmSync(profile, { recursive: true, force: true }); } catch (e) { /* profile still closing */ }
})().catch(e => { console.error(e); proc.kill(); process.exit(1); });
