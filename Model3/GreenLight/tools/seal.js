/* GreenLight · tools/seal.js — encrypt the staged real-patient payloads (access rules of Model1 plan 13.6, carried over by Model2 plan 6.5 and Model3).
 * Reads Results/real/greenlight_stage (written by step G), keeps one AES-256-GCM data key and one random
 * file id per payload in Results/real/greenlight_admin/keystore.json (created on first use; --rotate makes
 * new keys and file ids for every payload), writes GreenLight/data/sealed/<fid>.js, re-wraps the
 * entitled data keys for every user (public keys; no password needed), removes sealed files and keystore
 * entries without a staged payload, and deletes the real-patient plaintext GreenLight/data/patients/P*.js.
 * Usage: node tools/seal.js [--rotate]
 */
const fs = require('fs'), path = require('path');
const A = require('./lib_access.js');

const rotateAll = process.argv.includes('--rotate');
const staged = A.stagePayloads();
if (!staged.length) { console.error('nothing staged in ' + A.P.stage + ' — run step G (g_export) first'); process.exit(1); }
const ks = A.loadKeystore();
if (rotateAll) ks.payloads = {};
for (const s of staged) A.sealPayload(ks, s.pid, fs.readFileSync(s.file));
const keep = new Set(staged.map(s => s.pid));
for (const pid of Object.keys(ks.payloads)) if (!keep.has(pid)) delete ks.payloads[pid];
const fids = new Set(Object.values(ks.payloads).map(e => e.fid));
let removed = 0;
for (const f of fs.readdirSync(A.P.sealed)) if (/\.js$/.test(f) && !fids.has(f.replace(/\.js$/, ''))) { fs.rmSync(path.join(A.P.sealed, f)); removed++; }
A.saveKeystore(ks);
const users = A.loadUsers();
A.rewrapAll(users, A.loadAssign(), ks);
A.saveUsers(users);
let plain = 0;
if (fs.existsSync(A.P.patients)) for (const f of fs.readdirSync(A.P.patients)) if (/^P\d+\.js$/.test(f)) { fs.rmSync(path.join(A.P.patients, f)); plain++; }
const n = (k) => staged.filter(s => s.pid.startsWith(k)).length;
console.log(`sealed ${staged.length} payloads (clinical ${n('clinical:')}, research ${n('research:')}, cohort ${n('cohort:')})${rotateAll ? ' with new keys' : ''}; ` +
  `${removed} stale sealed file(s) removed; ${plain} plaintext file(s) deleted from GreenLight; ${users.users.length} user(s) re-wrapped`);
