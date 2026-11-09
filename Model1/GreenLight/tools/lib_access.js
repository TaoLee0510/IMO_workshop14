/* GreenLight · tools/lib_access.js
 * Shared code of tools/seal.js and tools/users.js (Model1 plan section 13.6).
 *
 * Key hierarchy
 *   data key   : one random AES-256-GCM key per payload (clinical:<id>, research:<id>, cohort:research);
 *                kept by the data steward in Results/real/greenlight_admin/keystore.json (never in GreenLight)
 *   sealed file: GreenLight/data/sealed/<fid>.js = GLM1.addSealed(fid, iv, ciphertext+tag); fid is random
 *                (no patient id, no field name in the file)
 *   user       : RSA-OAEP-2048 key pair; the private key (PKCS#8) is encrypted with AES-256-GCM under a KEK
 *                derived from the password by PBKDF2-SHA256 (600 000 iterations, random 16-byte salt);
 *                every data key the user is entitled to is wrapped with the user's public key
 *                (RSA-OAEP-SHA256), listed under its opaque fid. The steward can therefore assign, revoke and
 *                re-seal without knowing any password; passwords are never stored.
 *   users.js   : GreenLight/data/users.js = window.GLM1_USERS (usernames, roles, salts, iterations,
 *                encrypted private keys, public keys, wrapped keys under fids)
 *   assignments: Results/real/greenlight_admin/assignments.json (clinician -> patient ids; steward only)
 * Entitlements: clinician -> clinical payloads of the assigned patients; researcher -> every research payload
 * and the research cohort aggregate; steward -> none (user list, audit export). Guests need no account.
 */
const fs = require('fs'), path = require('path'), crypto = require('crypto');

const GL = path.resolve(__dirname, '..');
const RES = path.resolve(GL, '..', 'Results', 'real');
const P = {
  gl: GL, stage: path.join(RES, 'greenlight_stage'), admin: path.join(RES, 'greenlight_admin'),
  sealed: path.join(GL, 'data', 'sealed'), users: path.join(GL, 'data', 'users.js'), patients: path.join(GL, 'data', 'patients')
};
P.keystore = path.join(P.admin, 'keystore.json');
P.assign = path.join(P.admin, 'assignments.json');
const ITER = 600000;
const ROLES = ['clinician', 'researcher', 'steward'];

const b64 = (b) => Buffer.from(b).toString('base64');
const unb64 = (s) => Buffer.from(s, 'base64');
const readJson = (f, dflt) => fs.existsSync(f) ? JSON.parse(fs.readFileSync(f, 'utf8')) : dflt;
const writeSecret = (f, obj) => { fs.mkdirSync(path.dirname(f), { recursive: true }); fs.writeFileSync(f, JSON.stringify(obj, null, 1), { mode: 0o600 }); fs.chmodSync(f, 0o600); };

// ---- AES-256-GCM (ciphertext || 16-byte tag, the WebCrypto layout) -------------------------------
function aesEncrypt(key, buf) {
  const iv = crypto.randomBytes(12), c = crypto.createCipheriv('aes-256-gcm', key, iv);
  const ct = Buffer.concat([c.update(buf), c.final(), c.getAuthTag()]);
  return { iv: b64(iv), ct: b64(ct) };
}
function aesDecrypt(key, ivB64, ctB64) {
  const ct = unb64(ctB64), d = crypto.createDecipheriv('aes-256-gcm', key, unb64(ivB64));
  d.setAuthTag(ct.subarray(ct.length - 16));
  return Buffer.concat([d.update(ct.subarray(0, ct.length - 16)), d.final()]);
}

// ---- stores ------------------------------------------------------------------------------------------
const loadKeystore = () => readJson(P.keystore, { schema: 'glm1-keystore-1', payloads: {} });
const saveKeystore = (ks) => { ks.updated = new Date().toISOString(); writeSecret(P.keystore, ks); };
const loadAssign = () => readJson(P.assign, {});
const saveAssign = (a) => writeSecret(P.assign, a);
function loadUsers() {
  if (!fs.existsSync(P.users)) return { schema: 'glm1-users-1', kdf: { name: 'PBKDF2', hash: 'SHA-256', iterations: ITER }, users: [] };
  const s = fs.readFileSync(P.users, 'utf8');
  return JSON.parse(s.slice(s.indexOf('=') + 1).trim().replace(/;\s*$/, ''));
}
function saveUsers(u) {
  u.updated = new Date().toISOString();
  fs.writeFileSync(P.users, '/* GreenLight · users (tools/users.js); no passwords, no patient ids */\nwindow.GLM1_USERS = ' + JSON.stringify(u) + ';\n');
}

// ---- payloads ------------------------------------------------------------------------------------------
// staged plaintext payloads (written by step G): { pid, file }
function stagePayloads() {
  const out = [];
  for (const kind of ['clinical', 'research']) {
    const dir = path.join(P.stage, kind);
    if (!fs.existsSync(dir)) continue;
    for (const f of fs.readdirSync(dir).filter(x => /\.json$/.test(x)).sort()) out.push({ pid: kind + ':' + f.replace(/\.json$/, ''), file: path.join(dir, f) });
  }
  const c = path.join(P.stage, 'cohort_research.json');
  if (fs.existsSync(c)) out.push({ pid: 'cohort:research', file: c });
  return out;
}
const newFid = () => crypto.randomBytes(12).toString('hex');
const sealedFile = (fid) => path.join(P.sealed, fid + '.js');
function writeSealed(fid, enc) {
  fs.mkdirSync(P.sealed, { recursive: true });
  fs.writeFileSync(sealedFile(fid), `GLM1.addSealed("${fid}","${enc.iv}","${enc.ct}");\n`);
}
function readSealed(fid) {
  const m = /GLM1\.addSealed\("([0-9a-f]+)","([^"]+)","([^"]+)"\)/.exec(fs.readFileSync(sealedFile(fid), 'utf8'));
  return { iv: m[2], ct: m[3] };
}
function sealPayload(ks, pid, buf) {
  if (!ks.payloads[pid]) ks.payloads[pid] = { fid: newFid(), key: b64(crypto.randomBytes(32)) };
  const e = ks.payloads[pid];
  writeSealed(e.fid, aesEncrypt(unb64(e.key), buf));
}
// new data key and file id for these payloads (revocation): decrypt with the old key, re-encrypt
function rotate(ks, pids) {
  for (const pid of pids) {
    const old = ks.payloads[pid];
    if (!old) continue;
    const s = readSealed(old.fid);
    const buf = aesDecrypt(unb64(old.key), s.iv, s.ct);
    ks.payloads[pid] = { fid: newFid(), key: b64(crypto.randomBytes(32)) };
    writeSealed(ks.payloads[pid].fid, aesEncrypt(unb64(ks.payloads[pid].key), buf));
    fs.rmSync(sealedFile(old.fid), { force: true });
  }
}

// ---- users -------------------------------------------------------------------------------------------------
function entitlements(user, assign, ks) {
  const ids = Object.keys(ks.payloads);
  // '*' = every real patient, including patients sealed later (granted with --all-patients / --patients all)
  if (user.role === 'clinician') return (assign[user.u] || []).includes('*') ? ids.filter(pid => pid.startsWith('clinical:')) :
    (assign[user.u] || []).map(id => 'clinical:' + id).filter(pid => ids.includes(pid));
  if (user.role === 'researcher') return ids.filter(pid => pid.startsWith('research:') || pid === 'cohort:research');
  return [];
}
const wrap = (pubB64, keyB64) => b64(crypto.publicEncrypt({ key: crypto.createPublicKey({ key: unb64(pubB64), format: 'der', type: 'spki' }),
  padding: crypto.constants.RSA_PKCS1_OAEP_PADDING, oaepHash: 'sha256' }, unb64(keyB64)));
function rewrapAll(users, assign, ks) {
  for (const u of users.users) u.files = entitlements(u, assign, ks).map(pid => ({ f: ks.payloads[pid].fid, k: wrap(u.pub, ks.payloads[pid].key) }));
}
function makeCredentials(password) {
  const { publicKey, privateKey } = crypto.generateKeyPairSync('rsa', { modulusLength: 2048, publicExponent: 0x10001,
    publicKeyEncoding: { type: 'spki', format: 'der' }, privateKeyEncoding: { type: 'pkcs8', format: 'der' } });
  const salt = crypto.randomBytes(16);
  const kek = crypto.pbkdf2Sync(Buffer.from(password, 'utf8'), salt, ITER, 32, 'sha256');
  return { salt: b64(salt), iter: ITER, ekey: aesEncrypt(kek, privateKey), pub: b64(publicKey) };
}
// password check (tests and CLI): decrypt the private key; throws on a wrong password
function openPrivateKey(user, password) {
  const kek = crypto.pbkdf2Sync(Buffer.from(password, 'utf8'), unb64(user.salt), user.iter, 32, 'sha256');
  return crypto.createPrivateKey({ key: aesDecrypt(kek, user.ekey.iv, user.ekey.ct), format: 'der', type: 'pkcs8' });
}
function openPayloads(user, password) {
  const pk = openPrivateKey(user, password);
  return user.files.map(e => {
    const key = crypto.privateDecrypt({ key: pk, padding: crypto.constants.RSA_PKCS1_OAEP_PADDING, oaepHash: 'sha256' }, unb64(e.k));
    const s = readSealed(e.f);
    return JSON.parse(aesDecrypt(key, s.iv, s.ct).toString('utf8'));
  });
}
const validName = (u) => /^[A-Za-z0-9._-]{3,32}$/.test(u);

module.exports = { P, ITER, ROLES, b64, unb64, aesEncrypt, aesDecrypt, loadKeystore, saveKeystore, loadAssign, saveAssign, loadUsers, saveUsers,
  stagePayloads, sealPayload, rotate, readSealed, sealedFile, entitlements, rewrapAll, makeCredentials, openPrivateKey, openPayloads, validName };
