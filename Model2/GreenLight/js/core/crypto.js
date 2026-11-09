/* GreenLight · core/crypto.js
 * WebCrypto side of the sealed bundle (tools/lib_access.js): PBKDF2-SHA256 KEK from the password,
 * AES-GCM decryption of the user's RSA-OAEP private key, RSA-OAEP unwrapping of the data keys and
 * AES-GCM decryption of the sealed payloads. Works on file:// in Chrome (secure context).
 */
(function (root) {
  const GL = root.GL = root.GL || {};
  const sub = () => { if (!(root.crypto && root.crypto.subtle)) throw new Error('WebCrypto unavailable: open GreenLight from file:// or http://localhost'); return root.crypto.subtle; };
  const unb64 = (s) => Uint8Array.from(atob(s), c => c.charCodeAt(0));
  GL.crypto = {
    available: () => !!(root.crypto && root.crypto.subtle),
    // returns the user's private key (CryptoKey, not extractable); throws on a wrong password
    async openPrivateKey(user, password) {
      const s = sub();
      const base = await s.importKey('raw', new TextEncoder().encode(password), 'PBKDF2', false, ['deriveKey']);
      const kek = await s.deriveKey({ name: 'PBKDF2', salt: unb64(user.salt), iterations: user.iter, hash: 'SHA-256' }, base, { name: 'AES-GCM', length: 256 }, false, ['decrypt']);
      const pkcs8 = await s.decrypt({ name: 'AES-GCM', iv: unb64(user.ekey.iv) }, kek, unb64(user.ekey.ct));
      return s.importKey('pkcs8', pkcs8, { name: 'RSA-OAEP', hash: 'SHA-256' }, false, ['decrypt']);
    },
    async openPayload(privateKey, wrappedB64, sealed) {
      const s = sub();
      const raw = await s.decrypt({ name: 'RSA-OAEP' }, privateKey, unb64(wrappedB64));
      const key = await s.importKey('raw', raw, { name: 'AES-GCM' }, false, ['decrypt']);
      const pt = await s.decrypt({ name: 'AES-GCM', iv: unb64(sealed.iv) }, key, unb64(sealed.ct));
      return JSON.parse(new TextDecoder().decode(pt));
    },
    async sha256hex(str) {
      const h = await sub().digest('SHA-256', new TextEncoder().encode(str));
      return Array.from(new Uint8Array(h), b => b.toString(16).padStart(2, '0')).join('');
    }
  };
})(typeof window !== 'undefined' ? window : global);
