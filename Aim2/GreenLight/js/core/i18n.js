/* GreenLight · core/i18n.js
 * Tiny i18n: GL.i18n.t(key, params) with fallback to English, then the key.
 * Dictionaries live in js/i18n/dict.js (GL.I18N[lang]).
 */
(function (root) {
  const GL = root.GL = root.GL || {};
  const I = GL.i18n = { lang: 'en', listeners: [] };
  I.languages = [
    { code: 'en', label: 'English' }, { code: 'ja', label: '日本語' }, { code: 'de', label: 'Deutsch' }, { code: 'fr', label: 'Français' },
    { code: 'es', label: 'Español' }, { code: 'zh-CN', label: '简体中文' }, { code: 'zh-TW', label: '繁體中文' }
  ];
  I.t = function (key, params) {
    const dicts = GL.I18N || {};
    let s = (dicts[I.lang] && dicts[I.lang][key]);
    if (s === undefined) s = dicts.en && dicts.en[key];
    if (s === undefined) s = key;
    if (params) for (const k in params) s = s.split('{' + k + '}').join(params[k]);
    return s;
  };
  I.has = (key) => !!(GL.I18N && ((GL.I18N[I.lang] && GL.I18N[I.lang][key] !== undefined) || (GL.I18N.en && GL.I18N.en[key] !== undefined)));
  I.set = function (lang) {
    if (!I.languages.some(l => l.code === lang)) lang = 'en';
    I.lang = lang;
    try { localStorage.setItem('gl.lang', lang); } catch (e) { /* ignore */ }
    if (typeof document !== 'undefined') document.documentElement.lang = lang;
    I.listeners.forEach(fn => fn(lang));
  };
  I.init = function () {
    let lang = 'en';
    try { lang = localStorage.getItem('gl.lang') || ((typeof navigator !== 'undefined' && navigator.language) || 'en'); } catch (e) { /* ignore */ }
    if (lang.startsWith('zh')) lang = /tw|hk|hant/i.test(lang) ? 'zh-TW' : 'zh-CN';
    else lang = lang.slice(0, 2);
    I.set(lang);
  };
  I.onChange = (fn) => I.listeners.push(fn);
  I.locale = () => ({ en: 'en-US', ja: 'ja-JP', de: 'de-DE', fr: 'fr-FR', es: 'es-ES', 'zh-CN': 'zh-CN', 'zh-TW': 'zh-TW' }[I.lang] || 'en-US');
  I.fmtDate = (iso) => { try { return new Intl.DateTimeFormat(I.locale(), { dateStyle: 'medium', timeStyle: 'short' }).format(new Date(iso)); } catch (e) { return iso; } };
  I.fmtNum = (x, d) => { if (x === null || x === undefined || Number.isNaN(x)) return '–'; try { return new Intl.NumberFormat(I.locale(), { maximumFractionDigits: d === undefined ? 2 : d, minimumFractionDigits: 0 }).format(x); } catch (e) { return String(x); } };
})(typeof window !== 'undefined' ? window : global);
