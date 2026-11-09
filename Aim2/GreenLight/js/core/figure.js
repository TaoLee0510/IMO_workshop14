/* GreenLight · core/figure.js — concept illustrations (inline SVG).
 * variants: 'pair' (patient + data-built counterpart), 'model' (counterpart only, risk-tinted),
 * 'planner' (counterpart with branching treatment paths), 'history' (versions sharpening over time),
 * 'integration' (counterpart wired to modules). */
(function (root) {
  const GL = root.GL = root.GL || {};
  const F = GL.figure = {};
  const SIL = 'M50 2a11 11 0 1 1-.1 0zM42 27c-8 2-15 6-19 13-5 9-7 24-8 40-1 8-1 15 1 19l7-1c1-11 3-24 6-34l1 31c0 12 1 20 3 30 1 18 2 31 2 45 0 6-1 12 0 17h11c0-8 0-17 0-25 0-16 2-30 4-42 2 12 4 26 4 42 0 8 0 17 0 25h11c1-5 0-11 0-17 0-14 1-27 2-45 2-10 3-18 3-30l1-31c3 10 5 23 6 34l7 1c2-4 2-11 1-19-1-16-3-31-8-40-4-7-11-11-19-13-3 3-7 4-10 4s-7-1-10-4z';
  const NODES = [[50, 13], [50, 44], [38, 66], [62, 66], [50, 84], [44, 112], [56, 112], [46, 150], [54, 150], [24, 60], [76, 60]];
  const LINKS = [[0, 1], [1, 2], [1, 3], [2, 4], [3, 4], [4, 5], [4, 6], [5, 7], [6, 8], [2, 9], [3, 10], [5, 6]];
  const RISK = { green: '#22d3ee', amber: '#fbbf24', red: '#f87171' };
  function defs(id, glow) {
    return `<defs>
  <path id="${id}-sil" d="${SIL}"/>
  <linearGradient id="${id}-body" x1="0" y1="0" x2="0" y2="1"><stop offset="0" stop-color="#4ade80" stop-opacity=".95"/><stop offset="1" stop-color="#15803d" stop-opacity=".55"/></linearGradient>
  <pattern id="${id}-dots" width="5" height="5" patternUnits="userSpaceOnUse"><circle cx="2.5" cy="2.5" r="1" class="glf-dot"/></pattern>
  <pattern id="${id}-dots-s" width="7" height="7" patternUnits="userSpaceOnUse"><circle cx="3.5" cy="3.5" r=".9" class="glf-dot"/></pattern>
  <radialGradient id="${id}-glow" cx=".5" cy=".5" r=".5"><stop offset="0" stop-color="${glow}" stop-opacity=".5"/><stop offset="1" stop-color="${glow}" stop-opacity="0"/></radialGradient>
  <clipPath id="${id}-clip"><use href="#${id}-sil"/></clipPath>
</defs>`;
  }
  const nodes = (scale) => NODES.map(([x, y], i) => `<circle class="glf-node" cx="${x}" cy="${y}" r="${(i === 0 || i === 4 ? 2.6 : 1.9) * (scale || 1)}"/>`).join('');
  const links = () => LINKS.map(([a, b]) => `<line class="glf-link" x1="${NODES[a][0]}" y1="${NODES[a][1]}" x2="${NODES[b][0]}" y2="${NODES[b][1]}"/>`).join('');
  const scan = () => Array.from({ length: 24 }, (_, i) => `<line class="glf-scan" x1="0" x2="100" y1="${8 + i * 7.8}" y2="${8 + i * 7.8}"/>`).join('');
  const human = (id, tx, ty, s) => `<g transform="translate(${tx},${ty}) scale(${s || 1})"><use href="#${id}-sil" fill="url(#${id}-body)"/><g clip-path="url(#${id}-clip)">${scan()}</g></g>`;
  const digital = (id, tx, ty, s, dense) => `<g transform="translate(${tx},${ty}) scale(${s || 1})"><use href="#${id}-sil" fill="url(#${id}-${dense === false ? 'dots-s' : 'dots'})"/><use href="#${id}-sil" class="glf-outline"/>${links()}${nodes()}</g>`;

  F.svg = function (opt) {
    opt = opt || {};
    const id = opt.id || 'glf', variant = opt.variant || 'pair', w = opt.width || 230, h = opt.height || 210;
    const glow = RISK[opt.risk] || '#22d3ee';
    let body = '';
    if (variant === 'pair') {
      body = `<ellipse cx="160" cy="110" rx="70" ry="95" fill="url(#${id}-glow)"/>${human(id, 14, 10)}${digital(id, 110, 10)}
<g class="glf-flows"><path d="M74 60 C 95 40, 125 40, 148 56"/><path d="M78 100 C 100 90, 120 90, 148 98"/><path d="M76 140 C 98 160, 126 160, 150 140"/></g>`;
    } else if (variant === 'model') {
      // the counterpart alone, centred; orbit rings = continuous updating; glow tinted by risk
      body = `<ellipse cx="115" cy="105" rx="80" ry="100" fill="url(#${id}-glow)"/>
<g class="glf-orbit"><ellipse cx="115" cy="105" rx="92" ry="34" transform="rotate(-18 115 105)"/><ellipse cx="115" cy="105" rx="92" ry="34" transform="rotate(18 115 105)"/></g>
${digital(id, 65, 10)}
<g class="glf-pulse" style="--glf-c:${glow}"><circle cx="115" cy="94" r="5"/><circle cx="115" cy="94" r="5" class="glf-pulse-ring"/></g>
${opt.completed ? `<g class="glf-seal"><circle cx="186" cy="34" r="14"/><path d="M179 34l5 5 9-10"/></g>` : ''}`;
    } else if (variant === 'planner') {
      body = `<ellipse cx="70" cy="105" rx="60" ry="95" fill="url(#${id}-glow)"/>${digital(id, 18, 10)}
<g class="glf-paths"><path d="M118 100 C 150 100, 160 40, 212 40" class="p1"/><path d="M118 100 C 150 100, 160 100, 212 100" class="p2"/><path d="M118 100 C 150 100, 160 160, 212 160" class="p3"/>
<circle cx="212" cy="40" r="5" class="n1"/><circle cx="212" cy="100" r="5" class="n2"/><circle cx="212" cy="160" r="5" class="n3"/><text x="196" y="30" class="glf-rank">#1</text><text x="196" y="90" class="glf-rank">#2</text><text x="196" y="150" class="glf-rank">#3</text></g>`;
    } else if (variant === 'history') {
      body = `<ellipse cx="180" cy="110" rx="60" ry="90" fill="url(#${id}-glow)"/>
<g opacity=".35">${digital(id, 4, 60, .55, false)}</g><g opacity=".6">${digital(id, 68, 40, .72, false)}</g>${digital(id, 150, 14, .95)}
<g class="glf-flows"><path d="M44 170 C 70 180, 90 178, 112 170"/><path d="M128 160 C 150 170, 160 168, 178 160"/></g>
<g class="glf-vers"><text x="30" y="200">v1</text><text x="96" y="200">v2</text><text x="180" y="200">v3</text></g>`;
    } else if (variant === 'integration') {
      body = `<ellipse cx="115" cy="105" rx="70" ry="95" fill="url(#${id}-glow)"/>${digital(id, 65, 10)}
<g class="glf-mods"><rect x="6" y="40" width="46" height="26" rx="6"/><text x="29" y="57">ODE</text><rect x="6" y="150" width="46" height="26" rx="6"/><text x="29" y="167">CSV</text><rect x="178" y="40" width="46" height="26" rx="6"/><text x="201" y="57">ODE</text><rect x="178" y="150" width="46" height="26" rx="6"/><text x="201" y="167">REST</text></g>
<g class="glf-flows"><path d="M52 53 C 70 53, 80 70, 92 76"/><path d="M52 163 C 70 163, 80 140, 92 130"/><path d="M178 53 C 160 53, 150 70, 138 76"/><path d="M178 163 C 160 163, 150 140, 138 130"/></g>`;
    }
    return `<svg class="glf glf-${variant} ${opt.cls || ''}" viewBox="0 0 230 210" width="${w}" height="${h}" aria-hidden="true">${defs(id, glow)}${body}</svg>`;
  };
})(typeof window !== 'undefined' ? window : global);
