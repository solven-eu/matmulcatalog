// Browser-side filter/display for catalog.json.
//
// Loads the manifest, wires up controls, re-renders on any change.
// No framework — vanilla DOM updates so this works on GitHub Pages
// without a build step.

const SCHEMES_GITHUB_BASE =
  "https://github.com/solven-eu/matmulcatalog/blob/master/src/main/resources/schemes/";
// Raw-content base for lazy-fetching an individual scheme's full JSON (its
// structured `lineage` DAG with hashes) when the entry modal opens — keeps
// catalog.json lean (the DAG is only needed for the one open scheme).
const SCHEMES_RAW_BASE =
  "https://raw.githubusercontent.com/solven-eu/matmulcatalog/master/src/main/resources/schemes/";

const REFERENCES_MD_URL =
  "https://github.com/solven-eu/matmulcatalog/blob/master/REFERENCES.md";

/**
 * Maps a registry source key → the anchor in REFERENCES.md. Each anchor
 * leads to the full BibTeX + annotation + local PDF link for that entry.
 * Returns null when no REFERENCES.md entry exists for this source.
 */
function referencesMdAnchor(sourceKey) {
  if (!sourceKey) return null;
  if (sourceKey.startsWith("Strassen 1969"))                       return "1-strassen69";
  if (sourceKey.startsWith("Hopcroft-Kerr 1971") ||
      sourceKey.startsWith("Hopcroft–Kerr 1971"))                  return "2-hk71";
  if (sourceKey.startsWith("Winograd 1971"))                       return "3-winograd71";
  if (sourceKey.startsWith("Laderman 1976"))                       return "4-lad76";
  if (sourceKey.startsWith("Pan 1978"))                            return "5-pan78";
  if (sourceKey.startsWith("Bini 1979"))                           return "6-bini79";
  if (sourceKey.startsWith("Schönhage 1981"))                      return "7-sch81";
  if (sourceKey.startsWith("Coppersmith") || sourceKey.startsWith("Winograd 1990")) return "8-cw90";
  if (sourceKey.startsWith("Bläser 2003"))                         return "9-blaser03";
  if (sourceKey.startsWith("Drevet") ||
      sourceKey.startsWith("DIS09") ||
      sourceKey === "DIS09 (Drevet–Islam–Schost 2009)")            return "10-drisc09";
  // Pan's trilinear-aggregation FRAMEWORK (Pan 1978) — when used as the
  // primary citation rather than DIS09's specific TA-derived ranks.
  if (sourceKey.startsWith("Pan ") && sourceKey.includes("trilinear")) return "5-pan78";
  if (sourceKey.startsWith("Smirnov 2013"))                        return "11-smirnov";
  if (sourceKey === "AlphaTensor 2022")                            return "12-alphatensor";
  if (sourceKey.startsWith("Williams 2024"))                       return "13-williams2024";
  if (sourceKey === "AlphaEvolve 2025")                            return "14-alphaevolve";
  if (sourceKey.startsWith("Wang 2026"))                           return "15-wang26";
  if (sourceKey.startsWith("Phials"))                              return "16-phialsbasement";
  if (sourceKey === "fmm-lille catalog")                           return "17-fmm-lille";
  if (sourceKey === "Perminov (FastMatrixMultiplication)")         return "18-perminov";
  if (sourceKey === "Perminov 2026 (serendipitous)")               return "86-perminov-2026-serendipitous";
  if (sourceKey === "Kaporin 2024")                                return "87-kaporin-2024-brent";
  if (sourceKey.startsWith("Schwartz") && sourceKey.includes("2025")) return "19-schwartz-zwecher25";
  if (sourceKey.startsWith("Drevet") && sourceKey.includes("2010")) return "20-drevet-schost-poster";
  if (sourceKey.startsWith("Sedoglavic 2017"))                     return "21-sedoglavic17";
  if (sourceKey.startsWith("Waksman 1970"))                        return "22-waksman70";
  if (sourceKey.startsWith("Makarov 1987"))                        return "23-makarov87";
  if (sourceKey.startsWith("Probert"))                             return "24-probert-fischer-80";
  if (sourceKey.startsWith("Smith 2002"))                          return "25-smith-2002";
  if (sourceKey.startsWith("Mezzarobba"))                          return "26-mezzarobba-2007";
  if (sourceKey.startsWith("Mårtensson"))                          return "27-mss25";
  if (sourceKey.startsWith("Rosowski 2019") ||
      sourceKey.startsWith("Rosowski 2020"))                       return "28-rosowski19";
  if (sourceKey.startsWith("Makarov 1986"))                        return "29-makarov86";
  return null;
}

/**
 * Expand a cited-bound family entry (symbolic ⟨N,2,P⟩ + rank_formula)
 * into concrete ⟨n,m,p⟩ entries. Each declared variable sweeps its
 * inclusive range; the formula is evaluated for every combination.
 *
 * Concrete (non-family) entries pass through unchanged.
 *
 * SECURITY: rank_formula is evaluated with `Function(...)` in a
 * sandboxed scope exposing only `max`, `min`, and the declared vars.
 * It runs over our own emitted bounds JSON — not user input — so the
 * eval surface is acceptable.
 */
function expandFamilyEntry(c) {
  if (!c.format_symbolic) return [c];
  const vars = c.vars || [];
  const ranges = c.var_ranges || {};
  // Build the formula evaluator once.
  const fnBody = "with({max:Math.max,min:Math.min}){return (" + c.rank_formula + ");}";
  let evalRank;
  try {
    evalRank = new Function(...vars, fnBody);
  } catch (e) {
    console.warn("[catalog] bad family formula:", c.rank_formula, e);
    return [];
  }
  // Resolve each symbolic-format slot to either a literal integer or
  // a variable name; non-numeric strings are treated as variable refs.
  const fmtSpec = c.format_symbolic.map(s => {
    const n = parseInt(s, 10);
    return Number.isFinite(n) && String(n) === s ? { lit: n } : { var: s };
  });
  // Cartesian product over var ranges.
  const varList = vars.slice();
  const out = [];
  const rec = (i, env) => {
    if (i === varList.length) {
      const fmt = fmtSpec.map(slot => slot.lit !== undefined ? slot.lit : env[slot.var]);
      // Canonicalise: sorted ascending (tensor symmetry — rank is permutation-invariant)
      // so duplicates from different orderings collapse.
      const sorted = fmt.slice().sort((a, b) => a - b);
      const args = varList.map(v => env[v]);
      let rk;
      try { rk = evalRank(...args); } catch { return; }
      if (!Number.isFinite(rk) || rk <= 0) return;
      out.push({
        ...c,
        format: sorted,
        rank: Math.round(rk),
        notes: (c.notes || "") + ` [family: ${c.rank_formula}]`,
      });
      return;
    }
    const v = varList[i];
    const r = ranges[v];
    if (!r) return;
    for (let x = r[0]; x <= r[1]; x++) {
      env[v] = x;
      rec(i + 1, env);
    }
  };
  rec(0, {});
  return out;
}

let allSchemes = [];
let lowerBounds = []; // entries from lower-bounds.json
let derivedBounds = []; // entries from derived-from-cited-bounds.json (formula-derived, no scheme file)
let citedBounds = []; // entries from cited-bounds.json (literature claims w/o scheme)
// Globally-stable reference registry: source-key → {index, source, year, url, citedFor: [LB...]}
// Indices are assigned ONCE at load time, ordered by first appearance, so [3] stays [3] across filter changes.
let refRegistry = new Map();

async function loadCatalog() {
  const [catResp, lbResp, dbResp, cbResp] = await Promise.all([
    fetch("catalog.json"),
    fetch("lower-bounds.json").catch(() => null),
    fetch("derived-from-cited-bounds.json").catch(() => null),
    fetch("cited-bounds.json").catch(() => null),
  ]);
  const data = await catResp.json();
  allSchemes = data.schemes;
  let baseSchemes = allSchemes.length;
  // Catalog schemes have no explicit year — derive it from the source string
  // via classifySource so sort-by-year works for ALL rows, not just cited bounds.
  for (const s of allSchemes) {
    if (s.year == null) s.year = classifySource(s.source).year;
  }
  if (lbResp && lbResp.ok) {
    lowerBounds = (await lbResp.json()).entries || [];
  }
  if (dbResp && dbResp.ok) {
    derivedBounds = (await dbResp.json()).entries || [];
    for (const d of derivedBounds) {
      // 'split' can be: undefined (Rosowski formula), [u, v] (cubic split),
      // or [[..], [..], [..]] (non-cubic per-axis allocs). Render compactly.
      let splitLabel = "";
      if (Array.isArray(d.split)) {
        if (Array.isArray(d.split[0])) {
          splitLabel = ` (allocs ${d.split.map(a => "[" + a.join(",") + "]").join("/")})`;
        } else {
          splitLabel = ` (split ${d.split.join("+")})`;
        }
      }
      const derivedSourceLabel = `derived: ${d.source}${splitLabel}`;
      allSchemes.push({
        format: d.format,
        max_dim: Math.max(...d.format),
        field: d.field || "derived",
        commutative: d.commutative === true,
        rank: d.rank,
        additions: null,
        source: derivedSourceLabel,
        file: null,
        derived: true,
        year: classifySource(derivedSourceLabel).year,
        breakdown: d.breakdown,
        construction: d.construction,
        verified: d.verified,
      });
    }
  }
  if (cbResp && cbResp.ok) {
    citedBounds = (await cbResp.json()).entries || [];
    // Expand FAMILY entries (cited bounds with a symbolic format like
    // ⟨N,2,P⟩ + rank_formula) into concrete (n,m,p,rank) rows on the
    // fly. Each var sweeps its declared range.
    citedBounds = citedBounds.flatMap(expandFamilyEntry);
    for (const c of citedBounds) {
      // A "border" cited bound carries border_rank (R̃ ≤ …), NOT an exact rank —
      // a fundamentally different, asymptotic quantity. Use it as the row's value
      // but FLAG it (border:true) so it renders tagged and is isolated from exact
      // rank comparisons (it must not dominate / shave a real constructive scheme).
      const isBorder = c.kind === "border";
      const ncRank = c.rank != null ? c.rank : (isBorder ? c.border_rank : null);
      // NC row (always emitted when c.commutative is not explicitly true).
      // Skip when there's no usable rank value (else the table shows ×undefined / NaN ω).
      if (c.commutative !== true && ncRank != null) {
        allSchemes.push({
          format: c.format,
          max_dim: Math.max(...c.format),
          field: c.field || "R",
          commutative: false,
          rank: ncRank,
          border: isBorder,
          additions: null,
          source: c.source,
          file: null,
          cited: true,
          year: c.year,
          notes: c.notes || "",
        });
      }
      // Commutative row (separate). Emitted when either c.commutative=true OR
      // the cited bound carries a distinct commutative_rank alongside NC.
      if (c.commutative === true || c.commutative_rank != null) {
        const cmtRank = c.commutative === true ? c.rank : c.commutative_rank;
        allSchemes.push({
          format: c.format,
          max_dim: Math.max(...c.format),
          field: c.field || "R",
          commutative: true,
          rank: cmtRank,
          additions: null,
          source: c.source,
          file: null,
          cited: true,
          year: c.year,
          notes: c.notes || "",
        });
      }
    }
  }
  // Drop cited rank-claims fully superseded by an on-disk constructive scheme.
  // A cited bound carries NO factor matrices — it is only informative when we
  // LACK the scheme (CLAUDE.md: cited bounds are for "only the rank claim is
  // published, no explicit matrices"). If, for EVERY field it covers, some
  // on-disk scheme of the same format + commutativity already achieves rank ≤
  // the cited rank, the claim is pure duplication (e.g. the ⟨2,2,2⟩
  // Strassen-1969 cited row vs the Strassen atom ×7) — remove it entirely so it
  // never shows in either the default or "show history" view.
  {
    const odKey = (fmt, comm, f) => `${fmt.join(",")}|${comm === true}|${f}`;
    const bestOnDisk = new Map(); // cell -> min rank among on-disk constructive schemes
    for (const s of allSchemes) {
      const onDisk = !s.cited && !s.derived && s.file != null && s.scheme_provided !== false;
      if (!onDisk) continue;
      for (const f of schemeFieldList(s)) {
        const k = odKey(s.format, s.commutative, f);
        const cur = bestOnDisk.get(k);
        if (cur == null || s.rank < cur) bestOnDisk.set(k, s.rank);
      }
    }
    const before = allSchemes.length;
    allSchemes = allSchemes.filter(s => {
      if (!s.cited) return true;
      if (s.border) return true; // border rank is a different quantity — never shave it against an exact-rank scheme
      const fields = schemeFieldList(s);
      if (!fields.length) return true; // unclassifiable → never silently drop
      const dominated = fields.every(f => {
        const best = bestOnDisk.get(odKey(s.format, s.commutative, f));
        return best != null && best <= s.rank;
      });
      return !dominated;
    });
    const dropped = before - allSchemes.length;
    if (dropped) console.info(`[catalog] dropped ${dropped} cited bound(s) superseded by an on-disk scheme`);
  }

  // Lower-bound rows (rank ≥ …): published theoretical floors with no factor
  // matrices, surfaced as their own row class behind the "Lower bound" selector
  // (default-excluded, like border rank). A LOWER bound is a fundamentally
  // different quantity from an achievable rank — it must never dedupe against,
  // shave, or be shaved by an exact-rank scheme, so it carries `lower: true` and
  // is isolated in every dedup/dominance key. Emitted AFTER the cited-supersession
  // drop above so an on-disk scheme can never delete a lower bound.
  // (Skip `kind: "border"` LB entries — those bound BORDER rank R̃, a different
  // quantity again; they stay in lower-bounds.json for the reference registry /
  // lookupLowerBound only.)
  for (const lb of lowerBounds) {
    if (lb.kind === "border") continue;
    if (!Array.isArray(lb.format) || lb.format.length !== 3 || lb.lb == null) continue;
    // "all" = field-agnostic floor → valid under every field selector; otherwise
    // keep the single declared field WITHOUT cross-field lifting (a lower bound
    // proven over R does NOT carry to C — extending the field can only lower rank).
    const fields = lbFieldList(lb.field);
    allSchemes.push({
      format: lb.format,
      max_dim: Math.max(...lb.format),
      field: lb.field === "all" ? "R" : lb.field,
      fields,
      commutative: false,
      rank: lb.lb,
      lower: true,
      lb_field: lb.field, // original published field token ("all" | "char0" | "F2" | "R" | "F2/Z/Q/R/C" | …)
      border: false,
      tight: lb.tight === true,
      additions: null,
      source: lb.source,
      file: null,
      cited: true,
      year: lb.year,
      url: lb.url,
      notes: lb.notes || "",
    });
  }

  buildRefRegistry();
  // First render of references uses the full registry; render() will re-filter on first call.
  renderReferences();
  document.getElementById("result-count").textContent =
    `${allSchemes.length} entries loaded (${baseSchemes} verified schemes + ${derivedBounds.length} derived + ${citedBounds.length} cited bounds; manifest ${data.generated}; ${lowerBounds.length} lower bounds; ${refRegistry.size} references)`;
  bindControls();
  // Restore filter state from URL hash (if any) BEFORE first render so a
  // shared link lands on the intended view, not the default.
  readHashState();
  render();
  renderOmegaTriangle();  // static view (filter-independent) — build after load
  // The heatmap has its own controls (colour mode + column order) AND now resyncs to
  // the Field filter (so the per-shape best ω reflects the selected algebra).
  for (const id of ["omega-heat-mode", "omega-heat-order", "f-field"]) {
    document.getElementById(id)?.addEventListener("change", renderOmegaTriangle);
  }
}

// Categorise the scheme that wins a heatmap cell, for the "source / composition
// type" colour mode. Composed/derived schemes bucket by their composition op
// (recursive / concat / kron / …); imported atoms bucket by author family.
// True iff the scheme is a derived/composed product (not a primary atom/import).
// Mirrors the derived branch of heatSourceCategory; used for heatmap tie-breaks.
function isHeatDerived(s) {
  const src = (s.source || "").toLowerCase();
  return s.atom === false || s.derived === true
      || src.startsWith("derived") || src.startsWith("composed") || src === "unknown";
}

function heatSourceCategory(s) {
  const src = (s.source || "").toLowerCase();
  // Composed/derived: detected by the atom/derived flags OR a derived_/composed_
  // source prefix (some materialiser outputs ship with a stale atom:true).
  if (s.atom === false || s.derived || src.startsWith("derived") || src.startsWith("composed")) {
    if (src.startsWith("derived:")) return "derived: formula";
    if (src.includes("recursive")) return "derived: recursive";
    if (src.includes("concat")) return "derived: concat";
    if (src.includes("kron")) return "derived: kron";
    if (src.includes("recombin")) return "derived: recombine";
    if (src.includes("strassen")) return "derived: strassen-rec";
    return "derived: other";
  }
  if (src.startsWith("perminov")) return "Perminov";
  if (src.startsWith("hopcroft")) return "Hopcroft";
  if (src.startsWith("alphatensor")) return "AlphaTensor";
  if (src.startsWith("alphaevolve")) return "AlphaEvolve";
  if (src.startsWith("smirnov")) return "Smirnov";
  if (src.startsWith("laderman")) return "Laderman";
  if (src.startsWith("strassen")) return "Strassen";
  if (src.startsWith("winograd")) return "Winograd";
  if (src.startsWith("waksman")) return "Waksman";
  if (src.startsWith("rosowski")) return "Rosowski";
  if (src.startsWith("dis09")) return "DIS09";
  if (src.startsWith("fmm")) return "FMM-Lille";
  if (src.includes("dumas") || src.includes("pernet") || src.includes("sedoglavic")) return "DPS";
  if (src.startsWith("moosbauer")) return "Moosbauer-Poole";
  return "other";
}

// Stable colours for the common categories; anything unmapped cycles HEAT_FALLBACK.
const HEAT_SRC_COLORS = {
  "Perminov": "#1f77b4", "AlphaTensor": "#ff7f0e", "AlphaEvolve": "#d62728",
  "Smirnov": "#2ca02c", "Laderman": "#9467bd", "Strassen": "#8c564b",
  "Winograd": "#e377c2", "Waksman": "#7f7f7f", "Rosowski": "#bcbd22",
  "DIS09": "#17becf", "FMM-Lille": "#aec7e8", "DPS": "#98df8a",
  "Moosbauer-Poole": "#c5b0d5", "Hopcroft": "#ffbb78",
  "derived: recursive": "#393b79", "derived: concat": "#637939",
  "derived: kron": "#8c6d31", "derived: recombine": "#843c39",
  "derived: strassen-rec": "#7b4173", "derived: formula": "#a55194",
  "derived: other": "#bbbbbb", "other": "#dddddd",
};
const HEAT_FALLBACK = ["#ff9896", "#c49c94", "#f7b6d2", "#dbdb8d", "#9edae5", "#c7c7c7"];

/**
 * ω heatmap (triangle): rows = maxDim p (2..D), right-anchored on the cube axis
 * ⟨p,p,p⟩; moving left visits the smaller ⟨n,m,p⟩ shapes (n≤m≤p). Two colour
 * modes: ω (self-calibrated gradient over the best non-commutative rank) or
 * source/composition type (categorical). Pure derived view over allSchemes —
 * no filter dependence.
 */
function renderOmegaTriangle() {
  const host = document.getElementById("omega-triangle");
  if (!host) return;
  const mode = document.getElementById("omega-heat-mode")?.value || "omega";
  // Resync the heatmap to the live Field filter: only schemes valid over the
  // selected field win their shape's cell (so e.g. F2 shows ⟨4,4,4⟩=47, not 48).
  const field = document.getElementById("f-field")?.value || "";
  const STRASSEN_OMEGA = Math.log(7) / Math.log(2);   // ω < this ⇒ better-than-Strassen
  const D = 32;
  // Best (lowest) non-commutative, non-border SCHEME per canonical shape "n,m,p"
  // (the whole row, so source mode can read its provenance — not just the rank).
  const best = new Map();
  for (const s of allSchemes) {
    if (!s.format || s.format.length !== 3 || s.rank == null) continue;
    if (s.commutative === true || s.border === true || s.lower === true) continue;
    if (field && !schemeValidForRequestedField(s.field, field, s)) continue;
    const f = [s.format[0], s.format[1], s.format[2]].sort((a, b) => a - b);
    if (f[0] < 2 || f[2] > D) continue;
    const key = f.join(",");
    const prev = best.get(key);
    // Lower rank always wins. On a rank TIE prefer a non-derived (atom/known)
    // scheme so the cell attributes to the primary source — otherwise an
    // arbitrary derived tie-winner (e.g. a trivial [1,1]³ self-recombination of
    // Strassen) would mislabel ⟨2,2,2⟩ as "derived".
    if (prev == null || s.rank < prev.rank
        || (s.rank === prev.rank && !isHeatDerived(s) && isHeatDerived(prev))) {
      best.set(key, s);
    }
  }
  const omega = (n, m, p, r) => 3 * Math.log(r) / Math.log(n * m * p);
  const isPrime = (x) => { if (x < 2) return false; for (let i = 2; i * i <= x; i++) if (x % i === 0) return false; return true; };

  // Column ordering of the (n,m) pairs within a row (toggle). "max" groups by the
  // larger small-axis m (2,2; 2,3; 3,3; 2,4; 3,4; 4,4; …); "min" groups by the
  // smaller axis n (2,2; 2,3; 2,4; …; 3,3; 3,4; …) so all small-n formats sit
  // together — handy for seeing the shared ω skew of thin (small-n) shapes.
  // Either way ⟨d,d⟩ is the max pair, so the cube stays each row's rightmost cell.
  const order = document.getElementById("omega-heat-order")?.value || "min";
  const cols = [];
  for (let m = 2; m <= D; m++) for (let n = 2; n <= m; n++) cols.push([n, m]);
  if (order === "min") cols.sort((a, b) => a[0] - b[0] || a[1] - b[1]);

  // First pass: compute every cell's ω and calibrate the colour scale to the
  // actual best (lowest) / worst (highest) ω in *this* dataset, not a fixed
  // [2,3] window (which would waste most of the gamut — real ω here clusters
  // well above 2). green = best held, red = worst held.
  const cells = [];   // {n, m, d, r, w, scheme}
  let wMin = Infinity, wMax = -Infinity;
  for (const [n, m] of cols) {
    for (let d = m; d <= D; d++) {
      const s = best.get([n, m, d].join(","));
      if (s == null) continue;
      const r = s.rank;
      const w = omega(n, m, d, r);
      cells.push({ n, m, d, r, w, scheme: s });
      if (w < wMin) wMin = w;
      if (w > wMax) wMax = w;
    }
  }
  const span = wMax - wMin;
  const omegaColor = (w) => {
    const t = span > 1e-9 ? (w - wMin) / span : 0;   // 0 = best → green, 1 = worst → red
    return `hsl(${120 * (1 - t)},70%,50%)`;
  };
  // Source/composition mode: assign a stable colour per category present, in
  // descending frequency so the dominant providers get the named palette and the
  // rare ones cycle the fallback. catColor(scheme) → fill; catList → legend rows.
  const catCount = new Map();
  for (const cell of cells) {
    const cat = heatSourceCategory(cell.scheme);
    catCount.set(cat, (catCount.get(cat) || 0) + 1);
  }
  const catList = [...catCount.entries()].sort((a, b) => b[1] - a[1]);
  const catColorMap = new Map();
  let fb = 0;
  for (const [cat] of catList) {
    catColorMap.set(cat, HEAT_SRC_COLORS[cat] || HEAT_FALLBACK[fb++ % HEAT_FALLBACK.length]);
  }
  // Layout: rows are the largest axis p = maxDim (2…D), one per row. Each row is
  // RIGHT-ANCHORED so the cube ⟨p,p,p⟩ — the rightmost shape of every row — lands
  // on a single vertical axis on the right (user 2026-06-07). Moving left within a
  // row visits the smaller shapes ⟨n,m,p⟩ (n≤m≤p), in the chosen column order; the
  // left edge is ragged (the short ⟨2,2,2⟩ row up top, the full ⟨…,…,32⟩ row at
  // the bottom). Cell WIDTH tapers with the offset k from the cube (k=0 = widest),
  // so the cube and near-cube columns read clearly while the long tail stays compact.
  const chh = 13, top = 16, leftPad = 4, labelPad = 30;
  const N = cols.length;
  const wOff = (k) => Math.max(3, 11 - 0.28 * k);   // block width at offset k from the cube
  const soff = [0];                                  // cumulative offset widths
  for (let k = 0; k < N; k++) soff.push(soff[k] + wOff(k));
  const xR = leftPad + soff[N];                       // shared right edge = the cube axis

  const W = xR + labelPad;
  const H = top + (D - 1) * chh + 8;
  // Scale to fill the container width: a viewBox + width:100%/height:auto lets the
  // browser stretch the intrinsic W×H layout to the full available width (height
  // follows proportionally), instead of rendering at a fixed pixel width and
  // leaving slack / a horizontal scrollbar.
  let svg = `<svg viewBox="0 0 ${W} ${H}" width="100%" height="auto" preserveAspectRatio="xMinYMin meet" style="width:100%;height:auto;display:block" font-family="monospace" font-size="9">`;
  // Place cells row by row. Row d holds every (n,m) with m ≤ d, in the chosen
  // column order; the cube ⟨d,d,d⟩ is its last (max) pair → offset 0, flush right.
  // Coloured cells carry data-* for the tooltip; no-scheme cells are grey; cube
  // cells get a thin outline marking the right-hand ⟨p,p,p⟩ axis.
  for (let d = 2; d <= D; d++) {
    const row = cols.filter(([, m]) => m <= d);
    const L = row.length - 1;
    const y = top + (d - 2) * chh;
    for (let i = 0; i <= L; i++) {
      const n = row[i][0], m = row[i][1];
      const k = L - i;                                // offset from the cube
      const x = xR - soff[k + 1];
      const bw = Math.max(1, wOff(k) - 1);
      const cube = n === m && m === d;
      const s = best.get([n, m, d].join(","));
      if (s == null) {
        svg += `<rect x="${x}" y="${y}" width="${bw}" height="${chh - 1}" fill="#eee"${cube ? ' stroke="#333" stroke-width="0.7"' : ''}/>`;
        continue;
      }
      const w = omega(n, m, d, s.rank);
      const fill = mode === "source"
        ? (catColorMap.get(heatSourceCategory(s)) || "#ccc")
        : omegaColor(w);
      const src = s.source || "?";
      svg += `<rect x="${x}" y="${y}" width="${bw}" height="${chh - 1}" fill="${fill}"`
        + (cube ? ' stroke="#333" stroke-width="0.7"' : '')
        + ` data-shape="${n},${m},${d}" data-rank="${s.rank}" data-omega="${w.toFixed(3)}"`
        + ` data-source="${escapeHtml(src)}"`
        + ` style="cursor:crosshair">`
        + `<title>⟨${n},${m},${d}⟩ = ${s.rank}   ω=${w.toFixed(3)}   ${src}</title></rect>`;
      // Special border: this shape's best ω is below Strassen's log₂7 ≈ 2.807 —
      // a "better-than-Strassen" shape. Overlaid so it coexists with the cube outline.
      if (w < STRASSEN_OMEGA) {
        svg += `<rect x="${x}" y="${y}" width="${bw}" height="${chh - 1}" fill="none"`
          + ` stroke="#0033cc" stroke-width="1.3" pointer-events="none"/>`;
      }
    }
  }
  // Right-side axis labels: the max dim p next to its cube, prime-flagged ◆.
  for (let d = 2; d <= D; d++) {
    const y = top + (d - 2) * chh;
    const pr = isPrime(d);
    svg += `<text x="${xR + 3}" y="${y + chh - 3}" fill="${pr ? '#b00' : '#333'}" font-weight="${pr ? 'bold' : 'normal'}">${d}${pr ? '◆' : ''}</text>`;
  }
  svg += `</svg>`;
  host.innerHTML = svg;

  // Reflect the live calibration range in the ω legend.
  const legBest = document.getElementById("omega-legend-best");
  const legWorst = document.getElementById("omega-legend-worst");
  if (legBest && legWorst && cells.length) {
    legBest.textContent = "ω=" + wMin.toFixed(3) + " (best held)";
    legWorst.textContent = "ω=" + wMax.toFixed(3) + " (worst held)";
  }

  // Swap legends to match the colour mode.
  const legOmega = document.getElementById("omega-legend-omega");
  const legSource = document.getElementById("omega-legend-source");
  if (legOmega) legOmega.style.display = mode === "source" ? "none" : "";
  if (legSource) {
    legSource.style.display = mode === "source" ? "" : "none";
    if (mode === "source") {
      const swatches = catList.map(([cat, count]) =>
        `<span style="display:inline-block;margin:0 8px 4px 0;white-space:nowrap">`
        + `<span style="background:${catColorMap.get(cat)};padding:0 7px;color:#000">&nbsp;</span> `
        + `${escapeHtml(cat)} <span style="color:#888">(${count})</span></span>`).join("");
      legSource.innerHTML = `Cell = the source of the best-rank scheme held for that shape `
        + `(imported atoms by author, composed schemes by composition type). `
        + `<span style="background:#eee;padding:0 6px">no scheme</span><br>${swatches}`;
    }
  }

  // A single floating tooltip, shared by all cells, following the cursor.
  let tip = document.getElementById("omega-tooltip");
  if (!tip) {
    tip = document.createElement("div");
    tip.id = "omega-tooltip";
    tip.style.cssText = "position:fixed;z-index:1000;pointer-events:none;display:none;"
      + "background:#222;color:#fff;font:11px/1.4 monospace;padding:4px 7px;"
      + "border-radius:4px;box-shadow:0 2px 6px rgba(0,0,0,.3);white-space:nowrap";
    document.body.appendChild(tip);
  }
  const place = (e) => {
    // Keep the tip on-screen (flip left of the cursor near the right edge).
    const pad = 14;
    let x = e.clientX + pad, y = e.clientY + pad;
    if (x + tip.offsetWidth + 8 > window.innerWidth) x = e.clientX - tip.offsetWidth - pad;
    if (y + tip.offsetHeight + 8 > window.innerHeight) y = e.clientY - tip.offsetHeight - pad;
    tip.style.left = x + "px";
    tip.style.top = y + "px";
  };
  host.onmousemove = (e) => {
    const rect = e.target.closest("rect[data-shape]");
    if (!rect) { tip.style.display = "none"; return; }
    const [n, m, d] = rect.getAttribute("data-shape").split(",");
    tip.innerHTML = `<strong>⟨${n},${m},${d}⟩</strong> = ${rect.getAttribute("data-rank")}`
      + `<br>ω = ${rect.getAttribute("data-omega")}`
      + `<br>${escapeHtml(rect.getAttribute("data-source") || "?")}`;
    tip.style.display = "block";
    place(e);
  };
  host.onmouseleave = () => { tip.style.display = "none"; };
}

/**
 * Classify a raw source string (from catalog.json / cited-bounds.json / derived-from-cited-bounds.json)
 * into a compact label + canonical registry key. Different raw strings that refer to the
 * same underlying work (e.g. "Dronperminov-cr60_cn97_ZT_reduced" and "Dronperminov-Z")
 * collapse to a single registry entry with one [N] number.
 */
/**
 * Source priority — lower = primary discovery, higher = aggregator/rehost.
 * Used to dedupe rows that share the same (format, field, rank, additions):
 * the lowest-priority row wins, the rest are hidden so we don't list a
 * rediscovery 50 years later (e.g. fmm-lille hosting Strassen 1969).
 */
function sourcePriority(rawSource) {
  if (!rawSource) return 99;
  // Primary 1969-era algorithms.
  if (/^Strassen$/i.test(rawSource)) return 1;
  if (/^Laderman$/i.test(rawSource)) return 1;
  // Primary recent algorithms.
  if (/^Alphatensor/i.test(rawSource)) return 2;
  if (/^Alphaevolve$/i.test(rawSource)) return 2;
  // Compositions we built ourselves.
  if (/^Composed-/i.test(rawSource)) return 4;
  // Cited bounds (literature claims with attributed author).
  if (rawSource.match(/\(via /)) return 3;  // "Strassen 1969 (via DIS09)" → primary
  // Derived bounds (formula).
  if (/^derived:/i.test(rawSource)) return 5;
  // Aggregator catalogs (republishing third-party schemes).
  if (/^(Perminov|Dronperminov)/i.test(rawSource)) return 8;
  if (/^Fmm-lille$/i.test(rawSource)) return 9;
  return 6;
}

function classifySource(rawSource) {
  if (!rawSource) return { label: "(unknown)", key: null, year: null, url: null };
  // The manifest's explicit placeholder for a scheme with no source on disk —
  // same rendering as a missing source, not a literal "unknown" author.
  if (/^unknown$/i.test(rawSource)) return { label: "(unknown)", key: null, year: null, url: null };
  // Catalog scheme sources (filename-derived).
  if (/^Strassen$/i.test(rawSource))       return { label: "Strassen 1969", key: "Strassen 1969", year: 1969 };
  if (/^Laderman$/i.test(rawSource))       return { label: "Laderman 1976", key: "Laderman 1976", year: 1976 };
  if (/^Alphatensor-F2$/i.test(rawSource)) return { label: "AlphaTensor 2022 (F2)", key: "AlphaTensor 2022", year: 2022 };
  if (/^Alphatensor-Z$/i.test(rawSource))  return { label: "AlphaTensor 2022 (Z)", key: "AlphaTensor 2022", year: 2022 };
  if (/^Alphaevolve$/i.test(rawSource))    return { label: "AlphaEvolve 2025", key: "AlphaEvolve 2025", year: 2025 };
  if (/^Kaporin 2024/i.test(rawSource))    return { label: "Kaporin 2024", key: "Kaporin 2024", year: 2024 };
  // Bare "AlphaTensor" and its published factorisation archive URL
  // (github.com/google-deepmind/alphatensor — the alphatensor_Q imports).
  if (/alphatensor/i.test(rawSource))      return { label: "AlphaTensor 2022", key: "AlphaTensor 2022", year: 2022 };
  // The serendipitous 17–32 band is a DISTINCT paper (arXiv:2606.02480) — keep it
  // off the generic repo catch so it links to its own REFERENCES.md entry [86].
  if (/^Perminov 2026 \(serendipitous\)/i.test(rawSource)) return { label: "Perminov 2026 (serendip.)", key: "Perminov 2026 (serendipitous)", year: 2026 };
  if (/^(Perminov|Dronperminov)/i.test(rawSource)) return { label: "Perminov (catalog)", key: "Perminov (FastMatrixMultiplication)", year: null };
  if (/^Fmm-lille$/i.test(rawSource))      return { label: "fmm-lille", key: "fmm-lille catalog", year: null };
  if (/^(Derived|Composed)[_-]/i.test(rawSource)) return { label: rawSource.replace(/^(Derived|Composed)[_-]/i, "derived "), key: null }; // no ref
  // Cited/derived sources like "Strassen 1969 (via DIS09)" or "derived: Sedoglavic 2017 (closed-form identity)".
  let stripped = rawSource.replace(/^derived:\s*/i, "").replace(/\s*\(via [^)]+\)\s*$/, "").trim();
  // Composed schemes (recursive split / Kronecker / concat) are derived, not a
  // dated upstream result — they carry NO citation (user 2026-06-03). Catch the
  // "derived: Composed (…)" form here (the "^Composed-" guard above only catches
  // the un-prefixed raw form).
  if (/^Composed\b/i.test(stripped)) return { label: stripped, key: null };
  const m = stripped.match(/^([A-Za-z][A-Za-zÀ-ſ\s\-–.,]*?)\s+(\d{4})\b/);
  if (m) {
    const name = m[1].replace(/\s+/g, " ").trim();
    const year = m[2];
    return { label: `${name} ${year}`, key: `${name} ${year}`, year: parseInt(year, 10) };
  }
  return { label: stripped, key: stripped, year: null };
}

function buildRefRegistry() {
  refRegistry = new Map();
  // LB sources first (highest-quality metadata: year + url + tightness).
  for (const lb of lowerBounds) {
    const k = lb.source;
    if (!refRegistry.has(k)) {
      refRegistry.set(k, { index: refRegistry.size + 1, source: k, year: lb.year, url: lb.url, citedFor: [] });
    }
    refRegistry.get(k).citedFor.push({ format: lb.format, kind: "LB", value: lb.lb, field: lb.field });
  }
  // Cited-bound and derived sources via classifySource.
  function ingest(rawSource, format, kind, value, field) {
    const c = classifySource(rawSource);
    if (!c.key) return;
    if (!refRegistry.has(c.key)) {
      refRegistry.set(c.key, { index: refRegistry.size + 1, source: c.key, year: c.year, url: null, citedFor: [] });
    }
    refRegistry.get(c.key).citedFor.push({ format, kind, value, field });
  }
  for (const c of citedBounds) ingest(c.source, c.format, "UB", c.rank, c.field || "R");
  // Catalog scheme sources — register so e.g. "Strassen 69 [N]" link works on verified rows too.
  for (const s of allSchemes) {
    if (s.derived || s.cited) continue;
    ingest(s.source, s.format, "UB", s.rank, s.field);
  }
  // Derived bounds — register too (otherwise their compact labels can't link).
  for (const d of derivedBounds) ingest(d.source, d.format, "UB", d.rank, d.field || "R");
}

/**
 * Renders the References section. If `usedKeys` is provided, only refs whose
 * registry key is in that set are shown (so the footer shrinks to match the
 * currently filtered table).
 */
function renderReferences(usedKeys) {
  const ol = document.getElementById("references-list");
  if (!ol) return;
  ol.innerHTML = "";
  const sorted = [...refRegistry.values()].sort((a, b) => a.index - b.index);
  let shown = 0;
  for (const r of sorted) {
    if (usedKeys && !usedKeys.has(r.source)) continue;
    const li = document.createElement("li");
    li.id = "ref-" + r.index;
    // Primary link: REFERENCES.md entry (full BibTeX + PDF link + annotation).
    // Secondary link: external URL (when known, e.g. arXiv direct).
    const anchor = referencesMdAnchor(r.source);
    let sourceText;
    if (anchor) {
      sourceText = `<a href="${REFERENCES_MD_URL}#${anchor}" target="_blank" rel="noopener" title="Full BibTeX + local PDF + annotation in REFERENCES.md">${escapeHtml(r.source)}</a>`;
      if (r.url) {
        sourceText += ` <a href="${escapeHtml(r.url)}" target="_blank" rel="noopener" title="External URL">↗</a>`;
      }
    } else if (r.url) {
      sourceText = `<a href="${escapeHtml(r.url)}" target="_blank" rel="noopener">${escapeHtml(r.source)}</a>`;
    } else {
      sourceText = escapeHtml(r.source);
    }
    // Compact summary instead of dumping every cited format: heavy citers
    // (Perminov, FMM-Lille, AlphaTensor) generate hundreds of ⟨n,m,p⟩ entries
    // which makes the footer an unreadable block. Show counts + a small
    // sample, and keep the full list reachable via a tooltip on the count.
    const seen = new Set();
    const dedup = r.citedFor.filter(c => {
      const k = c.format.join(",") + "|" + (c.field || "") + "|" + c.kind;
      return !seen.has(k) && seen.add(k);
    });
    const ubCount = dedup.filter(c => c.kind !== "LB").length;
    const lbCount = dedup.length - ubCount;
    const parts = [];
    if (ubCount) parts.push(`${ubCount} upper-bound${ubCount === 1 ? "" : "s"}`);
    if (lbCount) parts.push(`${lbCount} lower-bound${lbCount === 1 ? "" : "s"}`);
    const summary = parts.join(" + ");
    // Tooltip carries the full per-shape list for the curious; the visible
    // text stays one line.
    const fullList = dedup
      .map(c => `⟨${c.format.join(",")}⟩ ${c.kind === "LB" ? "LB" : "UB"}`
          + (c.field && c.field !== "all" ? ` (${c.field})` : ""))
      .join(", ");
    // Prefix with the stable [N] so filtering hidden references doesn't
    // renumber the visible ones — match the source-column "[N]" exactly.
    li.innerHTML = `<span class="refnum-anchor">[${r.index}]</span> ${sourceText} — <span title="${escapeHtml(fullList)}">${summary}</span>.`;
    ol.appendChild(li);
    shown++;
  }
  // If filter wiped everything, show a placeholder.
  if (shown === 0) {
    const li = document.createElement("li");
    li.style.color = "var(--muted)";
    li.style.listStyleType = "none";
    li.textContent = "(no references for current filter)";
    ol.appendChild(li);
  }
}

// Expand a lower-bound entry's `field` token into the explicit list of fields it
// is valid over. Tokens: "all" = field-agnostic (every field); "char0" = the
// characteristic-0 fields only (Z/Q/R/C, NOT F2/F3 — a char-0 floor does not
// transfer to positive characteristic); otherwise a "/"-separated explicit set
// (e.g. "F2/Z/Q/R/C" = a floor known over F2 and char 0 but weaker over F3).
function lbFieldList(field) {
  if (field === "all") return ["Z", "Q", "R", "C", "F2", "F3"];
  if (field === "char0") return ["Z", "Q", "R", "C"];
  return field.split("/").map(x => x.trim()).filter(Boolean);
}

/**
 * Lookup the best (largest) known LB matching this scheme's (format, field).
 * Matches when the scheme is valid over any field the LB covers.
 */
function lookupLowerBound(scheme) {
  const fmt = scheme.format;
  // Match any LB whose field the scheme is valid over. The singular `field` is
  // gone — `fields` membership is authoritative.
  const fset = new Set(schemeFieldList(scheme));
  let best = null;
  for (const lb of lowerBounds) {
    if (lb.format[0] !== fmt[0] || lb.format[1] !== fmt[1] || lb.format[2] !== fmt[2]) continue;
    if (!lbFieldList(lb.field).some(f => fset.has(f))) continue;
    if (!best || lb.lb > best.lb) best = lb;
  }
  return best;
}

// Keys the user has explicitly interacted with this session. A key only
// reaches the URL hash once the user has touched the corresponding
// control — so a fresh page has a clean URL, but checking-and-unchecking
// cubicOnly leaves the value (even when equal to the default) in the
// hash, signalling intent.
const CTRL_TO_HASH_KEY = {
  "f-field": "field",
  "f-cubic": "cubic",
  "f-commutative": "commutative",
  "f-composed": "composed",
  "f-border": "border",
  "f-lower": "lower",
  "f-max-dim": "max",
  "f-min-dim": "min",
  "f-source": "source",
  "f-shape": "shape",
  "f-sort": "sort",
  "f-secondary": "secondary",
  "f-omega-beats-strassen": "betterthan",
  "f-show-history": "history",
};
const _touchedKeys = new Set();

function bindControls() {
  Object.entries(CTRL_TO_HASH_KEY).forEach(([id, key]) => {
    const el = document.getElementById(id);
    const handler = () => { _touchedKeys.add(key); writeHashState(); render(); };
    el.addEventListener("input", handler);
    el.addEventListener("change", handler);
  });
  // External hash changes (Back/Forward, manual edit, paste-in-new-tab) — reflect into controls.
  window.addEventListener("hashchange", () => { readHashState(); render(); });
}

// ── URL ↔ filter-state bijection ────────────────────────────────────────────
// All filter controls round-trip through location.hash so users can share a
// URL that restores the exact view (GitHub Pages = no server-side routing,
// hence the hash fragment rather than query string).
//
// Hash format: #key=value&key=value (URI-encoded). Keys are the short
// readable form of each control:
//   field       → f-field  (e.g. "Q", "R/Q/Z", "F2", "" for "all")
//   cubic       → f-cubic  ("1" / "0")
//   commutative → f-commutative  ("include" / "exclude" / "require")
//   min         → f-min-dim (integer)
//   max         → f-max-dim (integer)
//   source      → f-source (substring search)
//   shape       → f-shape (positional ⟨n,m,p⟩ pattern, `*` = wildcard, e.g. 2.*.16)
//   sort        → f-sort
//
// Missing keys → keep the DOM default (i.e. whatever the HTML <select>/<input>
// declares).  Unrecognised keys → silently ignored. Empty / default values
// are stripped from the written hash to keep URLs compact.
const HASH_KEY_TO_CTRL = {
  field: "f-field",
  cubic: "f-cubic",
  commutative: "f-commutative",
  composed: "f-composed",
  border: "f-border",
  lower: "f-lower",
  min: "f-min-dim",
  max: "f-max-dim",
  source: "f-source",
  shape: "f-shape",
  sort: "f-sort",
  secondary: "f-secondary",
  betterthan: "f-omega-beats-strassen",
  history: "f-show-history",
};

function parseHash() {
  const h = (location.hash || "").replace(/^#/, "");
  const out = {};
  if (!h) return out;
  for (const pair of h.split("&")) {
    if (!pair) continue;
    const eq = pair.indexOf("=");
    const k = eq >= 0 ? pair.slice(0, eq) : pair;
    const v = eq >= 0 ? pair.slice(eq + 1) : "";
    try { out[decodeURIComponent(k)] = decodeURIComponent(v); }
    catch { /* malformed segment — skip */ }
  }
  return out;
}

let _suppressHashWrite = false;

function readHashState() {
  const state = parseHash();
  _suppressHashWrite = true;
  try {
    for (const [k, ctrlId] of Object.entries(HASH_KEY_TO_CTRL)) {
      if (!(k in state)) continue; // missing → keep DOM default
      // Any key present in the URL is considered touched — restoring a
      // shared link preserves its sticky-overlay-the-default semantics.
      _touchedKeys.add(k);
      const el = document.getElementById(ctrlId);
      if (!el) continue;
      let v = state[k];
      // Legacy links used the empty value for the commutative/composed
      // tri-state selectors to mean "include (any)"; those options now carry
      // the explicit value "include". Map old "" → "include" so shared URLs
      // minted before the rename keep their semantics.
      if (v === "" && (k === "commutative" || k === "composed")) v = "include";
      if (el.type === "checkbox") {
        el.checked = (v === "1" || v === "true");
      } else if (el.tagName === "SELECT") {
        // Accept only values present as <option> — silently ignore unknown.
        const opts = [...el.options].map(o => o.value);
        if (opts.includes(v)) el.value = v;
      } else {
        el.value = v;
      }
    }
  } finally {
    _suppressHashWrite = false;
  }
}

function writeHashState() {
  if (_suppressHashWrite) return;
  const parts = [];
  for (const [k, ctrlId] of Object.entries(HASH_KEY_TO_CTRL)) {
    const el = document.getElementById(ctrlId);
    if (!el) continue;
    let v;
    if (el.type === "checkbox") {
      v = el.checked ? "1" : "0";
    } else {
      v = el.value;
    }
    // Strip values that match the HTML default UNLESS the user has
    // explicitly touched this control this session. The "touched" bit lets
    // someone toggle and untoggle cubicOnly (and similar) while preserving
    // the value in the URL — signalling that the default was intentional,
    // not just untouched. A fresh page has no touched keys, so the URL
    // stays clean.
    if (v === HASH_DEFAULTS[k] && !_touchedKeys.has(k)) continue;
    parts.push(encodeURIComponent(k) + "=" + encodeURIComponent(v));
  }
  const newHash = parts.length ? "#" + parts.join("&") : "";
  // Use history.replaceState (rather than location.hash = …) so we don't
  // spam the browser history on every keystroke in the source-search input,
  // and we don't fire our own hashchange listener.
  const url = location.pathname + location.search + newHash;
  history.replaceState(null, "", url);
}

// Defaults captured from index.html — kept in sync there. If you change a
// <select>/<input>'s `selected`/`value` attribute, update the matching entry
// here so the URL strip-default logic stays correct.
const HASH_DEFAULTS = {
  field: "Q",
  cubic: "1",
  commutative: "exclude",
  composed: "include",
  border: "exclude",
  lower: "exclude",
  min: "2",
  max: "32",
  source: "",
  sort: "format",
  secondary: "additions",
  betterthan: "0",
  history: "0",
};

// The "real" addition count for a scheme: the CSE-scheduled / min count when
// known (e.g. Strassen-Winograd's 15), NOT the flat structural additionCount(U,V,W)
// — that flat number (24 for Winograd) over-counts whenever common subexpressions
// are shared, so it must never be the headline figure. Falls back to the flat
// count when no schedule is recorded.
function effectiveAdditions(s) {
  const cands = [s.scheduled_additions, s.min_additions, s.additions].filter(v => v != null);
  return cands.length ? Math.min(...cands) : null;
}

// Additions table cell: headline the effective (scheduled) count, annotate the
// flat structural count when it is strictly larger (so the saving is visible).
function additionsCell(s) {
  const eff = effectiveAdditions(s);
  if (eff == null) return "—";
  const flat = s.additions;
  // When the flat structural count exceeds the CSE-scheduled (effective) count,
  // prefix it as "24 flat +15". The headline "+N" stays LAST so that — the column
  // being right-aligned — the relevant (effective) addition count lines up
  // vertically across every row, whether or not a flat annotation is present.
  const note = (flat != null && flat > eff)
    ? `<small class="field-narrow" title="flat structural additions (additionCount over U/V/W); the headline ${eff} is the CSE-scheduled count after common subexpressions are reused — e.g. Strassen-Winograd is 15 scheduled vs ${flat} flat">${flat} flat</small> `
    : "";
  return `${note}+${eff}`;
}

// Secondary ranking key (the tiebreak once rank is equal). Returns a score
// where LOWER = better, so the same comparator works for all keys:
//   additions          — fewer is better (use the cheaper of flat / scheduled)
//   buds               — more is better (higher bud_score), so negate
//   projection_margin  — more is better (stronger downward-projection parent), negate
// Missing values sort last (worst) in every key.
function secondaryRankScore(s, key) {
  switch (key) {
    case "buds":
      return -(s.bud_score || 0);
    case "projection_margin":
      return -(s.projection_margin ?? -Infinity);
    case "year":
      // Earliest scheme wins a tied cell — credits the rank-ORIGIN over a later
      // re-discovery or rationalisation (e.g. ℂ⟨4,4,4⟩=48: Kaporin 2024 wins over
      // AlphaEvolve 2025 and DPS 2025, which the additions tie-break would hide).
      // Lower year = better; an undated scheme ranks last (Infinity) — an unknown
      // date must never claim "earliest".
      return s.year ?? Infinity;
    case "additions":
    default:
      return effectiveAdditions(s) ?? Infinity;
  }
}

// Positional shape-pattern match, e.g. "2.*.16" against the DISPLAYED ⟨n,m,p⟩.
// Per-axis token grammar:
//   N    exact    (axis === N)
//   N+   at-least (axis >= N)   — e.g. "30+" means "30 or bigger"
//   N-   at-most  (axis <= N)   — the trailing `-` is a BOUND only when it is
//                                 NOT followed by a digit
//   *    wildcard (any value)
// The separator is LAX: tokens are the digit-runs (optionally bounded by `+` or a
// suffix `-`) and `*`s; ANY other character(s) between them are a separator. The
// `-` disambiguation: a `-` BETWEEN digits is a separator, so "2-3-16" is the
// exact ⟨2,3,16⟩, while "16-" (dash before a separator/end) is the ≤16 bound.
// Thus "2.*.16", "2x16", "2 16", "2*16", "2-3-16" parse positionally, and
// "30+.2.16" means ⟨≥30, 2, 16⟩. Sub-matching: fewer tokens than axes → trailing
// axes are wildcards ("3.3" matches every ⟨3,3,·⟩, "2x16" every ⟨2,16,·⟩).
// Requested specifically for Perminov, whose catalog is browsed one format at a
// time — a wildcard (or an `N+`/`N-` band) on an axis pulls a whole family at once.
function matchShapePattern(format, pattern) {
  // `-(?!\d)` = a `-` suffix only when not followed by a digit (else separator).
  const toks = pattern.match(/\d+(?:\+|-(?!\d))?|\*/g);
  if (!toks) return true; // no digit/`*` token → no constraint
  for (let i = 0; i < toks.length && i < format.length; i++) {
    const t = toks[i];
    if (t === "*") continue;
    const n = parseInt(t, 10); // parseInt("30+") === parseInt("30-") === 30
    if (t.endsWith("+")) {
      if (format[i] < n) return false; // "N+" → axis ≥ N
    } else if (t.endsWith("-")) {
      if (format[i] > n) return false; // "N-" → axis ≤ N
    } else if (format[i] !== n) {
      return false;
    }
  }
  return true;
}

function render() {
  const field = document.getElementById("f-field").value;
  const cubic = document.getElementById("f-cubic").checked;
  const commutative = document.getElementById("f-commutative").value;
  const composed = document.getElementById("f-composed").value;
  const border = document.getElementById("f-border")?.value || "exclude";
  const lower = document.getElementById("f-lower")?.value || "exclude";
  const maxDim = +document.getElementById("f-max-dim").value || 32;
  const minDim = +document.getElementById("f-min-dim").value || 2;
  const sourceQ = document.getElementById("f-source").value.toLowerCase();
  const shapeQ = document.getElementById("f-shape").value.trim();
  const sort = document.getElementById("f-sort").value;
  const secondary = document.getElementById("f-secondary").value;
  const omegaBeatsStrassen = document.getElementById("f-omega-beats-strassen").checked;
  const showHistory = document.getElementById("f-show-history").checked;
  // Strassen's exponent — log₂(7) ≈ 2.80735... The catalog-progress metric
  // is "how many schemes achieve ω strictly below Strassen?", per Perminov
  // 2026 arXiv:2603.02398. We use the catalog's impliedOmega for the test.
  const STRASSEN_OMEGA = Math.log(7) / Math.log(2);

  let filtered = allSchemes.filter(s => {
    if (field && !schemeValidForRequestedField(s.field, field, s)) return false;
    if (cubic && !(s.format[0] === s.format[1] && s.format[1] === s.format[2])) return false;
    if (commutative === "exclude" && s.commutative === true) return false;
    if (commutative === "require" && s.commutative !== true) return false;
    // Border-rank rows (R̃ ≤ …) are a different quantity from exact rank.
    // Default-excluded; "only" isolates them; "include" mixes them in.
    if (border === "exclude" && s.border === true) return false;
    if (border === "only" && s.border !== true) return false;
    // Lower-bound rows (rank ≥ …) are a different quantity from achievable rank.
    // Default-excluded; "only" isolates them; "include" mixes them in.
    if (lower === "exclude" && s.lower === true) return false;
    if (lower === "only" && s.lower !== true) return false;
    if (composed === "exclude" && isComposedScheme(s)) return false;
    if (composed === "require" && !isComposedScheme(s)) return false;
    // min-dim filter: hide rows where ANY axis is below the threshold
    // max-dim filter: hide rows where ANY axis exceeds the threshold
    const fmtMin = Math.min(s.format[0], s.format[1], s.format[2]);
    const fmtMax = Math.max(s.format[0], s.format[1], s.format[2]);
    if (fmtMax > maxDim) return false;
    if (fmtMin < minDim) return false;
    if (sourceQ && !s.source.toLowerCase().includes(sourceQ)) return false;
    if (shapeQ && !matchShapePattern(s.format, shapeQ)) return false;
    if (omegaBeatsStrassen) {
      const w = impliedOmega(s.format, s.rank);
      if (w == null || w >= STRASSEN_OMEGA) return false;
    }
    return true;
  });

  // When a single CONCRETE field is selected, collapse every row's cell
  // membership to THAT field, so the dominance/dedup passes below compare only
  // WITHIN the viewed field. Without this, a row dominated in the selected field
  // (e.g. Bläser ⟨3,3,3⟩≥19 vs Wang ≥20 over F₂) still wins an UNVIEWED field's
  // cell (R/Z/…) and so wrongly appears in the filtered view. Set/sub-class
  // selectors ("" = all, "ZT", "R/Q/Z") keep the full membership list.
  const SINGLE_FIELDS = new Set(["F2", "F3", "Z", "Q", "R", "C"]);
  const cellFieldsOf = (s) =>
    field === "ZT" ? ["ZT"]                                  // ZT view → one cell per format
      : (field && SINGLE_FIELDS.has(field)) ? [field]
      : schemeFieldList(s);

  // Dedupe in two passes:
  //
  // Pass 1: collapse rows with the SAME SOURCE + same (format, field, rank,
  // commutative). These represent the same algorithm; differing only in
  // whether 'additions' is populated. Prefer the row with additions != null
  // (more informative); if tie, lower source-priority wins.
  // Example: catalog 'Strassen' row (a=18) + cited 'Strassen 1969 (via DIS09)'
  // row (a=null) → collapse to the catalog one.
  //
  // Pass 2: collapse rows with DIFFERENT sources but same (format, field,
  // rank, additions, commutative). Lower source-priority (= more primary)
  // wins. Example: 'Strassen' (a=18) and 'fmm-lille' (a=18) → Strassen wins.
  {
    function keyFmt(s) {
      // `fields` membership signature replaces the old singular `field` (gone
      // 2026-06-04). Two rows only dedupe when they're valid over the same set.
      return `${s.format.join(",")}|${schemeFieldList(s).join("/")}|${s.rank}|${s.commutative === true}|${s.border === true}|${s.lower === true}`;
    }
    function priority(s) {
      const adds = effectiveAdditions(s) == null ? 1 : 0; // -1 step penalty for null adds
      return sourcePriority(s.source) * 10 + adds;
    }
    // Pass 1: same source-key
    {
      const byKey = new Map();
      for (const s of filtered) {
        const srcKey = classifySource(s.source).key || s.source;
        const key = `${keyFmt(s)}|src=${srcKey}`;
        const prio = priority(s);
        const prev = byKey.get(key);
        if (!prev || prio < prev.prio) byKey.set(key, { row: s, prio });
      }
      filtered = [...byKey.values()].map(v => v.row);
    }
    // Pass 2: cross-source dedupe by (format, field, rank, additions, commutative)
    {
      const byKey = new Map();
      for (const s of filtered) {
        const eff = effectiveAdditions(s);
        const adds = eff == null ? "-" : eff;
        const key = `${keyFmt(s)}|adds=${adds}`;
        // PREFER DERIVED: when collapsing a cross-source duplicate (same rank+additions+field),
        // keep our OWN derivation (composed) over an import — otherwise it'd be dropped here
        // before the dominance tie-break below ever sees it. [prefer-derived]
        const prio = sourcePriority(s.source) - (isComposedScheme(s) ? 10000 : 0);
        const prev = byKey.get(key);
        if (!prev || prio < prev.prio) byKey.set(key, { row: s, prio });
      }
      filtered = [...byKey.values()].map(v => v.row);
    }

    // Pass 3: drop derived-bounds rows that are STRICTLY WORSE than any
    // catalog/cited row in EVERY field they share. Per user feedback 2026-05-28:
    // it's noise to show "solven-strassen 2026: ⟨5,5,5⟩ ×98" when the catalog
    // already has AlphaEvolve at 93. Keep derived rows when they're a strict
    // improvement (lower rank) in at least one of their fields.
    //
    // Cells are per (format, commutative, FIELD) now that the singular `field`
    // cluster is gone (2026-06-04): a scheme participates in one cell per field
    // it is valid over (its `fields` membership).
    {
      const cellKey = (s, f) => `${s.format.join(",")}|${s.commutative === true}|${f}`;
      // First pass: per cell, find the min non-derived rank.
      const bestNonDerived = new Map();
      for (const s of filtered) {
        if (s.derived || s.border || s.lower) continue; // border / lower bound ≠ achievable rank — don't let them shave real derived rows
        for (const f of cellFieldsOf(s)) {
          const k = cellKey(s, f);
          const cur = bestNonDerived.get(k);
          if (cur == null || s.rank < cur) bestNonDerived.set(k, s.rank);
        }
      }
      // Second pass: keep a derived row if it strictly improves the best
      // non-derived rank in at least one of its fields (or no non-derived
      // competitor exists there).
      filtered = filtered.filter(s => {
        if (!s.derived) return true;
        const fields = cellFieldsOf(s);
        if (!fields.length) return true; // unclassifiable → never silently drop
        return fields.some(f => {
          const best = bestNonDerived.get(cellKey(s, f));
          return best == null || s.rank < best;
        });
      });
    }

    // Pass 4: dominance collapse (default ON; disabled by "show history").
    // Per user 2026-06-03/-06-04: by default show only the best representative
    // per (format, commutative, FIELD) cell, where a scheme participates in one
    // cell per field in its `fields` membership. A row survives iff it is the
    // best representative for at least one of its fields; a row beaten in every
    // field it covers is "dominated" and hidden until "Show history".
    //
    // This is what hides the AlphaTensor ⟨2,2,2⟩=7,+24 F₂ row: Strassen's entry
    // lists F₂ in its `fields`, so it competes in — and wins — the F₂ cell
    // (7,+18), leaving the AlphaTensor row the representative of no field.
    // Tie-break on equal rank: prefer scheme-on-disk (not derived/cited), then
    // lower additions, then lower source priority.
    if (!showHistory) {
      // Among equal-rank schemes the representative is chosen by the
      // user-selected secondary key (default "additions"), then a preference
      // for an on-disk scheme over a derived/cited one, then source priority.
      const better = (a, b) => {
        // Lower-bound rows live in their own cells (cellKey carries the `lower`
        // flag), so a/b here are both lower bounds or both achievable. For two
        // lower bounds the STRONGER (higher) floor wins — ⟨5,5,5⟩ ≥48 beats ≥47;
        // for achievable ranks the lower rank wins as usual.
        if (a.rank !== b.rank) return (a.lower && b.lower) ? a.rank > b.rank : a.rank < b.rank;
        const as = secondaryRankScore(a, secondary);
        const bs = secondaryRankScore(b, secondary);            // secondary: user choice
        if (as !== bs) return as < bs;
        // A real scheme (factor matrices on disk) still beats a formula-only derived BOUND.
        const aDisk = (a.derived || a.scheme_provided === false) ? 1 : 0;
        const bDisk = (b.derived || b.scheme_provided === false) ? 1 : 0;
        if (aDisk !== bDisk) return aDisk < bDisk;
        // PREFER DERIVED (user 2026-06-28): among real schemes at equal rank+additions, our OWN
        // derivation (composed: Kron / concat / recombination from bases) represents the cell
        // over an import — a transparent, reproducible lineage rather than crediting an importer
        // for a value we derive ourselves (⟨4,4,4⟩=49 is Strassen², not an AlphaTensor discovery).
        const aComp = isComposedScheme(a) ? 0 : 1;
        const bComp = isComposedScheme(b) ? 0 : 1;
        if (aComp !== bComp) return aComp < bComp;
        return sourcePriority(a.source) < sourcePriority(b.source);
      };
      // Best representative per (format, commutative, field) cell.
      const repByCell = new Map();
      // Border-rank rows live in their own cells (a separate quantity) so a
      // theoretical R̃ bound never wins/hides an exact-rank constructive scheme.
      const cellKey = (s, f) => `${s.format.join(",")}|${s.commutative === true}|${s.border === true}|${s.lower === true}|${f}`;
      for (const s of filtered) {
        for (const f of cellFieldsOf(s)) {
          const k = cellKey(s, f);
          const cur = repByCell.get(k);
          if (cur == null || better(s, cur)) repByCell.set(k, s);
        }
      }
      const keep = new Set(repByCell.values());
      // Rows with no classifiable field can't win a cell — keep them rather
      // than silently dropping (defensive; catalog rows always carry fields).
      //
      // Default view is STRICTLY one best-rank row per (format, commutative,
      // field) cell (user 2026-06-07): a higher-rank scheme that is merely
      // (rank, buds) Pareto-best — a bud-richer building block — is NOT a
      // representative of any field and is hidden until "Show history". This
      // reverses the earlier §serendipitous carve-out that surfaced such rows
      // by default; the buds column is still inspectable via the history view.
      filtered = filtered.filter(s =>
        keep.has(s) || schemeFieldList(s).length === 0);
    }

    // Note: a further "shave noise" pass that drops rediscoveries with worse
    // adds counts is applied at CATALOG BUILD TIME inside
    // GenerateCatalogManifest.shaveByBestAdditions, so the manifest itself
    // is pre-filtered. The JS-side dedup above only handles the
    // catalog-source vs cited-bound and same-(rank,adds) cross-source cases.
  }

  // Field selector natural order — must match the <option> sequence in
  // index.html so the sort feels natural to someone reading the dropdown
  // top-to-bottom.
  // R/Q/Z is the unnarrowed catch-all; sort it LAST after the narrowed
  // characteristic-0 tags so the table reads as "narrowed first, catch-all
  // last". Matches the selector ordering.
  const FIELD_NATURAL_ORDER = ["F2", "Z", "Q", "R", "C", "R/Q/Z"];
  const fieldRank = f => {
    const i = FIELD_NATURAL_ORDER.indexOf(f);
    return i === -1 ? FIELD_NATURAL_ORDER.length : i;
  };
  // Sort key for a row now that the singular `field` is gone: rank it by the
  // lowest-natural-order field in its `fields` membership.
  const rowFieldRank = s => {
    let best = FIELD_NATURAL_ORDER.length;
    for (const f of schemeFieldList(s)) best = Math.min(best, fieldRank(f));
    return best;
  };
  // Order/group by the NORMALISED (sorted) shape. 1143 catalog entries carry a
  // genuinely unsorted format (e.g. ⟨11,17,2⟩ — the shape order is meaningful and
  // kept for display), so format[2] is NOT always the max axis. Precompute a
  // sorted [min,mid,max] per row so max-dim grouping + lex ordering are correct
  // regardless of stored axis order (fixes ⟨2,2,2⟩ → ⟨3,17,2⟩ → ⟨3,19,2⟩ jumps).
  for (const s of filtered) s.__s = [...s.format].sort((x, y) => x - y);
  filtered.sort((a, b) => {
    switch (sort) {
      // Back-compat alias: bookmarks with "sort=max_dim" still work but use
      // the new multi-key behaviour.
      case "max_dim":
      case "format":
        // Group by max(n,m,p) ascending, then lex within each group:
        //   ⟨2,2,2⟩
        //   ⟨2,2,3⟩, ⟨2,3,3⟩, ⟨3,3,3⟩
        //   ⟨2,2,4⟩, ⟨2,3,4⟩, ⟨2,4,4⟩, ⟨3,3,4⟩, ⟨3,4,4⟩, ⟨4,4,4⟩
        //   ⟨2,2,5⟩, … ⟨5,5,5⟩
        // Formats are stored sorted (n ≤ m ≤ p), so format[2] is the max.
        return a.__s[2] - b.__s[2]
            || a.__s[0] - b.__s[0]
            || a.__s[1] - b.__s[1]
            || rowFieldRank(a) - rowFieldRank(b)
            // NC (commutative=false) before commutative=true
            || (a.commutative === true ? 1 : 0) - (b.commutative === true ? 1 : 0)
            // Lower rank = better (current SOTA) appears first
            || a.rank - b.rank
            // Within tied rank: the user-selected secondary key (default adds).
            || (secondaryRankScore(a, secondary) - secondaryRankScore(b, secondary));
      case "format_lex":
        // Pure lexicographic order on the shape tuple ⟨n,m,p⟩ as displayed:
        //   ⟨2,2,2⟩, ⟨2,2,3⟩, …, ⟨2,2,16⟩, ⟨2,3,2⟩, ⟨2,3,3⟩, …, then ⟨3,…⟩, …
        // Unlike "format", this does NOT group by max(n,m,p): the FIRST axis is
        // the primary key, so every ⟨2,·,·⟩ shape precedes any ⟨3,·,·⟩ shape.
        // Uses the stored/displayed axis order (not the sorted [min,mid,max]),
        // so an unsorted shape like ⟨2,3,2⟩ sits where it reads.
        // NOTE: requested specifically for Perminov — it mirrors the per-format
        // lexicographic ordering of his FastMatrixMultiplication catalog, so the
        // two can be scanned in lockstep.
        return a.format[0] - b.format[0]
            || a.format[1] - b.format[1]
            || a.format[2] - b.format[2]
            || rowFieldRank(a) - rowFieldRank(b)
            || (a.commutative === true ? 1 : 0) - (b.commutative === true ? 1 : 0)
            || a.rank - b.rank
            || (secondaryRankScore(a, secondary) - secondaryRankScore(b, secondary));
      case "rank":
        return (a.rank - b.rank)
            || (secondaryRankScore(a, secondary) - secondaryRankScore(b, secondary));
      case "omega_asc": {
        // Ascending implied ω (lower = closer to the ω=2 ideal = "better").
        // Rows whose ω is undefined (e.g. n·m·p = 1) sort to the END.
        const aw = impliedOmega(a.format, a.rank);
        const bw = impliedOmega(b.format, b.rank);
        const aMissing = (aw == null);
        const bMissing = (bw == null);
        if (aMissing && !bMissing) return 1;
        if (!aMissing && bMissing) return -1;
        if (aMissing && bMissing) {
          return a.__s[2] - b.__s[2] || a.rank - b.rank;
        }
        return (aw - bw)
            || a.__s[2] - b.__s[2]
            || a.__s[0] - b.__s[0]
            || a.__s[1] - b.__s[1]
            || a.rank - b.rank;
      }
      case "buds_desc": {
        // Bud-richness descending (serendipitous potential). Schemes with no
        // buds (or none known) sort to the END; ties broken by format then rank.
        const ab = a.bud_score || 0;
        const bb = b.bud_score || 0;
        return (bb - ab)
            || a.__s[2] - b.__s[2]
            || a.__s[0] - b.__s[0]
            || a.__s[1] - b.__s[1]
            || a.rank - b.rank;
      }
      case "year":
      case "year_asc": {
        // Rows with unknown year always sort to the END (both directions).
        // For desc: known years come in descending order, nulls last.
        // For asc:  known years come in ascending order, nulls last.
        const aMissing = (a.year == null);
        const bMissing = (b.year == null);
        if (aMissing && !bMissing) return 1;
        if (!aMissing && bMissing) return -1;
        if (aMissing && bMissing) {
          // Stable tiebreak by format/rank when both years missing.
          return a.__s[0] - b.__s[0]
              || a.__s[1] - b.__s[1]
              || a.__s[2] - b.__s[2]
              || a.rank - b.rank;
        }
        const cmp = (sort === "year_asc") ? (a.year - b.year) : (b.year - a.year);
        if (cmp !== 0) return cmp;
        return a.__s[0] - b.__s[0]
            || a.__s[1] - b.__s[1]
            || a.__s[2] - b.__s[2]
            || a.rank - b.rank;
      }
    }
    return 0;
  });

  document.getElementById("result-count").textContent =
    `Showing ${filtered.length} of ${allSchemes.length} schemes`;
  // Refresh footer references to only those used by currently displayed rows.
  // (Built up as we iterate the filtered list below; re-render after the loop.)

  // Pre-pass: compute omega range for heatmap calibration on the current selection.
  let minOmega = null, maxOmega = null;
  for (const s of filtered) {
    const w = impliedOmega(s.format, s.rank);
    if (w == null) continue;
    if (minOmega == null || w < minOmega) minOmega = w;
    if (maxOmega == null || w > maxOmega) maxOmega = w;
  }

  const body = document.getElementById("results-body");
  body.innerHTML = "";
  const frag = document.createDocumentFragment();
  const usedRefKeys = new Set();
  for (const s of filtered) {
    const tr = document.createElement("tr");
    const cls = classifySource(s.source);
    // Lower-bound rows are registered in refRegistry under their RAW source
    // string (buildRefRegistry keys LB entries by lb.source verbatim), whereas
    // classifySource strips a trailing "(via …)" and requires a year — so a
    // source like "Bläser (via Landsberg 2008)" classified to "Bläser" missed
    // its registry entry and rendered no [N] (while "Bläser 2003" worked only
    // because classifySource left it unchanged). For LB rows, prefer the
    // raw-source key when the registry holds it.
    const refKey = (s.lower && refRegistry.has(s.source)) ? s.source : cls.key;
    if (refKey) usedRefKeys.add(refKey);
    const ref = refKey ? refRegistry.get(refKey) : null;
    const refLink = ref
      ? ` <sup class="refnum"><a href="#ref-${ref.index}">[${ref.index}]</a></sup>`
      : "";
    let fileCell, sourceCell;
    if (s.derived) {
      tr.classList.add("derived");
      sourceCell = `<span class="derived-tag" title="${escapeHtml(s.source)}">${escapeHtml(cls.label)}</span>${refLink}`;
      if (s.verified === false) {
        tr.classList.add("derived-unverified");
        fileCell = `<span class="muted" style="color:#a55" title="${escapeHtml(s.construction || s.source || "")}">⚠ rank only — no constructor yet</span>`;
      } else {
        fileCell = `<span class="muted" title="${escapeHtml(s.construction || "")}">↻ construct on demand</span>`;
      }
    } else if (s.cited) {
      tr.classList.add("cited");
      sourceCell = `<span class="cited-tag" title="${escapeHtml(s.source)}${s.notes ? " — " + escapeHtml(s.notes) : ""}">${escapeHtml(cls.label)}</span>${refLink}`;
      fileCell = `<span class="muted">📖 ${s.lower ? "lower bound (no scheme)" : "rank claim (no scheme)"}</span>` + sourceSchemeLink(s);
    } else {
      // For COMPOSED schemes show the lineage formula in place of the
      // "Solven-..." auto-attribution — the construction itself IS the
      // source. Atoms keep their historical attribution (Strassen 1969,
      // Laderman 1976, …).
      // Formula check FIRST: a projected formula ("DIS09Lemma4(n=30) ↓[…]")
      // also matches isComposedLineage (via ↓), but its head is a closed-form
      // constructor whose "disc." credit must survive the projection.
      if (isFormulaLineage(s.lineage_compact)) {
        // Formula-DERIVED scheme (e.g. DIS09Lemma4(n=20)): we materialise it
        // ourselves from a published closed-form, so the construction (lineage)
        // IS the source — but the formula/rank was DISCOVERED by the cited paper,
        // so we keep that credit as a "disc." prior-art highlight (user 2026-06-08).
        sourceCell = `<code class="lineage-formula" title="${escapeHtml(s.lineage_compact)}${s.source ? " — " + escapeHtml(s.source) : ""}">${escapeHtml(collapseKeepLists(s.lineage_compact))}</code>`
          + ` <span class="disc-tag" title="discovered by ${escapeHtml(cls.label)}">disc. ${escapeHtml(cls.label)}</span>${refLink}`;
      } else if (isComposedLineage(s.lineage_compact)) {
        // Show the SHAPE-ONLY composition (e.g. "3x5x5 ⊗ 5x5x7"); the base's
        // source/rank/adds are noise in this structural view — drop them, keeping
        // the full ref in the hover tooltip (user 2026-06-08).
        const shapes = collapseKeepLists(shapesOnlyLineage(s.lineage_compact));
        sourceCell = `<code class="lineage-formula" title="${escapeHtml(s.lineage_compact)}">${escapeHtml(shapes)}</code>`;
      } else {
        sourceCell = `<span title="${escapeHtml(s.source)}">${escapeHtml(cls.label)}</span>${refLink}`;
      }
      fileCell = `<a href="${SCHEMES_GITHUB_BASE}${s.file}" target="_blank" rel="noopener" title="${escapeHtml(s.file)}">JSON</a>` + sourceSchemeLink(s);
    }
    // Field cell — if the user filtered to C but this row is R-class, mark it as a fallback.
    // Also surface commutative status if set, and the narrowed-from-source field tag
    // (Q / Z) so users can see why a Q-only filter selected this row.
    // Display rewrite: the stored value "R/Q/Z" (legacy catch-all default
    // when narrowing hasn't been computed) is shown as "Z/Q/R" to reflect
    // the inclusion order Z ⊂ Q ⊂ R and to make the catch-all visually
    // distinct from a properly-narrowed Z / Q / R cell.
    // Field cell (#182): show the clean canonical field LIST the scheme is
    // valid over (e.g. "Z, Q, R") instead of the old combined "Z/Q/R ·Z"
    // string-plus-narrow-badge. The narrowed-from-source tag (if any) moves
    // into the tooltip so the cell stays readable.
    const flist = schemeFieldList(s);
    // Display only the NARROWEST field(s); the implied ones (e.g. a Z scheme's
    // Q/R/C/F2/F3) would be noise in the column. Full list → tooltip.
    const canon = canonicalFields(flist);
    // ZT is a sub-class marker on integer schemes (coefficients all in {−1,0,1}),
    // carried as the per-scheme `zt` boolean — NOT a member of fields[] (per the
    // field conventions: ZT is not a field). Surface it in the column by narrowing
    // the displayed "Z" to "ZT" when the flag is set; "ZT" is the more informative
    // descriptor for a ternary-integer scheme. (zt is only ever true when Z ∈ fields.)
    const canonDisp = (s.zt === true) ? canon.map(f => (f === "Z" ? "ZT" : f)) : canon;
    const ztNote = (s.zt === true) ? " · ZT: ternary integer (coeffs ∈ {−1,0,1})" : "";
    const fieldText = canonDisp.length ? canonDisp.join(", ") : "—";
    const fullText = (flist.length ? flist.join(", ") : "—") + ztNote;
    const impliedText = (canon.length && flist.length > canon.length)
      ? ` (implies ${flist.filter(f => !canon.includes(f)).join(", ")})` : "";
    // The "narrowed-from-source" note is a legacy aid for rows lacking an
    // explicit `fields` membership. With `fields` authoritative (2026-06-04) the
    // list IS the narrow info, so only annotate the legacy (no-fields) case.
    const narrow = narrowField(s);
    const narrowNote = (!flist.length && narrow && narrow !== "R/Q/Z")
      ? ` (narrowed-from-source: ${narrow})` : "";
    let fieldDisplay = `<span class="field" title="valid over: ${escapeHtml(fullText)}${impliedText}${narrowNote}">${escapeHtml(fieldText)}</span>`;
    // Lower bounds are published per field. A field-agnostic ("all") floor holds
    // in any field — show "any" rather than the canonical "Z", which would
    // understate it; a field-specific floor shows exactly the published field.
    if (s.lower && s.lb_field) {
      const lbFieldText = s.lb_field === "all" ? "any" : lbFieldList(s.lb_field).join(", ");
      const lbFieldTip = s.lb_field === "all" ? "any field (field-agnostic argument)"
        : s.lb_field === "char0" ? "characteristic 0 (Z, Q, R, C) — does NOT hold over F2/F3"
        : `valid over ${lbFieldList(s.lb_field).join(", ")}`;
      fieldDisplay = `<span class="field" title="lower bound published over: ${escapeHtml(lbFieldTip)}">${escapeHtml(lbFieldText)}</span>`;
    }
    if (field === "C" && !flist.includes("C") && schemeValidForRequestedField(s.field, "C", s)) {
      fieldDisplay = `<span class="field" title="R ⊂ C: R-class scheme is automatically valid over C${narrowNote}">${escapeHtml(fieldText)} <small class="field-fallback">↗ C</small></span>`;
    }
    if (s.commutative === true) {
      fieldDisplay += ` <span class="commutative-tag" title="commutative-only rank — does NOT lift to recursive matmul over non-commutative entries">cmt</span>`;
    }
    // Rank cell: show NC rank + commutative rank when meaningfully different.
    let rankCell;
    if (s.lower) {
      rankCell = `≥${s.rank} <span class="commutative-tag" title="lower bound on the bilinear rank — a published theoretical floor, NOT an achievable scheme${s.tight ? "; TIGHT: a matching construction exists in the catalog" : ""}; shown separately from achievable ranks">LB${s.tight ? " ✓tight" : ""}</span>`;
    } else if (s.border) {
      rankCell = `≤${s.rank} <span class="commutative-tag" title="border rank R̃ — an asymptotic/degenerate bound, NOT an exact-multiplication scheme; shown separately from achievable ranks">border</span>`;
    } else {
      rankCell = `×${s.rank}`;
    }
    // commutative_rank is now split into a separate row at load time —
    // no longer rendered inline.
    // Omega cell with heatmap (6-decimal precision).
    const omega = impliedOmega(s.format, s.rank);
    const omegaText = omega == null ? "—" : omega.toFixed(6);
    const omegaStyle = omega == null ? "" : `style="background:${omegaHeatColor(omega, minOmega, maxOmega)}"`;
    const lineageCell = classifyLineage(s);
    const budsCell = budStructureCell(s);
    tr.innerHTML = `
      <td class="format">⟨${s.format[0]},${s.format[1]},${s.format[2]}⟩</td>
      <td>${fieldDisplay}</td>
      <td class="rank" ${omegaStyle} title="ω = 3·ln(rank)/ln(n·m·p)">${omegaText}</td>
      <td class="rank">${rankCell}</td>
      <td class="rank">${additionsCell(s)}</td>
      <td class="rank">${budsCell}</td>
      <td title="${escapeHtml(lineageCell.tooltip)}">${lineageCell.label}</td>
      <td>${sourceCell}</td>
      <td>${fileCell}</td>
      <td class="copy-cell"><button class="copy-btn" title="Copy a one-line identifier of this row to the clipboard" aria-label="Copy row identifier">📋</button></td>
    `;
    // Copy-to-clipboard: a one-liner uniquely identifying this entry. stopPropagation
    // so it doesn't also trigger the row→inspect-modal handler below.
    const copyBtn = tr.querySelector(".copy-btn");
    if (copyBtn) {
      copyBtn.addEventListener("click", (e) => {
        e.stopPropagation();
        copyEntryLine(s, copyBtn);
      });
    }
    // Click a row (anywhere but a link) to inspect the raw catalog entry (#186).
    tr.style.cursor = "pointer";
    tr.title = "Click to inspect the raw catalog entry";
    tr.addEventListener("click", (e) => {
      if (e.target.closest("a")) return; // let JSON / reference links work normally
      openEntryModal(s);
    });
    frag.appendChild(tr);
  }
  body.appendChild(frag);
  renderReferences(usedRefKeys);
}

// A one-line, uniquely-identifying string for a catalog row — copied to the
// clipboard by the per-row 📋 button. Carries shape, field membership, rank,
// additions, commutativity, lineage kind, source, and the on-disk file (or
// lineage formula for derived/cited rows) so a reader can paste it into an
// issue / commit / chat and unambiguously refer to this exact entry.
function entryOneLiner(s) {
  const fmt = `⟨${s.format.join(",")}⟩`;
  const fields = schemeFieldList(s).join("/") || s.field || "?";
  const eff = effectiveAdditions(s);
  const adds = eff != null ? ` a=${eff}${s.additions != null && s.additions > eff ? `(${s.additions} flat)` : ""}` : "";
  const cmt = s.commutative === true ? " commutative" : "";
  const kind = s.cited ? "cited" : s.derived ? "derived" : "scheme";
  const tail = s.file || s.lineage_compact || "";
  return `${fmt} ${fields} r=${s.rank}${adds}${cmt} [${kind}] ${s.source || "?"}${tail ? " :: " + tail : ""}`;
}

function copyEntryLine(s, btn) {
  const text = entryOneLiner(s);
  const flash = () => {
    btn.textContent = "✓";
    btn.classList.add("copied");
    setTimeout(() => { btn.textContent = "📋"; btn.classList.remove("copied"); }, 1000);
  };
  if (navigator.clipboard && navigator.clipboard.writeText) {
    navigator.clipboard.writeText(text).then(flash).catch(() => fallbackCopy(text, flash));
  } else {
    fallbackCopy(text, flash);
  }
}

function fallbackCopy(text, done) {
  const ta = document.createElement("textarea");
  ta.value = text;
  ta.style.position = "fixed";
  ta.style.opacity = "0";
  document.body.appendChild(ta);
  ta.select();
  try { document.execCommand("copy"); } catch (e) { /* ignore */ }
  document.body.removeChild(ta);
  done();
}

// External community-catalog deep links for a given shape. FMM-Lille publishes
// a canonical per-format page (sorted n≤m≤p convention), so we link it directly.
// Perminov's repo path needs the field subdir (Z/Q/ZT) + exact rank token which
// we can't reliably reconstruct here, so we link a repo-scoped code search for
// the shape — robust ("as best as we can given the trickier syntax").
function externalCatalogLinks(format, rank) {
  if (!Array.isArray(format) || format.length !== 3) return "";
  const sorted = [...format].sort((a, b) => a - b);
  const [a, b, c] = sorted;
  const fmmUrl = `https://fmm.univ-lille.fr/${a}x${b}x${c}.html`;
  // Search the Perminov repo for files matching this shape (any field / rank).
  const perminovQuery = encodeURIComponent(
    `repo:dronperminov/FastMatrixMultiplication ${a}x${b}x${c}_`);
  const perminovUrl = `https://github.com/search?q=${perminovQuery}&type=code`;
  return (
    `<span style="opacity:.7">external:</span> ` +
    `<a href="${fmmUrl}" target="_blank" rel="noopener" title="FMM-Lille per-format page (n≤m≤p)">FMM-Lille ↗</a>` +
    ` &nbsp;·&nbsp; ` +
    `<a href="${perminovUrl}" target="_blank" rel="noopener" title="Search Perminov FastMatrixMultiplication for this shape">Perminov ↗</a>`
  );
}

// Bud-structure table cell: compact `bud_score` + certainty badge + a ★ when
// the scheme is on the (rank, buds) Pareto frontier. The full bud multiset
// (e.g. "576×U⟨1,1,2⟩ + 246×U⟨1,1,3⟩ + …") is shown on hover to keep the
// column narrow. "—" when the scheme has no buds (or none known yet).
function budStructureCell(s) {
  const star = s.pareto_rank_buds
    ? `<span title="on the (rank, buds) Pareto frontier — registered 'best' building block even if not rank-minimal">★</span> `
    : "";
  if (!s.has_buds || !s.buds) {
    // Still show the ★ for rank-minimal frontier members with zero buds.
    return star ? `${star}<span class="muted">—</span>` : "—";
  }
  const cert = s.buds_certainty === "structural-estimate"
    ? ` <small class="field-narrow" title="structural estimate — a lower bound inferred through a concat/estimate op (not certified)">est</small>`
    : ` <small class="field-fallback" title="exact — from a cancellation-free lineage (Kron/flip/transpose) or an expanded scheme">exact</small>`;
  return `${star}<span title="${escapeHtml(s.buds)}">${s.bud_score || 0}</span>${cert}`;
}

// Raw-catalog-entry modal (#186). Shows the precomputed catalog.json entry
// verbatim (no SPA-side computation — user principle). Opened by clicking a row.
function openEntryModal(s) {
  let ov = document.getElementById("entry-modal");
  if (ov) ov.remove();
  ov = document.createElement("div");
  ov.id = "entry-modal";
  ov.style.cssText =
    "position:fixed;inset:0;background:rgba(0,0,0,.5);display:flex;align-items:center;" +
    "justify-content:center;z-index:1000";
  ov.addEventListener("click", (e) => { if (e.target === ov) ov.remove(); });
  const box = document.createElement("div");
  box.style.cssText =
    "background:#fff;max-width:820px;width:90%;max-height:80vh;overflow:auto;padding:1em 1.2em;" +
    "border-radius:6px;font-family:monospace;font-size:12px;box-shadow:0 4px 24px rgba(0,0,0,.35)";
  const title = `⟨${s.format.join(",")}⟩  ×${s.rank}`;
  const fileLink = (s.file
    ? `<a href="${SCHEMES_GITHUB_BASE}${escapeHtml(s.file)}" target="_blank" rel="noopener">${escapeHtml(s.file)}</a>`
    : "<em>(no on-disk scheme — derived/cited)</em>")
    + (s.source_scheme_url
      ? ` &nbsp;·&nbsp; source file: <a href="${escapeHtml(s.source_scheme_url)}" target="_blank" rel="noopener">${escapeHtml(s.source_scheme_url.replace(/^https?:\/\//, ""))}</a>`
      : "")
    + (s.source_paper_url
      ? ` &nbsp;·&nbsp; paper: <a href="${escapeHtml(s.source_paper_url)}" target="_blank" rel="noopener">${escapeHtml(s.source_paper_url.replace(/^https?:\/\//, ""))}</a>`
      : "");
  // Composed schemes get a MermaidJS lineage tree above the raw JSON; atoms /
  // formula / cited rows have nothing to draw so the block is omitted.
  const hasLineage = isComposedLineage(s.lineage_compact);
  const lineageBlock = hasLineage
    ? `<div style="margin-bottom:.6em">` +
      `<div style="font-weight:bold;margin-bottom:.3em">Lineage ` +
      `<span style="font-weight:normal;color:#888;font-size:.85em">— root badge = this scheme's actual rank; leaf badge = the base's rank (≈ = best-known when the exact base isn't catalogued); click a node to open its scheme</span></div>` +
      `<code class="lineage-formula" style="display:block;margin-bottom:.4em">${escapeHtml(s.lineage_compact)}</code>` +
      `<div class="lineage-graph" id="lineage-graph-box">rendering lineage…</div></div>`
    : "";
  // Pan-TA highlight: when this recombination fuses cyclic-rotation product pairs,
  // show WHERE the rank was lowered (rank = unpaired leaves + fused-pair cost; TA
  // saved `saving`). Makes the final rank explainable, per the catalog's `ta_fusion`.
  const ta = s.ta_fusion;
  const taFusionBlock = ta
    ? `<div style="margin-bottom:.6em;padding:.5em .6em;background:#eef5ff;border-left:3px solid #3b6fb0;border-radius:3px">` +
      `<div style="font-weight:bold;margin-bottom:.3em">Pan trilinear aggregation (TA) ` +
      `<span style="font-weight:normal;color:#888;font-size:.85em">— a saving WITHIN this recombination's multiplications: each fused cyclic-rotation pair costs abc+ab+bc+ca instead of 2·R</span></div>` +
      `<div style="margin-bottom:.3em">${escapeHtml(ta.summary || "")}</div>` +
      `<div style="color:#555;margin-bottom:.3em">${ta.pairs} pair(s) fused · saved <strong>${ta.saving}</strong> multiplications · ${ta.unpaired_leaf_sum} unpaired leaves + ${ta.fused_cost} fused = ${s.rank}</div>` +
      (Array.isArray(ta.fused)
        ? `<ul style="margin:0;padding-left:1.2em">` + ta.fused.map(fp =>
            `<li>⟨${fp.shape.join(",")}⟩ &amp; its rot² : fused ${fp.fused_cost} vs naïve ${fp.naive_rank} <span style="color:#2a7">(save ${fp.saving})</span></li>`).join("") + `</ul>`
        : "") +
      `</div>`
    : "";
  box.innerHTML =
    `<div style="display:flex;justify-content:space-between;align-items:center;margin-bottom:.5em">` +
    `<strong>${title}</strong>` +
    `<button id="entry-modal-close" style="cursor:pointer;border:none;background:none;font-size:1.2em">✕</button></div>` +
    `<div style="margin-bottom:.6em">${fileLink}</div>` +
    `<div style="margin-bottom:.6em">${externalCatalogLinks(s.format, s.rank)}</div>` +
    lineageBlock +
    taFusionBlock +
    `<pre style="white-space:pre-wrap;word-break:break-word;margin:0">${escapeHtml(JSON.stringify(s, null, 2))}</pre>`;
  ov.appendChild(box);
  box.querySelector("#entry-modal-close").addEventListener("click", () => ov.remove());
  document.body.appendChild(ov);

  // Render the Mermaid diagram once the modal is in the DOM. renderLineageGraph
  // lazily fetches the scheme's own JSON for the structured DAG (exact root rank +
  // hash-resolved leaves) and falls back to the compact-string renderer on failure.
  if (hasLineage) {
    renderLineageGraph(box, s);
  }
}

// Field-implications modal (#field-implications-modal). Reminds the reader of the
// full inclusion map behind the Field filter — selecting a field includes every
// scheme valid over a NARROWER field. Replaces the old per-option "implies Q, R"
// hints (which listed only part of the chain) and the dropped "Z/Q/R" group.
function openFieldImplicationsModal() {
  let ov = document.getElementById("field-implications-modal");
  if (ov) ov.remove();
  ov = document.createElement("div");
  ov.id = "field-implications-modal";
  ov.style.cssText =
    "position:fixed;inset:0;background:rgba(0,0,0,.5);display:flex;align-items:center;" +
    "justify-content:center;z-index:1000";
  ov.addEventListener("click", (e) => { if (e.target === ov) ov.remove(); });
  const box = document.createElement("div");
  box.style.cssText =
    "background:#fff;max-width:680px;width:90%;max-height:80vh;overflow:auto;padding:1em 1.4em;" +
    "border-radius:6px;font-size:14px;line-height:1.5;box-shadow:0 4px 24px rgba(0,0,0,.35)";
  box.innerHTML =
    `<div style="display:flex;justify-content:space-between;align-items:center;margin-bottom:.5em">` +
    `<strong>Field implications</strong>` +
    `<button id="field-impl-close" style="cursor:pointer;border:none;background:none;font-size:1.2em">✕</button></div>` +
    `<p>Selecting a field returns every scheme valid over it — including schemes ` +
    `tagged with a <em>narrower</em> field, since a narrower field implies the wider ones.</p>` +
    `<ul style="margin:.4em 0 .8em;padding-left:1.2em">` +
    `<li><b>Characteristic 0:</b> <code>Z ⊂ Q ⊂ R ⊂ C</code> — an integer (Z) scheme is also ` +
    `valid over ℚ, ℝ, ℂ; a rational (Q) scheme over ℝ, ℂ; a real (R) scheme over ℂ.</li>` +
    `<li><b>Integer → finite fields (theorem):</b> a <code>Z</code> scheme reduces mod p, so it is ` +
    `also valid over <code>F2</code> and <code>F3</code>.</li>` +
    `<li><b>Rational → finite fields (conditional):</b> a <code>Q</code> scheme reduces mod p ` +
    `<em>only</em> when every denominator is coprime to p — so <code>F2</code>/<code>F3</code> ` +
    `validity is a per-scheme fact, <em>not</em> implied by <code>Q</code>.</li>` +
    `<li><b>F2, F3:</b> independent characteristic-p universes; neither implies the other, and the ` +
    `char-0 chain does not imply them (except from <code>Z</code> via the reduction theorem above).</li>` +
    `<li><b>ZT</b> is <em>not</em> a field — it is the sub-class of <code>Z</code> schemes whose ` +
    `coefficients are all in {−1, 0, +1} (“ternary integer”). Selecting it shows <code>Z</code> ` +
    `schemes with the ZT flag set.</li>` +
    `</ul>` +
    `<p class="muted">The <b>Field</b> column shows only the <b>narrowest</b> field(s) a scheme ` +
    `is valid over — the implied wider fields are omitted to reduce noise (hover the cell for the ` +
    `full list). See the “Fields covered” legend below the table for what each algebra means.</p>`;
  ov.appendChild(box);
  box.querySelector("#field-impl-close").addEventListener("click", () => ov.remove());
  document.body.appendChild(ov);
}

(function wireFieldImplicationsButton() {
  const btn = document.getElementById("fields-impl-btn");
  if (btn) btn.addEventListener("click", openFieldImplicationsModal);
})();

// Esc closes the entry modal (bound once).
document.addEventListener("keydown", (e) => {
  if (e.key === "Escape") {
    const m = document.getElementById("entry-modal");
    if (m) m.remove();
    const fm = document.getElementById("field-implications-modal");
    if (fm) fm.remove();
  }
});

/**
 * Classify the lineage of a scheme row for display.
 *
 * Inputs available per row:
 *  - `s.derived` (boolean): set for derived-from-cited-bounds.json entries (formula-only,
 *    no factor matrices on disk).
 *  - `s.cited` (boolean): set for cited-bounds.json entries (literature claim
 *    without an explicit scheme).
 *  - `s.source` (string): canonical source label. Catalog schemes whose source
 *    string starts with `Composed-` (e.g. `Composed-recursive`,
 *    `Composed-strassen-recombine`, `Composed-kron`) are derived from smaller
 *    schemes via the named composition mechanism.
 *
 * Output: { label: HTML, tooltip: text }. The label is short for the column;
 * the tooltip carries the full lineage description.
 *
 * TODO: once `GenerateCatalogManifest` emits the full `lineage` JSON sub-tree
 * into catalog.json, replace this string-pattern heuristic with a structural
 * walk of `s.lineage.op` and the parent ref(s).
 */
// A lineage_compact is "composed" if it contains any operator: ⊗ (Kronecker),
// ⊕ (disjoint sum), +p / +n / +m (concat), ↓[…] (projection / row-col drop),
// R[…] / Rta[…] / R*[…] (recombination), TA[…] (peel-via-TA), AS( (augment-square),
// AxisFlip, AxisPermute, DisjointSum.
// A bare ref (e.g. "Strassen<2,2,2>=7", "3x3x6_m40_a862", "winograd-1971-7mult")
// is an Atom. A bare re-orientation ("9x10x13→⟨10,9,13⟩") stays an Atom — the
// scheme is the same tensor, just re-oriented.
function isComposedLineage(compact) {
  if (!compact) return false;
  return /[⊗◊⊕↓]|[+]p\b|[+]n\b|[+]m\b|^R\[|^Rta\[|^R\*\[|^TA\[|^AS\(|^Kron\[|^Derived\[|^Composed\[|AxisFlip|AxisPermute|DisjointSum|->NPM|->PMN|->NMP/.test(compact);
}

// A formula-DERIVED lineage: a single parametric closed-form constructor call
// such as "DIS09Lemma4(n=20)" (we materialise it ourselves from a published
// formula). NOT a composition (Kron/concat — handled by isComposedLineage) and
// not a bare scheme reference. Rendered as a derived construction WITH a
// "discovered by <paper>" credit (user 2026-06-08).
// Reduce a composed lineage to SHAPES ONLY: each atom token (a source-prefixed
// canonical key like "alphatensor_Z-3x5x5_m58_a369") collapses to its bare shape
// "3x5x5"; operators (⊗, +p/+n, c=/r=, parens, spaces) are left intact because
// they contain no NxMxP shape and so never match. (user 2026-06-08)
function shapesOnlyLineage(compact) {
  if (!compact) return compact;
  return compact.replace(/[A-Za-z0-9_-]*\d+x\d+x\d+[A-Za-z0-9_-]*/g, (tok) => {
    const m = tok.match(/(\d+x\d+x\d+)/);
    return m ? m[1] : tok;
  });
}

// Strip trailing STRUCTURAL suffixes — projection "↓[keepN|keepM|keepP]" and
// re-orientation "→⟨n,m,p⟩" — so a projected/re-oriented formula (e.g.
// "DIS09Lemma4(n=30) ↓[…]", the ⟨29,30,30⟩ row) still classifies by its
// constructive head. Suffixes can stack, so strip until fixpoint.
function lineageHead(compact) {
  let c = compact.trim(), prev;
  do {
    prev = c;
    c = c.replace(/\s*↓\[[^\]]*\]$/, "").replace(/\s*→⟨\d+,\d+,\d+⟩$/, "").trim();
  } while (c !== prev);
  return c;
}

function isFormulaLineage(compact) {
  if (!compact) return false;
  const c = lineageHead(compact);
  if (isComposedLineage(c)) return false;
  // A single parametric closed-form constructor call. Restricted to known
  // published formula families (not composition ops, which render infix and are
  // caught by isComposedLineage) so a future "Op(...)" form can't slip through.
  return /^(DIS09Lemma4|PanTrilinear\w*|Waksman\w*|Rosowski\w*|Makarov\w*|Sedoglavic\w*|Pan\w*)\([^)]*\)$/.test(c);
}

// Collapse each verbose projection keep-list "↓[0,1,2,…|…|…]" to a bare "↓"
// for the source cell — the target shape is already the row's format column,
// and the full keep-lists stay available in the hover tooltip.
function collapseKeepLists(compact) {
  return compact.replace(/\s*↓\[[^\]]*\]/g, " ↓").trim();
}

// The clean field LIST for a row (#182). The manifest (GenerateCatalogManifest)
// ALWAYS emits a canonically-ordered fields[] now — the SPA must NOT recompute
// it (user 2026-06-03: "the SPA should barely ever compute"). We render fields[]
// verbatim; the s.field split is a last-resort guard for a stale manifest only.
function schemeFieldList(s) {
  if (Array.isArray(s.fields) && s.fields.length) return s.fields;
  if (s.field) {
    return s.field === "R/Q/Z" ? ["Z", "Q", "R"] : s.field.split("/").map(x => x.trim()).filter(Boolean);
  }
  return [];
}

// Field-inclusion map: "X implies Y" = a scheme valid over X is also valid over Y.
// Char-0 chain Z ⊂ Q ⊂ R ⊂ C; plus the integer-reduction theorem (a Z scheme
// reduces mod p, so it is also valid over F2 and F3). Q does NOT imply F2/F3 in
// general (only when every denominator is coprime to p — a per-scheme fact), so
// F2/F3 are deliberately absent from Q's implications.
const FIELD_IMPLIES = {
  Z: ["Q", "R", "C", "F2", "F3"],
  Q: ["R", "C"],
  R: ["C"],
  C: [],
  F2: [],
  F3: [],
};

// Reduce a field list to its minimal generators: drop any field that another
// present field already implies. So an integer scheme [Z,Q,R,C,F2,F3] → [Z]
// (Z implies all the rest), while a mod-reducible rational [Q,R,C,F2,F3] → [Q,F2,F3]
// (Q implies R,C but not F2,F3). The full list stays available for the tooltip.
function canonicalFields(flist) {
  const present = new Set(flist);
  const implied = new Set();
  for (const g of flist) {
    for (const x of (FIELD_IMPLIES[g] || [])) {
      if (present.has(x)) implied.add(x);
    }
  }
  return flist.filter(f => !implied.has(f));
}

// Predicate for the "Composed" filter (#Q11): a scheme is composed if its
// lineage carries a composition operator OR its source is tagged "Composed"
// (the derived recursive-split / Kronecker rows). Atoms + imported schemes
// (Strassen, Perminov, AlphaTensor, FMM, …) are NOT composed.
function isComposedScheme(s) {
  if (isComposedLineage(s.lineage_compact)) return true;
  // Source provenance now reads "Derived_*" (was "Composed_*"); accept both so
  // older bookmarks / un-regenerated manifests still classify correctly.
  return /\b(derived|composed)\b/i.test(s.source || "");
}

function classifyLineage(s) {
  if (s.derived) {
    return { label: "<em>Formula</em>", tooltip: "Derived bound — no factor matrices on disk", composed: false };
  }
  if (s.cited) {
    return { label: "<em>Cited</em>", tooltip: "Upper-bound claim cited from literature without an explicit scheme", composed: false };
  }
  const compact = s.lineage_compact || "";
  if (isComposedLineage(compact)) {
    return {
      label: `<span class="lineage-atom-tag" title="Built (derived) by composing smaller schemes">Derived</span>`,
      tooltip: compact,
      composed: true,
    };
  }
  // Atom: either has a bare lineage_compact (canonical key or short label),
  // or has none and isn't a Composed-* source.
  return {
    label: `<span class="lineage-atom-tag" title="Primary scheme — not derived from others">Atom</span>`,
    tooltip: compact || s.source || "Primary scheme",
    composed: false,
  };
}

// ─── Lineage → MermaidJS (entry modal, #lineage-graph-box) ──────────────────
//
// Turn a `lineage_compact` string into a Mermaid flowchart so a composed
// scheme reads as a tree instead of a one-line operator soup. The grammar is
// the same one `isComposedLineage` recognises:
//
//   expr   := concat                              (lowest precedence)
//   concat := kron ( (' +p '|' +n '|' +m ') kron )*
//   kron   := factor ( ' ⊗ ' factor )*            (binds tighter than concat)
//   factor := 'R[' base ';' alloc ']'             (recombination)
//           | 'Composed[' mech ';' child ']'
//           | atom                                 (bare ref, e.g. perminov_ZT-4x5x6_m90_a1023)
//
// Atoms (Sedoglavic formulas, DIS09 lemmas, …) and Formula/Cited rows have no
// composition operator → lineageToMermaid returns null and the modal skips the
// graph. The string heuristic mirrors the TODO on classifyLineage: once the
// manifest emits a structured lineage sub-tree we can walk that directly.

// Split `s` on any separator in `seps` (longest-first), but only at top level
// — i.e. outside any '[' … ']' so an R[…] / Composed[…] payload stays intact.
// Returns { parts: trimmed operands, seps: the separators actually used }.
function splitTopLevel(s, seps) {
  const parts = [];
  const used = [];
  let depth = 0;
  let last = 0;
  let i = 0;
  while (i < s.length) {
    const ch = s[i];
    if (ch === "[") { depth++; i++; continue; }
    if (ch === "]") { depth--; i++; continue; }
    if (depth === 0) {
      const sep = seps.find((sp) => s.startsWith(sp, i));
      if (sep) {
        parts.push(s.slice(last, i));
        used.push(sep.trim());
        i += sep.length;
        last = i;
        continue;
      }
    }
    i++;
  }
  parts.push(s.slice(last));
  return { parts: parts.map((p) => p.trim()), seps: used };
}

// Escape text for a Mermaid quoted node/edge label: ASCII arrows and angle
// brackets become their unicode forms (Mermaid would otherwise mangle them or
// read them as markup even under securityLevel:loose), quotes are downgraded,
// and newlines become <br/> (htmlLabels is on).
function mmEsc(t) {
  return String(t)
    .replace(/->/g, "→")
    .replace(/</g, "⟨")
    .replace(/>/g, "⟩")
    .replace(/"/g, "'")
    .replace(/\n/g, "<br/>");
}

// ── Per-node metadata: shape + multiplication count ─────────────────────────
//
// To annotate every lineage node with the underlying scheme's metadata (the
// user-facing ask: "each node should indicate its number of multiplications")
// we resolve a shape and a multiplication count for each node:
//   - leaf  → shape parsed from the ref; rank from the ref's explicit token
//             (-r23-, _m90_) when present, else the best-known rank for the
//             shape across the catalog (orientation-agnostic fallback).
//   - ⊗     → shape = elementwise product of factors; rank = product of factors.
//   - +p/n/m→ shape = sum on the joined axis; rank = sum of operands.
//   - R[…]  → recombination keeps the matmul SHAPE but reduces the product
//             count by an amount not readable from the string, so we assert NO
//             rank for the recombine node (optimality-discipline: never print a
//             number we haven't derived).

// shape → best-known multiplication count, built once from the loaded catalog.
// Keyed both by exact orientation ("2x13x9") and by sorted shape ("2x9x13"):
// tensor rank is invariant under the 6 axis permutations, so the sorted key is a
// sound fallback when the exact orientation isn't catalogued.
let SHAPE_RANK_INDEX = null;
function shapeRankIndex() {
  if (SHAPE_RANK_INDEX) return SHAPE_RANK_INDEX;
  const exact = new Map();
  const sorted = new Map();
  for (const s of allSchemes || []) {
    if (!s.format || s.format.length !== 3 || s.rank == null) continue;
    // Skip non-achievable rows: a lower bound (rank ≥ …) or a border rank (R̃ ≤ …)
    // is NOT a multiplication count we can realise, so it must not become the
    // shape's "best known" rank used by lineage / omega annotations.
    if (s.lower === true || s.border === true) continue;
    const ek = s.format.join("x");
    if (!exact.has(ek) || s.rank < exact.get(ek)) exact.set(ek, s.rank);
    const sk = [...s.format].sort((a, b) => a - b).join("x");
    if (!sorted.has(sk) || s.rank < sorted.get(sk)) sorted.set(sk, s.rank);
  }
  SHAPE_RANK_INDEX = { exact, sorted };
  return SHAPE_RANK_INDEX;
}

// Parse ⟨n,m,p⟩ from a leaf ref ("5x5x6@f2ad9c", "16x16x16-direct", "⟨2,3,3⟩",
// "perminov_ZT-2x6x11_m103_a869") → [n,m,p], or null.
function shapeOfRef(ref) {
  let m = ref.match(/(\d+)x(\d+)x(\d+)/);
  if (m) return [+m[1], +m[2], +m[3]];
  m = ref.match(/⟨(\d+),(\d+),(\d+)⟩/);
  if (m) return [+m[1], +m[2], +m[3]];
  return null;
}

// Explicit, scheme-specific multiplication count from the ref token (-r23-,
// _m90_, _r45_), or null. Preferred over the shape lookup: it is the rank of
// THIS scheme, not merely the best known for the shape.
function rankOfRef(ref) {
  const m = ref.match(/[-_](?:m|r)(\d+)[-_]/);
  return m ? +m[1] : null;
}

// Best-known multiplication count for a shape, or null. {rank, oriented}:
// oriented=false means it came from the permutation-agnostic fallback.
function bestRankForShape(shape) {
  if (!shape) return null;
  const idx = shapeRankIndex();
  const ek = shape.join("x");
  if (idx.exact.has(ek)) return { rank: idx.exact.get(ek), oriented: true };
  const sk = [...shape].sort((a, b) => a - b).join("x");
  if (idx.sorted.has(sk)) return { rank: idx.sorted.get(sk), oriented: false };
  return null;
}

function sumOrNull(xs) { return xs.some((x) => x == null) ? null : xs.reduce((a, b) => a + b, 0); }
function productOrNull(xs) { return xs.some((x) => x == null) ? null : xs.reduce((a, b) => a * b, 1); }
function kronShape(shapes) {
  if (shapes.some((s) => !s)) return null;
  return shapes.reduce((a, s) => [a[0] * s[0], a[1] * s[1], a[2] * s[2]]);
}
// Resulting shape of a top-level concat: +p joins the p axis, +n the n axis,
// +m the m axis (the other two must match). null if any operand shape or
// separator is unknown.
function concatShape(shapes, seps) {
  if (shapes.some((s) => !s)) return null;
  const axisOf = { "+p": 2, "+n": 0, "+m": 1 };
  let acc = shapes[0].slice();
  for (let i = 1; i < shapes.length; i++) {
    const ax = axisOf[seps[i - 1]];
    if (ax == null) return null;
    acc = acc.slice();
    acc[ax] += shapes[i][ax];
  }
  return acc;
}

// Inline-HTML badges appended under a node label (htmlLabels:true renders them
// as real HTML; `title` becomes a native tooltip). Self-contained inline styles
// so no external CSS is required inside the Mermaid SVG.
function multBadge(rank, approx) {
  const bg = approx ? "#8a93a0" : "#1f6feb";
  const tip = approx ? "best-known multiplications for this shape" : "multiplications";
  return `<span style='display:inline-block;margin-top:3px;padding:0 6px;border-radius:8px;` +
    `background:${bg};color:#fff;font-size:.78em;font-weight:600' title='${tip}'>` +
    `${approx ? "≈×" : "×"}${rank}</span>`;
}
function shapeBadge(shape) {
  return `<span style='color:#555;font-size:.82em'>⟨${shape.join(",")}⟩</span>`;
}

// Open the catalog entry for a clicked lineage node. The shape is encoded in
// the Mermaid `call` directive ("2x4x3"); we resolve it to the best-known
// scheme for that shape and open its modal, letting the reader walk down the
// lineage. (The @hash is stripped by the generator, so a leaf no longer pins a
// specific base here — best-by-shape is the right navigation target anyway.)
function lineageNodeClick(shapeStr) {
  const shape = String(shapeStr).split("x").map(Number);
  if (shape.length !== 3 || shape.some((x) => !Number.isFinite(x))) return;
  const s = bestSchemeForShape(shape);
  if (s) openEntryModal(s);
}

// Best-known scheme for a shape: lowest-rank entry at that exact orientation,
// else (orientation-agnostic, since tensor rank is permutation-invariant) the
// lowest-rank entry of any axis permutation. Returns the scheme object or null.
function bestSchemeForShape(shape) {
  const exact = (s) => s.format[0] === shape[0] && s.format[1] === shape[1] && s.format[2] === shape[2];
  const sortKey = [...shape].sort((a, b) => a - b).join(",");
  const permOf = (s) => [...s.format].sort((a, b) => a - b).join(",") === sortKey;
  let best = null;
  for (const pass of [exact, permOf]) {
    for (const s of allSchemes || []) {
      if (!s.format || s.format.length !== 3 || s.rank == null) continue;
      if (s.lower === true || s.border === true) continue; // not an achievable scheme
      if (pass(s) && (!best || s.rank < best.rank)) best = s;
    }
    if (best) break;
  }
  return best;
}

// Build a Mermaid `graph TD` definition for a composed lineage, or null when
// the lineage is an atom / formula / cited row (nothing to draw). Each `parse`
// returns { id, rank, shape } so operator nodes can derive their own
// multiplication count from their children (see the metadata note above).
function lineageToMermaid(lc) {
  if (!isComposedLineage(lc)) return null;
  const lines = ["graph TD"];
  let seq = 0;
  const id = () => "L" + seq++;
  // baseLabel is escaped as text; `extra` are raw inline-HTML lines (badges),
  // stacked under it with <br/>. When `shape` is known the node becomes
  // clickable: it opens the best-known catalog entry for that shape, so the
  // reader can walk down into the involved schemes (securityLevel:loose lets
  // the `call` directive reach the global lineageNodeClick).
  const node = (baseLabel, cls, extra, shape) => {
    const me = id();
    const html = [mmEsc(baseLabel), ...(extra || [])].join("<br/>");
    lines.push(`  ${me}["${html}"]:::${cls}`);
    if (shape) lines.push(`  click ${me} call lineageNodeClick("${shape.join("x")}")`);
    return me;
  };
  const edge = (a, b, lbl) => {
    lines.push(lbl ? `  ${a} -->|${mmEsc(lbl)}| ${b}` : `  ${a} --> ${b}`);
  };
  // Badges for an operator node (shows resulting shape + derived rank).
  const opBadges = (shape, rank) => {
    const out = [];
    if (shape) out.push(shapeBadge(shape));
    if (rank != null) out.push(multBadge(rank, false));
    return out;
  };
  function parse(raw) {
    const str = raw.trim();
    // Concat first (lowest precedence) so `a ⊗ b +p c ⊗ d` groups as
    // concat(kron(a,b), kron(c,d)).
    const cc = splitTopLevel(str, [" +p ", " +n ", " +m "]);
    if (cc.parts.length > 1) {
      const kids = cc.parts.map(parse);
      const shape = concatShape(kids.map((k) => k.shape), cc.seps);
      const rank = sumOrNull(kids.map((k) => k.rank));
      const me = node("Concat", "opConcat", opBadges(shape, rank), shape);
      kids.forEach((k, idx) => edge(me, k.id, idx === 0 ? "" : cc.seps[idx - 1]));
      return { id: me, rank, shape };
    }
    // Kronecker: " ⊗ " (space both sides) and " ⊗ˢ" (symmetric variant, operand
    // is a bare ⟨n,m,p⟩ with no trailing space).
    const kk = splitTopLevel(str, [" ⊗ˢ", " ⊗ "]);
    if (kk.parts.length > 1) {
      const kids = kk.parts.map(parse);
      const shape = kronShape(kids.map((k) => k.shape));
      const rank = productOrNull(kids.map((k) => k.rank));
      const me = node("⊗ Kronecker", "opKron", opBadges(shape, rank), shape);
      kids.forEach((k) => edge(me, k.id, ""));
      return { id: me, rank, shape };
    }
    // Projection suffix (row/col drop): "child ↓[keepN|keepM|keepP]". The
    // target shape is the keep-list lengths; the rank is the child's (a
    // projection only zeroes rows/cols, it never adds products).
    const pj = str.match(/^([\s\S]*\S)\s*↓\[([^\]]*)\]$/);
    if (pj) {
      const kid = parse(pj[1]);
      const keeps = pj[2].split("|").map((x) => x.split(",").filter(Boolean).length);
      const shape = keeps.length === 3 && keeps.every((k) => k > 0) ? keeps : null;
      const me = node("Project ↓", "opComposed", shape ? [shapeBadge(shape)] : [], shape);
      edge(me, kid.id, "");
      return { id: me, rank: kid.rank, shape };
    }
    // Re-orientation suffix: "child →⟨n,m,p⟩" — same tensor, permuted axes.
    const or = str.match(/^([\s\S]*\S)\s*→⟨(\d+),(\d+),(\d+)⟩$/);
    if (or) {
      const kid = parse(or[1]);
      const shape = [+or[2], +or[3], +or[4]];
      const me = node("Orient →", "opComposed", [shapeBadge(shape)], shape);
      edge(me, kid.id, "");
      return { id: me, rank: kid.rank, shape };
    }
    if (str.startsWith("R[") && str.endsWith("]")) {
      const inner = str.slice(2, -1);
      const semi = inner.indexOf(";");
      const base = (semi >= 0 ? inner.slice(0, semi) : inner).trim();
      const alloc = semi >= 0 ? inner.slice(semi + 1).trim() : "";
      const baseShape = shapeOfRef(base);
      // Recombination keeps the shape but reduces the count by an amount not in
      // the string → show the shape, assert no rank for the recombine node.
      const me = node("Recombine" + (alloc ? "\n" + alloc : ""), "opRecomb",
        baseShape ? [shapeBadge(baseShape)] : [], baseShape);
      edge(me, parse(base).id, "base");
      return { id: me, rank: null, shape: baseShape };
    }
    if ((str.startsWith("Derived[") || str.startsWith("Composed[")) && str.endsWith("]")) {
      const inner = str.slice(str.indexOf("[") + 1, -1);
      const semi = inner.indexOf(";");
      const mech = (semi >= 0 ? inner.slice(0, semi) : inner).trim();
      const child = semi >= 0 ? inner.slice(semi + 1).trim() : "";
      const childRes = child ? parse(child) : null;
      const me = node("Derived\n" + mech, "opComposed",
        childRes && childRes.shape ? [shapeBadge(childRes.shape)] : [],
        childRes ? childRes.shape : null);
      if (childRes) edge(me, childRes.id, "");
      return { id: me, rank: null, shape: childRes ? childRes.shape : null };
    }
    // Leaf: explicit ref rank wins; else best-known for the shape (marked ≈).
    // (The lineage_compact is already hash-free — stripped by the generator.)
    const shape = shapeOfRef(str);
    const explicit = rankOfRef(str);
    const best = explicit == null ? bestRankForShape(shape) : null;
    const rank = explicit != null ? explicit : best ? best.rank : null;
    const approx = explicit == null && !!best;
    const me = node(str, "leaf", rank != null ? [multBadge(rank, approx)] : [], shape);
    return { id: me, rank, shape };
  }
  parse(lc);
  lines.push("  classDef opKron fill:#e3f2ff,stroke:#4a90d9,color:#123;");
  lines.push("  classDef opConcat fill:#fff0e0,stroke:#d98a4a,color:#321;");
  lines.push("  classDef opRecomb fill:#efe6ff,stroke:#8a6ad9,color:#212;");
  lines.push("  classDef opComposed fill:#e6ffe6,stroke:#5aaa5a,color:#121;");
  lines.push("  classDef leaf fill:#f6f6f6,stroke:#bbb,color:#333;");
  return lines.join("\n");
}

// ─── Structured lineage → Mermaid (exact: real root rank + hash-resolved leaves) ──
//
// Unlike lineageToMermaid (which parses the now-hash-free human string and can
// only ESTIMATE ranks), this walks the structured `lineage` DAG fetched from the
// scheme's own JSON. The ROOT node shows the scheme's ACTUAL stored rank; Atom
// leaves are resolved by `shape@hash` to the EXACT sub-scheme (≈ best-known only
// when the pinned base isn't in the catalog). Internal composition nodes carry no
// rank badge — their rank is NOT a simple function of the children (serendipitous
// / recombination land below the naive product).

// "{shape}@{hash7}" → rank, lazily built from catalog file paths (the last
// filename token is the content hash7).
let _schemeHashIndex = null;
function schemeRankByHash(shape, fullHash) {
  if (_schemeHashIndex == null) {
    _schemeHashIndex = new Map();
    for (const e of allSchemes) {
      if (!e.file || e.rank == null || !e.format) continue;
      const m = e.file.match(/-([0-9a-f]{4,})\.json$/);
      if (m) _schemeHashIndex.set(e.format.join("x") + "@" + m[1].slice(0, 7), e.rank);
    }
  }
  return _schemeHashIndex.get(shape + "@" + String(fullHash).slice(0, 7)) ?? null;
}

// Children + display metadata for any structured lineage op.
function lineageChildren(node) {
  switch (node.op) {
    case "KronProduct":  return { label: "⊗ Kronecker", cls: "opKron", kids: [{ n: node.outer }, { n: node.inner }] };
    case "KronChain":    return { label: "⊗ Kron chain", cls: "opKron", kids: (node.factors || []).map((f) => ({ n: f })) };
    case "ConcatCols":   return { label: "Concat ∥p", cls: "opConcat", kids: [{ n: node.left }, { n: node.right }] };
    case "ConcatRows":   return { label: "Concat ∥n", cls: "opConcat", kids: [{ n: node.top }, { n: node.bottom }] };
    case "ConcatInner":  return { label: "Concat ∥m", cls: "opConcat", kids: [{ n: node.left }, { n: node.right }] };
    case "SumInner":     return { label: "SumInner", cls: "opConcat", kids: [{ n: node.left }, { n: node.right }] };
    case "Recombination":
    case "RecombinationN": {
      const alloc = `${(node.allocA || []).join(",")} | ${(node.allocB || []).join(",")} | ${(node.allocC || []).join(",")}`;
      const kids = [{ n: node.base, lbl: "base" }].concat((node.leaves || []).map((l) => ({ n: l, lbl: "leaf" })));
      return { label: "Recombine\n" + alloc, cls: "opRecomb", kids };
    }
    case "SerendipitousProduct":
      return { label: "⊗ˢ Serendipitous", cls: "opKron",
        kids: [{ n: node.base }, { n: { op: "Atom", ref: `⟨${node.n2},${node.m2},${node.p2}⟩` } }] };
    case "AxisFlip":   return { label: "AxisFlip" + (node.mask != null ? " m" + node.mask : ""), cls: "opComposed", kids: [{ n: node.child }] };
    case "OrientAs":   return { label: "OrientAs", cls: "opComposed", kids: [{ n: node.child }] };
    case "Transpose":  return { label: "Transpose " + (node.perm || ""), cls: "opComposed", kids: [{ n: node.child }] };
    case "Project":    return { label: "Project", cls: "opComposed", kids: [{ n: node.child }] };
    default:           return { label: node.op || "?", cls: "opComposed",
        kids: ["child", "base", "left", "right", "outer", "inner", "top", "bottom"].filter((k) => node[k]).map((k) => ({ n: node[k] })) };
  }
}

// Best-effort shape of a descendant node (the root shape is passed in explicitly).
function lineageNodeShape(node) {
  if (!node) return null;
  switch (node.op) {
    case "Atom":      return shapeOfRef(node.ref || "");
    case "OrientAs":  return [node.n, node.m, node.p];
    case "SerendipitousProduct": {
      const b = lineageNodeShape(node.base);
      return b ? [b[0] * node.n2, b[1] * node.m2, b[2] * node.p2] : null;
    }
    case "KronProduct": return kronShape([lineageNodeShape(node.outer), lineageNodeShape(node.inner)]);
    case "KronChain":   return kronShape((node.factors || []).map(lineageNodeShape));
    case "AxisFlip":
    case "Transpose":
    case "Project":     return lineageNodeShape(node.child);
    default:            return null;
  }
}

function structuredLineageToMermaid(root, rootFormat, rootRank) {
  if (!root || !root.op) return null;
  const lines = ["graph TD"];
  let seq = 0;
  const id = () => "L" + seq++;
  const mk = (label, cls, badges, shape) => {
    const me = id();
    const html = [mmEsc(label), ...(badges || [])].join("<br/>");
    lines.push(`  ${me}["${html}"]:::${cls}`);
    if (shape) lines.push(`  click ${me} call lineageNodeClick("${shape.join("x")}")`);
    return me;
  };
  const edge = (a, b, lbl) => lines.push(lbl ? `  ${a} -->|${mmEsc(lbl)}| ${b}` : `  ${a} --> ${b}`);

  function walk(node, isRoot, forcedShape, forcedRank) {
    if (node.op === "Atom") {
      const ref = node.ref || "";
      const at = ref.indexOf("@");
      const shape = forcedShape || shapeOfRef(ref);
      let rank = null, approx = false, danglingPin = false;
      if (ref === "naive-1x1x1") rank = 1;
      else if (at >= 0 && shape) {
        rank = schemeRankByHash(shape.join("x"), ref.slice(at + 1));
        if (rank == null) danglingPin = true; // pinned base has no catalog entry
      }
      if (rank == null && shape) { const b = bestRankForShape(shape); if (b) { rank = b.rank; approx = true; } }
      // Leaf label already IS the shape (the ref), so no shapeBadge — just the rank.
      const badges = [];
      if (isRoot && forcedRank != null) badges.push(multBadge(forcedRank, false));
      else if (rank != null) badges.push(multBadge(rank, approx));
      if (danglingPin) {
        badges.push(`<span title="the pinned base @${ref.slice(at + 1, at + 8)} is not in the catalog — this is the best-known rank for the shape, not the exact base used" style="color:#b54708;font-size:.72em">⚠ base not catalogued</span>`);
      }
      return mk(at >= 0 ? ref.slice(0, at) : ref, danglingPin ? "leafDangling" : "leaf", badges, shape);
    }
    const { label, cls, kids } = lineageChildren(node);
    const shape = forcedShape || lineageNodeShape(node);
    const badges = [];
    if (shape) badges.push(shapeBadge(shape));
    // Only the root carries a rank badge — the real, stored rank. Internal nodes
    // get none (their rank is not the naive product/sum of the children).
    if (isRoot && forcedRank != null) badges.push(multBadge(forcedRank, false));
    const me = mk(label, cls, badges, shape);
    for (const k of kids) {
      if (k.n) edge(me, walk(k.n, false, null, null), k.lbl || "");
    }
    return me;
  }
  walk(root, true, rootFormat, rootRank);
  lines.push("  classDef opKron fill:#e3f2ff,stroke:#4a90d9,color:#123;");
  lines.push("  classDef opConcat fill:#fff0e0,stroke:#d98a4a,color:#321;");
  lines.push("  classDef opRecomb fill:#efe6ff,stroke:#8a6ad9,color:#212;");
  lines.push("  classDef opComposed fill:#e6ffe6,stroke:#5aaa5a,color:#121;");
  lines.push("  classDef leaf fill:#f6f6f6,stroke:#bbb,color:#333;");
  lines.push("  classDef leafDangling fill:#fff4e5,stroke:#f0a020,stroke-dasharray:4 3,color:#7a4a10;");
  return lines.join("\n");
}

// Render the lineage graph into the modal's box: prefer the structured DAG from
// the scheme's own JSON (exact root rank + hash-resolved leaves), fall back to
// the hash-free compact string when the fetch fails or the file has no lineage.
async function renderLineageGraph(box, s) {
  const target = box.querySelector("#lineage-graph-box");
  if (!target || !window.mermaid) return;
  let graphDef = null;
  if (s.file) {
    try {
      const resp = await fetch(SCHEMES_RAW_BASE + s.file);
      if (resp.ok) {
        const js = await resp.json();
        if (js && js.lineage) graphDef = structuredLineageToMermaid(js.lineage, s.format, s.rank);
      }
    } catch (e) { /* fall through to the compact-string renderer */ }
  }
  if (!graphDef) graphDef = lineageToMermaid(s.lineage_compact);
  if (!graphDef) { target.textContent = ""; return; }
  try {
    const { svg } = await mermaid.render("lineage-svg-" + lineageGraphSeq++, graphDef);
    target.innerHTML = svg;
  } catch (err) {
    target.textContent = "(could not render lineage graph: " + err.message + ")";
  }
}

// Monotonic id for the throwaway element mermaid.render() needs (must be a
// valid, unique DOM id per call).
let lineageGraphSeq = 0;

/**
 * Derive the narrowest known field tag for a scheme by inspecting
 * `file_source_label` and `file` (which encode `_Q_`, `_Z_`, `_ZT_`,
 * `_F2_`, `_C_`, `_R_` suffixes per the project convention) and
 * falling back to the coarse `field` attribute.
 *
 * Catalog `field` is "R/Q/Z" for most characteristic-0 entries — too
 * coarse for "Q only" or "Z only" filtering. We extract the narrower
 * tag from the source string where available.
 *
 * Returns one of "Z", "Q", "ZT", "R", "C", "F2", or "R/Q/Z" (when
 * nothing narrower is known) or the raw `field` for unusual values
 * like "derived" / "mixed".
 */
function narrowField(scheme) {
  const tag = (scheme.field || "").trim();
  // F2 / C / mixed / derived → trust the coarse tag, no narrowing from source.
  if (tag === "F2" || tag === "C" || tag === "mixed" || tag === "derived") return tag;
  const src = String(scheme.file_source_label || "") + "|" + String(scheme.file || "") + "|" + String(scheme.source || "");
  // Order matters: ZT and F2 share a 2-char prefix with single-letter tags;
  // check the longer variants first so "ZT" doesn't get mis-classified as "Z".
  // ZT is the ternary ({-1,0,1}) sub-class of Z — for the coarse FIELD axis it
  // IS just Z (an integer scheme), so it narrows to "Z". The ternary-ness is a
  // separate boolean (`scheme.ZT`) handled by the dedicated "ZT" filter value
  // in schemeValidForRequestedField — it is NOT, and never was, an F₂ thing.
  if (/[-_]ZT(?:[_/]|$)/.test(src)) return "Z";
  if (/[-_]F2(?:[_/]|$)/.test(src)) return "F2";
  if (/[-_]Q(?:[_/]|$)/.test(src)) return "Q";
  if (/[-_]Z(?:[_/]|$)/.test(src)) return "Z";
  if (/[-_]C(?:[_/]|$)/.test(src)) return "C";
  if (/[-_]R(?:[_/]|$)/.test(src)) return "R";
  return tag || "R/Q/Z";
}

/**
 * Field-inclusion rule (task #94 — strict equivalence classes per user spec):
 *
 * Each narrowed tag (Q, Z, ZT, R, C, F2) is treated as its own equivalence
 * class — we do NOT automatically promote Q-tagged schemes under an R-only
 * filter, even though Q ⊂ R mathematically. The reason: in this catalog "Q"
 * often means "the construction uses rationals in a way that exploits Q's
 * special properties" (small denominators, Gauss-like cancellation) and a
 * user looking for R-only schemes does not want to be flooded with Q-only
 * results that may not generalise the way they expect.
 *
 * The historic exception is the legacy "R/Q/Z" coarse bucket: a scheme
 * tagged R/Q/Z (no narrower source-side hint, e.g. Strassen, Laderman,
 * AlphaTensor-Z imported under that label) surfaces under any of R, Q,
 * or Z filters — these are the truly field-agnostic algorithms. R/Q/Z
 * also lifts to C (since R ⊂ C); this preserves the existing C-fallback
 * badge behaviour.
 *
 * F₂ is in its own characteristic-2 universe — no cross-acceptance either way.
 */
function schemeValidForRequestedField(schemeField, requestedField, scheme) {
  if (!requestedField) return true; // "" = all
  // ZT is NOT a field — it is the sub-class of Z (integer) schemes whose
  // coefficients are all in {-1,0,1}. Each scheme JSON carries a stamped `zt`
  // boolean (present only when Z ∈ fields), surfaced verbatim into catalog.json.
  // "ZT" as a filter value means "show Z schemes whose zt flag is true".
  // (Historically ZT was wrongly bucketed with F₂/Z₂ — it has nothing to do
  // with characteristic 2; it is ternary INTEGER, not ternary modular.)
  if (requestedField === "ZT") {
    return scheme ? scheme.zt === true : false;
  }
  // Authoritative path: explicit per-scheme fields[] (task #175/#181) — the
  // verified set of fields this scheme is valid over. Match by membership and
  // skip the legacy cluster/lift heuristics entirely.
  if (scheme && Array.isArray(scheme.fields) && scheme.fields.length) {
    // "R/Q/Z" is a SET filter (not a literal tag, which no longer exists after
    // the field→fields refactor): match any characteristic-0 scheme valid over
    // at least one of Z, Q, or R.
    if (requestedField === "R/Q/Z") {
      return scheme.fields.some(f => f === "Z" || f === "Q" || f === "R");
    }
    return scheme.fields.includes(requestedField);
  }
  // Legacy fallback (stub / non-bilinear / pre-migration entries without fields[]).
  // Prefer the narrowed-from-filename tag when we have the scheme; otherwise
  // fall back to the raw schemeField string (some callers only have the raw field).
  const narrow = scheme ? narrowField(scheme) : schemeField;

  if (narrow === requestedField) return true;

  // Legacy "R/Q/Z" bucket: surface under any of R, Q, Z requests; also lifts to C.
  if (narrow === "R/Q/Z") {
    return requestedField === "R" || requestedField === "Q" || requestedField === "Z"
        || requestedField === "C";
  }

  // R-class lifts to C (R ⊂ C). This is the only cross-class lift we keep —
  // R-valid schemes are unconditionally valid over C.
  if (narrow === "R" && requestedField === "C") return true;
  // Z and Q narrow-tagged schemes do NOT auto-promote to R (treat as own
  // equivalence class — user spec). They DO lift to C as integer/rational
  // data is trivially valid over C.
  if ((narrow === "Z" || narrow === "Q") && requestedField === "C") return true;

  // Everything else (F2 / ZT / C / mismatch) → no cross-acceptance.
  return false;
}

/**
 * Schönhage-style omega bound implied by a rank R for ⟨n,m,p⟩ matmul:
 * geometric-mean size = (n·m·p)^(1/3), ω = 3·ln(R) / ln(n·m·p).
 * For cubic ⟨n,n,n⟩=R this reduces to log_n(R) (e.g., Strassen 7 → ≈2.807).
 * Returns null if undefined (n=m=p=1 or rank 1).
 */
function impliedOmega(format, rank) {
  const product = format[0] * format[1] * format[2];
  if (product <= 1 || rank <= 1) return null;
  return 3 * Math.log(rank) / Math.log(product);
}

/**
 * Linear-interpolation heatmap color from low (green, hsl 120°) to high
 * (red, hsl 0°). Returns a CSS background string.
 */
function omegaHeatColor(omega, minOmega, maxOmega) {
  if (omega == null || minOmega == null || maxOmega == null || minOmega === maxOmega) {
    return "transparent";
  }
  const t = (omega - minOmega) / (maxOmega - minOmega);
  const hue = 120 * (1 - t); // 120 = green (low ω, good); 0 = red (high ω, bad)
  return `hsl(${hue}, 70%, 88%)`;
}

function escapeHtml(s) {
  return String(s).replace(/[&<>"']/g, c =>
    ({ "&": "&amp;", "<": "&lt;", ">": "&gt;", "\"": "&quot;", "'": "&#39;" }[c]));
}

/**
 * Link to the scheme's file in the originating author's OWN repository
 * (`source_scheme_url`, e.g. Perminov's FastMatrixMultiplication) — distinct
 * from `source_paper_url` (the publication) and from our own JSON copy. Empty
 * string when the scheme has no such pointer.
 */
function sourceSchemeLink(s) {
  if (!s || !s.source_scheme_url) return "";
  return ` <a href="${escapeHtml(s.source_scheme_url)}" target="_blank" rel="noopener"`
    + ` title="Scheme file in the source author's own repository">src ↗</a>`;
}

loadCatalog().catch(err => {
  document.getElementById("result-count").textContent =
    "ERROR loading catalog.json: " + err.message;
  console.error(err);
});

// ── ω-history section ───────────────────────────────────────────────────────
// Renders docs/omega-history.json (NC ω) + docs/omega-history-commutative.json
// (commutative ω_c, single-shape upper bounds) as a combined table + SVG
// timeline. The two series are visually distinct (blue = NC ω with a
// connecting path, green = commutative ω_c as scatter — NO path, because
// commutative single-shape bounds don't aggregate into an asymptotic
// sequence the same way NC ω does).
async function loadOmegaHistory() {
  const [resp, respC] = await Promise.all([
    fetch("omega-history.json").catch(() => null),
    fetch("omega-history-commutative.json").catch(() => null),
  ]);
  if (!resp || !resp.ok) return;
  const data = await resp.json();
  const ncEntries = data.entries || [];
  const cEntries = (respC && respC.ok) ? (await respC.json()).entries || [] : [];
  if (!ncEntries.length && !cEntries.length) return;

  const tbody = document.getElementById("omega-history-body");
  if (tbody) {
    const ncRows = ncEntries.map(e => ({
      kind: "NC ω", year: e.year, omega: e.omega_upper,
      method: e.method, source: e.source, notes: e.notes,
    }));
    const cRows = cEntries.map(e => ({
      kind: "commutative ω<sub>c</sub>", year: e.year, omega: e.omega_c_upper,
      method: e.method, source: e.source,
      notes: (e.notes || "") + (e.scheme_file ? ` (${e.scheme_file})` : ""),
    }));
    const all = ncRows.concat(cRows).sort((a, b) => a.year - b.year || a.kind.localeCompare(b.kind));
    tbody.innerHTML = all.map(r => `
      <tr>
        <td>${r.kind}</td>
        <td>${r.year}</td>
        <td><code>${r.omega}</code></td>
        <td>${escapeHtml(r.method || "")}</td>
        <td>${escapeHtml(r.source || "")}</td>
        <td class="muted" style="font-size:0.9em">${escapeHtml(r.notes || "")}</td>
      </tr>
    `).join("");
  }

  const svg = document.getElementById("omega-timeline");
  const renderWithMode = () => {
    if (svg) {
      const logScale = document.getElementById("omega-log-scale")?.checked ?? false;
      renderOmegaTimeline(svg, ncEntries, cEntries, data.lower_bound, logScale);
    }
  };
  document.getElementById("omega-log-scale")?.addEventListener("change", renderWithMode);
  renderWithMode();
}

function renderOmegaTimeline(svg, ncEntries, cEntries, lowerBound, logScale = false) {
  const W = 800, H = 240, padL = 60, padR = 20, padT = 20, padB = 50;
  const allYears = ncEntries.map(e => e.year).concat(cEntries.map(e => e.year));
  const allOmegas = ncEntries.map(e => e.omega_upper).concat(cEntries.map(e => e.omega_c_upper));
  // Floor the visible axis at 2.35 (task #153): every plotted ω is ≥ 2.37, so the
  // 2.0–2.35 band is dead space. The true information-theoretic bound ω ≥ 2 is
  // noted in the floor label rather than drawn (it would waste a third of the axis).
  const AXIS_FLOOR = 2.35;
  const yMin = Math.max(AXIS_FLOOR, lowerBound != null ? lowerBound : AXIS_FLOOR);
  const yMax = Math.max(...allOmegas) + 0.05;
  const xMin = Math.min(...allYears) - 2;
  const xMax = Math.max(...allYears) + 2;
  const xScale = y => padL + (y - xMin) / (xMax - xMin) * (W - padL - padR);
  // y-axis scaling — linear by default, optional log on (ω − yMin + ε)
  // so the dense band near ω ≈ 2.37 is expanded. Both modes pin yMin to
  // the bottom of the plot.
  const EPS = 0.01;
  const linearY = w => padT + (yMax - w) / (yMax - yMin) * (H - padT - padB);
  const logY = w => {
    const lo = Math.log(EPS);
    const hi = Math.log(yMax - yMin + EPS);
    const lw = Math.log(Math.max(EPS, w - yMin + EPS));
    return padT + (hi - lw) / (hi - lo) * (H - padT - padB);
  };
  const yScale = logScale ? logY : linearY;

  // Background gridlines.
  let svgBody = "";
  for (let w = Math.ceil(yMin * 10) / 10; w <= yMax; w += 0.1) {
    const y = yScale(w);
    svgBody += `<line x1="${padL}" x2="${W - padR}" y1="${y}" y2="${y}" stroke="#eee"/>`;
    if (Math.abs(w - Math.round(w * 10) / 10) < 0.001 && Math.round(w * 10) % 1 === 0) {
      svgBody += `<text x="${padL - 5}" y="${y + 4}" text-anchor="end" font-size="11" fill="#666">${w.toFixed(1)}</text>`;
    }
  }
  // Axis-floor marker — red dashed horizontal at the 2.35 floor; the genuine
  // information-theoretic bound is ω ≥ 2 (open), noted in the label.
  const yLB = yScale(yMin);
  svgBody += `<line x1="${padL}" x2="${W - padR}" y1="${yLB}" y2="${yLB}" stroke="#c00" stroke-dasharray="5,3"/>`;
  svgBody += `<text x="${W - padR}" y="${yLB - 4}" text-anchor="end" font-size="11" fill="#c00">axis floor ${yMin.toFixed(2)} · true bound ω ≥ 2 (open)</text>`;

  // Year labels (x axis).
  for (const yr of [1969, 1978, 1987, 2010, 2024]) {
    if (yr < xMin || yr > xMax) continue;
    const x = xScale(yr);
    svgBody += `<text x="${x}" y="${H - padB + 16}" text-anchor="middle" font-size="11" fill="#666">${yr}</text>`;
    svgBody += `<line x1="${x}" x2="${x}" y1="${H - padB}" y2="${H - padB + 4}" stroke="#aaa"/>`;
  }
  // x-axis.
  svgBody += `<line x1="${padL}" x2="${W - padR}" y1="${H - padB}" y2="${H - padB}" stroke="#333"/>`;

  // NC ω: line + points (blue, connected — proper asymptotic sequence).
  if (ncEntries.length) {
    const pathD = ncEntries.map((e, i) =>
      `${i === 0 ? "M" : "L"}${xScale(e.year)},${yScale(e.omega_upper)}`).join(" ");
    svgBody += `<path d="${pathD}" fill="none" stroke="#0077cc" stroke-width="2"/>`;
    for (const e of ncEntries) {
      const cx = xScale(e.year), cy = yScale(e.omega_upper);
      svgBody += `<circle cx="${cx}" cy="${cy}" r="4" fill="#0077cc" stroke="white" stroke-width="1.5">
        <title>NC ω · ${escapeHtml(e.source)} (${e.year}): ω ≤ ${e.omega_upper}${e.method ? "\n" + e.method : ""}</title>
      </circle>`;
    }
  }

  // Commutative ω_c: points only (green, NO connecting path — these are
  // single-shape upper bounds, not an aggregated asymptotic sequence).
  if (cEntries.length) {
    for (const e of cEntries) {
      const cx = xScale(e.year), cy = yScale(e.omega_c_upper);
      svgBody += `<circle cx="${cx}" cy="${cy}" r="4" fill="#2a9d3a" stroke="white" stroke-width="1.5">
        <title>commutative ω_c · ${escapeHtml(e.source)} (${e.year}): ω_c ≤ ${e.omega_c_upper} via ⟨${e.n},${e.n},${e.n}⟩=${e.rank}${e.method ? "\n" + e.method : ""}</title>
      </circle>`;
    }
  }

  // Legend (top-right corner).
  const lx = W - padR - 200, ly = padT + 10;
  svgBody += `<rect x="${lx - 8}" y="${ly - 12}" width="200" height="48" fill="white" stroke="#ddd"/>`;
  svgBody += `<circle cx="${lx + 4}" cy="${ly + 4}" r="4" fill="#0077cc"/>`;
  svgBody += `<text x="${lx + 14}" y="${ly + 8}" font-size="11" fill="#333">NC ω (asymptotic, recursive)</text>`;
  svgBody += `<circle cx="${lx + 4}" cy="${ly + 22}" r="4" fill="#2a9d3a"/>`;
  svgBody += `<text x="${lx + 14}" y="${ly + 26}" font-size="11" fill="#333">commutative ω_c (single-shape)</text>`;

  // Axis labels.
  svgBody += `<text x="${padL / 2}" y="${padT - 6}" text-anchor="middle" font-size="11" fill="#333" font-weight="bold">ω ≤</text>`;
  svgBody += `<text x="${(W + padL) / 2}" y="${H - 6}" text-anchor="middle" font-size="11" fill="#333">year</text>`;

  svg.innerHTML = svgBody;
}

loadOmegaHistory().catch(err => console.error("ω-history load failed:", err));
