package eu.solven.matmul.docs.generate;

import lombok.extern.slf4j.Slf4j;

import java.io.BufferedWriter;
import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.io.PrintWriter;
import java.nio.file.Files;
import java.nio.file.Path;
import java.time.Instant;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.HashMap;
import java.util.List;
import java.util.Map;
import java.util.regex.Matcher;
import java.util.regex.Pattern;


import eu.solven.matmul.ComplexNonCubicBilinearAlgorithm;
import eu.solven.matmul.NonCubicBilinearAlgorithm;
import eu.solven.matmul.Verifier;
import eu.solven.matmul.catalog.FieldAwareLookup;
import eu.solven.matmul.catalog.HumanScheme;
import eu.solven.matmul.catalog.SchemeIO;
import tools.jackson.databind.JsonNode;
import tools.jackson.databind.json.JsonMapper;
import tools.jackson.databind.node.ArrayNode;
import tools.jackson.databind.node.ObjectNode;

/**
 * Generates {@code docs/catalog.json} — a flat manifest of every verified
 * scheme on disk, used by the GitHub Pages browser UI ({@code docs/index.html}).
 *
 * <p>Output schema:</p>
 * <pre>
 * {
 *   "generated": "ISO-8601 timestamp",
 *   "schemes": [
 *     { "format": [n, m, p], "max_dim": N, "fields": ["F2","Z","Q","R","C"],
 *       "rank": R, "additions": A|null, "source": "...", "verified": true/false,
 *       "complex": false, "file": "section{N}/{file}.json" },
 *     ...
 *   ]
 * }
 * </pre>
 */
@Slf4j
public class GenerateCatalogManifest {

	private static final File SCHEMES_DIR = new File("src/main/resources/schemes");
	private static final File OUTPUT = new File("docs/catalog.json");
	/** Characteristic-0 field tags: a scheme valid over any of these is comparable
	 *  to the (NC, char-0) FMM/Perminov reference catalogs. */
	private static final java.util.Set<String> CHAR0_FIELD_TAGS =
			java.util.Set.of("Z", "Q", "R", "C");
	/** Caps for emitting the human-friendly multiplications/elements (#188) into
	 *  catalog.json — small schemes only, where a textual form is readable and
	 *  the size cost is bounded. Larger schemes get the form on demand later. */
	private static final int HUMAN_MAX_RANK = 80;
	private static final int HUMAN_MAX_DIM = 12;
	public static void main(String[] args) throws IOException {
		List<ObjectNode> schemeNodes = new ArrayList<>();
		// Tally, across every scheme's lineage, how many times each base is USED as a
		// building block (by recombination / projection / kronecker / …). Stamped onto
		// each entry as `used_as_base` in a post-pass once all entries exist.
		eu.solven.matmul.catalog.BaseUsageStats baseUsage = new eu.solven.matmul.catalog.BaseUsageStats();
		// Every corrupted basename, captured BEFORE shape-dedup (so same-shape
		// siblings the dedup drops are still excluded by FieldAwareLookup's gating).
		java.util.Set<String> corruptedFiles = new java.util.TreeSet<>();
		JsonMapper mapper = JsonMapper.builder().build();

		// Optional fmm-lille catalog cross-reference: maps (n,m,p) sorted →
		// {best_rank, references}. If the JSON file doesn't exist, the
		// cross-ref columns are simply absent from the manifest.
		Map<String, JsonNode> fmmByCanonical = loadFmmCatalog(mapper);
		// Merged FMM ∪ Perminov best NC rank per format → drives solven_discovery.
		Map<String, ExternalBest> externalBest = loadExternalBest(mapper);

		File[] files;
		try (var s = Files.walk(SCHEMES_DIR.toPath())) {
			files = s.filter(p -> p.toString().endsWith(".json"))
					.map(Path::toFile)
					.sorted()
					.toArray(File[]::new);
		}
		Arrays.sort(files);

		// explicitable per file — computed ONCE over the whole dependency DAG (a precise
		// EXPLICIT scheme vs a CITED BOUND), then emitted per entry below.
		Map<Path, Boolean> explicitableByFile = eu.solven.matmul.docs.verify.ComputeExplicitable
				.computeAll(Arrays.stream(files).map(File::toPath).toList());

		// Pinned-ref resolution set: every file's "{shape}@{hash7}" (canonical
		// filenames encode both shape and content hash). A lineage Atom ref of form
		// "{shape}@{fullhash}" with no matching key here is a DANGLING base — the
		// sub-scheme it pins was never persisted (e.g. an un-saved serendipitous /
		// re-oriented base). We warn on these so the SPA's ≈ best-known fallback in
		// the lineage graph isn't a silent surprise, and so they can be fixed.
		java.util.regex.Pattern pinnedFn =
				java.util.regex.Pattern.compile("^(\\d+x\\d+x\\d+)-r\\d+-.+-([0-9a-f]{4,})$");
		java.util.Set<String> existingPinnedKeys = new java.util.HashSet<>();
		for (File af : files) {
			String stem = af.getName().replaceFirst("\\.json$", "");
			java.util.regex.Matcher m = pinnedFn.matcher(stem);
			if (m.matches()) {
				existingPinnedKeys.add(m.group(1) + "@" + m.group(2).substring(0, Math.min(7, m.group(2).length())));
			}
		}
		java.util.Set<String> danglingBaseFiles = new java.util.TreeSet<>();

		// Atom resolver for lineage-level bud inference: maps a lineage Atom ref
		// (canonical {NxMxP}_m{R}_a{A} key) to its small on-disk scheme, so a
		// composite's buds can be inferred WITHOUT expanding the composite — only
		// its atom leaves are read (and cached). See LineageBudInference.
		Map<String, File> atomFileByKey = new java.util.HashMap<>();
		for (File af : files) {
			String stem = af.getName().endsWith(".json")
					? af.getName().substring(0, af.getName().length() - 5) : af.getName();
			atomFileByKey.putIfAbsent(eu.solven.matmul.catalog.Lineage.canonicalKey(stem), af);
		}
		Map<String, NonCubicBilinearAlgorithm> atomCache = new java.util.HashMap<>();
		java.util.function.Function<String, NonCubicBilinearAlgorithm> atomResolver = ref -> {
			String key = eu.solven.matmul.catalog.Lineage.canonicalKey(ref);
			if (atomCache.containsKey(key)) return atomCache.get(key);
			NonCubicBilinearAlgorithm alg = null;
			File af = atomFileByKey.get(key);
			if (af != null) {
				try { alg = SchemeIO.read(af); } catch (Exception ignored) { alg = null; }
			}
			atomCache.put(key, alg);
			return alg;
		};

		// Leaf-by-shape resolver for the recombination bud estimate: maps a
		// sub-shape (n,m,p) to the bud profile of the smallest-rank catalog scheme
		// there, expanded once (cached). Capped to bounded dims so the manifest
		// never expands a huge leaf; over the cap → null → the recombination stays
		// UNKNOWN (honest). bud_score is invariant under leaf axis-orientation, so
		// the Pareto metric is robust even though a per-type split may be permuted.
		Map<String, File> schemeByShape = new java.util.HashMap<>();
		Map<String, Integer> rankByShape = new java.util.HashMap<>();
		for (File af : files) {
			int[] sh;
			int rk;
			try {
				tools.jackson.databind.JsonNode r = SchemeIO.parseJson(af);
				tools.jackson.databind.JsonNode nN = r.get("n");
				if (nN == null || !nN.isArray() || nN.size() != 3) continue;
				sh = new int[] { nN.get(0).asInt(), nN.get(1).asInt(), nN.get(2).asInt() };
				rk = r.has("m") ? r.get("m").asInt() : (r.has("rank") ? r.get("rank").asInt() : -1);
				if (rk < 0) continue;
			} catch (Exception e) {
				continue;
			}
			Arrays.sort(sh);
			String key = sh[0] + "x" + sh[1] + "x" + sh[2];
			Integer cur = rankByShape.get(key);
			// Pick the min-rank scheme per shape; break ties by filename so the
			// chosen leaf-bud representative is DETERMINISTIC (independent of
			// File.listFiles() order — otherwise renaming files perturbs which
			// same-rank sibling wins, shifting recombination bud estimates and the
			// shave count between regenerations → spurious manifest churn).
			boolean better = cur == null || rk < cur
					|| (rk == cur && af.getName().compareTo(schemeByShape.get(key).getName()) < 0);
			if (better) {
				rankByShape.put(key, rk);
				schemeByShape.put(key, af);
			}
		}
		Map<String, eu.solven.matmul.catalog.LineageBudInference.Profile> leafProfileCache =
				new java.util.HashMap<>();
		java.util.function.Function<int[], eu.solven.matmul.catalog.LineageBudInference.Profile>
				leafShapeResolver = sz -> {
			int[] sh = sz.clone();
			Arrays.sort(sh);
			if (sh[2] > 16) return null;  // bound expansion cost → recombination UNKNOWN
			String key = sh[0] + "x" + sh[1] + "x" + sh[2];
			if (leafProfileCache.containsKey(key)) return leafProfileCache.get(key);
			eu.solven.matmul.catalog.LineageBudInference.Profile prof = null;
			File af = schemeByShape.get(key);
			if (af != null) {
				try {
					prof = eu.solven.matmul.catalog.LineageBudInference.fromExpanded(SchemeIO.read(af));
				} catch (Exception ignored) {
					prof = null;
				}
			}
			leafProfileCache.put(key, prof);
			return prof;
		};

		for (File f : files) {
			// Skip files that don't verify standalone (KNOWN_BROKEN_FILES) — publishing
			// them would be a correctness bug; they remain on disk for investigation.
			if (FieldAwareLookup.isKnownBroken(f.getName())) continue;
			// CONTENT-DRIVEN: shape (n), rank (m), source — read from JSON, NOT the
			// filename. Filenames are pure labels now (post StampFields/Source/Additions).
			tools.jackson.databind.JsonNode root;
			try {
				root = SchemeIO.parseJson(f);
			} catch (IOException e) {
				continue;
			}
			tools.jackson.databind.JsonNode nNode = root.get("n");
			if (nNode == null || !nNode.isArray() || nNode.size() != 3) continue;
			int n = nNode.get(0).asInt();
			int mm = nNode.get(1).asInt();
			int p = nNode.get(2).asInt();
			int rank = root.has("m") ? root.get("m").asInt() : (root.has("rank") ? root.get("rank").asInt() : -1);
			if (rank < 0) continue;
			int maxDim = Math.max(n, Math.max(mm, p));
			String source = root.has("source") && !root.get("source").asString().isBlank()
					? root.get("source").asString() : "unknown";
			boolean isComplex = SchemeIO.isComplex(root);
			boolean isZ2 = SchemeIO.isZ2(root);

			String field = isComplex ? "C" : (isZ2 ? "F2" : "R/Q/Z");

			// Additions is a per-individual metric → READ, don't recompute. The
			// manifest must never expand a (possibly large) scheme just to count
			// additions; that is per-individual work owned by the validate/enrich
			// step (same rule already enforced for `min_additions` below).
			//
			// `additions` here is the FLAT/STRUCTURAL count = additionCount(U,V,W),
			// intrinsic to the matrices. Its canonical home is the filename `_a{N}`
			// token, which is part of the scheme's IDENTITY — lineage refs key on
			// (shape, rank, adds) (#155) — so it is authoritative, always present,
			// and stable. We therefore read the filename FIRST; an LSP/CSE
			// improvement never touches it (that lowers the *scheduled* count →
			// `min_additions`/`scheduled_additions`, below, not `_a`). Source order:
			// filename `_a{N}` → legacy JSON `additions` → expand-and-count last.
			// Content-driven (was: filename `_a{N}` first). StampAdditions backfilled it
			// from the token into JSON, so content is now authoritative.
			Integer additions = (root.has("additions") && root.get("additions").isInt())
					? root.get("additions").asInt() : null;

			// Expand the explicit scheme ONLY when a downstream projection needs the
			// matrices: the human-friendly render (small schemes, #188) or a
			// still-missing additions count. A large scheme whose additions came
			// from JSON/filename is never expanded here — the whole point of the
			// individual→catalog contract (user 2026-06-06).
			NonCubicBilinearAlgorithm realAlg = null;
			if (additions == null || maxDim <= HUMAN_MAX_DIM) {
				try {
					if (isComplex) {
						ComplexNonCubicBilinearAlgorithm alg = SchemeIO.readComplex(root);
						if (additions == null) additions = Verifier.additionCount(alg);
					} else if (SchemeIO.isReduced(root)) {
						// _reduced files have CSE-encoded structure; use the reduced reader.
						realAlg = SchemeIO.readReduced(root);
						if (additions == null) additions = Verifier.additionCount(realAlg);
					} else {
						realAlg = SchemeIO.read(root);
						if (additions == null) additions = Verifier.additionCount(realAlg);
					}
				} catch (Exception ignored) {
					// best-effort; large/stub schemes may not expand — additions then
					// stays as sourced from JSON/filename (or null).
					realAlg = null;
				}
			}

			// CSE-minimised additive complexity (#189/#190): READ ONLY here — the
			// minimal SLP is expensive to derive, so it is computed by a dedicated
			// band-by-band process (MaterialiseAdditionsSlp) that stamps
			// "min_additions" + the full "slp" into each scheme JSON. The manifest
			// merely surfaces the precomputed number (user 2026-06-04: "the SPA
			// should barely ever compute"; analyse() must not run in manifest gen).
			Integer minAdditions = null;
			if (root.has("min_additions") && root.get("min_additions").isInt()) {
				minAdditions = root.get("min_additions").asInt();
			}

			// The displayed Source is always the scheme file's OWN source (the
			// importer/discoverer encoded in the filename prefix), NOT the
			// historical rank origin. We deliberately ignore attribution_for_rank
			// here: e.g. AlphaTensor's ⟨2,2,2⟩=7 over F₂ shows as "AlphaTensor
			// 2022", not "Strassen 1969" (user 2026-06-03 — attribution_for_rank
			// may stay in the per-scheme JSON for provenance, but it does NOT drive
			// the catalog/SPA Source column).
			// Also propagate the commutative flag — needed so downstream consumers
			// (FrontierClosure, RankBasesByOmega, paper-table regen) can filter
			// without re-reading every per-scheme file.
			String displaySource = source;
			// Perminov is likewise an AGGREGATOR for his schemes/known/<sub> subtree:
			// those schemes ORIGINATE ELSEWHERE (he files them under "known"), so
			// crediting "Perminov 2023" misattributes the historical record — e.g.
			// ⟨2,7,7⟩=76 under known/meta_flip_graph is Kauers & Wood 2025, not Perminov.
			// AttributePerminovKnown re-stamps the on-disk `source`; this guard is the
			// durable safety net (mirrors the fmm-lille exception) so a straggler /
			// future re-import can never surface as a Perminov "discovery". Only
			// schemes/results/* stays Perminov. The Perminov-file link is preserved in
			// source_scheme_url / imported_via. (user 2026-06-29)
			if (source.toLowerCase(java.util.Locale.ROOT).contains("perminov")
					&& root.has("original_source_path")) {
				var attr = eu.solven.matmul.catalog.PerminovKnownAttribution
						.forPath(root.get("original_source_path").asString());
				if (attr.isPresent() && !attr.get().isPerminovOwn()) {
					displaySource = attr.get().source();
				}
			}
			// fmm-lille is an AGGREGATOR, not an originator: it republishes others'
			// schemes. If its digest cites the literature source for this (format,rank),
			// attribute to THAT, not to "fmm-lille" (user 2026-06-06; e.g. ⟨3,3,6⟩=40 is
			// Smirnov 2013, not fmm-lille). This is the deliberate exception to the
			// "Source = importer" rule (2026-06-03): fmm-lille doesn't discover. Checks
			// displaySource (not raw source) so a Perminov known/tensor scheme remapped
			// to "FMM-Lille" just above is further refined to its per-format literature.
			if (displaySource.toLowerCase(java.util.Locale.ROOT).startsWith("fmm")) {
				JsonNode fmmRow = fmmByCanonical.get(canonicalKey(n, mm, p));
				if (fmmRow != null && fmmRow.has("rank") && fmmRow.get("rank").asInt() == rank
						&& fmmRow.has("references") && fmmRow.get("references").isArray()
						&& !fmmRow.get("references").isEmpty()) {
					displaySource = formatFmmRef(fmmRow.get("references").get(0).asText());
				}
			}
			// A "Solven_*"-prefixed source is a DERIVED materialiser output, not a
			// Solven discovery — crediting "Solven" as the source overclaims (user
			// 2026-06-06: ⟨7,7,7⟩=249 is a Strassen-recursion derivation). Relabel to
			// the construction method; genuine "we beat both catalogs" credit lives in
			// the separate solven_discovery flag (and those are derived_recursive, not
			// Solven_*). #152/#121.
			if (displaySource.toLowerCase(java.util.Locale.ROOT).startsWith("solven")) {
				String s = displaySource.toLowerCase(java.util.Locale.ROOT);
				displaySource = s.contains("strassen") ? "Strassen (recursive)"
						: s.contains("closure") ? "derived (closure search)"
						: "derived";
			}
			boolean commutative = false;
			String lineageCompact = null;
			// atom = this scheme is a PRIMITIVE, not composed by us from other catalog
			// entries. Derived from the lineage root: absent lineage, or a single
			// Atom/Leaf node (an explicit import, or a formula-constructor ref such as
			// DIS09Lemma4(n=…)) ⇒ atom. Any composition op at the root (Kron / Concat /
			// SumInner / Recombination / DisjointSum / SerendipitousProduct / Project /
			// AxisFlip / …) ⇒ NOT an atom. Lets the catalog (and the FMM comparison)
			// segment "what we imported/derived-by-formula" from "what we composed".
			boolean atom = true;
			boolean hasDanglingBase = false;
			try {
				JsonNode raw = mapper.readTree(f);
				if (raw != null) {
					if (raw.has("commutative")) {
						commutative = raw.get("commutative").asBoolean(false);
					}
					if (raw.has("lineage_compact")) {
						lineageCompact = prettyLineageCompact(raw.get("lineage_compact").asText());
					}
					atom = isAtomLineage(raw.get("lineage"));
					// Flag any pinned base ref ("{shape}@{hash}") that does not resolve
					// to an on-disk scheme (dangling — see existingPinnedKeys above).
					java.util.List<String> pinned = new java.util.ArrayList<>();
					collectPinnedRefs(raw.get("lineage"), pinned);
					for (String ref : pinned) {
						int at = ref.indexOf('@');
						String key = ref.substring(0, at) + "@"
								+ ref.substring(at + 1, Math.min(at + 8, ref.length()));
						if (!existingPinnedKeys.contains(key)) {
							hasDanglingBase = true;
							danglingBaseFiles.add(f.getName() + "  ⟶  " + key);
						}
					}
				}
				// Safety net for derived schemes whose lineage field was NOT persisted
				// (e.g. some Solven_Strassen / Solven_Closure / derived_recursive
				// materialiser outputs): they are derived BY CONSTRUCTION, so they must
				// not read as atoms even with a missing lineage. Source-prefix override
				// (user 2026-06-06; the real fix is to backfill their lineage — #177/#152).
				if (atom && isComposerSource(source)) atom = false;
			} catch (RuntimeException ignored) {}
			// Fallback: the _commutative filename suffix is an older convention.
			if (!commutative && f.getName().contains("_commutative")) {
				commutative = true;
			}

			ObjectNode entry = mapper.createObjectNode();
			// Canonicalise the shape to sorted n ≤ m ≤ p. The matmul tensor is
			// S₃-symmetric, so ⟨n,m,p⟩ and every axis permutation share the same
			// rank / additions / bud_score; cataloguing the canonical (sorted)
			// shape — the FMM convention used elsewhere in this repo — collapses
			// equivalent-orientation duplicates (e.g. ⟨11,17,2⟩ vs ⟨2,11,17⟩,
			// emitted by our own recursive/Rosowski derivations) and lets a
			// better-rank variant found in a non-canonical orientation win its
			// shape-class. The on-disk file keeps its native axis order (it is an
			// equivalent witness); only this catalogue VIEW is canonicalised.
			int[] sortedShape = { n, mm, p };
			java.util.Arrays.sort(sortedShape);
			ArrayNode fmt = mapper.createArrayNode();
			fmt.add(sortedShape[0]); fmt.add(sortedShape[1]); fmt.add(sortedShape[2]);
			entry.set("format", fmt);
			entry.put("max_dim", maxDim);
			// Explicit per-scheme fields[] (task #175): the full set of fields the
			// scheme verifies over (e.g. [F2,F3,Z,Q,R,C] for an integer scheme,
			// [Q,R,C] for ½-coefficient, [F2] for AlphaTensor). Drives the SPA field
			// filter by membership rather than the coarse `field` cluster label.
			// ALWAYS emitted (canonically ordered) so the SPA never has to compute
			// the field list itself (user 2026-06-03: "the SPA should barely ever
			// compute"). When no verified tags exist we expand the coarse `field`
			// cluster label in Java here, not in catalog.js.
			java.util.List<String> fieldTags = SchemeIO.fieldTags(root);
			if (fieldTags.isEmpty()) {
				// NO FABRICATION: no fields[] on disk → no fields[] in the catalog
				// (user 2026-06-12, "no fields → no fields"). We do NOT invent a field
				// here — neither the coarse "R/Q/Z"→[Z,Q,R] cluster default nor a
				// lineage-inferred char-0 floor. Both are field-discipline lies: the
				// cluster default stamped a bogus Z on the derived ⟨3,3,13⟩=89 (and
				// ⟨3,3,8⟩=55) stubs even though their smirnov13 leaf holds 1/8
				// coefficients (rational, NOT integer), and the lineage floor silently
				// drops F₂/F₃ membership. The honest fix is to compute the real fields
				// ON DISK (materialise the stub → content truth, via
				// BackfillMissingFields); until a scheme carries a fields[], it is
				// emitted with an empty membership and the SPA renders it as
				// unclassifiable (never silently dropped — catalog.js guards on
				// `fields.length === 0`).
				log.warn("manifest: {} has no fields[]; emitting EMPTY (no fabrication — "
						+ "run BackfillMissingFields to compute it on disk)", f.getName());
			}
			ArrayNode fa = mapper.createArrayNode();
			orderFields(fieldTags).forEach(fa::add);
			// `fields` is the SINGLE source of truth for the algebra a scheme is
			// valid over (user 2026-06-04). We deliberately do NOT emit the coarse
			// singular `field` cluster string anymore: collapsing {F2,…,R,C} into
			// one "R/Q/Z" label hid F₂-membership, so e.g. Strassen (valid over F₂)
			// never competed with — and could not dominate — an F₂-only re-discovery
			// at the same shape (AlphaTensor ⟨2,2,2⟩=7,+24 vs Strassen 7,+18).
			// The SPA now groups/dominates/filters/sorts on this membership array.
			entry.set("fields", fa);
			// zt (#…): a SUB-CLASS marker on integer schemes, NOT a field. zt =
			// "Z scheme whose every U/V/W coefficient is in {-1,0,1}" (Perminov's
			// "ternary integer" Z-target — distinct from F₃'s ternary modular; it
			// has nothing to do with F₂/Z₂ / characteristic 2). The flag is
			// COMPUTED & STAMPED into each scheme JSON by the MaterialiseZT
			// sanitization procedure; here we PREFER the stored value (user
			// 2026-06-04: "the catalog would rely on the given field").
			//
			// CONTENT-DRIVEN FALLBACK: when no `zt` key is stamped (a folder the
			// stamp pass never covered — e.g. the whole `derived/` tree, which was
			// 0/7784 stamped — or a freshly written scheme) but the scheme IS an
			// integer scheme (Z ∈ fields) AND its matrices expanded above, compute
			// the flag straight from the coefficients via the SAME definition
			// (SchemeIO.isTernary, every U/V/W coeff ∈ {-1,0,1}). This keeps the
			// manifest correct regardless of whether MaterialiseZT has been re-run,
			// matching the "read metadata from CONTENT, never a side pass" rule. It
			// reuses the already-expanded realAlg, so it adds no read cost; stubs
			// (realAlg == null) stay omitted exactly as before.
			Boolean zt = SchemeIO.readZT(root);
			if (zt == null && realAlg != null && fieldTags.contains("Z")) {
				zt = SchemeIO.isTernary(realAlg);
			}
			if (zt != null) {
				entry.put("zt", zt);
			}
			entry.put("rank", rank);
			// NOTE: singular `field` intentionally NOT emitted — see the `fields`
			// set() above. Everything keys on the membership array now.
			if (additions != null) entry.put("additions", additions);
			if (minAdditions != null) entry.put("min_additions", minAdditions);
			// Scheduled (CSE-shared) addition count when the scheme carries one
			// (#146/#185 — e.g. Strassen-Winograd ⟨2,2,2⟩=7 is 15 scheduled vs 24
			// flat). Surfaced so the SPA can show the special scheduled count.
			if (root.has("scheduled_additions") && root.get("scheduled_additions").isInt()) {
				entry.put("scheduled_additions", root.get("scheduled_additions").asInt());
			}
			// Forward the straight-line program (the CSE schedule that realises the
			// scheduled additions) so the SPA modal can show it — provenance for a
			// record like Perminov 2025's 58-addition ⟨3,3,3⟩.
			if (root.has("slp") && root.get("slp").isArray()) {
				entry.set("slp", root.get("slp"));
			}
			// Human-friendly flat multiplications/elements (#188), Perminov-style.
			// Capped to small schemes (rank ≤ HUMAN_MAX_RANK, maxDim ≤ HUMAN_MAX_DIM)
			// to keep catalog.json lean — these are the schemes where a textual form
			// is actually readable; the SPA modal renders them verbatim.
			if (realAlg != null && realAlg.r <= HUMAN_MAX_RANK && maxDim <= HUMAN_MAX_DIM) {
				HumanScheme.Readable hr = HumanScheme.of(realAlg);
				ArrayNode mults = mapper.createArrayNode();
				hr.multiplications().forEach(mults::add);
				ArrayNode elems = mapper.createArrayNode();
				hr.elements().forEach(elems::add);
				entry.set("multiplications", mults);
				entry.set("elements", elems);
			}
			// Bud structure (#159): predicts serendipitous-product potential. The
			// SOURCE OF TRUTH is the individual JSON's stamped `buds`/`has_buds`
			// (written by the Phase-2 SchemeAnalysis.Buds validate step) — the
			// catalog merely FORWARDS it, so the manifest stays a metadata
			// projection rather than an (expensive) re-expander, and so bud-aware
			// search can read it ref-level from the catalog for ANY scheme,
			// including maxDim>16 stubs once validated (same lifecycle as
			// `verified`). For small schemes not yet stamped we fall back to an
			// inline O(r·dim) compute so nothing regresses before the backfill runs.
			if (root.has("buds") || root.has("has_buds")) {
				boolean hasBuds = root.has("buds") || root.path("has_buds").asBoolean(false);
				entry.put("has_buds", hasBuds);
				if (root.has("buds")) entry.put("buds", root.get("buds").asText());
			} else if (realAlg != null && maxDim <= HUMAN_MAX_DIM) {
				var bs = eu.solven.matmul.catalog.SerendipitousBudProduct.summarise(realAlg);
				if (bs.hasBuds()) {
					entry.put("has_buds", true);
					entry.put("buds", bs.summary());
				}
			}
			// Projection margin μ (paper §projmargin): forwarded from the stamped
			// JSON (Phase-2 SchemeAnalysis.ProjectionMargin — computed in the SAME
			// single expansion as buds), with an inline fallback for small schemes
			// not yet stamped. R−μ is the rank one index down; high μ ⇒ strong
			// downward (projection) parent even at higher rank.
			// Pan-TA highlight (user 2026-06-26): when this scheme's lineage root is a
			// naïve-grid recombination, surface WHERE Pan trilinear aggregation saved
			// multiplications, so the final rank is explainable (rank = unpaired leaves +
			// fused-pair cost; TA bought `saving`). TA is a saving WITHIN the recombination's
			// multiplications — not a separate strategy — so it is highlighted on the node,
			// recomputed from base+allocs exactly as build/replay does. Structured for the
			// SPA + a one-line summary for displays/logs.
			try {
				java.util.Optional<eu.solven.matmul.catalog.Lineage.Node> taLn = SchemeIO.readLineage(root);
					taLn.ifPresent(baseUsage::accumulate);   // base-usage stats by op
				if (taLn.isPresent()
						&& taLn.get() instanceof eu.solven.matmul.catalog.Lineage.RecombinationTaN taNode) {
					eu.solven.matmul.catalog.TaFusionExplainer
							.describe(taNode, (a, b, c) -> lazyLookup().findRank(a, b, c))
							.filter(eu.solven.matmul.recombination.Recombination.TaFusionBreakdown::hasFusion)
							.ifPresent(bd -> {
								ObjectNode ta = mapper.createObjectNode();
								ta.put("pairs", bd.fusedPairs().size());
								ta.put("saving", bd.taSaving());
								ta.put("fused_cost", bd.fusedCost());
								ta.put("unpaired_leaf_sum", bd.unpairedLeafSum());
								ArrayNode fused = mapper.createArrayNode();
								for (var fp : bd.fusedPairs()) {
									ObjectNode pr = mapper.createObjectNode();
									ArrayNode shape = mapper.createArrayNode();
									shape.add(fp.shapeA()[0]);
									shape.add(fp.shapeA()[1]);
									shape.add(fp.shapeA()[2]);
									pr.set("shape", shape);
									pr.put("fused_cost", fp.fusedCost());
									pr.put("naive_rank", fp.naiveRank());
									pr.put("saving", fp.saving());
									fused.add(pr);
								}
								ta.set("fused", fused);
								ta.put("summary", bd.summary());
								entry.set("ta_fusion", ta);
							});
				}
			} catch (Exception ignored) {
				// best-effort highlight: never block manifest generation
			}
			if (root.has("projection_margin")) {
				entry.put("projection_margin", root.get("projection_margin").asInt());
			} else if (realAlg != null && maxDim <= HUMAN_MAX_DIM) {
				entry.put("projection_margin",
						eu.solven.matmul.catalog.ProjectionSearch.projectionMargin(realAlg));
			}
			// Corrupted flag (#68): a lineage stub whose recipe can't be replayed.
			// Surfaced so the SPA can show it as "rank claimed, not reproducible" and
			// so FieldAwareLookup can treat it as absent for search gating.
			if (SchemeIO.isCorrupted(root)) {
				entry.put("corrupted", true);
				corruptedFiles.add(f.getName());
			}
			// Dangling pinned base (this run): a lineage ref pins a sub-scheme that
			// is not on disk → the SPA graph can only show ≈ best-known for it.
			if (hasDanglingBase) {
				entry.put("dangling_base", true);
			}
			// Bud certainty + derived-scheme inference. A scheme whose buds came
			// from an EXPANDED scheme (small / forwarded) is certified `exact`. A
			// big composite that wasn't expanded gets its buds INFERRED from the
			// lineage (no expansion — only atom leaves are read), tagged `exact`
			// for cancellation-free ops or `structural-estimate` otherwise. See
			// LineageBudInference / the buds-from-lineage contract.
			if (entry.has("has_buds") && entry.get("has_buds").asBoolean(false)) {
				entry.put("buds_certainty", "exact");
			} else if (realAlg == null) {
				try {
					java.util.Optional<eu.solven.matmul.catalog.Lineage.Node> ln =
							SchemeIO.readLineage(root);
					if (ln.isPresent()) {
						eu.solven.matmul.catalog.LineageBudInference.Profile prof =
								eu.solven.matmul.catalog.LineageBudInference.infer(
										ln.get(), atomResolver, leafShapeResolver);
						if (prof.known() && prof.hasBuds()) {
							entry.put("has_buds", true);
							entry.put("buds", prof.summary());
							entry.put("buds_certainty",
									prof.certainty() == eu.solven.matmul.catalog.LineageBudInference.Certainty.EXACT
											? "exact" : "structural-estimate");
						}
					}
				} catch (Exception ignored) {
					// best-effort: leave unbudded
				}
			}
			entry.put("source", displaySource);
			// Provenance links (content-driven, optional): source_scheme_url points at
			// the scheme's file in the originating author's own repo (e.g. Perminov's
			// FastMatrixMultiplication); source_paper_url / source_url point at the
			// publication. Forwarded verbatim so the SPA can render clickable links.
			for (String urlKey : new String[] { "source_scheme_url", "source_paper_url", "source_url" }) {
				if (root.has(urlKey) && !root.get(urlKey).asString().isBlank()) {
					entry.put(urlKey, root.get(urlKey).asString());
				}
			}
			// Publication year of the rank's origin (prefer the JSON `year` field
			// when present, else infer from the source/attribution label). 9999 =
			// "unknown / not a dated reference" (compositions) → emit null.
			int year = root.has("year") && root.get("year").isInt()
					? root.get("year").asInt()
					: yearOfSource(displaySource);
			if (year != 9999) entry.put("year", year);
			// Honesty flag (#159 pipeline): a scheme is verified once its expanded
			// matrices have passed the symbolic/spot-check. Read the per-scheme JSON
			// flag, defaulting to TRUE when absent — every legacy on-disk scheme was
			// built+verified before being written, so absence means verified. A
			// detect-only stub (search prediction, not yet machine-checked) writes
			// "verified": false and is confirmed (or pruned) by the Phase-2 batch.
			boolean verified = SchemeIO.isVerified(root);
			entry.put("verified", verified);
			entry.put("atom", atom);      // primitive (import/formula) vs composed (#198 analysis)
			entry.put("complex", isComplex);
			if (commutative) entry.put("commutative", true);
			if (lineageCompact != null) entry.put("lineage_compact", lineageCompact);

			// Relative path under schemes root, with section{N} prefix.
			Path rel = SCHEMES_DIR.toPath().relativize(f.toPath());
			entry.put("file", rel.toString().replace(File.separatorChar, '/'));
			// explicitable: is this a precise EXPLICIT scheme (matrices, or a lineage that
			// bottoms out entirely in @hash/@naive/named precise leaves) — vs a CITED BOUND
			// (a best-at-shape @sota/-direct/bare leaf, or a dangling pin). NOT "buildable":
			// an HK71 cited bound has a buildable formula yet is not (yet) an explicit scheme.
			entry.put("explicitable", explicitableByFile.getOrDefault(f.toPath(), Boolean.FALSE));

			// Cross-reference with fmm-lille catalog (sorted-format key).
			JsonNode fmm = fmmByCanonical.get(canonicalKey(n, mm, p));
			boolean fmmCited = false;
			if (fmm != null) {
				ObjectNode xref = mapper.createObjectNode();
				if (fmm.has("rank")) xref.put("best_rank", fmm.get("rank").asInt());
				if (fmm.has("details_url")) xref.put("details_url", fmm.get("details_url").asText());
				if (fmm.has("references") && fmm.get("references").isArray()
						&& !fmm.get("references").isEmpty()) {
					xref.set("references", fmm.get("references"));
					fmmCited = true;
				}
				entry.set("fmm_lille", xref);
			}

			// 3-catalog comparison (#199): flag where THIS scheme beats the best
			// known external (FMM ∪ Perminov) rank. Only comparable schemes count —
			// non-commutative AND valid over a char-0 field (the external catalogs
			// are NC char-0); commutative / F₂-only schemes are not comparable.
			ExternalBest ext = externalBest.get(canonicalKey(n, mm, p));
			boolean char0 = fieldTags.stream().anyMatch(CHAR0_FIELD_TAGS::contains);
			boolean comparable = !commutative && char0;
			if (ext != null && comparable) {
				entry.put("external_best_rank", ext.rank());
				entry.put("external_best_source", ext.source());
				// cited = the format has a literature reference in FMM (Perminov is a
				// search catalog → uncited). The honest target differs by case:
				//   cited   → match a published result (check our reference is consistent);
				//   uncited → FMM-computed bound we aim to reproduce by our OWN derivation
				//             (a derived/composed lineage, atom:false), not by importing it.
				entry.put("external_cited", fmmCited);
				// PROVISIONAL: dropped automatically on the next sync if FMM/Perminov
				// publish an equal-or-better rank (they don't import our results).
				// Requires verified=true: an unverified (detect-only) stub is NOT yet a
				// discovery — it must pass Phase-2 validation before we claim it beats
				// both catalogs.
				if (verified && rank < ext.rank()) entry.put("solven_discovery", true);
			}

			schemeNodes.add(entry);
		}

		// Filter worse-than-naive entries: any scheme whose rank exceeds
		// n·m·p brings no value as a standalone matmul algorithm — it's
		// strictly worse than the trivial outer-product baseline. These are
		// almost always internal building blocks (e.g. DIS09 Lemma 4 applied
		// at small n) used by composite searches but not meaningful as
		// catalog entries on their own. On-disk JSON stays for the search
		// pool; the manifest just doesn't surface them. See user feedback
		// 2026-06-02.
		int beforeNaiveFilter = schemeNodes.size();
		schemeNodes = schemeNodes.stream().filter(e -> {
			ArrayNode fmt = (ArrayNode) e.get("format");
			long naive = (long) fmt.get(0).asInt() * fmt.get(1).asInt() * fmt.get(2).asInt();
			return e.get("rank").asLong() <= naive;
		}).collect(java.util.stream.Collectors.toList());
		log.info(String.format("dropped %d worse-than-naive entries (rank > n·m·p)%n",
				beforeNaiveFilter - schemeNodes.size()));

		// Shave noise: within each (format, field, rank) group, drop entries
		// that don't strictly improve on the best-known additions count. A
		// later paper "rediscovering" a known rank with WORSE adds doesn't
		// belong in the catalog manifest. Underlying scheme JSON files are
		// kept on disk for provenance — this only filters the manifest view.
		// See user feedback 2026-05-28.
		int beforeShave = schemeNodes.size();
		schemeNodes = shaveByBestAdditions(schemeNodes);
		log.info(String.format("shaved %d redundant entries (rediscoveries with worse adds)%n",
				beforeShave - schemeNodes.size()));

		// Register the (rank, buds) Pareto frontier per (format, field): stamps
		// `bud_score` everywhere and `pareto_rank_buds:true` on entries not
		// dominated on both axes — so a bud-richer higher-rank scheme is kept as a
		// "best" building block, not hidden behind the rank-minimal one.
		int paretoBest = markRankBudPareto(schemeNodes);
		log.info(String.format("registered %d (rank,buds) Pareto-best entries%n", paretoBest));

			// Stamp per-base usage: how many times each scheme is USED as a building block
			// (by recombination / projection / kronecker / …) across the whole catalog.
			int usedBases = stampBaseUsage(schemeNodes, baseUsage, mapper);
			log.info(String.format("stamped used_as_base on %d entries%n", usedBases));

		ObjectNode root = mapper.createObjectNode();
		// Intentionally no "generated" timestamp — it triggers merge conflicts
		// every time the file is regenerated; the file's content already
		// reflects the schemes' on-disk state, which is what matters.
		root.put("total", schemeNodes.size());
		ArrayNode schemes = mapper.createArrayNode();
		schemeNodes.forEach(schemes::add);
		root.set("schemes", schemes);
		// Full list of corrupted (unreplayable) stub basenames, INCLUDING same-shape
		// siblings the dedup dropped — FieldAwareLookup reads this to treat them as
		// absent for search gating (so an unverifiable rank can't shadow a fresh one).
		if (!corruptedFiles.isEmpty()) {
			ArrayNode corr = mapper.createArrayNode();
			corruptedFiles.forEach(corr::add);
			root.set("corrupted_files", corr);
		}
		log.info(String.format("flagged %d corrupted (unreplayable) stub files%n", corruptedFiles.size()));

		// Dangling pinned bases: lineage refs that pin a sub-scheme not present on
		// disk. Replay can still fall back to shape (catalog-best), so these are NOT
		// "corrupted", but the pinned identity is lost — warn loudly + record them.
		if (!danglingBaseFiles.isEmpty()) {
			ArrayNode dang = mapper.createArrayNode();
			danglingBaseFiles.forEach(dang::add);
			root.set("dangling_base_files", dang);
			log.warn("{} schemes reference a base that does NOT resolve (dangling pinned ref); "
					+ "the SPA shows ≈ best-known for those leaves. Examples: {}",
					danglingBaseFiles.size(),
					danglingBaseFiles.stream().limit(5).collect(java.util.stream.Collectors.joining("; ")));
		} else {
			log.info("no dangling pinned bases — every lineage ref resolves on disk");
		}

		OUTPUT.getParentFile().mkdirs();
		try (PrintWriter pw = new PrintWriter(new BufferedWriter(new FileWriter(OUTPUT)))) {
			// MatrixJsonFormatter (not Jackson's default pretty printer, whose
			// array indenter is INLINE → multiplications/elements render as one
			// very wide line). Rule: number arrays inline (format, u-rows),
			// string arrays one-per-line (multiplications, elements, fields).
			pw.print(eu.solven.matmul.catalog.MatrixJsonFormatter.format(mapper.writeValueAsString(root)));
		}
		log.info(String.format("wrote %s (%d schemes)%n", OUTPUT, schemeNodes.size()));
	}

	/**
	 * Per (format, field, rank) group, keep only entries that are not
	 * strictly dominated by another entry with EARLIER year and SMALLER
	 * additions count. Returns the surviving list.
	 *
	 * <p>"Year" is inferred from source name (filename prefix) via a small
	 * lookup table; unknown sources are treated as "modern" (year 9999),
	 * so they only survive if they have a strictly better adds count than
	 * any known prior result for the same (format, field, rank).</p>
	 */
	/**
	 * A scheme is an <em>atom</em> (primitive) iff its lineage root is absent or a
	 * single {@code Atom}/{@code Leaf} node — i.e. it is an explicit import or a
	 * formula-constructor reference (e.g. {@code DIS09Lemma4(n=…)}), not something
	 * we built by composing other catalog entries. Any composition op at the root
	 * (Kron/Concat/SumInner/Recombination/DisjointSum/SerendipitousProduct/Project/
	 * AxisFlip/AxisPermute/Transpose/Dce/AugmentSquareDiscard) ⇒ NOT an atom.
	 */
	/** Turn an FMM-digest reference token like {@code "smirnov:2013a"} or
	 *  {@code "hopcroft:1971"} into a display label {@code "Smirnov 2013"}. */
	static String formatFmmRef(String token) {
		if (token == null || token.isBlank()) return "fmm-lille";
		String author = token;
		String year = "";
		int colon = token.indexOf(':');
		if (colon >= 0) {
			author = token.substring(0, colon);
			java.util.regex.Matcher y = java.util.regex.Pattern.compile("(\\d{4})").matcher(token.substring(colon + 1));
			if (y.find()) year = " " + y.group(1);
		}
		if (author.isEmpty()) return token;
		String cap = Character.toUpperCase(author.charAt(0)) + author.substring(1);
		return cap + year;
	}

	/** True if the source label denotes one of OUR composer/materialiser outputs
	 *  (a derived scheme), never an imported primitive. Used as a fallback when
	 *  the lineage field is missing so derived schemes aren't mislabelled atoms. */
	private static boolean isComposerSource(String source) {
		if (source == null) return false;
		String s = source.toLowerCase(java.util.Locale.ROOT);
		return s.startsWith("solven") || s.startsWith("composed")
				|| s.contains("closure") || s.contains("recursive");
	}

	/** A real pinned base ref: {@code NxMxP@hash}. (An internal {@code @ref?:L0} lineage-id
	 *  fallback is NOT a catalog base — excluding it avoids a false dangling.) */
	private static final java.util.regex.Pattern PINNED_REF =
			java.util.regex.Pattern.compile("\\d+x\\d+x\\d+@[0-9a-fA-F]{4,}");

	/** Collect every pinned Atom ref ("{shape}@{hash}") anywhere in a lineage DAG. */
	private static void collectPinnedRefs(JsonNode node, List<String> out) {
		if (node == null) {
			return;
		}
		if (node.isObject()) {
			JsonNode ref = node.get("ref");
			// A pinned BASE ref is "{shape}@{hash}". Skip internal lineage-id references
			// ("@ref?:L0" — a dedup/cycle fallback that points at a lineage node, not a catalog
			// scheme); counting those as missing bases is a false dangling.
			if (ref != null && ref.isString() && PINNED_REF.matcher(ref.asString()).find()) {
				out.add(ref.asString());
			}
			node.properties().forEach(e -> collectPinnedRefs(e.getValue(), out));
		} else if (node.isArray()) {
			node.forEach(c -> collectPinnedRefs(c, out));
		}
	}

	private static boolean isAtomLineage(JsonNode lineage) {
		if (lineage == null || lineage.isNull()) return true; // no lineage = primitive import
		String op = lineage.path("op").asText("");
		return op.isEmpty() || op.equals("Atom") || op.equals("Leaf");
	}

	// ── lineage_compact display reformatting (#184) ─────────────────────────
	// The JSON `lineage` tree stays the authoritative, replayable record; this
	// only prettifies the COMPACT DISPLAY string copied into catalog.json.
	private static final java.util.regex.Pattern KRON =
			java.util.regex.Pattern.compile("Kron\\[\\s*([^;\\]]+?)\\s*;\\s*([^\\]]+?)\\s*\\](?:\\s*=\\s*\\S+)?");
	private static final java.util.regex.Pattern RECOMBINE =
			java.util.regex.Pattern.compile(
				"R\\[\\s*[^;\\]<]*?<(\\d+,\\d+,\\d+)>=(\\d+)\\s*;\\s*([^|\\]]+?)\\s*\\|\\s*([^|\\]]+?)\\s*\\|\\s*([^\\]]+?)\\s*\\]");
	/** A {@code @hash7} pins a scheme for REPLAY (the structured lineage tree);
	 *  it is noise in the human-readable compact display string, so we strip it
	 *  here at generation — never in the SPA. */
	private static final java.util.regex.Pattern LINEAGE_HASH =
			java.util.regex.Pattern.compile("@[0-9a-fA-F]{6,}");

	/** Strip a {@code source-} prefix from a lineage ref, leaving the
	 *  {@code NxMxP_m{rank}} shape token (#155 / #184). */
	private static String stripLineageSource(String ref) {
		return ref.trim().replaceAll("^[^\\s|]*?-(?=\\d+x\\d+x\\d+)", "");
	}

	/**
	 * Reformat a stored {@code lineage_compact} for display (#184):
	 * <ul>
	 *   <li>{@code Kron[src-AxBxC_mR; src-DxExF_mS] = …} → {@code AxBxC_mR ⊗ DxExF_mS}
	 *       (drop source prefixes + the redundant result shape).</li>
	 *   <li>{@code R[Src<n,m,p>=r; t1 | t2 | t3]} → {@code R[<n,m,p>=r; c=t3 | r=t1]}
	 *       (drop source name; label the output column (axis-3) and row (axis-1)
	 *       partitions; the internal contraction split t2 is omitted from the
	 *       compact view — full splits remain in the lineage tree).</li>
	 * </ul>
	 * Applied globally so nested / chained ({@code … +n R[…] +p …}) strings are
	 * transformed in place. Unrecognised content is left untouched.
	 */
	static String prettyLineageCompact(String compact) {
		if (compact == null || compact.isEmpty()) return compact;
		// Replay-only hashes (@hash7) never belong in the human display string.
		String out = LINEAGE_HASH.matcher(compact).replaceAll("");
		java.util.regex.Matcher km = KRON.matcher(out);
		StringBuilder kb = new StringBuilder();
		while (km.find()) {
			String repl = stripLineageSource(km.group(1)) + " ⊗ " + stripLineageSource(km.group(2));
			km.appendReplacement(kb, java.util.regex.Matcher.quoteReplacement(repl));
		}
		km.appendTail(kb);
		out = kb.toString();

		java.util.regex.Matcher rm = RECOMBINE.matcher(out);
		StringBuilder rb = new StringBuilder();
		while (rm.find()) {
			String shape = rm.group(1), rank = rm.group(2);
			String t1 = rm.group(3).trim(), t3 = rm.group(5).trim();
			String repl = "R[<" + shape + ">=" + rank + "; c=" + t3 + " | r=" + t1 + "]";
			rm.appendReplacement(rb, java.util.regex.Matcher.quoteReplacement(repl));
		}
		rm.appendTail(rb);
		return rb.toString();
	}

	/** Effective additions for shaving: the cheaper of the flat count and the
	 *  scheduled (CSE-shared) count, when the latter is present (#185). */
	private static int effectiveAdds(ObjectNode e) {
		int flat = e.has("additions") ? e.get("additions").asInt() : Integer.MAX_VALUE;
		int sched = e.has("scheduled_additions") ? e.get("scheduled_additions").asInt() : Integer.MAX_VALUE;
		return Math.min(flat, sched);
	}

	/** Bud-richness score for shaving's third (maximised) axis: total rank-one
	 *  terms living in buds, parsed from the stamped {@code buds} summary. 0 when
	 *  buds are absent or uncomputed (so it never spuriously rescues an entry). */
	private static int budScoreOf(ObjectNode e) {
		String buds = e.has("buds") ? e.get("buds").asText() : null;
		return eu.solven.matmul.catalog.BudParetoSelection.budScore(buds);
	}

	/** Canonical signature of an entry's {@code fields} membership array, used as
	 *  a grouping key now that the singular {@code field} cluster label is gone. */
	private static String fieldsKey(ObjectNode e) {
		JsonNode fields = e.get("fields");
		if (fields == null || !fields.isArray()) return "";
		StringBuilder sb = new StringBuilder();
		for (JsonNode f : fields) sb.append(f.asText()).append(',');
		return sb.toString();
	}

	/**
	 * Register the (rank, buds) Pareto frontier within each (format, field,
	 * commutative) group: stamp {@code bud_score} on every entry, and
	 * {@code pareto_rank_buds:true} on the entries that are not dominated on both
	 * (rank ↓, bud-richness ↑) axes. Returns the count of frontier entries.
	 * See {@link BudParetoSelection}.
	 */
	/** The 7-char content hash a scheme filename ends with ({@code …-1a2b3c4.json}). */
	private static final Pattern FILE_HASH7 = Pattern.compile("-([0-9a-f]{7})\\.json$");

	/**
	 * Stamp each entry with {@code used_as_base}: how many OTHER schemes use it as a
	 * building block, by op ({@code by_recombination} / {@code by_projection} /
	 * {@code by_kronecker} / …) plus a {@code total}. Hash-pinned references attribute to
	 * the exact scheme ({@code shape@hash7}); bare-shape references attribute to the
	 * rank-best entry of that shape (what a bare ref resolves to at build time). Returns
	 * the number of entries that received a non-empty stamp.
	 */
	private static int stampBaseUsage(List<ObjectNode> entries,
			eu.solven.matmul.catalog.BaseUsageStats usage, JsonMapper mapper) {
		// Rank-best entry per shape — the target a bare-shape reference resolves to.
		Map<String, ObjectNode> bestByShape = new java.util.LinkedHashMap<>();
		for (ObjectNode e : entries) {
			String shape = shapeOf(e);
			ObjectNode cur = bestByShape.get(shape);
			if (cur == null || e.get("rank").asLong() < cur.get("rank").asLong()) {
				bestByShape.put(shape, e);
			}
		}
		int stamped = 0;
		for (ObjectNode e : entries) {
			String shape = shapeOf(e);
			String hash7 = hash7Of(e);
			Map<String, Integer> merged = new java.util.LinkedHashMap<>();
			if (hash7 != null) {
				mergeCounts(merged, usage.forKey(shape + "@" + hash7));
			}
			if (bestByShape.get(shape) == e) {
				mergeCounts(merged, usage.forKey(shape));
			}
			if (merged.isEmpty()) {
				continue;
			}
			ObjectNode node = mapper.createObjectNode();
			int total = 0;
			for (Map.Entry<String, Integer> en : merged.entrySet()) {
				node.put("by_" + en.getKey(), en.getValue());
				total += en.getValue();
			}
			node.put("total", total);
			e.set("used_as_base", node);
			stamped++;
		}
		return stamped;
	}

	private static void mergeCounts(Map<String, Integer> into, Map<String, Integer> from) {
		for (Map.Entry<String, Integer> en : from.entrySet()) {
			into.merge(en.getKey(), en.getValue(), Integer::sum);
		}
	}

	private static String shapeOf(ObjectNode e) {
		ArrayNode f = (ArrayNode) e.get("format");
		return f.get(0).asInt() + "x" + f.get(1).asInt() + "x" + f.get(2).asInt();
	}

	/** The 7-char content hash from the entry's file name, or {@code null}. */
	private static String hash7Of(ObjectNode e) {
		JsonNode fileNode = e.get("file");
		if (fileNode == null || !fileNode.isString()) {
			return null;
		}
		Matcher m = FILE_HASH7.matcher(fileNode.asString());
		return m.find() ? m.group(1) : null;
	}

	private static int markRankBudPareto(List<ObjectNode> entries) {
		Map<String, List<ObjectNode>> groups = new java.util.LinkedHashMap<>();
		for (ObjectNode e : entries) {
			ArrayNode fmt = (ArrayNode) e.get("format");
			boolean commutative = e.has("commutative") && e.get("commutative").asBoolean(false);
			String key = fmt.get(0).asInt() + "," + fmt.get(1).asInt() + "," + fmt.get(2).asInt()
					+ "|" + fieldsKey(e) + "|" + (commutative ? "c" : "nc");
			groups.computeIfAbsent(key, k -> new java.util.ArrayList<>()).add(e);
		}
		int best = 0;
		for (List<ObjectNode> g : groups.values()) {
			int[] ranks = new int[g.size()];
			int[] budScores = new int[g.size()];
			for (int i = 0; i < g.size(); i++) {
				ranks[i] = g.get(i).get("rank").asInt();
				String buds = g.get(i).has("buds") ? g.get(i).get("buds").asText() : null;
				budScores[i] = eu.solven.matmul.catalog.BudParetoSelection.budScore(buds);
				g.get(i).put("bud_score", budScores[i]);
			}
			boolean[] front = eu.solven.matmul.catalog.BudParetoSelection.frontierMask(ranks, budScores);
			for (int i = 0; i < g.size(); i++) {
				if (front[i]) {
					g.get(i).put("pareto_rank_buds", true);
					best++;
				}
			}
		}
		return best;
	}

	private static List<ObjectNode> shaveByBestAdditions(List<ObjectNode> entries) {
		Map<String, List<ObjectNode>> groups = new java.util.LinkedHashMap<>();
		for (ObjectNode e : entries) {
			ArrayNode fmt = (ArrayNode) e.get("format");
			// Group by (format, field, commutative, rank). Commutative is
			// load-bearing: a commutative-only scheme does NOT lift to
			// recursive matmul over a non-commutative ring, so it can't shave
			// an NC alternative at the same nominal (format, field, rank).
			// Pre-fix this was missing — Waksman 1970 commutative ⟨3,3,3⟩=23
			// silently shaved Laderman 1976 NC ⟨3,3,3⟩=23 because they fell
			// in the same group.
			boolean commutative = e.has("commutative") && e.get("commutative").asBoolean(false);
			// Also split by PROVENANCE CATEGORY (known / derived / curated /
			// constructed — the schemes/ top-level folders, with the atom flag as
			// fallback for legacy layouts): we keep ONE representative of EACH
			// methodology per (format,field,rank), since the "what can a derived vs
			// a flip-graph/imported vs a formula-constructed scheme reach"
			// comparison is itself catalog value (user 2026-06-06; extended
			// 2026-06-10: keep the best from derived AND from not-derived, same
			// rule for constructed — HK71 stays interesting for history; the SPA
			// is what merges/rejects entries). So a derived rediscovery of an
			// imported rank, or an HK71 constructed scheme tying an FMM import,
			// survives instead of being shaved as a worse-adds dup.
			String key = fmt.get(0).asInt() + "," + fmt.get(1).asInt() + "," + fmt.get(2).asInt()
					+ "|" + fieldsKey(e)
					+ "|" + (commutative ? "c" : "nc")
					+ "|" + provenanceCategory(e)
					+ "|" + e.get("rank").asInt();
			groups.computeIfAbsent(key, k -> new java.util.ArrayList<>()).add(e);
		}
		List<ObjectNode> kept = new java.util.ArrayList<>();
		for (List<ObjectNode> group : groups.values()) {
			if (group.size() == 1) { kept.addAll(group); continue; }
			// Sort by (year asc, EFFECTIVE additions asc nulls last). Effective =
			// min(flat additions, scheduled_additions): a scheme with a cheaper
			// SCHEDULE (CSE-shared) is a genuine improvement even if its flat count
			// is higher — e.g. Strassen-Winograd ⟨2,2,2⟩ flat 24 but scheduled 15
			// beats Strassen's 18, so it must NOT be shaved (#185).
			group.sort((a, b) -> {
				int yA = yearOfSource(a.get("source").asText());
				int yB = yearOfSource(b.get("source").asText());
				if (yA != yB) return Integer.compare(yA, yB);
				return Integer.compare(effectiveAdds(a), effectiveAdds(b));
			});
			// Two secondary axes within a (format,field,comm,category,rank) group:
			// additions (minimised) AND bud-richness (maximised). Buds is a genuine
			// third optimisation axis — a same-rank scheme with MORE reusable buds is
			// a better building block (paper §serendipitous), so it must NOT be shaved
			// just because a leaner-additions sibling exists. Keep the oldest, plus any
			// later entry that strictly improves EITHER axis (additions ↓ or buds ↑).
			int bestAdds = Integer.MAX_VALUE;
			int bestBuds = Integer.MIN_VALUE;
			for (ObjectNode e : group) {
				int adds = effectiveAdds(e);
				int buds = budScoreOf(e);
				if (kept.isEmpty() || e == group.get(0) || adds < bestAdds || buds > bestBuds) {
					kept.add(e);
					if (adds < bestAdds) bestAdds = adds;
					if (buds > bestBuds) bestBuds = buds;
				}
			}
		}
		return kept;
	}

	/**
	 * Provenance category of a manifest entry, from its on-disk top-level
	 * folder: {@code known} (imports), {@code derived} (our closure output),
	 * {@code curated} (preserve-only), {@code constructed} (formula
	 * constructors, e.g. HK71). Dedup keeps the best entry of EACH category per
	 * (format,field,comm,rank) — cross-category election happens in the SPA.
	 * Legacy layouts without a category folder fall back to the atom flag
	 * (atom=true → known-like import, atom=false → derived-like composite).
	 */
	private static String provenanceCategory(ObjectNode e) {
		String file = e.has("file") ? e.get("file").asText() : "";
		int slash = file.indexOf('/');
		if (slash > 0) {
			String top = file.substring(0, slash);
			switch (top) {
				case "known":
				case "derived":
				case "curated":
				case "constructed":
					return top;
				default:
					break;
			}
		}
		boolean atom = !e.has("atom") || e.get("atom").asBoolean(true);
		return atom ? "known" : "derived";
	}

	/**
	 * Publication year inferred from the source label. We look for an
	 * embedded 4-digit year (e.g. "laderman-1976", "Schwartz-Zwecher 2025",
	 * "AlphaTensor 2022") — that's the convention used in both
	 * attribution_for_rank strings and filename prefixes. Falls back to a
	 * small alias table for sources whose name doesn't embed the year
	 * (Strassen, Laderman bare, AlphaEvolve, FMM-Lille aggregator).
	 *
	 * <p>Returns 9999 when no year is recoverable — these entries lose
	 * sort priority to any known-dated competitor in shaveByBestAdditions.</p>
	 */
	/** Canonical field display order: prime fields first, then char-0 by
	 *  containment (Z ⊂ Q ⊂ R ⊂ C). */
	private static final java.util.List<String> FIELD_ORDER =
			java.util.List.of("F2", "F3", "Z", "Q", "R", "C");

	/** Expand the coarse `field` cluster label (the only non-tagged values the
	 *  manifest produces: "C", "F2", "R/Q/Z") into an explicit field list, so
	 *  every catalog.json entry carries a `fields[]` the SPA can render verbatim.
	 *  The "R/Q/Z" catch-all means "char-0, not narrowed" → [Z,Q,R] (we do NOT
	 *  claim F2/F3 without verification). */
	private static volatile FieldAwareLookup LAZY_LOOKUP;

	/** Shared lookup for lineage field-inference on the (now-exceptional) empty-fields
	 *  path; built once, lazily (index build is non-trivial). */
	private static FieldAwareLookup lazyLookup() {
		FieldAwareLookup l = LAZY_LOOKUP;
		if (l == null) {
			synchronized (GenerateCatalogManifest.class) {
				l = LAZY_LOOKUP;
				if (l == null) LAZY_LOOKUP = l = new FieldAwareLookup("Q");
			}
		}
		return l;
	}

	/** Narrowest field → inclusion-correct set (mirrors {@code StampFields.expand}):
	 *  Z grants F₂/F₃ (mod-p reduction); Q/R/C lift upward only. */
	private static java.util.List<String> expandField(eu.solven.matmul.algebra.Field f) {
		switch (f.name()) {
			case "F2": return java.util.List.of("F2");
			case "F3": return java.util.List.of("F3");
			// Z (integer): char-0 inclusion chain + mod-p reduction to every prime.
			case "Z": return java.util.List.of("F2", "F3", "Z", "Q", "R", "C");
			case "Q": return java.util.List.of("Q", "R", "C");
			case "R": return java.util.List.of("R", "C");
			case "C": return java.util.List.of("C");
			default: throw new IllegalArgumentException("Unknown field, refusing to silently widen: " + f);
		}
	}

	private static java.util.List<String> expandFieldCluster(String field) {
		if (field == null || field.isEmpty()) return java.util.List.of();
		switch (field) {
			case "R/Q/Z": return java.util.List.of("Z", "Q", "R");
			case "C": return java.util.List.of("C");
			case "F2": return java.util.List.of("F2");
			case "F3": return java.util.List.of("F3");
			default:
				java.util.List<String> out = new java.util.ArrayList<>();
				for (String t : field.split("/")) { t = t.trim(); if (!t.isEmpty()) out.add(t); }
				return out;
		}
	}

	/** Return the fields in canonical {@link #FIELD_ORDER}; unknown tags sort last. */
	private static java.util.List<String> orderFields(java.util.List<String> fields) {
		java.util.List<String> out = new java.util.ArrayList<>(fields);
		out.sort((a, b) -> {
			int ia = FIELD_ORDER.indexOf(a), ib = FIELD_ORDER.indexOf(b);
			return Integer.compare(ia < 0 ? 99 : ia, ib < 0 ? 99 : ib);
		});
		return out;
	}

	private static final java.util.regex.Pattern YEAR_PATTERN =
			java.util.regex.Pattern.compile("(?<![0-9])(1[89][0-9]{2}|20[0-2][0-9])(?![0-9])");
	private static int yearOfSource(String source) {
		if (source == null || source.isEmpty()) return 9999;
		// Extract embedded 4-digit year (1800-2029) first — the convention
		// used by attribution strings AND most filenames.
		java.util.regex.Matcher m = YEAR_PATTERN.matcher(source);
		if (m.find()) return Integer.parseInt(m.group(1));
		// Aliases for sources without embedded year.
		String s = source.toLowerCase();
		if (s.startsWith("strassen")) return 1969;
		if (s.startsWith("laderman")) return 1976;
		if (s.startsWith("waksman")) return 1970;
		if (s.startsWith("makarov")) return 1986;
		if (s.startsWith("smirnov")) return 2013;
		if (s.startsWith("rosowski")) return 2019;
		if (s.startsWith("islam")) return 2009;
		if (s.startsWith("hopcroft-kerr") || s.startsWith("hopcroftkerr")) return 1971;
		if (s.startsWith("pan-")) return 1980;
		if (s.startsWith("winograd")) return 1971;
		if (s.startsWith("dis09") || s.startsWith("dis2009")) return 2009;
		if (s.startsWith("alphatensor")) return 2022;
		if (s.startsWith("alphaevolve")) return 2025;
		if (s.startsWith("perminov") || s.startsWith("dronperminov")) return 2023;
		if (s.startsWith("fmm-lille")) return 2024;
		if (s.startsWith("derived-") || s.startsWith("solven-")) return 9999;
		return 9999;
	}

	private static String capitalize(String s) {
		if (s == null || s.isEmpty()) return s;
		return Character.toUpperCase(s.charAt(0)) + s.substring(1);
	}

	/**
	 * Canonicalises a source label. The raw value is either a filename source
	 * token (e.g. {@code dis09_Q}, {@code fmm_lille}) or a curated
	 * {@code attribution_for_rank} string (e.g. {@code "Moosbauer (symmetric
	 * flips)"}). We collapse the inconsistent variants the catalog accumulated
	 * (filename title-casing carried field cruft like {@code _Q}, and several
	 * spellings of the same source coexisted) into one canonical label per
	 * source.
	 *
	 * <ul>
	 *   <li>{@code dis09*} → {@code "DIS09"} (Dumas-Iliopoulos-Saunders 2009;
	 *       per-Lemma precision tracked separately).</li>
	 *   <li>{@code fmm[-_]lille*} → {@code "FMM-Lille"}.</li>
	 *   <li>{@code moosbauer*} → {@code "Moosbauer-Poole 2025"}.</li>
	 * </ul>
	 */
	public static String normalizeSource(String raw) {
		if (raw == null || raw.isEmpty()) return raw;
		// Strip a trailing field-tag token: the filename encodes the algebra as a
		// suffix (_F2 / _F3 / _Z / _Q / _R / _C / _ZT). The field lives in the
		// dedicated field/fields columns — it doesn't belong in the Source label
		// (e.g. "alphatensor_F2" → "AlphaTensor 2022", not "Alphatensor_F2").
		String stripped = raw.replaceAll("(?i)[_-](F2|F3|ZT|Z|Q|R|C)$", "");
		// "Derived_*" provenance: a recursion / Kronecker / concat / recombine
		// result is *derived* from smaller schemes. On-disk filenames now use the
		// `derived_` prefix (renamed from the historical `composed_` in 2026-06);
		// the legacy `composed_` token is still folded to `Derived_` here for any
		// stragglers or externally-supplied attribution strings.
		stripped = stripped.replaceFirst("(?i)^composed([_-])", "Derived$1");
		String low = stripped.toLowerCase();
		if (low.startsWith("dis09")) return "DIS09";
		if (low.startsWith("fmm_lille") || low.startsWith("fmm-lille")) return "FMM-Lille";
		if (low.startsWith("moosbauer")) return "Moosbauer-Poole 2025";
		if (low.startsWith("alphatensor")) return "AlphaTensor 2022";
		if (low.startsWith("alphaevolve")) return "AlphaEvolve 2025";
		if (low.startsWith("perminov_2025") || low.startsWith("perminov-2025")) return "Perminov 2025";
		if (low.startsWith("perminov")) return "Perminov 2023";
		// Default: Title-Case each '_'/'-'-separated segment, preserving separators
		// (e.g. dumas_pernet_sedoglavic_2025 → Dumas_Pernet_Sedoglavic_2025).
		return titleCaseSegments(stripped);
	}

	/** Upper-cases the first letter of each {@code _}/{@code -}/space-separated
	 * segment, keeping separators intact. */
	private static String titleCaseSegments(String s) {
		StringBuilder sb = new StringBuilder(s.length());
		boolean atSegStart = true;
		for (int i = 0; i < s.length(); i++) {
			char c = s.charAt(i);
			if (c == '_' || c == '-' || c == ' ') {
				sb.append(c);
				atSegStart = true;
			} else {
				sb.append(atSegStart ? Character.toUpperCase(c) : c);
				atSegStart = false;
			}
		}
		return sb.toString();
	}

	/** Best external (non-commutative, char-0) rank per canonical format across
	 *  the reference catalogs we track, with which source achieved it. */
	record ExternalBest(int rank, String source) {}

	/**
	 * Merge the external reference catalogs — FMM-Lille
	 * (`references/catalogs/fmm-lille-catalog.json`) and Perminov
	 * (`references/catalogs/perminov-catalog.json`) — into the best (lowest) known
	 * non-commutative rank per canonical format, recording which source has it.
	 * This is the baseline the manifest flags our schemes against: a
	 * {@code solven_discovery} is a format where our comparable rank beats the
	 * min of BOTH. The mark is provisional — re-derived on every catalog
	 * regeneration, so it drops automatically when a future FMM/Perminov sync
	 * publishes an equal-or-better rank (they do not import our results).
	 */
	private static Map<String, ExternalBest> loadExternalBest(JsonMapper mapper) throws IOException {
		Map<String, ExternalBest> out = new HashMap<>();
		String[][] sources = {
			{"references/catalogs/fmm-lille-catalog.json", "fmm-lille"},
			{"references/catalogs/perminov-catalog.json", "perminov"},
			// Perminov's serendipitous 17–32 band — a SEPARATE additional catalog
			// (his status.json stops at 16). Distinct source label so it attributes
			// to the June-2026 paper "Meta Flip Graph meets Serendipitous Product"
			// (arXiv:2606.02480), not folded into the generic status.json "perminov".
			{"references/catalogs/perminov-serendipitous-catalog.json", "perminov-serendipitous"},
		};
		for (String[] src : sources) {
			File f = new File(src[0]);
			if (!f.isFile()) continue;
			JsonNode entries = mapper.readTree(f).get("entries");
			if (entries == null || !entries.isArray()) continue;
			for (JsonNode row : entries) {
				JsonNode fmt = row.get("format");
				if (fmt == null || !fmt.isArray() || fmt.size() != 3 || !row.has("rank")) continue;
				String key = canonicalKey(fmt.get(0).asInt(), fmt.get(1).asInt(), fmt.get(2).asInt());
				int rank = row.get("rank").asInt();
				ExternalBest prev = out.get(key);
				if (prev == null || rank < prev.rank()) out.put(key, new ExternalBest(rank, src[1]));
			}
		}
		return out;
	}

	/**
	 * Load `references/catalogs/fmm-lille-catalog.json` (produced by
	 * `SyncReferenceCatalogs --fmm`) and index its rows by canonical
	 * (sorted-format) key, so manifest entries can attach fmm-lille's
	 * best-known rank and reference list per format.
	 */
	private static Map<String, JsonNode> loadFmmCatalog(JsonMapper mapper) throws IOException {
		File f = new File("references/catalogs/fmm-lille-catalog.json");
		if (!f.isFile()) return Map.of();
		JsonNode root = mapper.readTree(f);
		JsonNode entries = root.get("entries");
		Map<String, JsonNode> out = new HashMap<>();
		if (entries == null || !entries.isArray()) return out;
		for (JsonNode row : entries) {
			JsonNode fmt = row.get("format");
			if (fmt == null || !fmt.isArray() || fmt.size() != 3) continue;
			int a = fmt.get(0).asInt(), b = fmt.get(1).asInt(), c = fmt.get(2).asInt();
			String key = canonicalKey(a, b, c);
			// Keep the lowest-rank entry per canonical format.
			JsonNode prev = out.get(key);
			if (prev == null || row.get("rank").asInt() < prev.get("rank").asInt()) {
				out.put(key, row);
			}
		}
		return out;
	}

	private static String canonicalKey(int a, int b, int c) {
		int[] sorted = { a, b, c };
		java.util.Arrays.sort(sorted);
		return sorted[0] + "x" + sorted[1] + "x" + sorted[2];
	}
}
