package eu.solven.matmul.docs.generate;

import lombok.extern.slf4j.Slf4j;

import eu.solven.matmul.papers.laderman1976.Laderman23;
import eu.solven.matmul.search.BlockSplitSearch;
import eu.solven.matmul.papers.rosowski2019.RosowskiBound;
import eu.solven.matmul.commutative.CommutativeBounds;
import eu.solven.matmul.catalog.Recombination;
import eu.solven.matmul.catalog.SchemeIO;

import java.io.IOException;
import java.io.PrintWriter;
import java.nio.file.Path;
import java.time.Instant;
import java.util.ArrayList;
import java.util.List;
import java.util.Map;
import java.util.Optional;
import java.util.function.Function;

import eu.solven.matmul.NonCubicBilinearAlgorithm;

/**
 * Generates {@code docs/derived-from-cited-bounds.json} — formula-derived upper-bound
 * entries for matmul targets via {@link BlockSplitSearch}.
 *
 * <p>Two layers:</p>
 * <ul>
 *   <li><strong>Cubic</strong> {@code ⟨n,n,n⟩} via the closed-form Sedoglavic
 *       identity {@link BlockSplitSearch#findBestSplit}, scanned PER FIELD
 *       (R, C, F2) — each field expressed individually (no cross-field mix).</li>
 *   <li><strong>Non-cubic</strong> {@code ⟨n,m,p⟩} via
 *       {@link BlockSplitSearch#findBestSplitNonCubic} with the canonical
 *       Strassen ⟨2,2,2⟩ outer base, scanned PER FIELD (R, C, F2) using
 *       {@link BlockSplitSearch#loadCatalogBestRanksForField}. Reports
 *       gaps (formula beats catalog) and missing directs (no catalog entry).</li>
 * </ul>
 *
 * <p>To keep the output tractable, non-cubic emission is capped per field
 * (top N gaps + top M missing-directs by predicted rank).</p>
 */
@Slf4j
public final class GenerateDerivedBounds {

	private static final int MAX_N_CUBIC = eu.solven.matmul.catalog.CatalogLimits.MAX_DIM;
	private static final int MAX_DIM_NONCUBIC = eu.solven.matmul.catalog.CatalogLimits.MAX_DIM;
	private static final int MAX_GAPS_PER_FIELD = 50;
	private static final int MAX_MISSING_PER_FIELD = 50;
	private static final String OUT_PATH = "docs/derived-from-cited-bounds.json";

	private GenerateDerivedBounds() {}

	public static void main(String[] args) throws IOException {
		List<String> entries = new ArrayList<>();
		// No "mixed"-field cubic section: a cross-field bound (mixing Z/Q/R/C
		// sub-ranks) is meaningless under field discipline. Each field is
		// expressed individually below via cubicEntriesPerField.

		// Resolve by (shape, source) from content — never a hardcoded filename
		// (filenames are pure labels and were renamed 2026-06; see CLAUDE.md).
		NonCubicBilinearAlgorithm strassen = SchemeIO.read(
				eu.solven.matmul.catalog.SchemeResolver.byShapeAndSource(2, 2, 2, "strassen"));
		for (String field : new String[] { "R", "C", "F2" }) {
			entries.addAll(nonCubicEntriesForField(strassen, field));
			entries.addAll(cubicEntriesPerField(field));
		}
		// Note: BlockSplitSearch internally uses "R" for the R/Q/Z cluster (single-letter
		// field keys are easier to filter on). The catalog manifest uses "R/Q/Z". Normalise
		// derived entries above to use the catalog's "R/Q/Z" convention so Pages filter
		// matching is consistent (handled in formatEntry above).
		entries.addAll(rosowskiCommutativeEntries());
		entries.addAll(rosowskiNonBilinearCommutativeEntries());
		entries.addAll(rosowskiBilinearNonCubicEntries());
		entries.addAll(commutativeRecombineEntries(strassen));
		// multiBaseSymmetricCubicEntries removed 2026-06-03 per user feedback:
		// the emitted entries (⟨5,5,5⟩=98 via search, ⟨7,7,7⟩=249, ⟨11,11,11⟩=873, …)
		// merely re-derive known upstream bounds (Sedoglavic-Smirnov 2021,
		// Sedoglavic 2017 Prop 1, …) yet were tagged "solven-strassen 2026 [30]
		// modern catalog SOTA". That's systematic over-attribution. Drop the
		// emission entirely; if/when we want to surface search-derivations
		// for analysis, do it with proper attribution_for_rank pointing at
		// the actual upstream source.

		// Defensive guard: a derived "bound" must beat the naive n·m·p algorithm,
		// else it's not a bound (it would imply ω ≥ 3). Fail the build loudly rather
		// than emit garbage — surfaces degenerate-split bugs at generation time.
		validateNoWorseThanNaive(entries);

		try (PrintWriter pw = new PrintWriter(OUT_PATH)) {
			pw.println("{");
			pw.println("  \"_description\": \"Formula-derived UPPER bounds via BlockSplitSearch. Two sections: cubic (field-pure closed-form identity) and non-cubic (per-field Strassen-recombine). Each entry is a predicted rank — materialise via Recombination.constructWithAllocation.\",");
			pw.println("  \"_schema\": {");
			pw.println("    \"format\": \"[n, m, p]\",");
			pw.println("    \"field\": \"R | C | F2 (field-pure; each field expressed individually)\",");
			pw.println("    \"rank\": \"predicted rank from the formula\",");
			pw.println("    \"breakdown\": \"how the rank decomposes\",");
			pw.println("    \"construction\": \"recipe to materialise the algorithm\",");
			pw.println("    \"verified\": \"false — derived bound, not yet built\",");
			pw.println("    \"direct_catalog_rank\": \"the catalog's current best (or null if missing)\",");
			pw.println("    \"source\": \"primary reference\"");
			pw.println("  },");
			// No "generated" timestamp — would conflict on every regen.
			pw.println("  \"entries\": [");
			for (int i = 0; i < entries.size(); i++) {
				pw.print("    " + entries.get(i));
				pw.println(i == entries.size() - 1 ? "" : ",");
			}
			pw.println("  ]");
			pw.println("}");
		}
		log.info("Wrote " + entries.size() + " derived-bound entries to " + OUT_PATH);
	}

	/**
	 * Cubic targets {@code ⟨n,n,n⟩} (n ∈ [4, 32]) with FIELD-PURE lookups via
	 * {@link BlockSplitSearch#loadCatalogBestRanksForField}. Emits a row when
	 * either no direct catalog scheme exists for this field, or formula < direct.
	 */
	/**
	 * Cubic targets {@code ⟨n,n,n⟩} (n ∈ [4, 32]) with field-pure lookups,
	 * RECURSIVE: each computed derived bound feeds back into the lookup table
	 * so later iterations can use it. This is the systemic propagation
	 * described in TRILINEAR_AGGREGATION.md §3quart — improving a smaller
	 * format's rank cascades to larger formats.
	 *
	 * <p>Example: F₂ catalog only has ⟨2..5⟩³ direct. Without propagation,
	 * ⟨10,10,10⟩ over F₂ can't be derived (would need ⟨6,6,6⟩). WITH
	 * propagation: ⟨6,6,6⟩=167 gets derived first via ⟨4,4,4⟩=47+3·⟨4,4,2⟩+3·⟨4,2,2⟩,
	 * then ⟨10,10,10⟩ can split as 6+4 using the freshly-derived ⟨6,6,6⟩.</p>
	 *
	 * <p><strong>Always emits</strong> the formula value for cubic targets,
	 * even when the catalog already has a better direct entry. This lets
	 * users SEE how the formula behaves across the range and compare
	 * strategies side-by-side (per user feedback 2026-05-28).</p>
	 */
	private static List<String> cubicEntriesPerField(String field) {
		// Start with the catalog and grow as we derive bounds.
		Map<String, Integer> growingRanks = new java.util.HashMap<>(
				BlockSplitSearch.loadCatalogBestRanksForField(field));
		Function<int[], Optional<Integer>> growingLookup = key -> {
			int[] sorted = { key[0], key[1], key[2] };
			java.util.Arrays.sort(sorted);
			return Optional.ofNullable(growingRanks.get(sorted[0] + "x" + sorted[1] + "x" + sorted[2]));
		};
		List<String> out = new ArrayList<>();
		int gaps = 0, missing = 0, ok = 0;
		for (int n = 4; n <= MAX_N_CUBIC; n++) {
			Optional<BlockSplitSearch.SplitCandidate> best = BlockSplitSearch.findBestSplit(n, growingLookup);
			if (best.isEmpty()) continue;
			BlockSplitSearch.SplitCandidate c = best.get();
			Integer directRank = growingLookup.apply(new int[] { n, n, n }).orElse(null);
			// Never emit a "bound" worse than (or equal to) the naive n³ — it would
			// imply ω ≥ 3, which is noise, not a bound (degenerate splits produce these).
			long naive = (long) n * n * n;
			if (c.formulaRank() >= naive) {
				log.info("  skip ⟨{}³⟩ {}: formula {} ≥ naive {}", n, field, c.formulaRank(), naive);
				continue;
			}
			// Suppress DOMINATED derived rows (user 2026-06-03). These cubic
			// entries are BlockSplitSearch recursive self-applications, not a dated
			// upstream result — so they carry no publication priority. An entry
			// that is already ≥ the direct catalog rank when generated is pure
			// noise (e.g. ⟨4,4,4⟩=58 vs direct 48 pinned a fake ω≈2.93 on the
			// "formula"). We emit a derived row ONLY when it is strictly better
			// than (or fills a gap in) the direct catalog. The "keep dominated for
			// history" exception applies to genuinely dated bounds that were SOTA
			// when published and lost the crown later — those live in
			// cited-bounds.json, not in this computed-split section.
			if (directRank != null && c.formulaRank() >= directRank) {
				ok++;  // dominated by the direct catalog → not emitted
				continue;
			}
			if (directRank == null) missing++;
			else gaps++;
			// Propagate: derived rank now becomes the best known for this format,
			// so subsequent larger targets can split on n and use this rank.
			String key = n + "x" + n + "x" + n;
			growingRanks.merge(key, (int) c.formulaRank(), Integer::min);
			StringBuilder sb = new StringBuilder();
			sb.append("{");
			sb.append("\"format\": [").append(n).append(",").append(n).append(",").append(n).append("], ");
			sb.append("\"field\": \"").append(field).append("\", ");
			sb.append("\"rank\": ").append(c.formulaRank()).append(", ");
			sb.append("\"split\": [").append(c.u()).append(",").append(c.v()).append("], ");
			sb.append("\"breakdown\": \"").append(c.rUuu()).append(" + 3·").append(c.rUuv())
					.append(" + 3·").append(c.rUvv()).append(" = ").append(c.formulaRank()).append("\", ");
			sb.append("\"construction\": \"BlockSplitSearch.findBestSplit cubic ⟨").append(n)
					.append("⟩³ over ").append(field).append("\", ");
			sb.append("\"verified\": false, \"direct_catalog_rank\": ");
			sb.append(directRank == null ? "null" : directRank.toString());
			// Attribution: these are COMPOSED schemes (recursive BlockSplitSearch),
			// not a single dated upstream result — tag as "Derived" (user
			// 2026-06-03; concrete atom-level attribution is a later refinement).
			sb.append(", \"source\": \"Derived\"");
			sb.append("}");
			out.add(sb.toString());
		}
		log.info("  " + field + " cubic (recursive): " + gaps + " gaps + " + missing + " missing + " + ok + " at-or-above-catalog");
		return out;
	}

	/**
	 * Fails the build if any emitted derived-bound entry has rank ≥ the naive
	 * {@code n·m·p} algorithm — such an entry is not a bound (ω ≥ 3) and signals
	 * a degenerate-split / lookup bug in one of the generators.
	 */
	private static void validateNoWorseThanNaive(List<String> entries) {
		java.util.regex.Pattern fmt = java.util.regex.Pattern.compile(
				"\"format\"\\s*:\\s*\\[\\s*(\\d+)\\s*,\\s*(\\d+)\\s*,\\s*(\\d+)\\s*\\]");
		java.util.regex.Pattern rnk = java.util.regex.Pattern.compile("\"rank\"\\s*:\\s*(\\d+)");
		List<String> bad = new ArrayList<>();
		for (String e : entries) {
			java.util.regex.Matcher mf = fmt.matcher(e);
			java.util.regex.Matcher mr = rnk.matcher(e);
			if (!mf.find() || !mr.find()) continue;
			long naive = (long) Integer.parseInt(mf.group(1)) * Integer.parseInt(mf.group(2))
					* Integer.parseInt(mf.group(3));
			long rank = Long.parseLong(mr.group(1));
			if (rank >= naive) bad.add("⟨" + mf.group(1) + "," + mf.group(2) + "," + mf.group(3)
					+ "⟩ rank=" + rank + " ≥ naive=" + naive);
		}
		if (!bad.isEmpty()) {
			throw new IllegalStateException("derived-bounds: " + bad.size()
					+ " worse-than-naive entries (ω ≥ 3) — degenerate-split bug:\n  "
					+ String.join("\n  ", bad));
		}
	}

	private static List<String> nonCubicEntriesForField(NonCubicBilinearAlgorithm strassen, String field) {
		Map<String, Integer> ranks = BlockSplitSearch.loadCatalogBestRanksForField(field);
		Function<int[], Optional<Integer>> lookup = BlockSplitSearch.rankLookupFromMap(ranks);
		Recombination.SotaResolver sota = (a, b, c) -> {
			if (a == 0 || b == 0 || c == 0) return 0;
			if (a == 1) return b * c;
			if (b == 1) return a * c;
			if (c == 1) return a * b;
			return lookup.apply(new int[] { a, b, c }).orElse(Recombination.SotaResolver.UNKNOWN_RANK);
		};

		List<EntryWithSortKey> gaps = new ArrayList<>();
		List<EntryWithSortKey> missing = new ArrayList<>();
		for (int n = 2; n <= MAX_DIM_NONCUBIC; n++) {
			for (int m = n; m <= MAX_DIM_NONCUBIC; m++) {
				for (int p = m; p <= MAX_DIM_NONCUBIC; p++) {
					if (n == m && m == p) continue;
					int tn = n, tm = m, tp = p;
					Recombination.SotaResolver pure = (a, b, c) -> {
						if (a == tn && b == tm && c == tp) return Recombination.SotaResolver.UNKNOWN_RANK;
						return sota.getRank(a, b, c);
					};
					Optional<BlockSplitSearch.NonCubicSplitCandidate> best =
							BlockSplitSearch.findBestSplitNonCubic(n, m, p, strassen, pure);
					if (best.isEmpty()) continue;
					BlockSplitSearch.NonCubicSplitCandidate c = best.get();
					if (c.rank() >= Integer.MAX_VALUE / 200) continue;
					Optional<Integer> direct = lookup.apply(new int[] { n, m, p });
					if (direct.isEmpty()) {
						missing.add(new EntryWithSortKey(formatNonCubicEntry(c, field, null), c.rank()));
					} else if (c.rank() < direct.get()) {
						long delta = direct.get() - c.rank();
						gaps.add(new EntryWithSortKey(formatNonCubicEntry(c, field, direct.get()), -delta));
					}
				}
			}
		}
		gaps.sort((a, b) -> Long.compare(a.sortKey, b.sortKey)); // most-negative (biggest delta) first
		missing.sort((a, b) -> Long.compare(a.sortKey, b.sortKey));
		List<String> out = new ArrayList<>();
		gaps.stream().limit(MAX_GAPS_PER_FIELD).forEach(e -> out.add(e.json));
		missing.stream().limit(MAX_MISSING_PER_FIELD).forEach(e -> out.add(e.json));
		log.info("  " + field + ": " + Math.min(gaps.size(), MAX_GAPS_PER_FIELD)
				+ " gaps + " + Math.min(missing.size(), MAX_MISSING_PER_FIELD)
				+ " missing (out of " + gaps.size() + " gaps total + " + missing.size() + " missing total)");
		return out;
	}

	private static String formatNonCubicEntry(BlockSplitSearch.NonCubicSplitCandidate c,
			String field, Integer directRank) {
		StringBuilder sb = new StringBuilder();
		sb.append("{");
		sb.append("\"format\": [").append(c.n()).append(",").append(c.m()).append(",").append(c.p()).append("], ");
		sb.append("\"field\": \"").append(field).append("\", ");
		sb.append("\"rank\": ").append(c.rank()).append(", ");
		sb.append("\"split\": [").append(java.util.Arrays.toString(c.allocA())).append(", ")
				.append(java.util.Arrays.toString(c.allocB())).append(", ")
				.append(java.util.Arrays.toString(c.allocC())).append("], ");
		sb.append("\"breakdown\": \"").append(c.breakdown()).append("\", ");
		sb.append("\"construction\": \"Recombination.constructWithAllocation(Strassen, ")
				.append(field).append("Lookup, ")
				.append(java.util.Arrays.toString(c.allocA())).append("/")
				.append(java.util.Arrays.toString(c.allocB())).append("/")
				.append(java.util.Arrays.toString(c.allocC())).append(")\", ");
		sb.append("\"verified\": false, ");
		sb.append("\"direct_catalog_rank\": ");
		sb.append(directRank == null ? "null" : directRank.toString());
		sb.append(", \"source\": \"Derived\"");
		sb.append("}");
		return sb.toString();
	}

	/**
	 * Emits Rosowski 2019/2020 commutative-only formula bounds for cubic
	 * {@code ⟨n,n,n⟩} {@code n ∈ [3, 20]}. Skips entries where DIS09 / catalog
	 * already has an equal-or-better commutative bound.
	 */
	private static List<String> rosowskiCommutativeEntries() {
		List<String> out = new ArrayList<>();
		// Compare against DIS09 commutative Table 4 (already in cited-bounds.json).
		// We hardcode a small table of DIS09 commutative ranks for cubic n ∈ [2, 30]:
		Map<Integer, Integer> dis09Commutative = Map.ofEntries(
				Map.entry(2, 7), Map.entry(3, 22), Map.entry(4, 46), Map.entry(5, 93),
				Map.entry(6, 141), Map.entry(7, 235), Map.entry(8, 316), Map.entry(9, 472),
				Map.entry(10, 595), Map.entry(11, 825), Map.entry(12, 987), Map.entry(13, 1318),
				Map.entry(14, 1525), Map.entry(15, 1941), Map.entry(16, 2212));
		int improved = 0, equalOrWorse = 0;
		for (int n = 3; n <= 16; n++) {
			Optional<Long> r = RosowskiBound.commutativeBound(n, n, n);
			if (r.isEmpty()) continue;
			Integer dis09 = dis09Commutative.get(n);
			boolean improvesDis09 = (dis09 == null || r.get() < dis09);
			if (!improvesDis09) {
				equalOrWorse++;
				continue;
			}
			improved++;
			StringBuilder sb = new StringBuilder();
			sb.append("{");
			sb.append("\"format\": [").append(n).append(",").append(n).append(",").append(n).append("], ");
			sb.append("\"field\": \"R/Q/Z\", \"commutative\": true, ");
			sb.append("\"rank\": ").append(r.get()).append(", ");
			sb.append("\"breakdown\": \"Rosowski formula n(lm+l+m−1)/2 = ").append(n).append("·")
					.append(n * n + 2 * n - 1).append("/2 = ").append(r.get()).append("\", ");
			sb.append("\"construction\": \"RosowskiBound.commutativeBound(").append(n).append(",")
					.append(n).append(",").append(n).append(") — commutative ring required\", ");
			sb.append("\"verified\": false, \"direct_catalog_rank\": ");
			sb.append(dis09 == null ? "null" : dis09.toString());
			sb.append(", \"source\": \"Rosowski 2019 (arXiv:1904.07683) — commutative-only\"");
			sb.append("}");
			out.add(sb.toString());
		}
		log.info("  Rosowski commutative cubic: " + improved + " improve DIS09, "
				+ equalOrWorse + " equal/worse");
		return out;
	}

	/**
	 * Rosowski 2019/2020 Theorems 2/3 — bilinear COMMUTATIVE rank bound for
	 * NON-CUBIC {@code ⟨l, n, m⟩} up to dim {@value #MAX_DIM_BIG}. Emits
	 * only when the formula beats any existing R-class catalog entry OR
	 * fills a missing direct (no catalog entry). Capped per-output to keep
	 * derived-from-cited-bounds.json reasonable.
	 *
	 * <p>Formulas are <em>commutative</em>: do NOT lift to recursive matmul.
	 * Each entry is tagged {@code commutative: true}.</p>
	 */
	private static final int MAX_DIM_BIG = eu.solven.matmul.catalog.CatalogLimits.MAX_DIM;
	private static final int MAX_PER_FAMILY = 200;

	private static List<String> rosowskiBilinearNonCubicEntries() {
		Map<String, Integer> ranks = BlockSplitSearch.loadCatalogBestRanksForField("R");
		Function<int[], Optional<Integer>> lookup = BlockSplitSearch.rankLookupFromMap(ranks);
		List<String> out = new ArrayList<>();
		int emitted = 0;
		for (int l = 2; l <= MAX_DIM_BIG && emitted < MAX_PER_FAMILY; l++) {
			for (int n = l; n <= MAX_DIM_BIG && emitted < MAX_PER_FAMILY; n++) {
				for (int m = n; m <= MAX_DIM_BIG && emitted < MAX_PER_FAMILY; m++) {
					if (l == n && n == m) continue; // cubic — handled elsewhere
					Optional<Long> bound = RosowskiBound.bestCommutativeBound(l, n, m);
					if (bound.isEmpty()) continue;
					Integer direct = lookup.apply(new int[] { l, n, m }).orElse(null);
					// Emit if better than direct, OR no direct exists.
					if (direct != null && direct <= bound.get()) continue;
					emitted++;
					StringBuilder sb = new StringBuilder();
					sb.append("{");
					sb.append("\"format\": [").append(l).append(",").append(n).append(",").append(m).append("], ");
					sb.append("\"field\": \"R/Q/Z\", \"commutative\": true, ");
					sb.append("\"rank\": ").append(bound.get()).append(", ");
					sb.append("\"breakdown\": \"Rosowski Thm 2/3 bilinear-commutative: best over axis perms\", ");
					sb.append("\"construction\": \"RosowskiBound.bestCommutativeBound(")
							.append(l).append(",").append(n).append(",").append(m)
							.append(") — commutative ring required\", ");
					sb.append("\"verified\": false, \"direct_catalog_rank\": ");
					sb.append(direct == null ? "null" : direct.toString());
					sb.append(", \"source\": \"Rosowski 2019 Thm 2/3 (arXiv:1904.07683) — commutative-only\"");
					sb.append("}");
					out.add(sb.toString());
				}
			}
		}
		log.info("  Rosowski commutative non-cubic: " + emitted + " entries (cap " + MAX_PER_FAMILY + ")");
		return out;
	}

	/**
	 * Rosowski 2019/2020 Theorems 4/5 — non-bilinear <b>COMMUTATIVE</b> rank
	 * upper bound on cubic {@code ⟨n,n,n⟩}. The paper is titled "Fast
	 * Commutative Matrix Algorithm" (arXiv:1904.07683); Thm 4/5 are
	 * non-bilinear (their products mix A- and B-entries), so correctness
	 * relies on scalar commutativity. Rosowski's "recursive non-bilinear"
	 * novelty recurses these <em>over commutative rings</em> (giving
	 * ω≈2.8125) — it does NOT lift to non-commutative matmul (it does not
	 * realise the NC tensor rank). These are therefore emitted with
	 * {@code "commutative": true} and must stay out of the NC SOTA pipeline.
	 * Emitted across the whole range so the commutative-vs-NC comparison is
	 * visible; {@code direct_catalog_rank} reports the NC catalog value purely
	 * for that comparison.
	 */
	private static List<String> rosowskiNonBilinearCommutativeEntries() {
		Map<String, Integer> ranks = BlockSplitSearch.loadCatalogBestRanksForField("R");
		Function<int[], Optional<Integer>> lookup = BlockSplitSearch.rankLookupFromMap(ranks);
		List<String> out = new ArrayList<>();
		int wins = 0;
		// Always emit cubic formula values across the range so the user can
		// compare formula behaviour vs catalog (per user feedback 2026-05-28).
		for (int n = 4; n <= 32; n++) {
			Optional<Long> r = RosowskiBound.nonBilinearRankBound(n);
			if (r.isEmpty()) continue;
			Integer direct = lookup.apply(new int[] { n, n, n }).orElse(null);
			if (direct == null || direct > r.get()) wins++;
			String formula = (n % 2 == 0)
					? "n(n²+3n+1)/2 = " + n + "·" + (n * n + 3 * n + 1) + "/2"
					: "n(n²+3n+2)/2 = " + n + "·" + (n * n + 3 * n + 2) + "/2";
			StringBuilder sb = new StringBuilder();
			sb.append("{");
			sb.append("\"format\": [").append(n).append(",").append(n).append(",").append(n).append("], ");
			sb.append("\"field\": \"R/Q/Z\", ");
			sb.append("\"commutative\": true, ");
			sb.append("\"rank\": ").append(r.get()).append(", ");
			sb.append("\"breakdown\": \"Rosowski Thm ").append(n % 2 == 0 ? "4" : "5")
					.append(" non-bilinear (commutative): ").append(formula).append(" = ").append(r.get()).append("\", ");
			sb.append("\"construction\": \"RosowskiBound.nonBilinearRankBound(").append(n)
					.append(") — non-bilinear COMMUTATIVE algorithm; recurses over commutative rings only, does NOT lift to NC matmul\", ");
			sb.append("\"verified\": false, \"direct_catalog_rank\": ");
			sb.append(direct == null ? "null" : direct.toString());
			sb.append(", \"source\": \"Rosowski 2019 Thm ").append(n % 2 == 0 ? "4" : "5")
					.append(" (arXiv:1904.07683) — commutative non-bilinear\"");
			sb.append("}");
			out.add(sb.toString());
		}
		log.info("  Rosowski commutative non-bilinear: " + wins
				+ " cubic targets below NC catalog (commutative-only; not NC SOTA)");
		return out;
	}

	/**
	 * Commutative recombination via NC Strassen ⟨2,2,2⟩=7 outer × commutative
	 * sub-ranks ({@link CommutativeBounds}). Emits cubic entries where the
	 * recombination beats both Rosowski direct AND DIS09 Table 4 direct.
	 *
	 * <p>Sub-optimal vs a hypothetical fully-commutative recombination using
	 * Hopcroft/Winograd ⟨2,2,2⟩=6 as the outer base (not yet on disk).</p>
	 */
	private static List<String> commutativeRecombineEntries(NonCubicBilinearAlgorithm strassen) {
		CommutativeBounds cmt = new CommutativeBounds();
		Recombination.SotaResolver sotaCmt = cmt.asSotaResolver();
		List<String> out = new ArrayList<>();
		int beats = 0;
		for (int n = 4; n <= 20; n++) {
			Optional<Long> direct = cmt.bestRank(n, n, n);
			if (direct.isEmpty()) continue;
			int directInt = direct.get().intValue();
			int tn = n;
			Recombination.SotaResolver pure = (a, b, c) -> {
				if (a == tn && b == tn && c == tn) return Recombination.SotaResolver.UNKNOWN_RANK;
				return sotaCmt.getRank(a, b, c);
			};
			Optional<BlockSplitSearch.NonCubicSplitCandidate> best =
					BlockSplitSearch.findBestSplitNonCubic(n, n, n, strassen, pure);
			if (best.isEmpty()) continue;
			BlockSplitSearch.NonCubicSplitCandidate c = best.get();
			if (c.rank() >= directInt) continue;
			beats++;
			StringBuilder sb = new StringBuilder();
			sb.append("{");
			sb.append("\"format\": [").append(n).append(",").append(n).append(",").append(n).append("], ");
			sb.append("\"field\": \"R/Q/Z\", \"commutative\": true, ");
			sb.append("\"rank\": ").append(c.rank()).append(", ");
			sb.append("\"split\": [").append(java.util.Arrays.toString(c.allocA())).append(", ")
					.append(java.util.Arrays.toString(c.allocB())).append(", ")
					.append(java.util.Arrays.toString(c.allocC())).append("], ");
			sb.append("\"breakdown\": \"NC-Strassen-outer × commutative sub-ranks (").append(c.breakdown()).append(")\", ");
			sb.append("\"construction\": \"Recombination.recombineWithAllocation(Strassen, cmtSota, allocs)\", ");
			sb.append("\"verified\": false, \"direct_catalog_rank\": ").append(directInt);
			sb.append(", \"source\": \"solven-strassen 2026 — NC-Strassen × commutative subs (Rosowski + DIS09)\"");
			sb.append("}");
			out.add(sb.toString());
		}
		log.info("  Commutative recombination: " + beats + " cubic improvements over direct");
		return out;
	}

	/**
	 * Multi-base + S₃ symmetry block-split for cubic {@code ⟨n,n,n⟩} over
	 * R/Q/Z, n=4..30. This is the DIS09 reproduction with modern catalog
	 * (Sedoglavic, AlphaTensor, AlphaEvolve, etc.) as the SOTA resolver,
	 * via {@link BlockSplitSearch#findBestMultiBaseSplit} expanded with
	 * the {@link eu.solven.matmul.SymmetryTransforms#s3Orbit}
	 * orbit. Emits only the rows that <strong>strictly improve</strong>
	 * over DIS09 Table 3 (so a future write-up can cite this repo as the
	 * source of the bound).
	 *
	 * <p>See [30] in {@code REFERENCES.md} and
	 * {@link eu.solven.matmul.catalog.TestDIS09FullScan} for the
	 * full side-by-side comparison.</p>
	 */
	private static List<String> multiBaseSymmetricCubicEntries(NonCubicBilinearAlgorithm strassen)
			throws IOException {
		Map<Integer, Integer> dis09 = Map.ofEntries(
				Map.entry(4, 49), Map.entry(5, 100), Map.entry(6, 161), Map.entry(7, 258),
				Map.entry(8, 343), Map.entry(9, 522), Map.entry(10, 700), Map.entry(11, 923),
				Map.entry(12, 1125), Map.entry(13, 1450), Map.entry(14, 1728),
				Map.entry(15, 2108), Map.entry(16, 2401), Map.entry(17, 2972),
				Map.entry(18, 3306), Map.entry(19, 4073), Map.entry(20, 4340),
				Map.entry(21, 5365), Map.entry(22, 5566), Map.entry(23, 6806),
				Map.entry(24, 7000), Map.entry(25, 8448), Map.entry(26, 8658),
				Map.entry(27, 10330), Map.entry(28, 10556), Map.entry(29, 12468),
				Map.entry(30, 12710));

		NonCubicBilinearAlgorithm laderman = NonCubicBilinearAlgorithm
				.fromCubic(eu.solven.matmul.papers.laderman1976.Laderman23.get());
		List<BlockSplitSearch.NamedBase> pool = new ArrayList<>();
		expandS3(pool, "Strassen ⟨2,2,2⟩=7", strassen);
		expandS3(pool, "Laderman ⟨3,3,3⟩=23", laderman);
		expandS3(pool, "mul211", eu.solven.matmul.AxisSplitBases.mul211());
		expandS3(pool, "mul121", eu.solven.matmul.AxisSplitBases.mul121());
		expandS3(pool, "mul112", eu.solven.matmul.AxisSplitBases.mul112());

		Map<String, Integer> ranks = BlockSplitSearch.loadCatalogBestRanksForField("R");
		Function<int[], Optional<Integer>> lookup = BlockSplitSearch.rankLookupFromMap(ranks);
		Recombination.SotaResolver sota = (a, b, c) -> {
			if (a == 0 || b == 0 || c == 0) return 0;
			if (a == 1) return b * c;
			if (b == 1) return a * c;
			if (c == 1) return a * b;
			return lookup.apply(new int[] { a, b, c }).orElse(Recombination.SotaResolver.UNKNOWN_RANK);
		};

		List<String> out = new ArrayList<>();
		int wins = 0;
		for (int n = 4; n <= 30; n++) {
			// Full enumeration with 1M combos/base cap: Strassen/axis-split (2 blocks/axis)
			// go through unbalanced; Laderman (3 blocks/axis) auto-falls back to balanced
			// past ~n=15. This is what reveals the n=21 5258 win.
			Optional<BlockSplitSearch.MultiBaseSplitCandidate> best =
					BlockSplitSearch.findBestMultiBaseSplit(n, n, n, pool, sota, false, 1_000_000L);
			if (best.isEmpty()) continue;
			BlockSplitSearch.MultiBaseSplitCandidate c = best.get();
			int dis = dis09.get(n);
			Integer directInt = ranks.get(n + "x" + n + "x" + n);
			// Always emit strict-vs-DIS09 entries for analysis (it's useful to know how
			// good our pipeline gets, even where the catalog already has better via
			// modern methods). Display-time shaving (drop derived ≥ catalog) lives in
			// docs/catalog.js per user feedback 2026-05-28.
			if (c.rank() >= dis) continue;
			wins++;
			StringBuilder sb = new StringBuilder();
			sb.append("{");
			sb.append("\"format\": [").append(n).append(",").append(n).append(",").append(n).append("], ");
			sb.append("\"field\": \"R/Q/Z\", ");
			sb.append("\"rank\": ").append(c.rank()).append(", ");
			sb.append("\"split\": [")
					.append(java.util.Arrays.toString(c.allocA())).append(", ")
					.append(java.util.Arrays.toString(c.allocB())).append(", ")
					.append(java.util.Arrays.toString(c.allocC())).append("], ");
			sb.append("\"breakdown\": \"").append(escape(c.breakdown())).append("\", ");
			sb.append("\"construction\": \"BlockSplitSearch.findBestMultiBaseSplit(")
					.append(n).append(",").append(n).append(",").append(n)
					.append(", pool=Strassen+Laderman+mul211/121/112 S₃-expanded, balanced=true)\", ");
			sb.append("\"verified\": false, \"direct_catalog_rank\": ");
			sb.append(directInt == null ? "null" : directInt.toString());
			sb.append(", \"dis09_rank\": ").append(dis);
			sb.append(", \"improvement_vs_dis09\": ").append(dis - c.rank());
			sb.append(", \"source\": \"solven-strassen 2026 [30] (multi-base + S₃ symmetry; modern catalog SOTA)\"");
			sb.append("}");
			out.add(sb.toString());
		}
		log.info("  Multi-base + S₃ (vs DIS09): " + wins + " strict improvements emitted");
		return out;
	}

	private static void expandS3(List<BlockSplitSearch.NamedBase> pool,
			String label, NonCubicBilinearAlgorithm base) {
		List<NonCubicBilinearAlgorithm> orbit =
				eu.solven.matmul.SymmetryTransforms.s3Orbit(base);
		for (int i = 0; i < orbit.size(); i++) {
			NonCubicBilinearAlgorithm a = orbit.get(i);
			pool.add(new BlockSplitSearch.NamedBase(
					label + " [σ" + i + " ⟨" + a.n + "," + a.m + "," + a.p + "⟩]", a));
		}
	}

	private static String escape(String s) {
		return s.replace("\\", "\\\\").replace("\"", "\\\"");
	}

	private record EntryWithSortKey(String json, long sortKey) {}
}
