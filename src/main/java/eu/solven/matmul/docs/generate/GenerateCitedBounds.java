package eu.solven.matmul.docs.generate;

import lombok.extern.slf4j.Slf4j;

import eu.solven.matmul.papers.rosowski2019.RosowskiBound;

import java.io.IOException;
import java.io.PrintWriter;
import java.nio.file.Files;
import java.nio.file.Path;
import java.time.Instant;
import java.util.ArrayList;
import java.util.List;
import java.util.Map;

import tools.jackson.databind.JsonNode;
import tools.jackson.databind.json.JsonMapper;

/**
 * Generates {@code docs/cited-bounds.json} from literature data sources.
 *
 * <p>"Cited bounds" are upper-bound RANK CLAIMS published in the literature
 * WITHOUT (or with hard-to-extract) accompanying factor matrices. They
 * complement:</p>
 * <ul>
 *   <li>{@code catalog.json}: schemes with verified factor matrices on disk.</li>
 *   <li>{@code lower-bounds.json}: lower-bound claims (impossibility).</li>
 *   <li>{@code derived-from-cited-bounds.json}: formula-derived upper bounds (computable).</li>
 *   <li>{@code cited-bounds.json} <i>(this file)</i>: published claims we
 *       haven't materialised as schemes. Often superseded by later work that
 *       does provide schemes.</li>
 * </ul>
 *
 * <p>Sources mined here:</p>
 * <ul>
 *   <li>{@code references/dis09-cubic-tables.json} — DIS09 Tables 3 & 4
 *       (Drevet–Islam–Schost 2009): cubic ranks n ∈ [2, 30] in both
 *       non-commutative and commutative cases.</li>
 * </ul>
 */
@Slf4j
public final class GenerateCitedBounds {

	private static final String OUT_PATH = "docs/cited-bounds.json";

	public static void main(String[] args) throws IOException {
		List<String> entries = new ArrayList<>();
		entries.addAll(extractFromDis09());
		entries.addAll(extractPerminovSerendipitous());

		try (PrintWriter pw = new PrintWriter(OUT_PATH)) {
			pw.println("{");
			pw.println("  \"_description\": \"Upper-bound RANK CLAIMS from the literature without accompanying factor matrices. See GenerateCitedBounds.java for sources.\",");
			pw.println("  \"_schema\": {");
			pw.println("    \"format\": \"[n, m, p]\",");
			pw.println("    \"field\": \"R | C | F2 | all — claim applies in this field-class\",");
			pw.println("    \"commutative\": \"true if the bound only holds under commuting scalars (not useful for recursive matmul); false otherwise\",");
			pw.println("    \"rank\": \"upper bound on the bilinear rank\",");
			pw.println("    \"algorithm\": \"the technique/family that achieves this rank in the source\",");
			pw.println("    \"source\": \"citation key (DIS 2009, Strassen 1969, Laderman 1976, ...)\",");
			pw.println("    \"year\": \"publication year\",");
			pw.println("    \"scheme_provided\": \"whether the source provides explicit factor matrices (false = rank claim only)\",");
			pw.println("    \"source_paper_url\": \"OPTIONAL — the publication the rank claim comes from\",");
			pw.println("    \"source_scheme_url\": \"OPTIONAL — the scheme's file in the source author's own repo, too large to copy here; cite-don't-copy (e.g. Perminov's repo)\",");
			pw.println("    \"notes\": \"context, prior baselines, supersession notes\"");
			pw.println("  },");
			// No "generated" timestamp — would conflict on every regen; the
			// content reflects the source data which is what callers care about.
			pw.println("  \"entries\": [");
			for (int i = 0; i < entries.size(); i++) {
				pw.print("    " + entries.get(i));
				pw.println(i == entries.size() - 1 ? "" : ",");
			}
			pw.println("  ]");
			pw.println("}");
		}
		log.info("Wrote " + entries.size() + " cited-bound entries to " + OUT_PATH);
	}

	/**
	 * Merges DIS09 Table 3 (non-commutative) and Table 4 (commutative). Each format
	 * gets a SINGLE entry per source algorithm. If the algorithm's NC rank and
	 * commutative rank differ, both are recorded ({@code rank}, {@code commutative_rank}).
	 * If they're equal, the commutative_rank field is omitted (no redundant signal).
	 * Source attribution uses the ORIGINAL algorithm's author (e.g., "Strassen 1969",
	 * "Laderman 1976"), not the DIS09 survey.
	 */
	private static List<String> extractFromDis09() throws IOException {
		Path src = Path.of("references/dis09-cubic-tables.json");
		JsonMapper mapper = JsonMapper.builder().build();
		JsonNode root = mapper.readTree(Files.readString(src));

		// Per-format NC and commutative ranks keyed by "n,m,p|algorithm" → rank.
		// We merge into per-format rows when (NC algorithm) == (commutative algorithm).
		Map<String, RankRow> ncByFormat = new java.util.LinkedHashMap<>();
		JsonNode nc = root.get("non_commutative");
		if (nc != null && nc.isArray()) {
			for (JsonNode row : nc) {
				int[] fmt = new int[3];
				for (int i = 0; i < 3; i++) fmt[i] = row.get("format").get(i).asInt();
				String algorithm = row.get("algorithm").asText();
				int rank = row.get("rank").asInt();
				boolean newInDis09 = row.has("new_in_dis09") && row.get("new_in_dis09").asBoolean();
				RankRow rr = new RankRow(fmt, algorithm, rank, null, "");
				rr.newInDis09 = newInDis09;
				ncByFormat.put(fmtKey(fmt), rr);
			}
		}
		JsonNode cmt = root.get("commutative");
		if (cmt != null && cmt.isArray()) {
			for (JsonNode row : cmt) {
				int[] fmt = new int[3];
				for (int i = 0; i < 3; i++) fmt[i] = row.get("format").get(i).asInt();
				String algorithm = row.get("algorithm").asText();
				int rank = row.get("rank").asInt();
				String key = fmtKey(fmt);
				RankRow ncRow = ncByFormat.get(key);
				if (ncRow != null && ncRow.algorithm.equals(algorithm) && ncRow.rank == rank) {
					// Same algorithm and same rank in both tables — just one row, no commutative annotation.
					continue;
				}
				if (ncRow != null && ncRow.algorithm.equals(algorithm) && ncRow.rank != rank) {
					// Same algorithm, but different rank when commutativity is allowed.
					ncRow.commutativeRank = rank;
					continue;
				}
				// Distinct commutative-only algorithm (e.g., Makarov333 for ⟨3,3,3⟩=22).
				ncByFormat.put(key + "@" + algorithm,
						new RankRow(fmt, algorithm, rank, null, "commutative-only — does NOT lift to recursive matmul"));
				// Mark these as commutative-only via the rank being commutative.
				ncByFormat.get(key + "@" + algorithm).isCommutativeOnly = true;
			}
		}

		List<String> out = new ArrayList<>();
		for (RankRow r : ncByFormat.values()) out.add(r.toJson());

		// Also emit historical baselines (Probert-Fischer 1980, Smith 2002)
		// whose ranks are tabulated in dis09-cubic-tables.json. These fill
		// the 1976-2009 cubic gap that DIS09 itself improved on. Only emit
		// when their rank STRICTLY differs from the DIS09 rank (no point
		// duplicating DIS09 baseline value).
		if (nc != null && nc.isArray()) {
			for (JsonNode row : nc) {
				int[] fmt = new int[3];
				for (int i = 0; i < 3; i++) fmt[i] = row.get("format").get(i).asInt();
				int dis09Rank = row.get("rank").asInt();
				if (row.has("probert_fischer_1980")
						&& !row.get("probert_fischer_1980").isNull()) {
					int pfRank = row.get("probert_fischer_1980").asInt();
					if (pfRank != dis09Rank) {
						out.add(simpleCitedJson(fmt, pfRank, "Probert-Fischer 1980", 1980,
								"non-commutative cubic baseline (cited by DIS09 Table 3)", false));
					}
				}
				if (row.has("smith_2002") && !row.get("smith_2002").isNull()) {
					int smRank = row.get("smith_2002").asInt();
					if (smRank != dis09Rank) {
						out.add(simpleCitedJson(fmt, smRank, "Smith 2002", 2002,
								"non-commutative cubic baseline (cited by DIS09 Table 3)", false));
					}
				}
			}
		}

		// Schachtel 1978 — ⟨5,5,5⟩=103, non-commutative, predates DIS09.
		// Per user request 2026-05-29: register even though scheme not on disk.
		out.add(simpleCitedJson(new int[]{5, 5, 5}, 103, "Schachtel 1978", 1978,
				"non-commutative ⟨5,5,5⟩; Information Processing Letters 7 (1978) 180–182", false));

		// Sedoglavic-Smirnov 2021 — ⟨5,5,5⟩ rank ≤ 98, border rank ≤ 89.
		// Cited as reference [6] in Smirnov 2017 (⟨3,P,Q⟩ paper); paper title:
		// "The tensor rank of 5x5 matrices multiplication is bounded by 98
		// and its border rank by 89".
		// The exact-rank bound 98 is now dominated by AlphaEvolve 93 (2025)
		// and AlphaTensor 96 (2022), but border rank ≤ 89 remains the SOTA
		// border-rank upper bound for ⟨5,5,5⟩ and is registered separately.
		// Per user 2026-06-03.
		out.add(simpleCitedJson(new int[]{5, 5, 5}, 98, "Sedoglavic-Smirnov 2021", 2021,
				"non-commutative tensor rank ≤ 98; companion border rank ≤ 89 surfaced as a separate cited row (kind=border)",
				false));
		// Border-rank companion. We surface border-rank claims as cited
		// entries tagged kind=border so they don't compete with exact-rank
		// rows in the shave; SPA can display them under a "border rank" lane.
		out.add(borderCitedJson(new int[]{5, 5, 5}, 89, "Sedoglavic-Smirnov 2021", 2021,
				"R̃(⟨5,5,5⟩) ≤ 89 — best published border-rank upper bound at this shape"));

		// Johnson-McLoughlin 1986 — ⟨3,3,3⟩=23, non-commutative.
		// Independent re-derivation of Laderman 1976's bound (NOT a discovery).
		// SIAM J. Comput. 15(2) 1986: 595–603. Per "Distinguish discoveries
		// from re-discoveries" rule (CLAUDE.md): we register it for the
		// historical record but the rank is attributed to Laderman 1976.
		out.add(simpleCitedJson(new int[]{3, 3, 3}, 23, "Johnson-McLoughlin 1986", 1986,
				"non-commutative ⟨3,3,3⟩=23; re-derivation, NOT a discovery — rank established by Laderman 1976; SIAM J. Comput. 15(2) 1986: 595–603",
				false));

		// Schwartz-Zwecher 2025 — Table 1 of arXiv:2508.01748.
		// "TA-New25" family: trilinear aggregation + kin-row unification.
		// Non-commutative, over Q. For n ∈ {20, 22, 24, 26, 28, 30, 32} we
		// also ship the explicit factor matrices under section{n}/, so
		// hasSchemeFile auto-sets "scheme_provided". For n ∈ {34..50} the
		// .npz files live under references/schwartz-zwecher-2025/ but are
		// NOT imported into the catalog yet (file sizes 50–500 MB; awaiting
		// our own kin-row constructor — see ROADMAP entry "Schwartz-Zwecher
		// 2025 — n=34..50 import deferred"). All ranks taken from Table 1.
		// n=44 / rank=36110 is the headline (ω₀≈2.773203, fastest known
		// feasible algorithm). Discoveries: yes (strict improvements over
		// Pan 1982's same-table baseline, e.g. n=44 was 36133).
		int[][] szTable1 = {
				{28, 10550}, {30, 12688}, {32, 15096}, {34, 17790},
				{36, 20786}, {38, 24100}, {40, 27748}, {42, 31746},
				{44, 36110}, {46, 40856}, {48, 46000}, {50, 51558},
				{60, 86118},
		};
		for (int[] row : szTable1) {
			int n = row[0];
			int rank = row[1];
			String headlineNote = (n == 44)
					? " Headline result: optimal exponent ω₀≈2.773203 over all feasible-size algorithms."
					: "";
			out.add(simpleCitedJson(new int[]{n, n, n}, rank,
					"Schwartz-Zwecher 2025", 2025,
					"Q⟨" + n + "," + n + "," + n + "⟩:r=" + rank
							+ "; TA-New25 trilinear aggregation + kin-row unification"
							+ "; arXiv:2508.01748 Table 1; non-commutative."
							+ headlineNote,
					false));
		}

		// Hopcroft-Kerr 1971 ⟨a,2,c⟩ FAMILY — closed-form formula.
		// Single entry with symbolic N/P (not enumerated). The SPA renders
		// it once and evaluates on the fly for any concrete (N,P).
		int MAX_DIM = eu.solven.matmul.catalog.CatalogLimits.MAX_DIM;
		out.add(familyJson(
				new String[]{"N", "2", "P"},
				"(3*N*P + max(N,P))/2",
				new String[]{"N", "P"},
				new int[][]{ {2, MAX_DIM}, {2, MAX_DIM} },
				"Hopcroft-Kerr 1971", 1971,
				"R(⟨a,2,c⟩)=(3ac+max(a,c))/2; DIS09 Table 2 row 3; tensor-symmetric so also covers ⟨2,b,c⟩ / ⟨a,b,2⟩; NON-COMMUTATIVE — lifts to recursive matmul",
				false));

		// Islam 2009 (MSc thesis, Western Ontario, Schost supervised) —
		// Proposition 3 generalises Waksman 1970 to arbitrary RECTANGULAR
		// ⟨a,b,c⟩ COMMUTATIVE matmul. Numerically the Waksman formula on
		// the middle axis; the contribution is "first publication of the
		// rectangular case". Same numerical bound was independently
		// re-derived by Rosowski 2019 Thm 2/3 (already populated in
		// derived-bounds via RosowskiBound).
		//
		// We emit one concrete row per (A,B,C) ONLY where Islam's commutative
		// formula STRICTLY beats the best known R/Q/Z non-commutative bound.
		// A tie carries no information (non-commutative ⊆ commutative) and
		// would misattribute the rank to Islam 2009 — canonical example
		// suppressed by this filter: ⟨2,2,2⟩=7 (which is Strassen 1969).
		emitIslam2009Filtered(out, ncByFormat, MAX_DIM);

		return out;
	}

	/**
	 * Enumerate every (A,B,C) in {@code [2, MAX_DIM]³}, compute Islam 2009's
	 * commutative rank, compare against {@link #bestNonCommutativeRankRQZ}
	 * for the same shape, and emit a concrete cited-bound row only when
	 * commutative is strictly smaller. Ties (e.g. ⟨2,2,2⟩=7, ⟨3,3,3⟩=23)
	 * and losses are dropped — they're already covered by an earlier
	 * non-commutative bound and would misattribute the discovery.
	 *
	 * <p>The Islam formula is not symmetric in (A,B,C): B is the parity axis.
	 * Since the matmul tensor IS symmetric in its three axes, we evaluate
	 * Islam at all 6 permutations of the shape and take the minimum — that
	 * is the rank Islam 2009 actually achieves for the unordered shape.</p>
	 */
	private static void emitIslam2009Filtered(List<String> out,
			Map<String, RankRow> ncByFormat, int MAX_DIM) {
		int kept = 0, dropped = 0;
		for (int a = 2; a <= MAX_DIM; a++) {
			for (int b = a; b <= MAX_DIM; b++) {  // canonical a ≤ b ≤ c
				for (int c = b; c <= MAX_DIM; c++) {
					int islam = islamCommutativeMinOverPermutations(a, b, c);
					int nonCom = bestNonCommutativeRankRQZ(a, b, c, ncByFormat);
					if (islam < nonCom) {
						out.add(simpleCitedJson(new int[]{a, b, c}, islam,
								"Islam 2009", 2009,
								"Generalised Waksman; commutative-only; MSc thesis Ch. 6 Prop. 3"
										+ "; strictly improves R/Q/Z non-commutative best at this shape"
										+ "; numerically coincides with Rosowski 2019 Thm 2/3",
								true));
						kept++;
					} else {
						dropped++;
					}
				}
			}
		}
		log.info("Islam 2009 shape filter: kept {} (strict commutative win), dropped {} (tied or dominated by non-com)",
				kept, dropped);
	}

	/** Islam 2009 formula on (A,B,C), B is the parity-driving "middle" axis. */
	private static int islamRank(int A, int B, int C) {
		if (B % 2 == 0) return B * (A * C + A + C - 1) / 2;
		return (B - 1) * (A * C + A + C - 1) / 2 + A * C;
	}

	/**
	 * Minimum Islam rank over all 6 axis permutations of (a,b,c). Matmul
	 * tensor is symmetric in n,m,p, so the rank Islam achieves at the
	 * unordered shape is the best of any axis-assignment for B.
	 */
	private static int islamCommutativeMinOverPermutations(int a, int b, int c) {
		int best = Integer.MAX_VALUE;
		int[][] perms = { {a,b,c}, {a,c,b}, {b,a,c}, {b,c,a}, {c,a,b}, {c,b,a} };
		for (int[] p : perms) {
			best = Math.min(best, islamRank(p[0], p[1], p[2]));
		}
		return best;
	}

	/**
	 * Best known upper bound on the R/Q/Z NON-commutative rank of ⟨n,m,p⟩,
	 * used to filter dominated commutative claims. Consults:
	 * <ol>
	 *   <li>DIS09 Table 3 non-commutative entry, if {@code n=m=p} and in range.</li>
	 *   <li>Hopcroft-Kerr 1971 formula {@code (3ac+max(a,c))/2}, when any
	 *       axis is 2 (rank ceiling if the numerator is odd).</li>
	 *   <li>Trivial bound {@code n·m·p} as the fallback (always valid).</li>
	 * </ol>
	 * Returns the minimum over all applicable bounds. Considered an UPPER
	 * bound — we use it as a "non-com is at least this good" baseline, so a
	 * commutative claim must beat this strictly to carry information.
	 */
	private static int bestNonCommutativeRankRQZ(int n, int m, int p,
			Map<String, RankRow> ncByFormat) {
		int best = n * m * p;  // trivial upper bound, always valid
		// DIS09 cubic.
		if (n == m && m == p) {
			RankRow r = ncByFormat.get(fmtKey(new int[]{n, n, n}));
			if (r != null && !r.isCommutativeOnly && r.rank < best) {
				best = r.rank;
			}
		}
		// Hopcroft-Kerr ⟨a,2,c⟩ family — applies to any shape with a 2-axis.
		int[] sorted = { n, m, p };
		java.util.Arrays.sort(sorted);
		if (sorted[0] == 2) {
			int a = sorted[1], c = sorted[2];
			int num = 3 * a * c + Math.max(a, c);
			int hk = (num + 1) / 2;  // ceiling: HK is an upper bound; ceiling stays valid
			if (hk < best) best = hk;
		}
		return best;
	}

	/**
	 * Emit a family entry with a symbolic formula instead of a single
	 * concrete {@code (n,m,p,rank)} tuple. The SPA evaluates the
	 * formula on-demand for any concrete shape that matches the
	 * symbolic format pattern.
	 *
	 * @param symbolicFormat e.g. {@code {"N","2","P"}} — strings are
	 *                       evaluated as variables, integer strings as constants
	 * @param rankFormula    e.g. {@code "(3*N*P + max(N,P))/2"}
	 * @param vars           the symbol names that appear in the formula
	 * @param ranges         per-var inclusive {@code [lo, hi]}
	 */
	private static String familyJson(String[] symbolicFormat, String rankFormula,
			String[] vars, int[][] ranges, String source, int year,
			String notes, boolean commutative) {
		StringBuilder sb = new StringBuilder();
		sb.append("{");
		sb.append("\"format_symbolic\": [")
				.append(json(symbolicFormat[0])).append(",")
				.append(json(symbolicFormat[1])).append(",")
				.append(json(symbolicFormat[2])).append("], ");
		sb.append("\"field\": \"R/Q/Z\", ");
		if (commutative) sb.append("\"commutative\": true, ");
		sb.append("\"rank_formula\": ").append(json(rankFormula)).append(", ");
		sb.append("\"vars\": [");
		for (int i = 0; i < vars.length; i++) {
			if (i > 0) sb.append(",");
			sb.append(json(vars[i]));
		}
		sb.append("], ");
		sb.append("\"var_ranges\": {");
		for (int i = 0; i < vars.length; i++) {
			if (i > 0) sb.append(", ");
			sb.append(json(vars[i])).append(": [")
					.append(ranges[i][0]).append(",").append(ranges[i][1]).append("]");
		}
		sb.append("}, ");
		sb.append("\"algorithm\": ").append(json(source)).append(", ");
		sb.append("\"source\": ").append(json(source)).append(", ");
		sb.append("\"year\": ").append(year).append(", ");
		sb.append("\"scheme_provided\": false");
		if (notes != null && !notes.isEmpty()) sb.append(", \"notes\": ").append(json(notes));
		sb.append("}");
		return sb.toString();
	}

	private static String simpleCitedJson(int[] fmt, int rank, String source, int year,
			String notes, boolean commutative) {
		StringBuilder sb = new StringBuilder();
		sb.append("{");
		sb.append("\"format\": [").append(fmt[0]).append(",").append(fmt[1]).append(",").append(fmt[2]).append("], ");
		sb.append("\"field\": \"R/Q/Z\", ");
		if (commutative) sb.append("\"commutative\": true, ");
		sb.append("\"rank\": ").append(rank).append(", ");
		sb.append("\"algorithm\": ").append(json(source)).append(", ");
		sb.append("\"source\": ").append(json(source)).append(", ");
		sb.append("\"year\": ").append(year).append(", ");
		sb.append("\"scheme_provided\": ").append(hasSchemeFile(fmt, source));
		if (notes != null && !notes.isEmpty()) sb.append(", \"notes\": ").append(json(notes));
		sb.append("}");
		return sb.toString();
	}

	/**
	 * Border-rank companion of {@link #simpleCitedJson}. Same shape, but
	 * marked {@code kind="border"} so the SPA can lane-separate border-rank
	 * upper bounds from exact-rank claims (they're not directly comparable —
	 * border rank is the limit of rational-coefficient algorithms as the
	 * parameter {@code x → 0}, and only converts to exact rank via the
	 * Schönhage τ-theorem at the asymptotic level).
	 */
	private static String borderCitedJson(int[] fmt, int borderRank, String source, int year,
			String notes) {
		StringBuilder sb = new StringBuilder();
		sb.append("{");
		sb.append("\"format\": [").append(fmt[0]).append(",").append(fmt[1]).append(",").append(fmt[2]).append("], ");
		sb.append("\"field\": \"R/Q/Z\", ");
		sb.append("\"kind\": \"border\", ");
		sb.append("\"border_rank\": ").append(borderRank).append(", ");
		sb.append("\"algorithm\": ").append(json(source)).append(", ");
		sb.append("\"source\": ").append(json(source)).append(", ");
		sb.append("\"year\": ").append(year).append(", ");
		sb.append("\"scheme_provided\": false");
		if (notes != null && !notes.isEmpty()) sb.append(", \"notes\": ").append(json(notes));
		sb.append("}");
		return sb.toString();
	}

	private static String fmtKey(int[] fmt) {
		return fmt[0] + "," + fmt[1] + "," + fmt[2];
	}

	/**
	 * Produce the full attribution for a DIS09 row. Two-part logic:
	 *
	 * <ol>
	 * <li><strong>Base-case algorithms</strong> (the row's rank IS the
	 *     original author's algorithm at its base format): attribute to
	 *     the original author. E.g. Strassen at ⟨2,2,2⟩ = Strassen 1969;
	 *     Laderman at ⟨3,3,3⟩ = Laderman 1976.</li>
	 * <li><strong>Composed / derived results</strong> (DIS09 used the
	 *     named algorithm as a recursive base or composition technique):
	 *     attribute to DIS09 with the base named. E.g. ⟨7,7,7⟩=258
	 *     "Winograd" → "DIS 2009 (Winograd 1971 base)".</li>
	 * </ol>
	 *
	 * Family algorithms (Waksman 1970's closed-form for any cubic)
	 * always attribute to the family author regardless of n.
	 */
	private static String attributionLabel(String algorithm, int[] fmt) {
		int n = fmt[0], m = fmt[1], p = fmt[2];
		switch (algorithm) {
			// FAMILY formulas — apply to any cubic, attribute to family author.
			case "Waksman": return "Waksman 1970";

			// BASE-CASE algorithms — original author at the base format,
			// otherwise DIS09's composition with that algorithm as base.
			case "Strassen":
				return (n == 2 && m == 2 && p == 2)
						? "Strassen 1969"
						: "DIS 2009 (Strassen 1969 base)";
			case "Laderman":
				return (n == 3 && m == 3 && p == 3)
						? "Laderman 1976"
						: "DIS 2009 (Laderman 1976 base)";
			case "Makarov":  // ⟨5,5,5⟩=100, non-commutative (Makarov 1987)
				return (n == 5 && m == 5 && p == 5)
						? "Makarov 1987"
						: "DIS 2009 (Makarov 1987 base)";
			case "Makarov333":  // ⟨3,3,3⟩=22, commutative (Makarov 1986)
				return "Makarov 1986";
			case "Winograd":  // Winograd's ⟨2,2,2⟩=7 alternative to Strassen
				return (n == 2 && m == 2 && p == 2)
						? "Winograd 1971"
						: "DIS 2009 (Winograd 1971 base)";
			case "Winograd2":
				return "DIS 2009 (Winograd 1971 symmetric variant)";

			// DIS09-internal techniques — always DIS09.
			case "Hopcroft332": return "DIS 2009 (Hopcroft-Kerr 1971 base)";
			case "TA":          return "DIS 2009 (Pan trilinear aggregation framework)";
			case "mul121":      return "DIS 2009 (mul121 composition technique)";

			default: return algorithm + " (per DIS09)";
		}
	}

	/** Year for a given attribution label (used as the cited-bound year field). */
	private static String yearForLabel(String label) {
		if (label.startsWith("Strassen 1969")) return "1969";
		if (label.startsWith("Winograd 1971")) return "1971";
		if (label.startsWith("Waksman 1970")) return "1970";
		if (label.startsWith("Laderman 1976")) return "1976";
		if (label.startsWith("Makarov 1987")) return "1987";
		if (label.startsWith("Makarov 1986")) return "1986";
		if (label.startsWith("DIS 2009")) return "2009";
		return "";
	}

	private static final String PERMINOV_REPO_BLOB =
			"https://github.com/dronperminov/FastMatrixMultiplication/blob/master/";

	/** Perminov's June-2026 paper "Meta Flip Graph meets Serendipitous Product:
	 *  new Fast Matrix Multiplication results" — the publication these 17–32
	 *  serendipitous-product bounds come from (NOT status.json / the 2025 paper). */
	private static final String PERMINOV_SERENDIPITOUS_PAPER = "https://arxiv.org/abs/2606.02480";

	/**
	 * Perminov's serendipitous-product improvements over the 17–32 band
	 * ({@code references/perminov-serendipitous-17-32.json}; 971 formats). Each is
	 * a rank claim {@code serendipitous_rank = s1 ⊗ˢ s2} that betters his prior
	 * {@code curr_rank}. We CITE rather than copy — these schemes are large and
	 * already published — so {@code scheme_provided:false}, a {@code source_paper_url}
	 * to the June-2026 paper, and a {@code source_scheme_url} pointing back to
	 * Perminov's repo (the base {@code s1}; the product is rebuilt from
	 * {@code s1, s2}). Field-class from the upstream path tag: {@code ZT/Z} ⇒
	 * integer (Z, reduces everywhere), {@code Q} ⇒ rational.
	 */
	private static List<String> extractPerminovSerendipitous() throws IOException {
		Path src = Path.of("references/perminov-serendipitous-17-32.json");
		List<String> out = new ArrayList<>();
		if (!Files.isRegularFile(src)) {
			return out;
		}
		JsonNode root = JsonMapper.builder().build().readTree(Files.readString(src));
		for (Map.Entry<String, JsonNode> e : root.properties()) {
			String[] dims = e.getKey().split("x");
			if (dims.length != 3) {
				continue;
			}
			JsonNode v = e.getValue();
			int rank = v.get("serendipitous_rank").asInt();
			String path = v.has("path") ? v.get("path").asString() : "";
			String tag = path.isEmpty() ? "ZT"
					: path.substring(path.lastIndexOf('/') + 1).replace(".json", "").replaceAll(".*_", "");
			String field = tag.equals("Q") ? "Q" : "Z";
			int curr = v.has("curr_rank") ? v.get("curr_rank").asInt() : rank;
			StringBuilder sb = new StringBuilder();
			sb.append("{");
			sb.append("\"format\": [").append(dims[0]).append(",").append(dims[1]).append(",").append(dims[2]).append("], ");
			sb.append("\"field\": ").append(json(field)).append(", ");
			sb.append("\"commutative\": false, ");
			sb.append("\"rank\": ").append(rank).append(", ");
			sb.append("\"algorithm\": ").append(json("serendipitous product "
					+ dimStr(v.path("s1").path("dimension")) + " ⊗ˢ " + dimStr(v.path("s2").path("dimension")))).append(", ");
			sb.append("\"source\": \"Perminov 2026 (serendipitous)\", ");
			sb.append("\"year\": 2026, ");
			// The serendipitous 17–32 results are published in Perminov's June-2026
			// paper "Meta Flip Graph meets Serendipitous Product" (arXiv:2606.02480),
			// NOT in status.json or the 2025 ternary-meta-flip-graph paper.
			sb.append("\"source_paper_url\": ").append(json(PERMINOV_SERENDIPITOUS_PAPER)).append(", ");
			sb.append("\"scheme_provided\": false, ");
			if (!path.isEmpty()) {
				sb.append("\"source_scheme_url\": ").append(json(PERMINOV_REPO_BLOB + path)).append(", ");
			}
			sb.append("\"notes\": ").append(json("improves Perminov curr_rank " + curr
					+ "; large scheme not copied — reconstruct from base (s1) + s2 at source_scheme_url"));
			sb.append("}");
			out.add(sb.toString());
		}
		log.info("Perminov serendipitous 17-32: {} cited-bound entries", out.size());
		return out;
	}

	private static String dimStr(JsonNode dim) {
		if (dim == null || !dim.isArray() || dim.size() != 3) {
			return "⟨?⟩";
		}
		return "⟨" + dim.get(0).asInt() + "," + dim.get(1).asInt() + "," + dim.get(2).asInt() + "⟩";
	}

	private static String json(String s) {
		return "\"" + s.replace("\\", "\\\\").replace("\"", "\\\"") + "\"";
	}

	private static class RankRow {
		final int[] fmt;
		final String algorithm;
		int rank;
		Integer commutativeRank;
		String notes;
		boolean isCommutativeOnly = false;
		boolean newInDis09 = false;

		RankRow(int[] fmt, String algorithm, int rank, Integer commutativeRank, String notes) {
			this.fmt = fmt;
			this.algorithm = algorithm;
			this.rank = rank;
			this.commutativeRank = commutativeRank;
			this.notes = notes == null ? "" : notes;
		}

		String toJson() {
			String label = attributionLabel(algorithm, fmt);
			String year = yearForLabel(label);

			StringBuilder sb = new StringBuilder();
			sb.append("{");
			sb.append("\"format\": [").append(fmt[0]).append(",").append(fmt[1]).append(",").append(fmt[2]).append("], ");
			sb.append("\"field\": \"R/Q/Z\", ");
			if (isCommutativeOnly) sb.append("\"commutative\": true, ");
			sb.append("\"rank\": ").append(rank).append(", ");
			if (commutativeRank != null) {
				sb.append("\"commutative_rank\": ").append(commutativeRank).append(", ");
			}
			sb.append("\"algorithm\": ").append(json(algorithm)).append(", ");
			sb.append("\"source\": ").append(json(label)).append(", ");
			if (!year.isEmpty()) sb.append("\"year\": ").append(year).append(", ");
			// Pass `label` (publication-style name) before `algorithm` (internal
			// tag): the on-disk filename prefix tracks the publication
			// (`makarov-1986_*`) not the internal tag (`Makarov333`). Falling
			// back to `algorithm` handles cases where the on-disk name uses
			// the internal tag directly.
			sb.append("\"scheme_provided\": ").append(hasSchemeFile(fmt, label) || hasSchemeFile(fmt, algorithm));
			if (!notes.isEmpty()) sb.append(", \"notes\": ").append(json(notes));
			sb.append("}");
			return sb.toString();
		}
	}

	/**
	 * Filesystem check: does {@code src/main/resources/schemes/{known,derived,curated}/section{maxDim}/}
	 * contain a JSON file matching this {@code (algorithm, fmt)} tuple?
	 * Used to dynamically set {@code "scheme_provided"} so the regenerated
	 * cited-bounds reflects on-disk reality instead of a stale hardcoded
	 * {@code false}. Match is loose: filename starts with the lowercased
	 * algorithm name (spaces → hyphens) and contains the exact shape.
	 */
	private static boolean hasSchemeFile(int[] fmt, String algorithm) {
		int max = Math.max(fmt[0], Math.max(fmt[1], fmt[2]));
		String prefix = algorithm.toLowerCase(java.util.Locale.ROOT).replace(' ', '-');
		String shape = "_" + fmt[0] + "x" + fmt[1] + "x" + fmt[2] + "_";
		// sectionN dirs live under known/derived/curated since the 2026-06-09 split,
		// not at the schemes root — check each subtree.
		for (String sub : new String[] { "known", "derived", "curated" }) {
			java.nio.file.Path dir = java.nio.file.Path.of("src/main/resources/schemes", sub, "section" + max);
			if (!java.nio.file.Files.isDirectory(dir)) continue;
			try (java.util.stream.Stream<java.nio.file.Path> walk = java.nio.file.Files.list(dir)) {
				boolean hit = walk.anyMatch(p -> {
					String name = p.getFileName().toString().toLowerCase(java.util.Locale.ROOT);
					return name.endsWith(".json") && name.contains(shape) && name.startsWith(prefix);
				});
				if (hit) return true;
			} catch (java.io.IOException e) {
				// skip this subtree
			}
		}
		return false;
	}
}
