package eu.solven.matmul.docs.explore;

import java.io.IOException;
import java.io.PrintWriter;
import java.nio.file.Files;
import java.nio.file.Path;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Comparator;
import java.util.HashMap;
import java.util.List;
import java.util.Map;
import java.util.Optional;

import eu.solven.matmul.catalog.FieldAwareLookup;
import eu.solven.matmul.recombination.BlockSplitSearch;
import eu.solven.matmul.search.CitedBound;
import lombok.extern.slf4j.Slf4j;
import tools.jackson.databind.JsonNode;
import tools.jackson.databind.json.JsonMapper;

/**
 * Compare what our search algorithm <em>can find</em> for every shape
 * against fmm-lille's published best ranks. Surfaces three categories
 * of shapes:
 *
 * <ul>
 *   <li><strong>WIN</strong>: our predicted rank is strictly lower than
 *       fmm-lille's. Worth materialising + verifying + registering.</li>
 *   <li><strong>TIE</strong>: same rank as fmm-lille but our search
 *       produced a strategy → we have a constructive lineage
 *       (Strassen / Pan pair / concat / recombination tree) where
 *       fmm-lille may have only the rank claim.</li>
 *   <li><strong>GAP</strong>: fmm-lille has a strictly lower rank than
 *       our search can reach. Indicates a missing technique in our
 *       pipeline. The most likely culprit is full Pan TA aggregation
 *       (Islam 2009 §5.4-5.6, not yet implemented as a constructor —
 *       currently only the bound formula is exposed).</li>
 * </ul>
 *
 * <p>This is distinct from {@code references/fmm-lille-discrepancies.md},
 * which compares fmm-lille against schemes <em>materialised on disk</em>.
 * That report tells you which leaf you could update. This report tells
 * you which <em>technique</em> our search pipeline is missing.</p>
 *
 * <p>Output:
 * {@code references/fmm-lille-search-comparison.md} (markdown report
 * with per-category tables).</p>
 *
 * <p>Run:</p>
 * <pre>
 *   mvn -q exec:java -Dexec.mainClass=eu.solven.matmul.docs.explore.FmmSearchComparison
 * </pre>
 */
@Slf4j
public final class FmmSearchComparison {

	private static final int MAX_DIM = eu.solven.matmul.catalog.CatalogLimits.MAX_DIM;
	/**
	 * Non-cubic shapes only walked up to this dim; cubic ⟨n,n,n⟩ goes to {@link #MAX_DIM}.
	 * Slightly wider than {@link eu.solven.matmul.catalog.CatalogLimits#MAX_NONCUBIC_DIM} = 8
	 * because this sweep targets FMM-Lille coverage which carries a few non-trivial entries
	 * up to dim 12.
	 */
	private static final int MAX_NONCUBIC_DIM = 12;
	private static final Path FMM_CATALOG = Path.of("references/catalogs/fmm-lille-catalog.json");
	private static final Path OUT = Path.of("references/fmm-lille-search-comparison.md");

	private FmmSearchComparison() {}

	private record Row(int n, int m, int p, long fmmRank, long ourRank, String ourStrategy) {
		int delta() { return (int) (ourRank - fmmRank); }
	}

	public static void main(String[] args) throws IOException {
		log.info("Loading fmm-lille catalog from {}", FMM_CATALOG);
		Map<String, Long> fmmRanks = loadFmm();
		log.info("  loaded {} fmm-lille entries", fmmRanks.size());

		FieldAwareLookup lookup = new FieldAwareLookup("Q");
		List<BlockSplitSearch.NamedBase> pool = BlockSplitSearch.defaultPool();
		CitedBound sota = new CitedBound(lookup);

		List<Row> wins = new ArrayList<>();
		List<Row> ties = new ArrayList<>();
		List<Row> gaps = new ArrayList<>();
		List<Row> noFmm = new ArrayList<>();   // we found a strategy; fmm has no entry
		List<Row> noOurs = new ArrayList<>();  // fmm has an entry; we found nothing

		int processed = 0;
		long t0 = System.nanoTime();

		for (int n = 2; n <= MAX_DIM; n++) {
			for (int m = n; m <= MAX_DIM; m++) {
				for (int p = m; p <= MAX_DIM; p++) {
					// Scope: cubic up to MAX_DIM, non-cubic only up to
					// MAX_NONCUBIC_DIM (the non-cubic shape space is dense
					// and dominated by larger shapes' kron decomposition;
					// the small slice is the actionable one).
					boolean cubic = (n == m && m == p);
					if (!cubic && p > MAX_NONCUBIC_DIM) continue;
					processed++;
					String key = canon(n, m, p);
					Long fmmRank = fmmRanks.get(key);

					Optional<BlockSplitSearch.NonCubicStrategy> best =
							BlockSplitSearch.findBestStrategy(n, m, p, pool, sota, true);

					// Also consider direct SOTA lookup (catalog leaf + Pan
					// TA formula via CitedBound). For cubic shapes
					// where Pan TA dominates every decomposition, this is
					// the constructive optimum that findBestStrategy alone
					// would miss.
					long sotaDirect = sota.getRank(n, m, p);
					boolean haveSotaLeaf = sotaDirect > 0 && sotaDirect < (long) n * m * p;
					long ourRank;
					String label;
					if (best.isEmpty()) {
						if (!haveSotaLeaf) {
							if (fmmRank != null) {
								noOurs.add(new Row(n, m, p, fmmRank, -1, "—"));
							}
							continue;
						}
						ourRank = sotaDirect;
						label = "sota-leaf";
					} else if (haveSotaLeaf && sotaDirect < best.get().rank()) {
						ourRank = sotaDirect;
						label = "sota-leaf";
					} else {
						ourRank = best.get().rank();
						label = best.get().label();
					}

					if (fmmRank == null) {
						noFmm.add(new Row(n, m, p, -1, ourRank, label));
					} else if (ourRank < fmmRank) {
						wins.add(new Row(n, m, p, fmmRank, ourRank, label));
					} else if (ourRank == fmmRank) {
						ties.add(new Row(n, m, p, fmmRank, ourRank, label));
					} else {
						gaps.add(new Row(n, m, p, fmmRank, ourRank, label));
					}
				}
			}
			if (n % 4 == 0) {
				long elapsedMs = (System.nanoTime() - t0) / 1_000_000L;
				log.info("[progress] processed up to n={}, {} shapes, {}ms",
						n, processed, elapsedMs);
			}
		}
		long totalMs = (System.nanoTime() - t0) / 1_000_000L;
		log.info("Search complete: {} shapes in {}ms", processed, totalMs);
		log.info("  WIN  : {}  (our search beats fmm-lille)", wins.size());
		log.info("  TIE  : {}  (match fmm-lille, we add lineage)", ties.size());
		log.info("  GAP  : {}  (fmm-lille has a better rank — missing technique)", gaps.size());
		log.info("  noFMM: {}  (we have it, fmm-lille doesn't list it)", noFmm.size());
		log.info("  noOUR: {}  (fmm-lille has it, our search found nothing)", noOurs.size());

		emitReport(wins, ties, gaps, noFmm, noOurs);
		log.info("Wrote {}", OUT);
	}

	private static Map<String, Long> loadFmm() throws IOException {
		JsonNode root = JsonMapper.builder().build().readTree(FMM_CATALOG.toFile());
		JsonNode entries = root.get("entries");
		Map<String, Long> out = new HashMap<>();
		if (entries == null || !entries.isArray()) return out;
		for (JsonNode e : entries) {
			JsonNode fmt = e.get("format");
			JsonNode rank = e.get("rank");
			if (fmt == null || rank == null) continue;
			int n = fmt.get(0).asInt();
			int m = fmt.get(1).asInt();
			int p = fmt.get(2).asInt();
			out.put(canon(n, m, p), rank.asLong());
		}
		return out;
	}

	private static String canon(int n, int m, int p) {
		int[] s = { n, m, p };
		Arrays.sort(s);
		return s[0] + "x" + s[1] + "x" + s[2];
	}

	private static void emitReport(List<Row> wins, List<Row> ties,
			List<Row> gaps, List<Row> noFmm, List<Row> noOurs) throws IOException {
		try (PrintWriter pw = new PrintWriter(Files.newBufferedWriter(OUT))) {
			pw.println("# FMM-Lille vs our search pipeline");
			pw.println();
			pw.println("Generated by `eu.solven.matmul.docs.explore.FmmSearchComparison`.");
			pw.println();
			pw.println("**Data sources (no network during the comparison):**");
			pw.println();
			pw.println("- FMM-Lille snapshot: `references/catalogs/fmm-lille-catalog.json`");
			pw.println("  (refresh via `SyncReferenceCatalogs --fmm` or wait for");
			pw.println("  the periodic `.github/workflows/sync-reference-catalogs.yml` run).");
			pw.println("- Our prediction: `BlockSplitSearch.findBestStrategy` (recombination /");
			pw.println("  concat / Kronecker / Pan pair-fused) + a direct `sota-leaf` candidate");
			pw.println("  from `CitedBound` (catalog lookup + Pan TA formula + HK).");
			pw.println();
			pw.println("Compares the best rank our search pipeline (Strassen + concat +");
			pw.println("Kronecker + recombination + Pan pair-fused, with a formula-aware");
			pw.println("SOTA resolver knowing Pan TA / HK formulas) can predict against");
			pw.println("fmm-lille's published best rank.");
			pw.println();
			pw.println("**Scope**:");
			pw.println();
			pw.println("- Cubic ⟨n,n,n⟩: n ∈ [2, " + MAX_DIM + "].");
			pw.println("- Non-cubic ⟨n,m,p⟩ with n ≤ m ≤ p and not-all-equal:");
			pw.println("  p ≤ " + MAX_NONCUBIC_DIM + " (the small-non-cubic slice that");
			pw.println("  isn't already covered by direct Kronecker decomposition of");
			pw.println("  larger cubic shapes).");
			pw.println();
			pw.println("**This is the search-level diff** — distinct from");
			pw.println("`fmm-lille-discrepancies.md` which compares schemes already");
			pw.println("materialised on disk. A GAP here is a *missing technique* in our");
			pw.println("pipeline; a TIE means our search would produce a constructive");
			pw.println("lineage that matches fmm-lille's rank claim, even if no scheme is");
			pw.println("on disk yet.");
			pw.println();
			pw.println("## Headline counts");
			pw.println();
			pw.println("| category | count | meaning |");
			pw.println("|---|---:|---|");
			pw.println("| **WIN** | " + wins.size() + " | our predicted rank `<` fmm-lille — materialise + register |");
			pw.println("| **TIE** | " + ties.size() + " | match fmm-lille, we have a constructive strategy |");
			pw.println("| **GAP** | " + gaps.size() + " | fmm-lille beats us — missing technique |");
			pw.println("| no FMM  | " + noFmm.size() + " | we have a strategy, fmm-lille lists nothing |");
			pw.println("| no ours | " + noOurs.size() + " | fmm-lille has it, our search returned empty |");

			emitSection(pw, "WINS — our search beats fmm-lille", wins, false);
			emitSection(pw, "GAPS — fmm-lille beats our search (missing technique)", gaps, true);
			emitSection(pw, "TIES — same rank, we provide lineage", ties, false);
			emitSection(pw, "fmm-lille has no entry (we have a strategy)", noFmm, false);
			emitSection(pw, "Our search returned nothing (fmm-lille has it)", noOurs, false);
		}
	}

	private static void emitSection(PrintWriter pw, String title, List<Row> rows,
			boolean sortByDeltaDesc) {
		pw.println();
		pw.println("## " + title + " (" + rows.size() + ")");
		pw.println();
		if (rows.isEmpty()) {
			pw.println("_(none)_");
			return;
		}
		// Order: by Δ desc for GAPs (worst-first), else by shape ascending.
		List<Row> sorted = new ArrayList<>(rows);
		if (sortByDeltaDesc) {
			sorted.sort(Comparator.comparingInt(Row::delta).reversed());
		} else {
			sorted.sort(Comparator
					.<Row>comparingInt(r -> r.n)
					.thenComparingInt(r -> r.m)
					.thenComparingInt(r -> r.p));
		}
		pw.println("| shape | fmm rank | our rank | Δ | our strategy |");
		pw.println("|---|---:|---:|---:|---|");
		for (Row r : sorted) {
			String shape = "⟨" + r.n + "," + r.m + "," + r.p + "⟩";
			String fmm = r.fmmRank < 0 ? "—" : Long.toString(r.fmmRank);
			String ours = r.ourRank < 0 ? "—" : Long.toString(r.ourRank);
			String delta;
			if (r.fmmRank < 0 || r.ourRank < 0) {
				delta = "—";
			} else {
				int d = r.delta();
				delta = (d > 0 ? "+" : "") + d;
			}
			pw.println("| " + shape + " | " + fmm + " | " + ours + " | " + delta
					+ " | " + r.ourStrategy + " |");
		}
	}
}
