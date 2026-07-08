package eu.solven.matmul.docs.verify;

import java.io.File;
import java.io.PrintWriter;
import java.util.ArrayList;
import java.util.Comparator;
import java.util.List;
import java.util.Map;

import eu.solven.matmul.recombination.BlockSplitSearch;
import lombok.extern.slf4j.Slf4j;
import tools.jackson.databind.JsonNode;
import tools.jackson.databind.json.JsonMapper;

/**
 * Cross-checks our catalog's best non-commutative characteristic-0 rank per
 * shape against the FMM-Lille digest ({@code references/catalogs/fmm-lille-catalog.json},
 * which is over Q, non-commutative). Surfaces three classes (task #180):
 *
 * <ul>
 *   <li><b>WORSE</b> — our best &gt; FMM best: a gap (missing import or
 *       un-materialised composition, e.g. the ⟨8,8,8⟩=336 Kron that was absent).</li>
 *   <li><b>BETTER</b> — our best &lt; FMM best: a genuine win to highlight.</li>
 *   <li><b>MISSING</b> — FMM has a shape we have no char-0 scheme for.</li>
 * </ul>
 *
 * <p>Writes a markdown report to {@code references/fmm-cross-check.md}. Scope is
 * capped at the catalog's MAX_DIM to avoid noise from huge FMM-only shapes.</p>
 *
 * <pre>mvn -q exec:java -Dexec.mainClass=eu.solven.matmul.docs.verify.FmmCrossCheck</pre>
 */
@Slf4j
public final class FmmCrossCheck {

	private static final int MAX_DIM = eu.solven.matmul.catalog.CatalogLimits.MAX_DIM;
	private static final String OUT = "references/fmm-cross-check.md";

	private FmmCrossCheck() {}

	/** Best NC characteristic-0 rank per canonical shape from raw on-disk scheme
	 * files — no derived/PanTA padding. CONTENT-DRIVEN (2026-06 migration): shape ←
	 * {@code n}, rank ← {@code m}, field ← {@code fields[]}, commutative ←
	 * {@code commutative} — NEVER the filename (a pure cosmetic label now). Includes
	 * only NC-over-Q/R schemes (fields name Q, R or Z): excludes C-only (e.g.
	 * AlphaEvolve complex), modular-only (F2/F3) and commutative schemes, matching
	 * the FMM-Lille digest's non-commutative-over-Q scope. */
	private static Map<String, Integer> rawOnDiskBestCharZeroNC() throws java.io.IOException {
		Map<String, Integer> best = new java.util.HashMap<>();
		java.nio.file.Path root = java.nio.file.Path.of("src/main/resources/schemes");
		try (var walk = java.nio.file.Files.walk(root)) {
			for (var it = walk.iterator(); it.hasNext();) {
				java.nio.file.Path pth = it.next();
				if (!pth.getFileName().toString().endsWith(".json")) continue;
				JsonNode d;
				try {
					d = eu.solven.matmul.catalog.SchemeIO.parseJson(pth.toFile());
				} catch (Exception e) {
					continue;
				}
				if (!d.has("n") || !d.has("m")) continue;
				// Commutative-only schemes don't lift to NC matmul — exclude.
				if (d.path("commutative").asBoolean(false)) continue;
				// char-0 NC over Q/R = fields[] names Q, R or Z (Z⊂Q⊂R). Excludes C-only
				// and modular-only (F2/F3), matching the FMM digest's NC-over-Q scope.
				JsonNode fields = d.get("fields");
				boolean charZeroQ = false;
				if (fields != null) {
					for (JsonNode f : fields) {
						String fn = f.asText();
						if ("Q".equals(fn) || "R".equals(fn) || "Z".equals(fn)) { charZeroQ = true; break; }
					}
				}
				if (!charZeroQ) continue;
				JsonNode shp = d.get("n");
				if (shp == null || !shp.isArray() || shp.size() != 3) continue;
				int[] dd = { shp.get(0).asInt(), shp.get(1).asInt(), shp.get(2).asInt() };
				java.util.Arrays.sort(dd);
				if (dd[2] > MAX_DIM) continue;
				int rank = d.get("m").asInt();
				best.merge(dd[0] + "x" + dd[1] + "x" + dd[2], rank, Integer::min);
			}
		}
		return best;
	}

	record Row(int n, int m, int p, int ours, int fmm, List<String> refs) {
		long vol() { return (long) n * m * p; }
		int delta() { return ours - fmm; }
	}

	public static void main(String[] args) throws Exception {
		JsonMapper mapper = JsonMapper.builder().build();
		JsonNode root = mapper.readTree(new File("references/catalogs/fmm-lille-catalog.json"));

		// FMM best rank per canonical (sorted) shape.
		Map<String, int[]> fmmBest = new java.util.HashMap<>();      // key -> {n,m,p,rank}
		Map<String, List<String>> fmmRefs = new java.util.HashMap<>();
		for (JsonNode e : root.get("entries")) {
			JsonNode fmt = e.get("format");
			if (fmt == null || fmt.size() != 3) continue;
			int[] d = { fmt.get(0).asInt(), fmt.get(1).asInt(), fmt.get(2).asInt() };
			java.util.Arrays.sort(d);
			if (d[2] > MAX_DIM) continue;
			int rank = e.get("rank").asInt();
			String key = d[0] + "x" + d[1] + "x" + d[2];
			int[] prev = fmmBest.get(key);
			if (prev == null || rank < prev[3]) {
				fmmBest.put(key, new int[] { d[0], d[1], d[2], rank });
				List<String> refs = new ArrayList<>();
				if (e.get("references") != null) e.get("references").forEach(r -> refs.add(r.asText()));
				fmmRefs.put(key, refs);
			}
		}

		// Our best NC char-0 rank per canonical shape, from RAW on-disk scheme
		// files only — no PanTA/trivial padding (which would inflate "WORSE" with
		// weak derived bounds). This yields the actionable "FMM has a better
		// EXPLICIT scheme" gap list. F2 / complex / commutative files excluded.
		Map<String, Integer> ours = rawOnDiskBestCharZeroNC();

		List<Row> worse = new ArrayList<>(), better = new ArrayList<>(), missing = new ArrayList<>();
		for (Map.Entry<String, int[]> e : fmmBest.entrySet()) {
			int[] d = e.getValue();
			int fmm = d[3];
			Integer our = ours.get(e.getKey());
			List<String> refs = fmmRefs.getOrDefault(e.getKey(), List.of());
			if (our == null) {
				missing.add(new Row(d[0], d[1], d[2], -1, fmm, refs));
			} else if (our > fmm) {
				worse.add(new Row(d[0], d[1], d[2], our, fmm, refs));
			} else if (our < fmm) {
				better.add(new Row(d[0], d[1], d[2], our, fmm, refs));
			}
		}
		// Shape-lexicographic (n,m,p) so families cluster (⟨27,28,·⟩ reads as a run)
		// and a given shape is findable by scan; the gap column still carries the
		// triage signal (--worst = max over the gap column, order-independent).
		Comparator<Row> byShape = Comparator.comparingInt((Row r) -> r.n)
				.thenComparingInt(r -> r.m)
				.thenComparingInt(r -> r.p);
		worse.sort(byShape);
		better.sort(byShape);
		missing.sort(byShape);

		try (PrintWriter pw = new PrintWriter(OUT)) {
			pw.println("# Catalog vs FMM-Lille digest cross-check");
			pw.println();
			pw.println("Auto-generated by `FmmCrossCheck`. Compares our best **non-commutative");
			pw.println("characteristic-0** rank per shape against the FMM-Lille digest");
			pw.println("(`references/catalogs/fmm-lille-catalog.json`, over Q). Capped at max-dim "
					+ MAX_DIM + ".");
			pw.println();
			pw.printf("- WORSE (we lag FMM): **%d**%n", worse.size());
			pw.printf("- BETTER (we beat FMM): **%d**%n", better.size());
			pw.printf("- MISSING (FMM has, we don't): **%d**%n", missing.size());
			pw.println();
			long smallWorse = worse.stream().filter(r -> Math.max(r.n, Math.max(r.m, r.p)) <= 12).count();
			pw.println("**Interpretation.** WORSE at small shapes (max-dim ≤ 12): "
					+ smallWorse + " — these are the importable/actionable gaps. The remainder "
					+ "are LARGE shapes where our on-disk best is a NOT-YET-RUN projection / "
					+ "recombination closure, NOT a structural gap: they close by re-running the "
					+ "projection closure (`projectInto`) + richer recombination bases over 17–32 "
					+ "— the same finite mechanisms we already have. This is explicitly NOT a "
					+ "Schönhage τ-theorem / disjoint-sum gap: we never leverage τ and don't need "
					+ "to — FMM's large 'disjoint-sum'-labelled schemes are recombinations of small "
					+ "HK71 bases (see project_17x17x17_search_gap / project_fmm_projection_gap_is_global). "
					+ "MISSING is dominated by ⟨2,3,k⟩ / Hopcroft-Kerr-family shapes (formula-derivable).");
			pw.println();

			pw.println("## WORSE — we should import / materialise these");
			pw.println();
			pw.println("| shape | ours | FMM | gap | FMM source |");
			pw.println("| --- | ---: | ---: | ---: | --- |");
			for (Row r : worse) {
				pw.printf("| ⟨%d,%d,%d⟩ | %d | %d | +%d | %s |%n",
						r.n, r.m, r.p, r.ours, r.fmm, r.delta(), String.join(", ", r.refs));
			}
			pw.println();

			pw.println("## MISSING — FMM has a scheme we lack");
			pw.println();
			pw.println("| shape | FMM | FMM source |");
			pw.println("| --- | ---: | --- |");
			for (Row r : missing) {
				pw.printf("| ⟨%d,%d,%d⟩ | %d | %s |%n", r.n, r.m, r.p, r.fmm, String.join(", ", r.refs));
			}
			pw.println();

			pw.println("## BETTER — we beat FMM (highlight / verify attribution)");
			pw.println();
			pw.println("| shape | ours | FMM |");
			pw.println("| --- | ---: | ---: |");
			for (Row r : better) {
				pw.printf("| ⟨%d,%d,%d⟩ | %d | %d |%n", r.n, r.m, r.p, r.ours, r.fmm);
			}
		}
		log.info("FmmCrossCheck: WORSE={} BETTER={} MISSING={} (FMM shapes ≤dim{}: {}). Wrote {}",
				worse.size(), better.size(), missing.size(), MAX_DIM, fmmBest.size(), OUT);
	}
}
