package eu.solven.matmul.docs.verify;

import java.io.BufferedWriter;
import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.io.PrintWriter;
import java.util.LinkedHashMap;
import java.util.List;
import java.util.Map;
import java.util.TreeMap;

import tools.jackson.databind.JsonNode;
import tools.jackson.databind.json.JsonMapper;
import lombok.extern.slf4j.Slf4j;

/**
 * 3-catalog comparison report: OURS vs the UNION of the external reference
 * catalogs (FMM-Lille ∪ Perminov). Java replacement for the field-naive Python
 * {@code tools/compare_fmm_lille.py} (which was FMM-only).
 *
 * <p>It reads {@code docs/catalog.json}, where {@link GenerateCatalogManifest}
 * has already done the field-correct, commutative-excluded, vs-both comparison
 * per scheme — emitting {@code external_best_rank} / {@code external_best_source}
 * / {@code external_cited} and {@code solven_discovery}. So this tool is a pure
 * aggregation/segmentation over those facts; the comparison itself lives in the
 * manifest generator (single source of truth, regenerated on every sync).</p>
 *
 * <h2>Per-format classification (user 2026-06-06)</h2>
 * Take our best <em>comparable</em> rank per format (non-commutative, char-0)
 * and compare to {@code external_best = min(FMM, Perminov)}:
 * <ul>
 *   <li><b>better</b> → {@code solven_discovery} (provisional; drops when a sync
 *       publishes an equal/better rank).</li>
 *   <li><b>tie + FMM-cited</b> → we match a <em>published</em> result.</li>
 *   <li><b>tie + uncited + derived</b> ({@code atom:false}) → we reproduced FMM's
 *       <em>computed</em> bound by our OWN derivation — the good case.</li>
 *   <li><b>tie + uncited + atom</b> → we hold the computed number only as an
 *       import; a <em>derivation opportunity</em>.</li>
 *   <li><b>worse + FMM-cited</b> → <em>import gap</em> (import the cited scheme).</li>
 *   <li><b>worse + uncited</b> → <em>derivation target</em> (we aim to derive it).</li>
 * </ul>
 *
 * <p>Run: {@code mvn -q -ntp exec:java
 * -Dexec.mainClass=eu.solven.matmul.docs.verify.CompareReferenceCatalogs}</p>
 */
@Slf4j
public final class CompareReferenceCatalogs {

	private CompareReferenceCatalogs() {}

	private static final File CATALOG = new File("docs/catalog.json");
	private static final File OUT = new File("references/reference-comparison.md");

	/** Our best comparable entry for one format, plus the external baseline. */
	private static final class Best {
		int ourRank = Integer.MAX_VALUE;
		boolean ourAtom = true;       // atom-ness of the current best
		int externalRank;
		String externalSource;
		boolean cited;
		int maxDim;
	}

	public static void main(String[] args) throws IOException {
		JsonMapper mapper = JsonMapper.builder().build();
		JsonNode root = mapper.readTree(CATALOG);
		JsonNode schemes = root.get("schemes");
		if (schemes == null || !schemes.isArray()) {
			throw new IOException("docs/catalog.json has no schemes[] — run GenerateCatalogManifest first");
		}

		// Per canonical format, the best comparable (has external_best_rank) entry.
		Map<String, Best> byFormat = new TreeMap<>();
		for (JsonNode s : schemes) {
			if (!s.has("external_best_rank")) continue; // not comparable (commutative / F2-only / no external)
			JsonNode fmt = s.get("format");
			int n = fmt.get(0).asInt(), m = fmt.get(1).asInt(), p = fmt.get(2).asInt();
			int[] sorted = { n, m, p };
			java.util.Arrays.sort(sorted);
			String key = sorted[0] + "x" + sorted[1] + "x" + sorted[2];
			int rank = s.get("rank").asInt();
			Best b = byFormat.computeIfAbsent(key, k -> new Best());
			b.externalRank = s.get("external_best_rank").asInt();
			b.externalSource = s.path("external_best_source").asText("?");
			b.cited = s.path("external_cited").asBoolean(false);
			b.maxDim = sorted[2];
			if (rank < b.ourRank) {
				b.ourRank = rank;
				b.ourAtom = s.path("atom").asBoolean(true);
			}
		}

		// Tally: classification × band.
		Map<String, int[]> bucket = new LinkedHashMap<>(); // label -> [≤16, 17-32, >32]
		String[] labels = { "discovery (better than both)", "tie: matches-cited",
				"tie: derived-match (reproduced computed)", "tie: import-of-computed (derive!)",
				"worse: import-gap (cited)", "worse: derivation-target (uncited)" };
		for (String l : labels) bucket.put(l, new int[3]);
		int tie = 0, better = 0, worse = 0;
		List<String> discoveries = new java.util.ArrayList<>();
		for (Map.Entry<String, Best> e : byFormat.entrySet()) {
			Best b = e.getValue();
			int bandIdx = b.maxDim <= 16 ? 0 : b.maxDim <= 32 ? 1 : 2;
			String label;
			if (b.ourRank < b.externalRank) {
				label = labels[0]; better++;
				discoveries.add(String.format("%s ours=%d ext=%d (%s)",
						e.getKey(), b.ourRank, b.externalRank, b.externalSource));
			} else if (b.ourRank == b.externalRank) {
				tie++;
				label = b.cited ? labels[1] : (b.ourAtom ? labels[3] : labels[2]);
			} else {
				worse++;
				label = b.cited ? labels[4] : labels[5];
			}
			bucket.get(label)[bandIdx]++;
		}

		StringBuilder sb = new StringBuilder();
		sb.append("# Our catalog vs min(FMM-Lille, Perminov)\n\n");
		sb.append("Field-clean (non-commutative, char-0), 3-catalog comparison. Computed by\n");
		sb.append("`GenerateCatalogManifest` (per-scheme `external_best_*` / `solven_discovery`);\n");
		sb.append("this report aggregates `docs/catalog.json`. Supersedes `tools/compare_fmm_lille.py`.\n\n");
		sb.append(String.format("**Formats compared:** %d  —  tie %d, better %d, worse %d%n%n",
				byFormat.size(), tie, better, worse));
		sb.append("## Classification × max-dim band\n\n");
		sb.append("| classification | ≤16 | 17–32 | >32 | total |\n|---|--:|--:|--:|--:|\n");
		for (String l : labels) {
			int[] c = bucket.get(l);
			int t = c[0] + c[1] + c[2];
			if (t == 0) continue;
			sb.append(String.format("| %s | %d | %d | %d | %d |%n", l, c[0], c[1], c[2], t));
		}
		sb.append(String.format("%n## Discoveries (better than BOTH catalogs) — %d, provisional%n%n", discoveries.size()));
		sb.append("Dropped automatically when a future FMM/Perminov sync publishes an equal/better rank.\n\n");
		discoveries.sort(null);
		for (String d : discoveries) sb.append("- ").append(d).append('\n');

		OUT.getParentFile().mkdirs();
		try (PrintWriter pw = new PrintWriter(new BufferedWriter(new FileWriter(OUT)))) {
			pw.print(sb);
		}
		log.info("compared {} formats: tie={} better={} worse={} → wrote {}",
				byFormat.size(), tie, better, worse, OUT);
		System.out.print(sb);
	}
}
