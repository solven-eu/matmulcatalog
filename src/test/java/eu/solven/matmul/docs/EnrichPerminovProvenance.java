package eu.solven.matmul.docs;

import java.io.FileReader;
import java.io.Reader;
import java.nio.file.Files;
import java.nio.file.Path;
import java.util.HashMap;
import java.util.List;
import java.util.Map;
import java.util.stream.Stream;

import tools.jackson.databind.JsonNode;
import tools.jackson.databind.json.JsonMapper;
import tools.jackson.databind.node.ObjectNode;

/**
 * One-shot tool: reads
 * {@code dronperminov/FastMatrixMultiplication/schemes/status.json}
 * and enriches each {@code perminov-*.json} scheme file in our catalog
 * with the true {@code attribution_for_rank} extracted from
 * status.json's {@code source} field.
 *
 * <p>Use the published Perminov status index as the source of truth
 * for provenance: the {@code source} path encodes the originating
 * sub-catalog ({@code classic/}, {@code alpha_tensor/},
 * {@code alpha_evolve/}, {@code a_60_addition/}, {@code tensor/},
 * etc.), which gives proper attribution instead of tagging everything
 * as "Perminov".</p>
 *
 * <p>Matches each catalog file by {@code (format, rank, additions)}
 * against status.json entries. Adds two JSON fields:</p>
 * <ul>
 *   <li>{@code attribution_for_rank}: the simplified label
 *       ("AlphaTensor 2022", "AlphaEvolve 2025", "Stapleton 2025
 *       (a=60)", "Moosbauer", "classic", "tensor decomposition")</li>
 *   <li>{@code original_source_path}: the raw {@code source} string
 *       from status.json for audit purposes</li>
 * </ul>
 *
 * <p>Run modes:</p>
 * <pre>
 *   mvn -q -o exec:java -Dexec.mainClass=eu.solven.matmul.catalog.EnrichPerminovProvenance \
 *       -Dexec.args="dry-run"     # print proposed changes, modify nothing
 *   mvn -q -o exec:java -Dexec.mainClass=eu.solven.matmul.catalog.EnrichPerminovProvenance \
 *       -Dexec.args="apply"       # write changes in place
 * </pre>
 *
 * <p>Fetch status.json yourself first:</p>
 * <pre>
 *   curl -fL "https://raw.githubusercontent.com/dronperminov/FastMatrixMultiplication/master/schemes/status.json" \
 *        -o /tmp/perminov-status.json
 * </pre>
 */
public final class EnrichPerminovProvenance {

	private EnrichPerminovProvenance() {}

	private static final String STATUS_JSON = "/tmp/perminov-status.json";
	private static final String SCHEMES_DIR = "src/main/resources/schemes";

	public static void main(String[] args) throws Exception {
		boolean apply = args.length > 0 && "apply".equalsIgnoreCase(args[0]);
		System.out.println("Mode: " + (apply ? "APPLY (writing JSON files)" : "DRY-RUN (no writes)"));

		// (1) Build lookup: (n×m×p, rank, additions) → source path
		Map<String, String> lookup = buildLookup();
		System.out.println("Loaded " + lookup.size() + " unique (format, rank, additions) → source mappings from "
				+ STATUS_JSON);

		// (2) Walk perminov-*.json catalog files, match, propose changes
		JsonMapper mapper = new JsonMapper();
		int audited = 0, matched = 0, missing = 0, alreadyTagged = 0;
		Map<String, Integer> attributionHist = new HashMap<>();

		try (Stream<Path> walk = Files.walk(Path.of(SCHEMES_DIR))) {
			List<Path> files = walk
					.filter(p -> p.toString().endsWith(".json"))
					.filter(p -> p.getFileName().toString().toLowerCase().startsWith("perminov-"))
					.sorted()
					.toList();
			for (Path p : files) {
				audited++;
				JsonNode node = mapper.readTree(p.toFile());
				if (!(node instanceof ObjectNode root)) continue;
				if (root.has("attribution_for_rank")) {
					alreadyTagged++;
					continue;
				}
				int n = root.path("n").get(0).asInt();
				int m = root.path("n").get(1).asInt();
				int p3 = root.path("n").get(2).asInt();
				int rank = parseRank(p.getFileName().toString());
				String key = n + "x" + m + "x" + p3 + ":r" + rank;
				String origSrc = lookup.get(key);
				if (origSrc == null) {
					missing++;
					if (missing <= 10) {
						System.out.println("  no match: " + p.getFileName() + " (key=" + key + ")");
					}
					continue;
				}
				matched++;
				String attribLabel = simplifyAttribution(origSrc);
				attributionHist.merge(attribLabel, 1, Integer::sum);
				if (apply) {
					root.put("attribution_for_rank", attribLabel);
					root.put("original_source_path", origSrc);
					// Canonical formatter, NOT Jackson's default pretty-printer,
					// so we don't re-introduce a divergent style on edit.
					eu.solven.matmul.catalog.MatrixJsonFormatter.write(p.toFile(), root);
				}
			}
		}

		System.out.println();
		System.out.println("== Summary ==");
		System.out.println("Files audited:        " + audited);
		System.out.println("Already had tag:      " + alreadyTagged);
		System.out.println("Matched in status:    " + matched);
		System.out.println("No match:             " + missing);
		System.out.println();
		System.out.println("Attribution distribution (proposed):");
		attributionHist.entrySet().stream()
				.sorted((a, b) -> Integer.compare(b.getValue(), a.getValue()))
				.forEach(e -> System.out.printf("  %-30s %d%n", e.getKey(), e.getValue()));
		if (!apply) {
			System.out.println();
			System.out.println("Re-run with -Dexec.args=apply to write changes.");
		}
	}

	/** Build {@code "(n×m×p, rank, additions)" → source path} from status.json. */
	private static Map<String, String> buildLookup() throws Exception {
		JsonMapper mapper = new JsonMapper();
		JsonNode root;
		try (Reader r = new FileReader(STATUS_JSON)) {
			root = mapper.readTree(r);
		}
		// Collect all candidates per (fmt, rank, adds), then pick best by priority.
		Map<String, java.util.List<String>> candidates = new HashMap<>();
		for (String fmt : root.propertyNames()) {
			JsonNode schemesByField = root.get(fmt).get("schemes");
			if (schemesByField == null) continue;
			for (String field : schemesByField.propertyNames()) {
				for (JsonNode entry : schemesByField.get(field)) {
					int rank = entry.path("rank").asInt();
					int adds = entry.path("complexity").asInt();
					String src = entry.path("source").asText(null);
					if (src == null) continue;
					// Match on (format, rank) only — Perminov's addition-reduction often gives a different
					// 'complexity' than the originating algorithm, but the RANK provenance is what
					// we want to attribute (the original author is responsible for the rank;
					// Perminov for the addition-count improvement).
					String key = fmt + ":r" + rank;
					candidates.computeIfAbsent(key, k -> new java.util.ArrayList<>()).add(src);
					// Status.json also lists `duplicates` — same scheme via alternative sources.
					JsonNode dups = entry.get("duplicates");
					if (dups != null) {
						for (JsonNode d : dups) candidates.get(key).add(d.asText());
					}
				}
			}
		}
		Map<String, String> out = new HashMap<>();
		for (var e : candidates.entrySet()) {
			// Pick the source with the LOWEST priority value (most-authoritative).
			out.put(e.getKey(), e.getValue().stream()
					.min(java.util.Comparator.comparingInt(EnrichPerminovProvenance::sourcePriority))
					.orElseThrow());
		}
		return out;
	}

	/** Lower = more authoritative (original author). Higher = Perminov-derived. */
	private static int sourcePriority(String sourcePath) {
		if (sourcePath.contains("/classic/")) return 1;
		if (sourcePath.contains("/alpha_tensor/")) return 2;
		if (sourcePath.contains("/alpha_evolve/")) return 2;
		if (sourcePath.contains("/a_60_addition/")) return 2;
		if (sourcePath.contains("/jakobmoosbauer_symmetric_flips/")) return 3;
		if (sourcePath.contains("/jakobmoosbauer_flips/")) return 3;
		if (sourcePath.contains("/fmm_add_reduction/")) return 4;
		if (sourcePath.contains("/meta_flip_graph/")) return 5;
		if (sourcePath.contains("/tensor/")) return 6;
		// schemes/results/* and anything else: Perminov-derived (lowest priority)
		return 9;
	}

	/** Map a status.json source path to a human label. */
	private static String simplifyAttribution(String sourcePath) {
		if (sourcePath.contains("/alpha_tensor/")) return "AlphaTensor 2022";
		if (sourcePath.contains("/alpha_evolve/")) return "AlphaEvolve 2025";
		if (sourcePath.contains("/a_60_addition/")) return "Stapleton 2025 (a=60)";
		if (sourcePath.contains("/jakobmoosbauer_symmetric_flips/")
				|| sourcePath.contains("/jakobmoosbauer_flips/")) return "Moosbauer (symmetric flips)";
		if (sourcePath.contains("/fmm_add_reduction/")) return "fmm reduction (Moosbauer derivative)";
		if (sourcePath.contains("/meta_flip_graph/")) return "meta-flip-graph search";
		if (sourcePath.contains("/classic/")) {
			// e.g. Strassen-222-7-18.m or Laderman-333-23.m
			String fn = sourcePath.substring(sourcePath.lastIndexOf('/') + 1);
			if (fn.toLowerCase().startsWith("strassen")) return "Strassen 1969";
			if (fn.toLowerCase().startsWith("laderman")) return "Laderman 1976";
			if (fn.toLowerCase().startsWith("smirnov")) return "Smirnov 2013";
			return "classic (" + fn + ")";
		}
		if (sourcePath.contains("/tensor/")) return "Perminov (tensor decomposition)";
		// schemes/results/* — Perminov's own derivations (no known/* primary source for this rank)
		return "Perminov (FastMatrixMultiplication)";
	}

	private static int parseRank(String filename) {
		var m = java.util.regex.Pattern.compile("_r(\\d+)_").matcher(filename);
		return m.find() ? Integer.parseInt(m.group(1)) : -1;
	}

	private static int countAdditions(String filename) {
		var m = java.util.regex.Pattern.compile("_a(\\d+)\\.json$").matcher(filename);
		return m.find() ? Integer.parseInt(m.group(1)) : -1;
	}
}
