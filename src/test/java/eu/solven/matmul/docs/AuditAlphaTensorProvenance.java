package eu.solven.matmul.docs;

import java.nio.file.Files;
import java.nio.file.Path;
import java.util.HashMap;
import java.util.HashSet;
import java.util.List;
import java.util.Map;
import java.util.Set;
import java.util.stream.Stream;

import tools.jackson.databind.JsonNode;
import tools.jackson.databind.json.JsonMapper;
import tools.jackson.databind.node.ObjectNode;

/**
 * AT 2022 provenance audit. Distinguishes:
 *
 * <ul>
 *   <li>AT discoveries (the rank is a genuine AlphaTensor improvement over
 *       prior SOTA) — tags {@code "discovery": true,
 *       "attribution_for_rank": "AlphaTensor 2022"}.</li>
 *   <li>AT rediscoveries (the rank matched a previously known bound, AT
 *       just recomputed it) — tags {@code "discovery": false,
 *       "attribution_for_rank": "<prior author>"}.</li>
 * </ul>
 *
 * <p>Source of truth: AT 2022 (Fawzi et al., Nature) Figure 3 prior-art
 * column + Extended Data Table 1 (composed improvements). F₂-specific
 * improvements are also in the discoveries set.</p>
 *
 * <p>Run modes (output is to stdout for review; pass {@code apply} to
 * write back to the JSON files):</p>
 * <pre>
 *   mvn -q -o exec:java -Dexec.classpathScope=test \
 *       -Dexec.mainClass=eu.solven.matmul.catalog.AuditAlphaTensorProvenance \
 *       -Dexec.args="dry-run"
 *   mvn -q -o exec:java -Dexec.classpathScope=test \
 *       -Dexec.mainClass=eu.solven.matmul.catalog.AuditAlphaTensorProvenance \
 *       -Dexec.args="apply"
 * </pre>
 */
public final class AuditAlphaTensorProvenance {

	private AuditAlphaTensorProvenance() {}

	private static final String SCHEMES_DIR = "src/main/resources/schemes";

	/**
	 * (sorted-format, rank) for entries AT 2022 reports as
	 * <strong>discoveries</strong> in standard arithmetic (Fawzi 2022,
	 * Extended Data Table 1 + Figure 3 "AlphaTensor rank Standard").
	 * Format is sorted so {@code 4×4×5} matches {@code ⟨4,4,5⟩} /
	 * {@code ⟨4,5,4⟩} / {@code ⟨5,4,4⟩}.
	 */
	private static final Set<String> AT_DISCOVERIES_STANDARD = new HashSet<>(List.of(
			// Direct discoveries (marked * in Table 1)
			"3x4x5:r47", "4x4x5:r63", "4x5x5:r76",
			// Composed improvements over prior SOTA (Extended Data Table 1)
			"3x4x11:r103", "3x5x9:r105", "3x9x11:r225",
			"4x5x9:r139", "4x5x10:r152", "4x5x11:r169",
			"4x9x10:r255", "4x9x11:r280", "4x11x11:r343", "4x11x12:r366",
			"5x5x7:r134", "5x7x9:r234", "5x7x10:r257", "5x7x11:r280",
			"5x8x9:r262", "5x8x10:r287", "5x8x11:r317",
			"5x9x9:r296", "5x9x10:r323", "5x9x11:r358", "5x9x12:r381",
			"6x7x9:r270", "6x7x10:r296", "6x7x11:r322",
			"6x8x10:r329", "6x8x11:r365",
			"6x9x9:r342", "6x9x10:r373", "6x9x11:r411",
			"7x7x9:r318", "7x7x10:r350", "7x7x11:r384",
			"7x8x9:r354", "7x8x10:r393", "7x8x11:r432", "7x8x12:r462",
			"7x9x9:r399", "7x9x10:r441", "7x9x11:r481", "7x9x12:r510",
			"7x10x10:r478", "7x11x11:r582",
			"8x8x10:r441", "8x8x11:r489",
			"8x9x10:r489", "8x9x11:r533", "8x9x12:r560",
			"8x10x10:r532", "8x10x11:r596", "8x10x12:r636",
			"8x11x11:r649", "8x11x12:r691",
			"9x9x9:r498", "9x9x10:r534", "9x9x11:r576",
			"9x10x10:r606", "9x10x11:r657", "9x10x12:r696",
			"9x11x11:r725", "9x11x12:r760",
			"10x10x10:r682", "10x10x11:r746", "10x10x12:r798",
			"10x11x11:r821", "10x11x12:r874", "10x12x12:r928",
			"11x11x11:r896", "11x11x12:r941", "11x12x12:r990"
	));

	/**
	 * AT 2022 discoveries in F₂ (modular Z₂) arithmetic. These hold ONLY
	 * over F₂; the standard-arithmetic rank may be higher.
	 */
	private static final Set<String> AT_DISCOVERIES_F2 = new HashSet<>(List.of(
			"4x4x4:r47",     // 47/F₂ vs 49 standard (Strassen²)
			"5x5x5:r96"      // 96/F₂ vs 98 standard
	));

	/** Coarse prior-author attribution for (sorted-format, rank) when AT did NOT discover the bound. */
	private static String priorAuthor(String formatKey, int rank, boolean isF2) {
		// Direct lookup table from Fawzi 2022 Figure 3 prior-art column
		Map<String, String> priors = new HashMap<>();
		priors.put("2x2x2:r7", "Strassen 1969");
		priors.put("3x3x3:r23", "Laderman 1976");
		priors.put("4x4x4:r49", "Strassen 1969");        // ⟨2,2,2⟩² recursion
		priors.put("5x5x5:r98", "Sedoglavic-Smirnov 2021"); // composition
		priors.put("2x2x3:r11", "Hopcroft-Kerr 1971");
		priors.put("2x2x4:r14", "Hopcroft-Kerr 1971");
		priors.put("2x2x5:r18", "Hopcroft-Kerr 1971");
		priors.put("2x2x6:r21", "Hopcroft-Kerr 1971");
		priors.put("2x2x7:r25", "Hopcroft-Kerr 1971");
		priors.put("2x2x8:r28", "Hopcroft-Kerr 1971");
		priors.put("2x3x3:r15", "Hopcroft-Kerr 1971");
		priors.put("2x3x4:r20", "Hopcroft-Kerr 1971");
		priors.put("2x3x5:r25", "Hopcroft-Kerr 1971");
		priors.put("2x4x4:r26", "Hopcroft-Kerr 1971");
		priors.put("2x4x5:r33", "Hopcroft-Kerr 1971");
		priors.put("2x5x5:r40", "Hopcroft-Kerr 1971");
		priors.put("3x3x4:r29", "Smirnov 2013");
		priors.put("3x3x5:r36", "Smirnov 2013");
		priors.put("3x4x4:r38", "Smirnov 2013");
		priors.put("3x4x5:r48", "Smirnov 2013");           // AT improved to 47
		priors.put("3x5x5:r58", "Sedoglavic-Smirnov 2021");
		priors.put("4x4x5:r64", "Smirnov 2013");           // AT improved to 63
		priors.put("4x5x5:r80", "Smirnov 2013");           // AT improved to 76
		// F₂-specific prior bounds (mirror standard up to AT's known modular improvements)
		// For F₂ we conservatively cite the standard-arithmetic prior.
		return priors.getOrDefault(formatKey + ":r" + rank,
				isF2 ? "Strassen 1969 (recursive F₂)" : "Strassen 1969 (recursive)");
	}

	public static void main(String[] args) throws Exception {
		boolean apply = args.length > 0 && "apply".equalsIgnoreCase(args[0]);
		System.out.println("Mode: " + (apply ? "APPLY" : "DRY-RUN"));

		JsonMapper mapper = new JsonMapper();
		int audited = 0, discoveries = 0, rediscoveries = 0;
		Map<String, Integer> attributionHist = new HashMap<>();

		try (Stream<Path> walk = Files.walk(Path.of(SCHEMES_DIR))) {
			List<Path> files = walk
					.filter(p -> p.toString().endsWith(".json"))
					.filter(p -> p.getFileName().toString().toLowerCase().startsWith("alphatensor"))
					.sorted()
					.toList();
			for (Path p : files) {
				audited++;
				JsonNode node = mapper.readTree(p.toFile());
				if (!(node instanceof ObjectNode root)) continue;
				int n = root.path("n").get(0).asInt();
				int m = root.path("n").get(1).asInt();
				int p3 = root.path("n").get(2).asInt();
				int rank = parseRank(p.getFileName().toString());
				boolean isF2 = p.getFileName().toString().contains("F2");
				String key = sortedKey(n, m, p3) + ":r" + rank;
				boolean isDiscovery = AT_DISCOVERIES_STANDARD.contains(key)
						|| (isF2 && AT_DISCOVERIES_F2.contains(key));
				String attribution;
				if (isDiscovery) {
					discoveries++;
					attribution = "AlphaTensor 2022";
				} else {
					rediscoveries++;
					attribution = priorAuthor(sortedKey(n, m, p3), rank, isF2);
				}
				attributionHist.merge(attribution, 1, Integer::sum);
				if (apply) {
					root.put("discovery", isDiscovery);
					root.put("attribution_for_rank", attribution);
					// Canonical formatter, NOT Jackson's default pretty-printer,
					// so we don't re-introduce a divergent style on edit.
					eu.solven.matmul.catalog.MatrixJsonFormatter.write(p.toFile(), root);
				}
			}
		}

		System.out.println();
		System.out.println("== Summary ==");
		System.out.println("Files audited:    " + audited);
		System.out.println("AT discoveries:   " + discoveries);
		System.out.println("AT rediscoveries: " + rediscoveries);
		System.out.println();
		System.out.println("Attribution distribution:");
		attributionHist.entrySet().stream()
				.sorted((a, b) -> Integer.compare(b.getValue(), a.getValue()))
				.forEach(e -> System.out.printf("  %-40s %d%n", e.getKey(), e.getValue()));
		if (!apply) {
			System.out.println();
			System.out.println("Re-run with -Dexec.args=apply to write changes.");
		}
	}

	private static String sortedKey(int n, int m, int p) {
		int[] s = { n, m, p };
		java.util.Arrays.sort(s);
		return s[0] + "x" + s[1] + "x" + s[2];
	}

	private static int parseRank(String filename) {
		var mt = java.util.regex.Pattern.compile("_r(\\d+)_").matcher(filename);
		return mt.find() ? Integer.parseInt(mt.group(1)) : -1;
	}
}
