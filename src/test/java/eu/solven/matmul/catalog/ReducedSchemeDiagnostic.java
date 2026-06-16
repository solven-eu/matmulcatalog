package eu.solven.matmul.catalog;

import java.io.File;
import java.util.Arrays;
import java.util.Map;
import java.util.TreeMap;

import eu.solven.matmul.NonCubicBilinearAlgorithm;
import eu.solven.matmul.Verifier;
import eu.solven.matmul.catalog.SchemeIO;

/**
 * Diagnostic: scans only the {@code _reduced} schemes, reports each as
 * `OK / load-err / verify-err`. Groups by error class for triage.
 */
public class ReducedSchemeDiagnostic {

	public static void main(String[] args) throws Exception {
		File dir = new File("src/main/resources/schemes");
		File[] files = dir.listFiles((d, n) -> n.contains("_reduced") && n.endsWith(".json"));
		if (files == null) { System.out.println("no schemes/"); return; }
		Arrays.sort(files);

		int ok = 0;
		Map<String, Integer> errCounts = new TreeMap<>();
		Map<String, String> firstExample = new TreeMap<>();

		for (File f : files) {
			try {
				tools.jackson.databind.JsonNode root = SchemeIO.parseJson(f);
				NonCubicBilinearAlgorithm alg = SchemeIO.readReduced(root);
				boolean z2 = SchemeIO.isZ2(root);
				int wrong = z2 ? Verifier.residualNonCubicF2(alg)
						: (int) Math.round(Verifier.residualNonCubic(alg) * 1000);
				if (z2 ? wrong == 0 : wrong == 0) {
					ok++;
				} else {
					String key = z2 ? "verify-fail-F2" : "verify-fail-real";
					errCounts.merge(key, 1, Integer::sum);
					firstExample.putIfAbsent(key, f.getName() + " (wrong=" + wrong + ")");
				}
			} catch (Exception e) {
				String key = "load-" + e.getClass().getSimpleName();
				errCounts.merge(key, 1, Integer::sum);
				firstExample.putIfAbsent(key, f.getName() + ": " + e.getMessage());
			}
		}

		System.out.printf("Reduced schemes: total=%d ok=%d failed=%d%n",
				files.length, ok, files.length - ok);
		System.out.println();
		System.out.println("Error category counts:");
		for (Map.Entry<String, Integer> e : errCounts.entrySet()) {
			System.out.printf("  %-30s %4d   example: %s%n",
					e.getKey(), e.getValue(), firstExample.get(e.getKey()));
		}
	}
}
