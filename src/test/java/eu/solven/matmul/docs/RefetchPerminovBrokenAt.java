package eu.solven.matmul.docs;

import java.io.BufferedReader;
import java.io.FileReader;
import java.io.InputStreamReader;
import java.io.Reader;
import java.net.URI;
import java.net.http.HttpClient;
import java.net.http.HttpRequest;
import java.net.http.HttpResponse;
import java.nio.file.Files;
import java.nio.file.Path;
import java.util.regex.Matcher;
import java.util.regex.Pattern;

import eu.solven.matmul.NonCubicBilinearAlgorithm;
import eu.solven.matmul.Verifier;
import eu.solven.matmul.catalog.FieldAwareLookup;
import eu.solven.matmul.catalog.SchemeIO;
import tools.jackson.databind.JsonNode;
import tools.jackson.databind.json.JsonMapper;

/**
 * For each {@link FieldAwareLookup#KNOWN_BROKEN_FILES} entry that
 * {@link AuditBrokenAtSchemes} found IN Perminov's status.json, fetch
 * the actual factor matrices from the Perminov GitHub repo, verify
 * them, and write a CORRECTLY-tagged scheme JSON to disk.
 *
 * <p>Output file naming reflects TRUE provenance per Perminov:</p>
 * <ul>
 *   <li>{@code /known/alpha_tensor/} source → {@code alphatensor-{field}_{shape}_r{r}_a{a}.json}</li>
 *   <li>{@code /known/tensor/} source → {@code perminov-tensor-{field}_{shape}_r{r}_a{a}.json}</li>
 *   <li>{@code /results/} source → {@code perminov-results-{field}_{shape}_r{r}_a{a}.json}</li>
 * </ul>
 *
 * <p>Output field is sourced from Perminov's section (Z / ZT / Q),
 * matching the actual coefficient ring the scheme uses. This fixes
 * the original mis-attribution: many "alphatensor-Z" files were
 * really Q-coefficient schemes by other authors.</p>
 */
public final class RefetchPerminovBrokenAt {

	private static final String STATUS_JSON = "/tmp/perminov-status.json";
	private static final String PERMINOV_BASE =
			"https://raw.githubusercontent.com/dronperminov/FastMatrixMultiplication/master/";

	private static final Pattern FILENAME = Pattern.compile(
			"alphatensor-Z_(\\d+)x(\\d+)x(\\d+)_r(\\d+)_a\\d+\\.json");

	public static void main(String[] args) throws Exception {
		boolean apply = args.length > 0 && "apply".equalsIgnoreCase(args[0]);
		System.out.println("Mode: " + (apply ? "APPLY (writing files)" : "DRY-RUN (no writes)"));

		JsonMapper mapper = new JsonMapper();
		JsonNode statusRoot;
		try (Reader r = new FileReader(STATUS_JSON)) {
			statusRoot = mapper.readTree(r);
		}

		HttpClient http = HttpClient.newHttpClient();
		int refetchedAt = 0, refetchedOther = 0, skipped = 0, failed = 0;

		for (String name : getKnownBrokenAtFiles()) {
			Matcher m = FILENAME.matcher(name);
			if (!m.matches()) continue;
			int n = Integer.parseInt(m.group(1));
			int mm = Integer.parseInt(m.group(2));
			int p = Integer.parseInt(m.group(3));
			int rank = Integer.parseInt(m.group(4));
			String fmtKey = n + "x" + mm + "x" + p;

			JsonNode fmt = statusRoot.get(fmtKey);
			if (fmt == null || fmt.get("schemes") == null) {
				System.out.printf("SKIP (no status entry): %s%n", name);
				skipped++;
				continue;
			}
			JsonNode schemes = fmt.get("schemes");

			// Look across all fields, pick best (lowest-priority) source per Perminov.
			// We prefer the `source` field (actual scheme location) over `path` (status view).
			// `source` may point to .mpl (Maple) files we can't parse — skip those.
			String bestField = null;
			String bestPath = null;
			String bestSourceTag = null;  // "alphatensor", "perminov-tensor", "perminov-results"
			int bestPriority = Integer.MAX_VALUE;
			for (String field : new String[] { "Z", "ZT", "Q" }) {
				JsonNode arr = schemes.get(field);
				if (arr == null) continue;
				for (JsonNode entry : arr) {
					if (entry.path("rank").asInt() != rank) continue;
					String src = entry.path("source").asText("");
					String path = entry.path("path").asText("");
					// Prefer source when it's a JSON file; fall back to path otherwise.
					String fetchPath = src.endsWith(".json") ? src : path;
					int prio;
					String tag;
					if (src.contains("/alpha_tensor/")) {
						prio = 1; tag = "alphatensor";
					} else if (src.contains("/known/tensor/")) {
						prio = 5; tag = "perminov-tensor";
					} else if (src.contains("/results/")) {
						prio = 6; tag = "perminov-results";
					} else {
						prio = 9; tag = "perminov-other";
					}
					if (prio < bestPriority) {
						bestPriority = prio;
						bestField = field;
						bestPath = fetchPath;
						bestSourceTag = tag;
					}
				}
			}

			if (bestPath == null) {
				System.out.printf("SKIP (no matching rank in status): %s%n", name);
				skipped++;
				continue;
			}

			String url = PERMINOV_BASE + bestPath;
			System.out.printf("FETCH %s → %s (field=%s, source=%s)%n",
					name, bestPath, bestField, bestSourceTag);

			String body;
			try {
				HttpResponse<String> resp = http.send(
						HttpRequest.newBuilder(URI.create(url)).GET().build(),
						HttpResponse.BodyHandlers.ofString());
				if (resp.statusCode() != 200) {
					System.err.printf("  FAIL HTTP %d%n", resp.statusCode());
					failed++;
					continue;
				}
				body = resp.body();
			} catch (Exception e) {
				System.err.printf("  FAIL fetch: %s%n", e.getMessage());
				failed++;
				continue;
			}

			// Try to parse + verify via SchemeIO. Perminov's JSON format may need
			// translation; for now try the direct path and report.
			NonCubicBilinearAlgorithm alg;
			try {
				alg = SchemeIO.read(body);
			} catch (Exception e) {
				System.err.printf("  FAIL parse: %s%n", e.getMessage());
				failed++;
				continue;
			}

			boolean ok = Verifier.passesRandomMatmulSpotCheck(alg);
			System.out.printf("  Verifier: %s (rank=%d)%n", ok ? "PASS" : "FAIL", alg.r);
			if (!ok) {
				failed++;
				continue;
			}

			if (apply) {
				int adds = Verifier.additionCount(alg);
				Path dir = Path.of("src/main/resources/schemes/section" + Math.max(n, Math.max(mm, p)));
				Files.createDirectories(dir);
				String fieldTag = bestField;  // Q / Z / ZT
				String outName = String.format("%s-%s_%dx%dx%d_r%d_a%d.json",
						bestSourceTag, fieldTag, n, mm, p, alg.r, adds);
				Path out = dir.resolve(outName);
				SchemeIO.write(alg, out.toFile());
				System.out.printf("  WROTE %s%n", out);
			}

			if ("alphatensor".equals(bestSourceTag)) refetchedAt++;
			else refetchedOther++;
		}

		System.out.println();
		System.out.printf("== Summary == True-AT refetched: %d  | Other refetched: %d  | Skipped: %d  | Failed: %d%n",
				refetchedAt, refetchedOther, skipped, failed);
	}

	private static java.util.Set<String> getKnownBrokenAtFiles() throws Exception {
		var f = FieldAwareLookup.class.getDeclaredField("KNOWN_BROKEN_FILES");
		f.setAccessible(true);
		@SuppressWarnings("unchecked")
		var set = (java.util.Set<String>) f.get(null);
		return set;
	}
}
