package eu.solven.matmul.docs;

import java.io.FileReader;
import java.io.Reader;
import java.util.ArrayList;
import java.util.List;
import java.util.regex.Matcher;
import java.util.regex.Pattern;

import eu.solven.matmul.catalog.FieldAwareLookup;
import tools.jackson.databind.JsonNode;
import tools.jackson.databind.json.JsonMapper;

/**
 * For each file in {@link FieldAwareLookup} known-broken denylist, look up
 * Perminov's {@code status.json} and report:
 *
 * <ul>
 *   <li>What FIELD (Q, Z, ZT) the {@code (format, rank)} tuple is actually
 *       attainable in.</li>
 *   <li>What SOURCE Perminov attributes that rank to (alpha_tensor, tensor,
 *       fmm_add_reduction, etc.).</li>
 * </ul>
 *
 * <p>This tells us, for each broken AT-Z file in our local catalog,
 * whether the rank IS in fact an AlphaTensor discovery (so we should
 * re-fetch the actual factor matrices) or if our local file is just
 * mis-attributed (in which case the denylist entry is correct and
 * permanent).</p>
 *
 * <p>Requires {@code /tmp/perminov-status.json} to be present (fetched
 * via {@code curl} from the Perminov GitHub repo).</p>
 */
public final class AuditBrokenAtSchemes {

	private static final String STATUS_JSON = "/tmp/perminov-status.json";

	private static final Pattern FILENAME = Pattern.compile(
			"alphatensor-Z_(\\d+)x(\\d+)x(\\d+)_r(\\d+)_a\\d+\\.json");

	public static void main(String[] args) throws Exception {
		JsonMapper mapper = new JsonMapper();
		JsonNode root;
		try (Reader r = new FileReader(STATUS_JSON)) {
			root = mapper.readTree(r);
		}

		List<String> brokenFiles = new ArrayList<>(getKnownBrokenAtFiles());
		int trueAt = 0, notAt = 0, missing = 0;
		System.out.printf("%-60s %-8s %-8s %s%n", "filename", "field", "complexity", "Perminov-source");
		System.out.println("-".repeat(120));
		for (String name : brokenFiles) {
			Matcher m = FILENAME.matcher(name);
			if (!m.matches()) {
				System.out.println("SKIP (filename parse): " + name);
				continue;
			}
			int n = Integer.parseInt(m.group(1));
			int mm = Integer.parseInt(m.group(2));
			int p = Integer.parseInt(m.group(3));
			int rank = Integer.parseInt(m.group(4));
			String fmtKey = n + "x" + mm + "x" + p;

			JsonNode formatNode = root.get(fmtKey);
			if (formatNode == null) {
				System.out.printf("%-60s MISSING in status.json%n", name);
				missing++;
				continue;
			}

			JsonNode schemes = formatNode.get("schemes");
			if (schemes == null) {
				System.out.printf("%-60s NO schemes in status.json%n", name);
				missing++;
				continue;
			}

			boolean found = false;
			boolean isAt = false;
			for (String field : new String[] { "Z", "ZT", "Q" }) {
				JsonNode arr = schemes.get(field);
				if (arr == null) continue;
				for (JsonNode entry : arr) {
					if (entry.path("rank").asInt() != rank) continue;
					String source = entry.path("source").asText("");
					int complexity = entry.path("complexity").asInt(-1);
					boolean hereIsAt = source.contains("/alpha_tensor/");
					System.out.printf("%-60s %-8s %-8d %s%n",
							name, field, complexity, source);
					found = true;
					if (hereIsAt) isAt = true;
					break;
				}
			}
			if (!found) {
				System.out.printf("%-60s NOT FOUND at rank=%d%n", name, rank);
				missing++;
			} else if (isAt) {
				trueAt++;
			} else {
				notAt++;
			}
		}
		System.out.println();
		System.out.printf("== Summary == AT-attributable: %d  | NOT-AT (mis-attribution): %d  | MISSING in Perminov: %d%n",
				trueAt, notAt, missing);
	}

	/** Read the denylist from FieldAwareLookup via reflection. */
	private static java.util.Set<String> getKnownBrokenAtFiles() throws Exception {
		var f = FieldAwareLookup.class.getDeclaredField("KNOWN_BROKEN_FILES");
		f.setAccessible(true);
		@SuppressWarnings("unchecked")
		var set = (java.util.Set<String>) f.get(null);
		return set;
	}
}
