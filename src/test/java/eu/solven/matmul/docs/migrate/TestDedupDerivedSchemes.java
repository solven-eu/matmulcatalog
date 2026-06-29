package eu.solven.matmul.docs.migrate;

import static org.assertj.core.api.Assertions.assertThat;

import java.nio.file.Files;
import java.nio.file.Path;
import java.util.ArrayList;
import java.util.List;
import java.util.Map;
import java.util.TreeMap;

import org.junit.jupiter.api.Test;

import tools.jackson.databind.JsonNode;

import eu.solven.matmul.catalog.SchemeIO;

/**
 * Regression guard for the stale-duplicate-clutter bug that {@link DedupDerivedSchemes}
 * cleaned: {@code derived/} held SEVERAL dense schemes for the same ⟨n,m,p⟩ at the same
 * rank (the canonical ⟨2,2,7⟩=25 had a clean {@code ConcatCols(…2x2x2…)} twin <em>and</em> a
 * convoluted {@code Project(2x2x3→2x2x2)×3} one). A re-derivation that re-introduces a
 * dominated duplicate must fail here rather than silently bloat the catalog again.
 *
 * <p>Scope matches the cleaner: DENSE derived files only (lineage-only stubs carry no
 * {@code u} and are replay-on-demand, out of scope). Fast — JSON parse only, no verify.
 */
public class TestDedupDerivedSchemes {

	private static final Path DERIVED = Path.of("src/main/resources/schemes/derived");

	@Test
	public void no_dense_derived_shape_rank_bucket_has_a_duplicate() throws Exception {
		// (shape|rank) → dense derived filenames at that exact shape+rank.
		Map<String, List<String>> buckets = new TreeMap<>();
		try (var s = Files.walk(DERIVED)) {
			for (Path p : s.filter(x -> x.toString().endsWith(".json")).sorted().toList()) {
				JsonNode root = SchemeIO.parseJson(p.toFile());
				JsonNode u = root.get("u");
				if (u == null || !u.isArray()) {
					// lineage-only stub — out of scope (the cleaner skips these too).
					continue;
				}
				JsonNode n = root.get("n");
				JsonNode m = root.get("m");
				if (n == null || !n.isArray() || n.size() != 3 || m == null) {
					continue;
				}
				String key = n.get(0).asInt() + "x" + n.get(1).asInt() + "x" + n.get(2).asInt()
						+ "|r" + m.asInt();
				buckets.computeIfAbsent(key, k -> new ArrayList<>()).add(p.getFileName().toString());
			}
		}
		List<String> dups = new ArrayList<>();
		for (var e : buckets.entrySet()) {
			if (e.getValue().size() > 1) {
				dups.add(e.getKey() + " → " + e.getValue());
			}
		}
		assertThat(dups)
				.as("each ⟨n,m,p⟩+rank should keep ONE dense derived witness "
						+ "(run DedupDerivedSchemes --apply); duplicates: %s", dups)
				.isEmpty();
	}
}
