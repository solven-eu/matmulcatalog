package eu.solven.matmul.catalog;

import static org.assertj.core.api.Assertions.assertThat;

import java.io.IOException;
import java.nio.file.Path;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.List;
import java.util.Map;
import java.util.regex.Matcher;
import java.util.regex.Pattern;

import org.junit.jupiter.api.Test;

import tools.jackson.databind.JsonNode;
import tools.jackson.databind.json.JsonMapper;

/**
 * Field-discipline guard for Kronecker-derived schemes.
 *
 * <p>A scheme whose lineage is literally {@code Kron(A, B, …)} is the Kronecker
 * product of its factors, so it is valid over a field {@code K} <b>only if every
 * factor is valid over {@code K}</b> — and then its rank is the product of the
 * factor ranks. Therefore a Kron scheme claiming {@code K} must satisfy</p>
 *
 * <pre>  rank ≥ ∏ bestRank_K(factor shape)</pre>
 *
 * <p>over the {@code K}-valid catalog. Violating it means the scheme was built
 * from a factor that is <em>not</em> valid over {@code K} (typically a
 * ½-symmetrisation / {@code 1/8}-rational leaf that is F₃/Q-only, not Z), so the
 * {@code K} tag is an over-claim.</p>
 *
 * <p>The canonical example (caught in the wild, 2026-06): {@code Z⟨6,6,16⟩=385}
 * = {@code ⟨2,2,2⟩=7 ⊗ ⟨3,3,8⟩=55}, but the only {@code ⟨3,3,8⟩=55} is F₃/Q-only
 * (the Z best is 56). So {@code 385 < 7·56 = 392}: a Z scheme of rank 385 cannot
 * exist here, and the genuine {@code Z⟨6,6,16⟩} SOTA is Perminov's 392. The same
 * shape over F₃/Q legitimately holds 385.</p>
 *
 * <p>This is intentionally cheap: a few catalog reads + multiplications, no
 * matrix replay. It guards against the lineage-field-inference over-grant
 * silently re-appearing when {@code derived/} is regenerated.</p>
 */
public class TestKronStubFieldConsistency {

	/** A bare shape factor, e.g. {@code 3x3x8} (the Kron lineage uses U+2297 "⊗"). */
	private static final Pattern SHAPE = Pattern.compile("^(\\d+)x(\\d+)x(\\d+)$");

	/** Fields with a containment lattice where the rank-product bound is sound. */
	private static final List<String> CHECKED_FIELDS = List.of("F2", "F3", "Z");

	@Test
	public void kron_schemes_do_not_overclaim_a_field_their_rank_cannot_support() throws IOException {
		JsonNode root = new JsonMapper().readTree(Path.of("docs/catalog.json").toFile());
		JsonNode schemes = root.get("schemes");

		// bestRank[field][sortedShapeKey] = lowest catalog rank valid over that field.
		Map<String, Map<String, Integer>> bestRank = new HashMap<>();
		for (JsonNode s : schemes) {
			if (s.path("commutative").asBoolean(false)) continue;
			int[] fmt = format(s);
			if (fmt == null || !s.has("rank")) continue;
			int rank = s.get("rank").asInt();
			String key = shapeKey(fmt);
			for (JsonNode f : s.path("fields")) {
				bestRank.computeIfAbsent(f.asString(), k -> new HashMap<>())
						.merge(key, rank, Math::min);
			}
		}

		List<String> violations = new ArrayList<>();
		for (JsonNode s : schemes) {
			if (s.path("commutative").asBoolean(false)) continue;
			int[] fmt = format(s);
			if (fmt == null || !s.has("rank")) continue;
			List<int[]> leaves = cleanKronLeaves(s.path("lineage_compact").asString(""), fmt);
			if (leaves == null) continue; // not a clean Kron chain — nothing to check
			int rank = s.get("rank").asInt();
			for (JsonNode fNode : s.path("fields")) {
				String K = fNode.asString();
				if (!CHECKED_FIELDS.contains(K)) continue;
				Map<String, Integer> best = bestRank.getOrDefault(K, Map.of());
				long product = 1;
				boolean complete = true;
				for (int[] leaf : leaves) {
					Integer b = best.get(shapeKey(leaf));
					if (b == null) { complete = false; break; }
					product *= b;
				}
				if (complete && rank < product) {
					violations.add(String.format("⟨%d,%d,%d⟩=%d claims %s but best-%s factor product = %d (%s)",
							fmt[0], fmt[1], fmt[2], rank, K, K, product, s.path("lineage_compact").asString("")));
				}
			}
		}

		assertThat(violations)
				.as("Kron-derived schemes whose claimed field is impossible at their rank "
						+ "(rank < product of best K-valid factor ranks) — a field over-claim. "
						+ "Correct the scheme's fields[] to the intersection of its factors' fields.")
				.isEmpty();
	}

	private static int[] format(JsonNode s) {
		JsonNode f = s.get("format");
		if (f == null || !f.isArray() || f.size() != 3) return null;
		return new int[] { f.get(0).asInt(), f.get(1).asInt(), f.get(2).asInt() };
	}

	/** Canonical n≤m≤p key so transposed shapes share a best-rank bucket. */
	private static String shapeKey(int[] f) {
		int[] d = { f[0], f[1], f[2] };
		java.util.Arrays.sort(d);
		return d[0] + "x" + d[1] + "x" + d[2];
	}

	/**
	 * Parse a {@code lineage_compact} of the form {@code "AxBxC ⊗ DxExF ⊗ …"} into
	 * its factor shapes, but ONLY when the factors' axiswise product equals
	 * {@code fmt} — i.e. it is a genuine, fully-parsed Kron chain. Returns
	 * {@code null} for anything else (pinned {@code shape@hash} refs, nested
	 * non-Kron ops, or a product that does not reconstruct the shape), so the
	 * guard never fires on a lineage it did not fully understand.
	 */
	private static List<int[]> cleanKronLeaves(String lineageCompact, int[] fmt) {
		if (lineageCompact == null || !lineageCompact.contains("⊗")) return null;
		String[] parts = lineageCompact.split("⊗");
		List<int[]> leaves = new ArrayList<>();
		int[] ax = { 1, 1, 1 };
		for (String p : parts) {
			Matcher m = SHAPE.matcher(p.trim());
			if (!m.matches()) return null;
			int[] leaf = { Integer.parseInt(m.group(1)), Integer.parseInt(m.group(2)),
					Integer.parseInt(m.group(3)) };
			leaves.add(leaf);
			for (int i = 0; i < 3; i++) ax[i] *= leaf[i];
		}
		if (leaves.size() < 2) return null;
		if (ax[0] != fmt[0] || ax[1] != fmt[1] || ax[2] != fmt[2]) return null; // not a clean reconstruction
		return leaves;
	}
}
