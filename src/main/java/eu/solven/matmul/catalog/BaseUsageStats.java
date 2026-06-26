package eu.solven.matmul.catalog;

import java.util.LinkedHashMap;
import java.util.Map;
import java.util.regex.Matcher;
import java.util.regex.Pattern;

/**
 * Accumulates, across the whole catalog, how many times each base scheme is actually
 * USED as a building block by other schemes — broken down by the using lineage op
 * (recombination / projection / kronecker / concat / …). Answers "which bases earn
 * their keep": a base referenced by 200 recombinations is load-bearing; one referenced
 * by none is a dead end.
 *
 * <p>Usage is read from the lineage DAG: each {@link Lineage.Atom} reference is
 * attributed to the NEAREST ENCLOSING COMPOSITION op (transform wrappers —
 * Transpose / OrientAs / AxisFlip / AxisPermute / Dce — pass the context through, since
 * they only re-orient a base before a real composition consumes it). Keys are the
 * referenced base's {@code shape@hash7} (pinned ref) or bare {@code shape} (catalog-best
 * ref); primitive {@code naive-}/{@code direct-} leaves are not catalog bases and are
 * skipped.</p>
 */
public final class BaseUsageStats {

	/** key (shape@hash7 | shape) → op category → count. Insertion-ordered for stable JSON. */
	private final Map<String, Map<String, Integer>> counts = new LinkedHashMap<>();

	/** {@code "12x13x14@1a2b3c4"} or a filename stem ending in {@code -1a2b3c4}. */
	private static final Pattern HASH_REF =
			Pattern.compile("(\\d+x\\d+x\\d+)[@-].*?([0-9a-f]{7})[0-9a-f]*$");
	/** A bare {@code "12x13x14"} shape ref (no hash / source decoration). */
	private static final Pattern BARE_SHAPE = Pattern.compile("^(\\d+x\\d+x\\d+)$");

	/** The op category a node attributes its base references to, or {@code null} for a
	 *  transform wrapper / leaf (which carries no usage of its own). */
	public static String category(Lineage.Node n) {
		return switch (n) {
			case Lineage.RecombinationN x -> "recombination";
			case Lineage.RecombinationTaN x -> "recombination";
			case Lineage.RecombinationWithPairN x -> "recombination";
			case Lineage.Project x -> "projection";
			case Lineage.KronProduct x -> "kronecker";
			case Lineage.KronChain x -> "kronecker";
			case Lineage.ConcatCols x -> "concat";
			case Lineage.ConcatRows x -> "concat";
			case Lineage.SumInner x -> "sum";
			case Lineage.DisjointSum x -> "disjoint_sum";
			case Lineage.SerendipitousProduct x -> "serendipitous";
			case Lineage.AugmentSquareDiscard x -> "augment";
			case Lineage.PeeledViaTa x -> "recombination";
			default -> null;   // Atom / Transpose / OrientAs / AxisFlip / AxisPermute / Dce
		};
	}

	/** Walk one scheme's lineage, attributing every base reference it uses. */
	public void accumulate(Lineage.Node root) {
		walk(root, null);
	}

	private void walk(Lineage.Node n, String usingOp) {
		if (n instanceof Lineage.Atom a) {
			String key = baseKey(a.ref());
			if (key != null && usingOp != null) {
				counts.computeIfAbsent(key, k -> new LinkedHashMap<>())
						.merge(usingOp, 1, Integer::sum);
			}
			return;
		}
		String cat = category(n);
		String childOp = cat != null ? cat : usingOp;   // transforms pass the context through
		for (Lineage.Node c : Lineage.childrenOf(n)) {
			walk(c, childOp);
		}
	}

	/** Normalise an Atom ref to a usage key, or {@code null} if it is not a catalog base
	 *  (a {@code naive-}/{@code direct-} primitive, or an internal lineage-id ref). */
	static String baseKey(String ref) {
		if (ref == null || ref.startsWith("naive-") || ref.startsWith("direct-")
				|| ref.startsWith("@ref")) {
			return null;
		}
		Matcher h = HASH_REF.matcher(ref);
		if (h.matches()) {
			return h.group(1) + "@" + h.group(2);
		}
		Matcher b = BARE_SHAPE.matcher(ref);
		if (b.matches()) {
			return b.group(1);
		}
		return null;
	}

	/** Per-op usage counts for a base key, or empty if it is never used. */
	public Map<String, Integer> forKey(String key) {
		return counts.getOrDefault(key, Map.of());
	}

	/** All keys seen (for diagnostics). */
	public Map<String, Map<String, Integer>> all() {
		return counts;
	}
}
