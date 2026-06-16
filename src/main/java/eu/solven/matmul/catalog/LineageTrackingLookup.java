package eu.solven.matmul.catalog;

import java.nio.file.Path;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.HashSet;
import java.util.List;
import java.util.Optional;
import java.util.Set;

import eu.solven.matmul.NonCubicBilinearAlgorithm;

/**
 * Wraps a {@link Recombination.AlgorithmLookup} to record each
 * sub-shape consulted during a {@link Recombination#constructWithAllocation}
 * call, so the calling materialiser can attach the resulting
 * {@code List<Lineage.Node>} as the {@code leaves} field of a
 * {@link Lineage.RecombinationN} node.
 *
 * <p>Deduplicated by canonical (sorted) shape key — the same sub-shape
 * consulted multiple times during recombination contributes a single
 * leaf entry, in insertion order. Mirrors the convention used by
 * {@code RecursiveMaterialiser.RecursiveLookup}.</p>
 *
 * <p>Lifecycle: build one per materialiser run, hand it to
 * {@link Recombination#constructWithAllocation}, then read
 * {@link #leaves()} to assemble the lineage tree before writing the
 * scheme via {@link SchemeIO#write(NonCubicBilinearAlgorithm,
 * java.io.File, boolean, Lineage.Node)}.</p>
 */
public final class LineageTrackingLookup implements Recombination.AlgorithmLookup {

	private final FieldAwareLookup lookup;
	private final List<Lineage.Node> leaves = new ArrayList<>();
	private final Set<String> seen = new HashSet<>();

	public LineageTrackingLookup(FieldAwareLookup lookup) {
		this.lookup = lookup;
	}

	@Override
	public Optional<NonCubicBilinearAlgorithm> find(int n, int m, int p) {
		Optional<FieldAwareLookup.WithSource> hit = lookup.findWithSource(n, m, p);
		if (hit.isPresent()) {
			String key = canon(n, m, p);
			if (seen.add(key)) {
				// PIN the PRECISE scheme actually consulted (canon@hash), NOT "key-direct".
				// The old "-direct" marker recorded only "best at shape", a @sota / cited
				// bound — so every recombination leaf became non-bit-exactly-replayable
				// (the inherited bug behind ~14.7k stale "-direct" refs). The exact scheme
				// matters for downstream DCE/projection; "take the best is bad". Fail loud
				// rather than emit a bare/@sota leaf.
				leaves.add(preciseLeaf(key, hit.get().path()));
			}
		}
		return hit.map(FieldAwareLookup.WithSource::alg);
	}

	/** A bit-exact {@code canon@hash7} leaf for the scheme file actually used; throws if
	 *  the file has no stamped content hash (refusing to emit a bare/@sota leaf). */
	private static Lineage.Node preciseLeaf(String canonKey, Path src) {
		String hash = null;
		try {
			hash = SchemeIO.readHash(SchemeIO.parseJson(src.toFile()));
		} catch (Exception e) {
			// fall through to the filename fallback / fail-loud below
		}
		if (hash == null || hash.isBlank()) {
			// Older atoms (e.g. alphatensor imports) encode the hash only in the filename.
			java.util.regex.Matcher m = java.util.regex.Pattern
					.compile("-([0-9a-f]{4,})\\.json$").matcher(src.getFileName().toString());
			if (m.find()) hash = m.group(1);
		}
		if (hash == null || hash.isBlank()) {
			throw new IllegalStateException("LineageTrackingLookup: cannot pin a precise leaf for "
					+ canonKey + " from " + src + " (no content hash in JSON or filename) — refusing"
					+ " to emit a bare/@sota leaf in an explicit recombination.");
		}
		return new Lineage.Atom(canonKey + "@" + hash);
	}

	/** Sub-shape leaves consulted during recombination, in insertion order. */
	public List<Lineage.Node> leaves() {
		return List.copyOf(leaves);
	}

	private static String canon(int n, int m, int p) {
		int[] s = { n, m, p };
		Arrays.sort(s);
		return s[0] + "x" + s[1] + "x" + s[2];
	}
}
