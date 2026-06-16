package eu.solven.matmul.search;

import static org.assertj.core.api.Assertions.assertThat;

import org.junit.jupiter.api.Test;

import eu.solven.matmul.NonCubicBilinearAlgorithm;
import eu.solven.matmul.SymmetryTransforms;

/**
 * Guards the axis-flip dedup invariant: every entry in
 * {@link BlockSplitSearch#rootPool} must be in axis-flip canonical form (the
 * lex-min representative of its axis-flip orbit). Axis-flip variants are
 * expanded at search time via {@link BlockSplitSearch#defaultPool} or via
 * {@link eu.solven.matmul.search.AnalyticalMaskSearch}; keeping them in
 * rootPool would be cost-redundant.
 *
 * <p>Also verifies {@link BlockSplitSearch#defaultPool} restores the
 * coverage by expanding the axis-flip orbit (so callers retain the
 * ⟨17,17,17⟩=2930-path equivalents via Winograd's orbit).
 */
class TestRootPoolCanonicalAxisFlip {

	@Test
	void rootPool_entries_have_distinct_axisFlip_canonical_signatures() {
		// Dedup invariant: no two rootPool entries are axis-flip-equivalent.
		// (We DON'T require each entry to BE its own canonical — the on-disk
		// JSON files can be in any orbit member; what matters is that no two
		// entries share an orbit, because the search-time orbit expansion
		// would then duplicate the work.)
		java.util.Map<String, String> canonicalToLabel = new java.util.HashMap<>();
		for (BlockSplitSearch.NamedBase nb : BlockSplitSearch.rootPool()) {
			NonCubicBilinearAlgorithm alg = nb.base();
			String canonical = SymmetryTransforms.axisFlipCanonicalSignature(alg);
			String shapeKey = alg.n + "x" + alg.m + "x" + alg.p + ":r=" + alg.r + "|" + canonical;
			String prevLabel = canonicalToLabel.put(shapeKey, nb.label());
			assertThat(prevLabel)
					.as("rootPool entry '%s' is axis-flip-equivalent to earlier entry '%s' — "
							+ "one of them is redundant under defaultPool's orbit expansion",
							nb.label(), prevLabel)
					.isNull();
		}
	}

	@Test
	void defaultPool_expands_axisFlip_orbit_so_winograd_orbit_present() {
		// After dedup, rootPool has canonical Winograd-1971 (mask=0). The
		// defaultPool wrapper must expand the orbit so the mask=1 variant
		// (which reaches ⟨17,17,17⟩=2930 at (9,8)³) is reachable by
		// downstream search.
		long winogradOrbitSize = BlockSplitSearch.defaultPool().stream()
				.filter(nb -> nb.label().startsWith("Winograd<2,2,2>=7"))
				.count();
		assertThat(winogradOrbitSize)
				.as("defaultPool should contain Winograd's axis-flip orbit (canonical + flipped variants)")
				.isGreaterThanOrEqualTo(2);
	}

	@Test
	void defaultPool_strictly_larger_than_rootPool() {
		assertThat(BlockSplitSearch.defaultPool().size())
				.as("axis-flip expansion should grow defaultPool beyond rootPool")
				.isGreaterThan(BlockSplitSearch.rootPool().size());
	}

	private static String signatureOf(NonCubicBilinearAlgorithm a) {
		double[][] srcU = a.denseU();
		double[][] srcV = a.denseV();
		double[][] srcW = a.denseW();
		StringBuilder sb = new StringBuilder();
		sb.append(a.n).append(',').append(a.m).append(',').append(a.p).append('|');
		for (double[] row : srcU) { for (double v : row) sb.append(v).append(','); sb.append(';'); }
		sb.append('|');
		for (double[] row : srcV) { for (double v : row) sb.append(v).append(','); sb.append(';'); }
		sb.append('|');
		for (double[] row : srcW) { for (double v : row) sb.append(v).append(','); sb.append(';'); }
		return sb.toString();
	}
}
