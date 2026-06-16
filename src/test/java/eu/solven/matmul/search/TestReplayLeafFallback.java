package eu.solven.matmul.search;

import static org.assertj.core.api.Assertions.assertThat;

import org.junit.jupiter.api.Test;

import eu.solven.matmul.NonCubicBilinearAlgorithm;
import eu.solven.matmul.algebra.Field;
import eu.solven.matmul.catalog.FieldAwareLookup;
import eu.solven.matmul.catalog.Lineage;

/**
 * The resolveLeaf shape-extraction fallback: source-prefixed / canonical-key /
 * named leaf refs that the sweep surfaced as un-replayable must now resolve by
 * shape (so projection parents stop being skipped). Each resolves to a valid
 * scheme of the embedded shape.
 */
public class TestReplayLeafFallback {

	private final LineageReplayer replayer =
			LineageReplayer.withDefaultPool(new FieldAwareLookup(Field.Q));

	private NonCubicBilinearAlgorithm replay(String ref) {
		return replayer.replay(new Lineage.Atom(ref));
	}

	@Test
	public void source_prefixed_ref_resolves_by_shape() {
		NonCubicBilinearAlgorithm a = replay("alphatensor_Z-2x3x3_m15_a58");
		assertThat(new int[] { a.n, a.m, a.p }).containsExactly(2, 3, 3);
	}

	@Test
	public void perminov_canonical_key_ref_resolves() {
		NonCubicBilinearAlgorithm a = replay("perminov_Z-2x5x15_m118_a2913");
		assertThat(new int[] { a.n, a.m, a.p }).containsExactly(2, 5, 15);
	}

	@Test
	public void angle_named_ref_resolves() {
		NonCubicBilinearAlgorithm a = replay("Strassen<2,2,2>=7");
		assertThat(new int[] { a.n, a.m, a.p }).containsExactly(2, 2, 2);
	}
}
