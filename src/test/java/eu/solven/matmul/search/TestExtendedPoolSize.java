package eu.solven.matmul.search;

import eu.solven.matmul.recombination.BlockSplitSearch;

import static org.assertj.core.api.Assertions.assertThat;

import java.util.List;

import org.junit.jupiter.api.Test;

class TestExtendedPoolSize {

	@Test
	void capAt5() {
		List<BlockSplitSearch.NamedBase> pool = BlockSplitSearch.extendedPool(5);
		assertThat(pool).isNotEmpty();
		// Should include the 8 rootPool entries + every other Leaf NC Z
		// scheme with max(n,m,p) ≤ 5. Expect a meaningful expansion (≥ 20).
		assertThat(pool.size()).isGreaterThanOrEqualTo(20);
		// Every entry respects the cap.
		for (BlockSplitSearch.NamedBase nb : pool) {
			assertThat(Math.max(nb.base().n, Math.max(nb.base().m, nb.base().p)))
					.isLessThanOrEqualTo(5);
		}
		System.out.printf("extendedPool(5) → %d entries%n", pool.size());
		for (BlockSplitSearch.NamedBase nb : pool) {
			System.out.printf("   ⟨%d,%d,%d⟩=%d   %s%n",
					nb.base().n, nb.base().m, nb.base().p, nb.base().r, nb.label());
		}
	}
}
