package eu.solven.matmul.search;

import eu.solven.matmul.recombination.BlockSplitSearch;

import static org.assertj.core.api.Assertions.assertThat;

import java.util.List;

import org.junit.jupiter.api.Test;

class TestRootPoolContents {

	@Test
	void poolHasExpectedShapes() {
		List<BlockSplitSearch.NamedBase> pool = BlockSplitSearch.rootPool();
		assertThat(pool).isNotEmpty();
		// rootPool now uses fullCheapOrbit (S₃ shape × axis-flip).
		// Bound: 4 roots × up to 48 orbit variants = 192 entries upper-bound;
		// in practice 30-100 after content dedup. Lower bound is the
		// previous 8 (S₃ shape only); axis-flip adds more.
		// rootPool now spans Strassen + several AT-Z rectangular roots +
		// AT⟨4,4,4⟩=49 + AE⟨5,5,5⟩=93 + Sedoglavic⟨7,7,7⟩=250 plus their
		// shape orbits. Lower-bound asserted only; exact size is a moving
		// target as more historical roots get added.
		assertThat(pool.size()).isGreaterThanOrEqualTo(15);
		System.out.printf("rootPool() → %d entries%n", pool.size());
		// Sanity: every entry has a non-trivial scheme.
		for (BlockSplitSearch.NamedBase nb : pool) {
			assertThat(nb.base()).isNotNull();
			assertThat(nb.base().r).isPositive();
		}
	}
}
