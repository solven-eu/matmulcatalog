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

	/**
	 * Issue #8 (2026-09-30): ⟨2,4,4⟩=26 (Hopcroft–Kerr) must be a ROOT base — the
	 * FMM-Lille recipe ranks for the ⟨10,·,·⟩ / ⟨6,14,25⟩ family are ⟨2,4,4⟩-outer
	 * recombinations over unequal blocks, and before the registration the base was
	 * reachable only via the extended pool or a one-file {@code --base=2x4x4}. Two
	 * CONTENT-distinct representatives are required: their product supports differ
	 * and tile differently (Perminov-ZT rep wins the ⟨10,·,·⟩ family, AlphaTensor-Z
	 * rep wins ⟨6,14,25⟩); collapsing them to one silently loses ranks.
	 */
	@Test
	void poolHasTwoContentDistinct244Roots() {
		java.util.Set<String> hashes = new java.util.HashSet<>();
		for (BlockSplitSearch.NamedBase nb : BlockSplitSearch.rootPool()) {
			var b = nb.base();
			int[] dims = { b.n, b.m, b.p };
			java.util.Arrays.sort(dims);
			if (dims[0] == 2 && dims[1] == 4 && dims[2] == 4 && b.r == 26) {
				// canonical orientation so the two reps' orbit members don't inflate the count
				b.orientAs(2, 4, 4).ifPresent(o -> hashes.add(eu.solven.matmul.catalog.SchemeIO.contentHash(o)));
			}
		}
		assertThat(hashes)
				.as("rootPool must carry ≥ 2 content-distinct ⟨2,4,4⟩=26 bases (Perminov-ZT + AlphaTensor-Z reps)")
				.hasSizeGreaterThanOrEqualTo(2);
	}
}
