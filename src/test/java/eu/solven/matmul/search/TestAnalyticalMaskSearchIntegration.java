package eu.solven.matmul.search;

import eu.solven.matmul.recombination.BlockSplitSearch;

import static org.assertj.core.api.Assertions.assertThat;

import java.util.List;
import java.util.Optional;

import org.junit.jupiter.api.Test;

import eu.solven.matmul.catalog.FieldAwareLookup;

/**
 * Validates that task #105 (AnalyticalMaskSearch integration into
 * {@link BlockSplitSearch#findBestMultiBaseSplit}) lets a pool with NO
 * axis-flip orbit expansion (RecombinationPoolConfig.simple → CANONICAL orbit mode) still
 * discover the Winograd-mask=1 path to ⟨17,17,17⟩=2930 at (9,8)³.
 *
 * <p>Before #105: RecombinationPoolConfig.simple would stop at 2940 (canonical Winograd
 * mask=0 / Strassen). The analytical mask sweep inside findBestMultiBaseSplit
 * now explores all 8 axis-flip variants per 2×2×2 base at search time.
 */
class TestAnalyticalMaskSearchIntegration {

	@Test
	void poolConfig_simple_finds_2930_at_17x17x17() throws Exception {
		FieldAwareLookup lookup = new FieldAwareLookup("Q");
		CitedBound sota = new CitedBound(lookup);

		// RecombinationPoolConfig.simple uses InternalOrbitMode.CANONICAL — no axis-flip
		// expansion at pool-build time. Pre-#105 this would only have access
		// to canonical Strassen and canonical Winograd, both producing 2940 at
		// (9,8)³. The mask=1 variant of Winograd (= 2930) was unreachable
		// without orbit expansion.
		List<BlockSplitSearch.NamedBase> pool = BlockSplitSearch.buildPool(RecombinationPoolConfig.simple());

		// balancedOnly=false would enumerate the full Cartesian product of
		// allocations (slow); balancedOnly=true restricts to multisets
		// summing to the target which is enough to hit (9,8)³ — the
		// allocation Winograd mask=1 uses to reach 2930.
		Optional<BlockSplitSearch.MultiBaseSplitCandidate> best =
				BlockSplitSearch.findBestMultiBaseSplit(17, 17, 17, pool, sota::getRank, true);

		assertThat(best).isPresent();
		assertThat(best.get().rank())
				.as("RecombinationPoolConfig.simple should now reach 2930 at ⟨17,17,17⟩ via "
						+ "AnalyticalMaskSearch's in-loop mask sweep (task #105)")
				.isLessThanOrEqualTo(2930L);
	}
}
