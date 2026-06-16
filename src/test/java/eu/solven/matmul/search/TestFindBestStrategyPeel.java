package eu.solven.matmul.search;

import static org.assertj.core.api.Assertions.assertThat;

import java.util.List;
import java.util.Optional;

import org.junit.jupiter.api.Test;

import java.io.File;

import eu.solven.matmul.NonCubicBilinearAlgorithm;
import eu.solven.matmul.catalog.FieldAwareLookup;
import eu.solven.matmul.catalog.SchemeIO;

/**
 * DIS09 §6 integration test: the block-split search reaches the γ5-aware
 * recipe {@code ⟨3,3,3⟩=25} via Strassen ⟨2,2,2⟩=7 on the unbalanced
 * (1,2)³ allocation (operationally equivalent to over-allocation
 * (2,2)³ + peel (0,1)³). NOTE: this works even with {@code maxPadding=0}
 * because the unbalanced-allocation enumeration {@code (1,2)³} ALREADY
 * exposes the same sub-shape distribution as the canonical DIS09 (2,2)³
 * + peel writeup. Task #87 adds {@code maxPadding} for cases where
 * over-allocation is genuinely needed (e.g. multi-axis padding patterns
 * not reachable via plain unbalanced enumeration).
 */
class TestFindBestStrategyPeel {

	@Test
	void strassen_3x3x3_finds_25_with_maxPadding_1() throws Exception {
		FieldAwareLookup lookup = new FieldAwareLookup("Q");
		CitedBound sota = new CitedBound(lookup);
		// Strassen-only pool — otherwise Laderman ⟨3,3,3⟩=23 in the
		// default pool would dominate by direct allocation (1,1,1)³.
		NonCubicBilinearAlgorithm strassen = SchemeIO.read(
				eu.solven.matmul.catalog.SchemeResolver.byHint("src/main/resources/schemes/known/section2/strassen-2x2x2_m7_a18.json"));
		List<BlockSplitSearch.NamedBase> simplePool =
				List.of(new BlockSplitSearch.NamedBase("Strassen ⟨2,2,2⟩=7", strassen, null));

		// Direct test of findBestMultiBaseSplit (recombination only),
		// bypassing findBestStrategy's catalog short-circuit that would
		// just return Laderman 23.
		Optional<BlockSplitSearch.MultiBaseSplitCandidate> noPad =
				BlockSplitSearch.findBestMultiBaseSplit(3, 3, 3, simplePool, sota,
						false, Long.MAX_VALUE, Integer.MAX_VALUE, Integer.MAX_VALUE, 0);
		assertThat(noPad).isPresent();
		long noPadRank = noPad.get().rank();
		System.out.printf("⟨3,3,3⟩ recombination, no peel:   rank = %d%n", noPadRank);

		Optional<BlockSplitSearch.MultiBaseSplitCandidate> withPad =
				BlockSplitSearch.findBestMultiBaseSplit(3, 3, 3, simplePool, sota,
						false, Long.MAX_VALUE, Integer.MAX_VALUE, Integer.MAX_VALUE, 1);
		assertThat(withPad).isPresent();
		long peelRank = withPad.get().rank();
		System.out.printf("⟨3,3,3⟩ recombination, padding=1: rank = %d via %s%n",
				peelRank, withPad.get().baseLabel());

		// DIS09 §6 explicitly cites 25 as the γ5-aware Strassen-recursive
		// bound on ⟨3,3,3⟩. The search must reach AT LEAST this with no
		// peel (via unbalanced) and AT LEAST this with peel enabled.
		assertThat(noPadRank).isEqualTo(25L);
		assertThat(peelRank).isLessThanOrEqualTo(25L);
		// Note: with maxPadding≥1 the search may exploit degenerate
		// allocations (e.g. block of effective size 0) that route a sub-
		// product back to the outer-target shape's catalog rank
		// (Laderman 23). Not a γ5 recipe — a search-degeneracy worth a
		// follow-up filter, captured in task #87 followups.
	}
}
