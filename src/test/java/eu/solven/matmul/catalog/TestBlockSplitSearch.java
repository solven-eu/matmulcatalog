package eu.solven.matmul.catalog;

import static org.assertj.core.api.Assertions.assertThat;

import java.util.Map;
import java.util.Optional;
import java.util.function.Function;

import org.junit.jupiter.api.Test;

import eu.solven.matmul.search.BlockSplitSearch;
import eu.solven.matmul.search.BlockSplitSearch.SplitCandidate;

/**
 * Phase 1 search validation. The {@code ⟨7,7,7⟩} target should pick
 * the {@code 4+3} split because it minimises the formula RHS.
 */
public class TestBlockSplitSearch {

	@Test
	public void best_split_777_is_4_plus_3() {
		Map<String, Integer> ranks = BlockSplitSearch.loadCatalogBestRanks();
		Function<int[], Optional<Integer>> lookup = BlockSplitSearch.rankLookupFromMap(ranks);

		Optional<SplitCandidate> best = BlockSplitSearch.findBestSplit(7, lookup);
		assertThat(best).isPresent();
		SplitCandidate c = best.get();
		assertThat(c.u()).isEqualTo(4);
		assertThat(c.v()).isEqualTo(3);
		// 49 (Strassen²) + 3·29 (⟨3,3,4⟩) + 3·38 (⟨3,4,4⟩) = 250 — matches Sedoglavic 2017.
		// Catalog gained better leaves (DPS 2025 ⟨4,4,4⟩=48 / AlphaTensor F₂=47,
		// Perminov refinements) so the formula now lands as low as 242. Track
		// monotone-non-increasing with a loose bound.
		assertThat(c.formulaRank()).isLessThanOrEqualTo(250);
	}

	@Test
	public void best_split_555_uses_a_valid_split() {
		Map<String, Integer> ranks = BlockSplitSearch.loadCatalogBestRanks();
		Function<int[], Optional<Integer>> lookup = BlockSplitSearch.rankLookupFromMap(ranks);

		Optional<SplitCandidate> best = BlockSplitSearch.findBestSplit(5, lookup);
		assertThat(best).isPresent();
		// Whatever the catalog picks, the predicted rank should be > direct (AlphaEvolve 93)
		// because formula doesn't yet apply the Sedoglavic algebraic saving and ⟨5,5,5⟩ has
		// a strong direct scheme.
		Optional<Integer> direct = lookup.apply(new int[] { 5, 5, 5 });
		assertThat(direct).isPresent();
		assertThat(best.get().formulaRank()).isGreaterThan(direct.get());
	}
}
