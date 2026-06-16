package eu.solven.matmul.catalog;

import eu.solven.matmul.search.ConcatSplitSearch;

import static org.assertj.core.api.Assertions.assertThat;

import java.util.Optional;

import org.junit.jupiter.api.Test;

import eu.solven.matmul.NonCubicBilinearAlgorithm;
import eu.solven.matmul.Verifier;

public class TestConcatSplitSearch {

	@Test
	public void search_picks_p10x6_plus_p10x10_for_target_2x10x16() {
		// FMM's own decomposition for ⟨2,10,16⟩=248:
		//     ⟨2,10,6⟩ + ⟨2,10,10⟩
		// Verify our search finds the same (or better) split when fed the
		// catalog's actual SOTA via FieldAwareLookup. Predicted total rank
		// from our catalog: 94 + 155 = 249 (1 above FMM's 248 since we
		// can't reproduce FMM's single-shot HK construction yet — that's
		// task #45).
		FieldAwareLookup lookup = new FieldAwareLookup("Q");
		Recombination.SotaResolver sota = (a, b, c) ->
				lookup.find(a, b, c).map(alg -> alg.r).orElse(Integer.MAX_VALUE / 100);

		Optional<ConcatSplitSearch.ConcatSplit> best =
				ConcatSplitSearch.findBest(2, 10, 16, sota);
		assertThat(best).isPresent();
		// Best should be a p-axis split (axis=2). Total ≤ 249.
		assertThat(best.get().axis()).isEqualTo(2);
		assertThat(best.get().totalRank()).isLessThanOrEqualTo(249);

		// Materialise + spot-check.
		NonCubicBilinearAlgorithm scheme = ConcatSplitSearch.materialise(best.get(), lookup);
		assertThat(scheme.n).isEqualTo(2);
		assertThat(scheme.m).isEqualTo(10);
		assertThat(scheme.p).isEqualTo(16);
		assertThat((long) scheme.r).isEqualTo(best.get().totalRank());
		assertThat(Verifier.passesRandomMatmulSpotCheck(scheme)).isTrue();

		System.out.printf("⟨2,10,16⟩ via concat: split=[%d, %d] on axis %d, rank=%d%n",
				best.get().leftSize(), best.get().rightSize(),
				best.get().axis(), best.get().totalRank());
	}
}
