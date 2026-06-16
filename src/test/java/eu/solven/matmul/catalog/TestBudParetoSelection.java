package eu.solven.matmul.catalog;

import static org.assertj.core.api.Assertions.assertThat;

import org.junit.jupiter.api.Test;

/**
 * The (rank, buds) Pareto registration: a bud-richer higher-rank scheme must
 * stay on the frontier (it's a better serendipitous building block), while a
 * higher-rank scheme that is no bud-richer is dominated.
 */
public class TestBudParetoSelection {

	@Test
	public void budScore_parses_summary_strings() {
		assertThat(BudParetoSelection.budScore("4×U⟨1,1,2⟩")).isEqualTo(8);
		assertThat(BudParetoSelection.budScore("2×U⟨1,1,3⟩ + 1×W⟨1,2,1⟩")).isEqualTo(6 + 2);
		assertThat(BudParetoSelection.budScore("")).isZero();
		assertThat(BudParetoSelection.budScore(null)).isZero();
	}

	@Test
	public void strassen7_and_naive8_are_both_on_frontier() {
		// Strassen ⟨2,2,2⟩=7 (0 buds) vs naive ⟨2,2,2⟩=8 (4 U-buds, score 8):
		// neither dominates — lower rank vs richer buds → both registered.
		boolean[] f = BudParetoSelection.frontierMask(new int[] { 7, 8 }, new int[] { 0, 8 });
		assertThat(f).containsExactly(true, true);
	}

	@Test
	public void higher_rank_without_more_buds_is_dominated() {
		// rank 9 / 0 buds is strictly worse than rank 7 / 0 buds → off frontier.
		boolean[] f = BudParetoSelection.frontierMask(new int[] { 7, 8, 9 }, new int[] { 0, 8, 0 });
		assertThat(f).containsExactly(true, true, false);
	}

	@Test
	public void richer_buds_at_each_higher_rank_all_survive() {
		boolean[] f = BudParetoSelection.frontierMask(new int[] { 7, 8, 9 }, new int[] { 0, 8, 12 });
		assertThat(f).containsExactly(true, true, true);
	}

	@Test
	public void at_equal_rank_the_bud_richer_dominates() {
		boolean[] f = BudParetoSelection.frontierMask(new int[] { 7, 7 }, new int[] { 0, 5 });
		assertThat(f).containsExactly(false, true);
	}

	@Test
	public void co_optimal_ties_both_kept() {
		boolean[] f = BudParetoSelection.frontierMask(new int[] { 7, 7 }, new int[] { 5, 5 });
		assertThat(f).containsExactly(true, true);
	}
}
