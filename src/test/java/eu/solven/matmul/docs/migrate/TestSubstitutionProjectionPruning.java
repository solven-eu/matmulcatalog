package eu.solven.matmul.docs.migrate;

import static org.assertj.core.api.Assertions.assertThat;

import org.junit.jupiter.api.Test;

/**
 * The parent-pruning rule of {@link MaterialiseSubstitutionProjections} (fast, no catalog).
 *
 * <p>The substitution lemma saves a <em>specific</em> dimension per dropped axis — an index of
 * {@code p} kills {@code m} products, of {@code m} kills {@code n}, of {@code n} kills {@code p}.
 * Until 2026-11-25 the rule estimated the saving by the <em>smaller</em> of the two other
 * dimensions and skipped a parent that could not win with twice that: every parent whose real
 * saving exceeded twice the minimum was dropped without a log line. ⟨6,24,30⟩:2520 was never
 * tried for ⟨6,24,29⟩ (ours 2506; the lemma gives 2520 − 24 = 2496, FMM's rank), which left
 * four rows of the 2026-11-24 gap list open after a closure run that reported "no parent drops
 * below it".</p>
 */
public class TestSubstitutionProjectionPruning {

	@Test
	public void a_parent_whose_lemma_saving_wins_is_never_pruned() {
		// ⟨6,24,30⟩:2520 → ⟨6,24,29⟩: p-drop kills m = 24 → 2496 < 2506.
		assertThat(MaterialiseSubstitutionProjections.hopeless(new int[] { 6, 24, 30 }, 0, 0, 1, 2520, 2506))
				.as("the p-drop saves m = 24, not min(n, m) = 6").isFalse();
		// ⟨3,20,20⟩-like: an n-drop kills p products, an m-drop kills n.
		assertThat(MaterialiseSubstitutionProjections.hopeless(new int[] { 4, 20, 20 }, 1, 0, 0, 1000, 985)).isFalse();
		assertThat(MaterialiseSubstitutionProjections.hopeless(new int[] { 20, 5, 4 }, 0, 1, 0, 1000, 985)).isFalse();
	}

	/** Exhaustive over small shapes: whenever the exact lemma saving of a single drop wins, the parent is kept. */
	@Test
	public void pruning_is_never_tighter_than_the_lemma_for_any_single_drop() {
		for (int n = 1; n <= 12; n++) {
			for (int m = 1; m <= 12; m++) {
				for (int p = 1; p <= 12; p++) {
					int[] parent = { n, m, p };
					int pr = 10_000;
					// lemma: n-drop saves p, m-drop saves n, p-drop saves m
					int[][] drops = { { 1, 0, 0, p }, { 0, 1, 0, n }, { 0, 0, 1, m } };
					for (int[] d : drops) {
						long upper = pr - d[3] + 1L; // the lemma's result is strictly below upper
						assertThat(MaterialiseSubstitutionProjections.hopeless(parent, d[0], d[1], d[2], pr, upper))
								.as("parent <%d,%d,%d>, drop %d%d%d", n, m, p, d[0], d[1], d[2]).isFalse();
					}
				}
			}
		}
	}

	@Test
	public void a_parent_far_above_the_target_is_still_pruned_and_no_target_rank_means_no_pruning() {
		assertThat(MaterialiseSubstitutionProjections.hopeless(new int[] { 6, 24, 30 }, 0, 0, 1, 2600, 2506)).isTrue();
		assertThat(MaterialiseSubstitutionProjections.hopeless(new int[] { 6, 24, 30 }, 0, 0, 1, 2600, Long.MAX_VALUE))
				.isFalse();
	}
}
