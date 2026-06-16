package eu.solven.matmul.docs.explore;

import static org.assertj.core.api.Assertions.assertThat;

import java.util.List;

import org.junit.jupiter.api.Test;

public class TestScanSerendipityFrontier {

	/** σ_W of ⟨3,2,3⟩-inner: {2:1, 3:5, 4:5}. Needing 5 must prefer one triple
	 *  (burden 2) over five pairs (burden 5). */
	@Test
	public void cheapest_profile_prefers_low_burden() {
		long[] sigma = { 0, 0, 1, 5, 5 };
		List<int[]> out = ScanSerendipityFrontier.cheapestProfiles(sigma, 20, 7, 5, 2);
		assertThat(out).isNotEmpty();
		assertThat(out.get(0)).containsExactly(3);
	}

	/** Class-count constraint: at rank 8 with minClasses 7, only one pair fits
	 *  (burden 1 → 7 classes); a triple would leave 6 classes. */
	@Test
	public void class_count_constraint_filters() {
		long[] sigma = { 0, 0, 1, 5 };
		List<int[]> out = ScanSerendipityFrontier.cheapestProfiles(sigma, 8, 7, 1, 5);
		assertThat(out).hasSize(1);
		assertThat(out.get(0)).containsExactly(2);
	}

	/** Unreachable needed → empty, never a silently weaker profile. */
	@Test
	public void unreachable_needed_returns_empty() {
		long[] sigma = { 0, 0, 1 };
		assertThat(ScanSerendipityFrontier.cheapestProfiles(sigma, 6, 1, 100, 3)).isEmpty();
	}
}
