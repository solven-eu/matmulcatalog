package eu.solven.matmul.catalog;

import static org.assertj.core.api.Assertions.assertThat;

import org.junit.jupiter.api.Test;

import eu.solven.matmul.algebra.Field;

/**
 * {@code findRank} returns the NAÏVE rank {@code a·b·c} for any shape with no catalog
 * scheme — never the old {@code MAX/100} sentinel. Naïve is always achievable as a
 * <em>formula</em> (the trivial scheme; no dense representation needed), so even shapes
 * far too big to ever be stored have a well-defined rank. The sentinel used to poison
 * sum-based lower bounds ({@code AllocationOptimizer} dropped good bases over a degenerate
 * relaxation block → the ⟨12,13,13⟩=1274 loss); this guard keeps {@code findRank}
 * constructive.
 */
public class TestFindRankNaive {

	@Test
	public void absent_huge_shapes_return_naive_not_sentinel() {
		FieldAwareLookup lookup = new FieldAwareLookup(Field.Q);
		// ⟨1,1234,1⟩ can never be on disk; naïve = 1·1234·1 = 1234.
		assertThat(lookup.findRank(1, 1234, 1)).isEqualTo(1234);
		// ⟨1,1234,5678⟩ likewise; naïve = 1·1234·5678 = 7,006,652 — and NOT the ~21M sentinel.
		assertThat(lookup.findRank(1, 1234, 5678)).isEqualTo(1 * 1234 * 5678);
		assertThat(lookup.findRank(1, 1234, 5678)).isLessThan(Integer.MAX_VALUE / 100);
	}

	@Test
	public void present_shape_returns_catalog_rank_below_naive() {
		FieldAwareLookup lookup = new FieldAwareLookup(Field.Q);
		// ⟨2,2,2⟩ = Strassen 7 (a real catalog scheme), strictly better than naïve 8.
		assertThat(lookup.findRank(2, 2, 2)).isEqualTo(7);
	}
}
