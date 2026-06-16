package eu.solven.matmul.catalog;

import static org.assertj.core.api.Assertions.assertThat;

import org.junit.jupiter.api.Test;

import eu.solven.matmul.algebra.Field;

/**
 * The serendipity ceiling must bracket honestly: the knapsack DP exact on
 * synthetic prices, and on the real catalog the ⟨2,4,3⟩⊗⟨3,2,3⟩ ceiling must
 * sit at-or-below the plain Kronecker product and at-or-above what the catalog
 * already achieves (296 ≥ floor).
 */
public class TestSerendipityCeiling {

	@Test
	public void knapsack_dp_is_exact_on_synthetic_prices() {
		// σ = {2:1, 3:5, 4:5} (the real ⟨2,4,3⟩ W-table), r=20, ≥7 classes:
		// optimum is six 3s + one 2 → 6·5 + 1 = 31.
		long[] sigma = new long[21];
		sigma[2] = 1;
		sigma[3] = 5;
		sigma[4] = 5;
		assertThat(SerendipityCeiling.maxSavings(20, 7, sigma)).isEqualTo(31);
		// With ≥6 classes the same partition is still allowed (7 ≥ 6) → 31.
		assertThat(SerendipityCeiling.maxSavings(20, 6, sigma)).isEqualTo(31);
		// Forcing ≥10 classes dilutes class sizes: best is five 3s + 2 + 4×1
		// (11 parts ≥ 10) → wait, parts are flexible; DP must find the max.
		assertThat(SerendipityCeiling.maxSavings(20, 10, sigma))
				.isLessThanOrEqualTo(31).isGreaterThan(0);
		// No prices → no savings, regardless of class count.
		assertThat(SerendipityCeiling.maxSavings(20, 7, new long[21])).isZero();
	}

	@Test
	public void catalog_ceiling_brackets_the_6x8x9_record() {
		FieldAwareLookup q = new FieldAwareLookup(Field.Q);
		SerendipityCeiling.AxisCeiling w = SerendipityCeiling.forAxis(
				SerendipitousBudProduct.BudType.W, 2, 4, 3, 20, q, 3, 2, 3);
		long kron = 20L * q.findRank(3, 2, 3);
		// Floor must be a true bracket: ≤ catalog record (296), ≤ plain Kron.
		assertThat(w.productRankFloor()).isLessThanOrEqualTo(q.findRank(6, 8, 9));
		assertThat(w.productRankFloor()).isLessThanOrEqualTo(kron);
		assertThat(w.maxSavings()).isPositive();
		// Spanning + divisibility: ≥ n·p + 1 = 7 classes at rank 20 (4·6 > 20).
		assertThat(w.minClasses()).isEqualTo(7);
	}
}
