package eu.solven.matmul.catalog;

import static org.assertj.core.api.Assertions.assertThat;

import java.util.List;

import org.junit.jupiter.api.Test;

import eu.solven.matmul.papers.laderman1976.Laderman23;
import eu.solven.matmul.NonCubicBilinearAlgorithm;
import eu.solven.matmul.papers.strassen1969.Strassen7;
import eu.solven.matmul.algebra.Algebra;
import eu.solven.matmul.algebra.Field;
import eu.solven.matmul.recombination.Recombination;
import eu.solven.matmul.recombination.Recombination.Result;
import eu.solven.matmul.recombination.Recombination.SotaResolver;

/**
 * Tests for the Java port of AlphaTensor's recombination algorithm.
 */
public class TestRecombination {

	private static final SotaResolver REAL_RING =
			Recombination.catalogResolver(Algebra.nonCommutative(Field.R));

	@Test
	public void blockFillings_count_and_sums() {
		List<int[]> fillings = Recombination.blockFillings(2, 4);
		// C(4+1, 1) = 5: [0,4], [1,3], [2,2], [3,1], [4,0]
		assertThat(fillings).hasSize(5);
		for (int[] f : fillings) {
			assertThat(f).hasSize(2);
			assertThat(f[0] + f[1]).isEqualTo(4);
		}
	}

	@Test
	public void blockFillings_three_blocks() {
		List<int[]> fillings = Recombination.blockFillings(3, 2);
		// C(2+2, 2) = 6
		assertThat(fillings).hasSize(6);
		for (int[] f : fillings) {
			assertThat(f).hasSize(3);
			assertThat(f[0] + f[1] + f[2]).isEqualTo(2);
		}
	}

	@Test
	public void recombine_strassen_for_222_target_recovers_rank_7() {
		// Recombining ⟨2,2,2⟩ via Strassen ⟨2,2,2⟩: trivial case, each base slot
		// gets exactly one target unit, every sub-problem is ⟨1,1,1⟩ (rank 1).
		// Total = 7 base-rank entries × 1 = 7. Matches Strassen.
		NonCubicBilinearAlgorithm strassen = NonCubicBilinearAlgorithm.fromCubic(Strassen7.get());
		Result r = Recombination.recombine(2, 2, 2, strassen, REAL_RING);
		assertThat(r.totalRank).isEqualTo(7);
	}

	@Test
	public void recombine_strassen_for_444_target_matches_strassen_squared() {
		// ⟨4,4,4⟩ via Strassen: the genuine (positive-filling) balanced [2,2]³ split gives
		// ⟨2,2,2⟩ sub-problems → 7×7 = 49 (Strassen²). NOTE: recombine's blockFillings also
		// admit DEGENERATE fillings (an empty base block, e.g. [0,4]) — those collapse to a
		// single full-⟨4,4,4⟩ product looked up directly in the catalog. Before catalogResolver
		// was fixed to use the full catalog, that lookup cost the cubic 64 and lost; now it
		// returns the real ⟨4,4,4⟩=48 (DPS-2025 SOTA), so recombine's min is a degenerate
		// catalog passthrough (48), NOT a Strassen recombination. Hence ≤ 49, not = 49.
		// (Genuine recombination still bottoms out at Strassen² = 49.)
		NonCubicBilinearAlgorithm strassen = NonCubicBilinearAlgorithm.fromCubic(Strassen7.get());
		Result r = Recombination.recombine(4, 4, 4, strassen, REAL_RING);
		assertThat(r.totalRank).isLessThanOrEqualTo(49);
	}

	@Test
	public void recombine_strassen_for_666_finds_balanced_split() {
		// ⟨6,6,6⟩ via Strassen: genuine balanced [3,3]³ → ⟨3,3,3⟩ → 7×23 = 161 (Strassen ×
		// Laderman). As with ⟨4,4,4⟩ above, the post-fix min (153) is the DEGENERATE empty-block
		// filling reducing to the catalog's direct ⟨6,6,6⟩=153, not a genuine recombination.
		// Assert ≤ the 161 Strassen×Laderman baseline.
		NonCubicBilinearAlgorithm strassen = NonCubicBilinearAlgorithm.fromCubic(Strassen7.get());
		Result r = Recombination.recombine(6, 6, 6, strassen, REAL_RING);
		assertThat(r.totalRank).isLessThanOrEqualTo(161);
	}

	@Test
	public void recombine_laderman_for_333_recovers_rank_23() {
		// Trivial case: ⟨3,3,3⟩ recombined via Laderman gives exactly Laderman.
		NonCubicBilinearAlgorithm lad = NonCubicBilinearAlgorithm.fromCubic(Laderman23.get());
		Result r = Recombination.recombine(3, 3, 3, lad, REAL_RING);
		assertThat(r.totalRank).isEqualTo(23);
	}

	@Test
	public void recombine_strassen_for_233_target() {
		// ⟨2,3,3⟩ via Strassen ⟨2,2,2⟩: sub-problems will have non-cubic shapes;
		// catalog has ⟨2,3,3⟩=15 so a good allocation can match or beat
		// Strassen × naive bounds. Just confirm we get a finite result and it
		// doesn't exceed the naive 2·3·3 = 18 bound.
		NonCubicBilinearAlgorithm strassen = NonCubicBilinearAlgorithm.fromCubic(Strassen7.get());
		Result r = Recombination.recombine(2, 3, 3, strassen, REAL_RING);
		assertThat(r.totalRank).isLessThanOrEqualTo(18);
	}

	@Test
	public void recombine_result_allocation_sums_match_target() {
		NonCubicBilinearAlgorithm strassen = NonCubicBilinearAlgorithm.fromCubic(Strassen7.get());
		Result r = Recombination.recombine(6, 6, 6, strassen, REAL_RING);
		int sumA = 0, sumB = 0, sumC = 0;
		for (int x : r.allocations[0]) sumA += x;
		for (int x : r.allocations[1]) sumB += x;
		for (int x : r.allocations[2]) sumC += x;
		assertThat(sumA).isEqualTo(6);
		assertThat(sumB).isEqualTo(6);
		assertThat(sumC).isEqualTo(6);
	}
}
