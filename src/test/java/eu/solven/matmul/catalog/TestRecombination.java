package eu.solven.matmul.catalog;

import static org.assertj.core.api.Assertions.assertThat;

import java.util.List;

import org.junit.jupiter.api.Test;

import eu.solven.matmul.papers.laderman1976.Laderman23;
import eu.solven.matmul.NonCubicBilinearAlgorithm;
import eu.solven.matmul.papers.strassen1969.Strassen7;
import eu.solven.matmul.algebra.Algebra;
import eu.solven.matmul.algebra.Field;
import eu.solven.matmul.catalog.Recombination;
import eu.solven.matmul.catalog.Recombination.Result;
import eu.solven.matmul.catalog.Recombination.SotaResolver;

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
		// ⟨4,4,4⟩ recombined via Strassen with each base slot getting 2 target
		// units. Each sub-problem becomes ⟨2,2,2⟩ (rank 7 from catalog).
		// Total = 7 × 7 = 49 — Strassen² baseline.
		NonCubicBilinearAlgorithm strassen = NonCubicBilinearAlgorithm.fromCubic(Strassen7.get());
		Result r = Recombination.recombine(4, 4, 4, strassen, REAL_RING);
		assertThat(r.totalRank).isEqualTo(49);
	}

	@Test
	public void recombine_strassen_for_666_finds_balanced_split() {
		// ⟨6,6,6⟩ via Strassen with each base slot getting 3 target units.
		// Each sub-problem is ⟨3,3,3⟩ — best catalog entry is Laderman = 23.
		// Total = 7 × 23 = 161 (matches Strassen × Laderman composition).
		NonCubicBilinearAlgorithm strassen = NonCubicBilinearAlgorithm.fromCubic(Strassen7.get());
		Result r = Recombination.recombine(6, 6, 6, strassen, REAL_RING);
		assertThat(r.totalRank).isEqualTo(161);
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
