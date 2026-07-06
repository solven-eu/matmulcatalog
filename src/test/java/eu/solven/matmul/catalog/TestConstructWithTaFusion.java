package eu.solven.matmul.catalog;

import eu.solven.matmul.recombination.Recombination;

import static org.assertj.core.api.Assertions.assertThat;

import org.junit.jupiter.api.BeforeAll;
import org.junit.jupiter.api.Test;

import eu.solven.matmul.NonCubicBilinearAlgorithm;
import eu.solven.matmul.verifiers.Verifier;

/**
 * Guards {@link Recombination#constructWithTaFusion} — generic Pan-TA fusion over an
 * ARBITRARY naïve-grid base (not just the ⟨1,2,2⟩ peel). Each disjoint cyclic-rotation
 * single-block product pair is fused at {@code fusedRank} instead of the two leaves'
 * summed rank, and the result is a CORRECT matmul.
 */
public class TestConstructWithTaFusion {

	private static FieldAwareLookup lk;
	private static Recombination.SotaResolver sota;

	@BeforeAll
	static void setUp() {
		lk = new FieldAwareLookup("R");
		sota = (a, b, c) -> lk.findRank(a, b, c);
	}

	/** find()-based sub resolver (materialisable leaves only). */
	private static NonCubicBilinearAlgorithm resolve(int[] sz) {
		return lk.find(sz[0], sz[1], sz[2]).orElse(null);
	}

	@Test
	public void unprofitable_peel_pair_is_refused_and_the_fallback_verifies() {
		// ⟨2,3,3⟩ via the ⟨1,2,2⟩ peel. Its cross-pair ⟨2,2,1⟩+⟨2,1,2⟩ is a cyclic
		// rotation GEOMETRICALLY, but fusing it is NEVER profitable: fusedRank(2,1,2)
		// = nrp+np+nr+rp = 12 ≥ 2·naive = 8 (a peel pair has shape ⟨a,k,a⟩-family,
		// whose correction terms always dominate). Historically this test expected 1
		// fused pair — an accident of the UNKNOWN_RANK sentinel (leaves priced as
		// infinite made ANY fusion look like a saving); the findRank naive fallback
		// exposed it. The honest behaviour: the economics gate refuses, and the
		// unfused fallback construction still computes matmul exactly.
		Recombination.TaFusedConstruction tc = Recombination.constructWithTaFusion(
				NonCubicBilinearAlgorithm.naive(1, 2, 2), TestConstructWithTaFusion::resolve, sota,
				new int[] { 2 }, new int[] { 2, 1 }, new int[] { 2, 1 });
		assertThat(tc.fusedPairs()).as("an unprofitable pair must NOT be fused").isEmpty();
		assertThat(tc.alg().n).isEqualTo(2);
		assertThat(tc.alg().m).isEqualTo(3);
		assertThat(tc.alg().p).isEqualTo(3);
		assertThat(Verifier.isExactNonCubic(tc.alg()))
				.as("the unfused fallback ⟨2,3,3⟩ must compute matmul exactly").isTrue();
	}

	@Test
	public void forced_fusion_is_bit_exact_on_the_small_peel() {
		// The bit-exactness guard the old test intended: exercise embedTaPair on the
		// smallest full-tensor-verifiable instance. Fusion only fires when priced
		// profitable, so steer the sota to price the cross leaves as if no fast
		// scheme existed (999 ≫ fusedRank 12) — the TA construction must be correct
		// whenever elected, whatever the economics that elected it.
		Recombination.SotaResolver expensiveLeaves = (a, b, c) ->
				(a * b * c == 4 && Math.min(a, Math.min(b, c)) == 1) ? 999 : sota.getRank(a, b, c);
		Recombination.TaFusedConstruction tc = Recombination.constructWithTaFusion(
				NonCubicBilinearAlgorithm.naive(1, 2, 2), TestConstructWithTaFusion::resolve,
				expensiveLeaves,
				new int[] { 2 }, new int[] { 2, 1 }, new int[] { 2, 1 });
		assertThat(tc.fusedPairs()).as("the steered peel fuses exactly one cross-pair").hasSize(1);
		assertThat(Verifier.isExactNonCubic(tc.alg()))
				.as("the TA-fused ⟨2,3,3⟩ must compute matmul exactly").isTrue();
	}

	@Test
	public void ta_fusion_works_over_a_2x3x3_grid_with_multiple_fusions() {
		// FMM's ⟨22,28,28⟩ recipe grid: 22→[12,10], 28→[10,9,9]². A ⟨2,3,3⟩ naïve grid
		// (NOT a ⟨1,2,2⟩ peel) with TWO disjoint cyclic-rotation cross-pairs — the generic
		// claim. Leaves are materialisable, so the find()-resolver suffices.
		Recombination.TaFusedConstruction tc = Recombination.constructWithTaFusion(
				NonCubicBilinearAlgorithm.naive(2, 3, 3), TestConstructWithTaFusion::resolve, sota,
				new int[] { 12, 10 }, new int[] { 10, 9, 9 }, new int[] { 10, 9, 9 });
		assertThat(tc.fusedPairs())
				.as("the ⟨2,3,3⟩ grid fuses two cyclic cross-pairs (FMM uses 2 TA legs)").hasSize(2);
		assertThat(tc.alg().n).isEqualTo(22);
		assertThat(tc.alg().m).isEqualTo(28);
		assertThat(tc.alg().p).isEqualTo(28);
		assertThat(Verifier.passesRandomMatmulSpotCheck(tc.alg()))
				.as("the TA-fused ⟨22,28,28⟩ must verify").isTrue();
	}

	@Test
	public void isNaiveGrid_distinguishes_grids_from_block_combining_bases() {
		assertThat(Recombination.isNaiveGrid(NonCubicBilinearAlgorithm.naive(1, 2, 2))).isTrue();
		assertThat(Recombination.isNaiveGrid(NonCubicBilinearAlgorithm.naive(2, 3, 3))).isTrue();
		// Strassen ⟨2,2,2⟩=7 combines blocks (M1=(A11+A22)(B11+B22)) — NOT a naïve grid.
		NonCubicBilinearAlgorithm strassen = lk.find(2, 2, 2).orElseThrow();
		assertThat(Recombination.isNaiveGrid(strassen)).isFalse();
	}
}
