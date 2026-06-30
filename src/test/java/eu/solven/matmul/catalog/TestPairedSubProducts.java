package eu.solven.matmul.catalog;

import static org.assertj.core.api.Assertions.assertThat;

import org.junit.jupiter.api.Test;

import eu.solven.matmul.isotropy.PairedSubProducts;
import eu.solven.matmul.recombination.Recombination;

public class TestPairedSubProducts {

	@Test
	public void cyclic_equivalence_detects_three_rotations() {
		int[] base = { 3, 4, 5 };
		assertThat(PairedSubProducts.cyclicallyEquivalent(base, new int[] { 3, 4, 5 })).isTrue();
		assertThat(PairedSubProducts.cyclicallyEquivalent(base, new int[] { 4, 5, 3 })).isTrue();
		assertThat(PairedSubProducts.cyclicallyEquivalent(base, new int[] { 5, 3, 4 })).isTrue();
		// Reflections (not cyclic) — Pan's formula requires cyclic, not full S₃
		assertThat(PairedSubProducts.cyclicallyEquivalent(base, new int[] { 3, 5, 4 })).isFalse();
		assertThat(PairedSubProducts.cyclicallyEquivalent(base, new int[] { 5, 4, 3 })).isFalse();
	}

	@Test
	public void panPairable_uses_transpose_so_two_equal_dims_fuse_but_all_distinct_dont() {
		// cubic: always pairable.
		assertThat(PairedSubProducts.panPairable(new int[] { 9, 9, 9 }, new int[] { 9, 9, 9 })).isTrue();
		// genuine cyclic rotation: pairable.
		assertThat(PairedSubProducts.panPairable(new int[] { 8, 9, 9 }, new int[] { 9, 9, 8 })).isTrue();
		// TWO IDENTICAL, exactly two equal dims: pairable VIA TRANSPOSE
		// (⟨8,9,9⟩ᵀ = ⟨9,9,8⟩ is the rotation).
		assertThat(PairedSubProducts.panPairable(new int[] { 8, 9, 9 }, new int[] { 8, 9, 9 })).isTrue();
		// TWO IDENTICAL, all dims distinct: NOT pairable even with transpose
		// (⟨8,9,10⟩ᵀ = ⟨10,9,8⟩, still not a rotation of ⟨8,9,10⟩).
		assertThat(PairedSubProducts.panPairable(new int[] { 8, 9, 10 }, new int[] { 8, 9, 10 })).isFalse();
		// two transposes of a distinct shape = two identical distinct → NOT pairable.
		assertThat(PairedSubProducts.panPairable(new int[] { 8, 9, 10 }, new int[] { 10, 9, 8 })).isFalse();
	}

	@Test
	public void pair_cost_reproduces_dis09_formula() {
		// Pan's pair-product cost: abc + ab + bc + ca
		assertThat(PairedSubProducts.pairCost(3, 4, 5))
				.isEqualTo(3 * 4 * 5L + 3 * 4L + 4 * 5L + 5 * 3L);
		assertThat(PairedSubProducts.pairCost(2, 2, 2)).isEqualTo(8L + 4L + 4L + 4L); // 20
	}

	@Test
	public void greedy_pairing_finds_savings_when_naive_loses() {
		// Two ⟨6,6,6⟩ sub-products: naive cost = 2·216 = 432, pair cost = 216+36+36+36 = 324
		// (savings = 108, assuming sota gives the naive rank n³).
		int[][] shapes = { { 6, 6, 6 }, { 6, 6, 6 } };
		Recombination.SotaResolver naive = (a, b, c) -> a * b * c; // worst-case
		long result = PairedSubProducts.applyPairing(shapes, naive);
		assertThat(result).isEqualTo(324L); // = pairCost(6,6,6)
	}

	@Test
	public void no_pairing_when_catalog_rank_already_below_pair_cost() {
		// Two ⟨2,2,2⟩ products. Pair cost = 8+4+4+4=20. Strassen says R=7 each → 14 total.
		// 14 < 20 so do NOT pair.
		int[][] shapes = { { 2, 2, 2 }, { 2, 2, 2 } };
		Recombination.SotaResolver strassen = (a, b, c) -> 7;
		long result = PairedSubProducts.applyPairing(shapes, strassen);
		assertThat(result).isEqualTo(14L);
	}

	@Test
	public void mixed_pool_pairs_only_cyclic_matches() {
		// ⟨3,4,5⟩, ⟨4,5,3⟩ (cyclic pair, naive 2·60=120, pair = 60+12+20+15=107) → SAVE 13
		// + a lonely ⟨2,2,2⟩ with no partner.
		int[][] shapes = { { 3, 4, 5 }, { 4, 5, 3 }, { 2, 2, 2 } };
		Recombination.SotaResolver naive = (a, b, c) -> a * b * c;
		long result = PairedSubProducts.applyPairing(shapes, naive);
		assertThat(result).isEqualTo(107L + 8L);
	}
}
