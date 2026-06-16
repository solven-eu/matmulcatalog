package eu.solven.matmul.catalog;

import eu.solven.matmul.papers.pan1978.PanPairProduct;

import static org.assertj.core.api.Assertions.assertThat;

import org.junit.jupiter.api.Test;

public class TestPanPairProduct {

	@Test
	public void rank_formula_matches_paired_sub_products_cost() {
		// Pan pair: abc + ab + bc + ca; same as PairedSubProducts.pairCost.
		assertThat(PanPairProduct.rank(2, 2, 2)).isEqualTo(20);
		assertThat(PanPairProduct.rank(3, 3, 3)).isEqualTo(54);
		assertThat(PanPairProduct.rank(7, 7, 7)).isEqualTo(490);
		assertThat(PanPairProduct.rank(2, 3, 4))
				.isEqualTo(2 * 3 * 4 + 2 * 3 + 3 * 4 + 2 * 4);

		for (int n : new int[] { 2, 3, 5, 7, 11 }) {
			assertThat((long) PanPairProduct.rank(n, n, n))
					.isEqualTo(PairedSubProducts.pairCost(n, n, n));
		}
	}

	@Test
	public void build_222_passes_random_spot_check() {
		PanPairProduct.PairScheme s = PanPairProduct.build(2, 2, 2);
		assertThat(s.rank()).isEqualTo(20);
		assertThat(PanPairProduct.spotCheck(s, 5, 1e-9, 0xC0FFEEL)).isTrue();
	}

	@Test
	public void build_333_passes_random_spot_check() {
		PanPairProduct.PairScheme s = PanPairProduct.build(3, 3, 3);
		assertThat(s.rank()).isEqualTo(54);
		assertThat(PanPairProduct.spotCheck(s, 5, 1e-9, 0xC0FFEEL)).isTrue();
	}

	@Test
	public void build_777_passes_random_spot_check() {
		// The use case driving task #42: pair two ⟨7,7,7⟩ sub-products
		// in Strassen[7,7]³ → ⟨14,14,14⟩ for a rank-1720 materialisation.
		PanPairProduct.PairScheme s = PanPairProduct.build(7, 7, 7);
		assertThat(s.rank()).isEqualTo(490);
		assertThat(PanPairProduct.spotCheck(s, 3, 1e-9, 0xC0FFEEL)).isTrue();
	}

	@Test
	public void build_234_non_cubic_passes_spot_check() {
		// Cyclic pair ⟨2,3,4⟩+⟨3,4,2⟩.
		PanPairProduct.PairScheme s = PanPairProduct.build(2, 3, 4);
		int expectedRank = 2 * 3 * 4 + 2 * 3 + 3 * 4 + 2 * 4;
		assertThat(s.rank()).isEqualTo(expectedRank);
		assertThat(PanPairProduct.spotCheck(s, 5, 1e-9, 0xC0FFEEL)).isTrue();
	}
}
