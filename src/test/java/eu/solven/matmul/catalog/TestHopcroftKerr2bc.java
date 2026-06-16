package eu.solven.matmul.catalog;

import eu.solven.matmul.papers.hopcroftkerr1971.HopcroftKerr2bc;

import static org.assertj.core.api.Assertions.assertThat;

import org.junit.jupiter.api.Test;

import eu.solven.matmul.NonCubicBilinearAlgorithm;
import eu.solven.matmul.Verifier;

public class TestHopcroftKerr2bc {

	@Test
	public void rank_formula_matches_dis09_table2() {
		// HK formula (3pn+max(p,n))/2 is for shape ⟨p, 2, n⟩ — outer dims p, n with middle 2.
		// In our shape convention this is ⟨p, 2, n⟩, so our gaps ⟨2,9,9⟩ correspond to HK(9, 9) = 126
		// via tensor symmetry ⟨2,9,9⟩ ≅ ⟨9,2,9⟩.
		assertThat(HopcroftKerr2bc.rank(2, 2)).isEqualTo(7);
		assertThat(HopcroftKerr2bc.rank(3, 3)).isEqualTo(15);
		assertThat(HopcroftKerr2bc.rank(9, 9)).isEqualTo(126);
		assertThat(HopcroftKerr2bc.rank(10, 15)).isEqualTo(233);
		assertThat(HopcroftKerr2bc.rank(10, 16)).isEqualTo(248);
		assertThat(HopcroftKerr2bc.rank(12, 16)).isEqualTo(296);
	}

	@Test
	public void squareEven_n_even_in_range_2_to_16() {
		for (int n = 2; n <= 16; n += 2) {
			NonCubicBilinearAlgorithm alg = HopcroftKerr2bc.buildSquareEven(n);
			assertThat(alg.n).isEqualTo(n);
			assertThat(alg.m).isEqualTo(2);
			assertThat(alg.p).isEqualTo(n);
			assertThat(alg.r).as("n=" + n).isEqualTo(HopcroftKerr2bc.rank(n, n));
			assertThat(Verifier.passesRandomMatmulSpotCheck(alg))
					.as("n=" + n + " spot-check").isTrue();
		}
	}

	@Test
	public void buildSquare_dispatches_correctly() {
		// Sanity: dispatcher returns the same result as the explicit odd/even calls.
		NonCubicBilinearAlgorithm odd = HopcroftKerr2bc.buildSquare(5);
		assertThat(odd.r).isEqualTo(40);
		NonCubicBilinearAlgorithm even = HopcroftKerr2bc.buildSquare(8);
		assertThat(even.r).isEqualTo(100);
	}

	@Test
	public void squareOdd_n_odd_in_range_3_to_15() {
		// The case-1 alternating coloring (1, 2, 1, 2, …, 1, 2, 3) places every
		// same-method pair as either (1,1) bridged by a method-2 neighbour, or
		// (2,2) bridged by a method-1 neighbour. Both cases are now implemented;
		// expect a valid scheme for every odd n in the range.
		for (int n = 3; n <= 15; n += 2) {
			NonCubicBilinearAlgorithm alg = HopcroftKerr2bc.buildSquareOdd(n);
			assertThat(alg.n).isEqualTo(n);
			assertThat(alg.m).isEqualTo(2);
			assertThat(alg.p).isEqualTo(n);
			assertThat(alg.r)
					.as("n=" + n)
					.isEqualTo(HopcroftKerr2bc.rank(n, n));
			assertThat(Verifier.passesRandomMatmulSpotCheck(alg))
					.as("n=" + n + " spot-check")
					.isTrue();
		}
	}

	@Test
	public void squareOdd_n5_with_paired_same_method_cases() {
		// n=5 with methods (1, 2, 1, 2, 3). Same-method pairs:
		//   (1,3) both method 1, bridge row 2 method 2 → (1,1,bridge2) ✓
		//   (2,4) both method 2, bridge row 3 method 1 → (2,2,bridge1) ✓
		// Both implemented; expect a valid 40-product scheme.
		NonCubicBilinearAlgorithm alg = HopcroftKerr2bc.buildSquareOdd(5);
		assertThat(alg.n).isEqualTo(5);
		assertThat(alg.m).isEqualTo(2);
		assertThat(alg.p).isEqualTo(5);
		assertThat(alg.r).isEqualTo(HopcroftKerr2bc.rank(5, 5));  // (3·5·5+5)/2 = 40
		assertThat(alg.r).isEqualTo(40);
		assertThat(Verifier.passesRandomMatmulSpotCheck(alg)).isTrue();
	}

	@Test
	public void squareOdd_n3_reproduces_paper_example() {
		// Page-4 example: ⟨3, 2, 3⟩ in 15 multiplications.
		NonCubicBilinearAlgorithm a = HopcroftKerr2bc.buildSquareOdd(3);
		assertThat(a.n).isEqualTo(3);
		assertThat(a.m).isEqualTo(2);
		assertThat(a.p).isEqualTo(3);
		assertThat(a.r).isEqualTo(15);
		assertThat(Verifier.passesRandomMatmulSpotCheck(a)).isTrue();
	}
}
