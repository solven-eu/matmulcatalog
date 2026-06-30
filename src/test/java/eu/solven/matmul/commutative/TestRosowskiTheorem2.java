package eu.solven.matmul.commutative;

import static org.assertj.core.api.Assertions.assertThat;
import static org.assertj.core.api.Assertions.assertThatThrownBy;

import org.junit.jupiter.api.Test;

import eu.solven.matmul.NonBilinearAlgorithm;
import eu.solven.matmul.verifiers.Verifier;
import eu.solven.matmul.papers.rosowski2019.RosowskiBound;
import eu.solven.matmul.papers.rosowski2019.RosowskiTheorem2;

/**
 * Validates Rosowski 2019 Theorem 2 (divisions-free non-bilinear commutative
 * {@code l×n · n×m}, even contraction {@code n}). The headline target is the
 * {@code ⟨2,2,p⟩ = 3p+1} family, but Theorem 2 is general so we exercise a
 * spread of {@code (l, evenContraction, m)} shapes.
 */
public class TestRosowskiTheorem2 {

	@Test
	public void two_two_p_family_verifies_and_has_rank_3p_plus_1() {
		for (int p = 1; p <= 16; p++) {
			NonBilinearAlgorithm alg = RosowskiTheorem2.build22p(p);
			assertThat(alg.n).isEqualTo(2);
			assertThat(alg.m).isEqualTo(2);
			assertThat(alg.p).isEqualTo(p);
			assertThat(alg.r).as("⟨2,2," + p + "⟩ rank").isEqualTo(3 * p + 1);
			assertThat(Verifier.residualNonBilinear(alg))
					.as("Rosowski Thm2 ⟨2,2," + p + "⟩=" + (3 * p + 1) + " residual")
					.isLessThan(1e-10);
		}
	}

	@Test
	public void rank_matches_rosowski_bilinear_bound() {
		// The construction's rank must equal the catalogued closed-form bound.
		for (int p = 1; p <= 16; p++) {
			long bound = RosowskiBound.commutativeBoundBilinear(2, 2, p).orElseThrow();
			assertThat((long) RosowskiTheorem2.build22p(p).r)
					.as("Thm2 rank == commutativeBoundBilinear(2,2," + p + ")")
					.isEqualTo(bound);
		}
	}

	@Test
	public void general_even_contraction_shapes_verify() {
		// (l, evenContraction, m): exercise non-square and contraction n=4.
		int[][] shapes = { { 3, 2, 3 }, { 2, 4, 3 }, { 3, 4, 5 }, { 4, 2, 4 }, { 2, 4, 4 } };
		for (int[] s : shapes) {
			int l = s[0], n = s[1], m = s[2];
			NonBilinearAlgorithm alg = RosowskiTheorem2.build(l, n, m);
			int expected = n * (l * m + l + m - 1) / 2;
			assertThat(alg.r)
					.as("rank ⟨" + l + "," + n + "," + m + "⟩")
					.isEqualTo(expected);
			assertThat(Verifier.residualNonBilinear(alg))
					.as("residual ⟨" + l + "," + n + "," + m + "⟩")
					.isLessThan(1e-10);
		}
	}

	@Test
	public void odd_contraction_is_rejected() {
		// Theorem 2 is the EVEN-contraction case; odd n is Theorem 3 (not this).
		assertThatThrownBy(() -> RosowskiTheorem2.build(2, 3, 4))
				.isInstanceOf(IllegalArgumentException.class);
	}
}
