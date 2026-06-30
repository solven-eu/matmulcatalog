package eu.solven.matmul.commutative;

import static org.assertj.core.api.Assertions.assertThat;

import org.junit.jupiter.api.Test;

import eu.solven.matmul.NonBilinearAlgorithm;
import eu.solven.matmul.papers.rosowski2019.Rosowski21;
import eu.solven.matmul.papers.rosowski2019.RosowskiAlgorithm1;
import eu.solven.matmul.verifiers.Verifier;

/**
 * Validates Rosowski 2019 Algorithm 1 ({@code ⟨n,3,3⟩ = 6n+3} non-bilinear
 * commutative) builds and verifies for {@code n = 1..32} — the full range
 * this project tracks.
 */
public class TestRosowskiAlgorithm1 {

	@Test
	public void rosowski_algorithm1_verifies_for_n_1_to_32() {
		for (int n = 1; n <= 32; n++) {
			NonBilinearAlgorithm alg = RosowskiAlgorithm1.build(n);
			assertThat(alg.n).isEqualTo(n);
			assertThat(alg.m).isEqualTo(3);
			assertThat(alg.p).isEqualTo(3);
			assertThat(alg.r).isEqualTo(6 * n + 3);
			double residual = Verifier.residualNonBilinear(alg);
			assertThat(residual)
					.as("Rosowski Algorithm 1 ⟨" + n + ",3,3⟩=" + (6 * n + 3) + " residual")
					.isLessThan(1e-10);
		}
	}

	@Test
	public void rosowski_algorithm1_n3_matches_corollary1_rank() {
		// At n=3 the rank is 21 — same as Rosowski Corollary 1 (Rosowski21).
		// The schemes aren't bit-identical (different index orderings of the
		// 21 products) but both compute ⟨3,3,3⟩ in 21 commutative mults.
		NonBilinearAlgorithm alg = RosowskiAlgorithm1.build(3);
		assertThat(alg.r).isEqualTo(21);
		assertThat(Verifier.isExactNonBilinear(alg)).isTrue();
		NonBilinearAlgorithm corollary = Rosowski21.build();
		assertThat(corollary.r).isEqualTo(alg.r);
	}
}
