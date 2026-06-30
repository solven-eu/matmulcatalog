package eu.solven.matmul.commutative;

import static org.assertj.core.api.Assertions.assertThat;

import org.junit.jupiter.api.Test;

import eu.solven.matmul.NonBilinearAlgorithm;
import eu.solven.matmul.papers.rosowski2019.Rosowski21;
import eu.solven.matmul.verifiers.Verifier;

/**
 * Verifies that the encoded Rosowski 2019 Corollary 1 (⟨3,3,3⟩=21
 * commutative non-bilinear) algorithm passes
 * {@link Verifier#isExactNonBilinear}.
 */
public class TestRosowski21 {

	@Test
	public void rosowski21_verifies_as_commutative_matmul() {
		NonBilinearAlgorithm alg = Rosowski21.build();
		assertThat(alg.n).isEqualTo(3);
		assertThat(alg.m).isEqualTo(3);
		assertThat(alg.p).isEqualTo(3);
		assertThat(alg.r).isEqualTo(21);
		assertThat(alg.isPurelyBilinear()).isFalse();
		double residual = Verifier.residualNonBilinear(alg);
		System.out.printf("Rosowski 21 residual: %.6e%n", residual);
		assertThat(Verifier.isExactNonBilinear(alg))
				.as("Rosowski 2019 Corollary 1 (⟨3,3,3⟩=21 commutative) should verify; residual=" + residual)
				.isTrue();
	}
}
