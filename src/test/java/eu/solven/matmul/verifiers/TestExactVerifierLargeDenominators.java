package eu.solven.matmul.verifiers;

import static org.assertj.core.api.Assertions.assertThat;

import org.junit.jupiter.api.Test;

import eu.solven.matmul.NonCubicBilinearAlgorithm;

/**
 * The exact (BigInteger) verifier must accept exact rational schemes whose
 * denominators are large, and must still reject inexact ones.
 *
 * <p>The silent failure (2026-09-30): denominators were recovered by scanning
 * {@code d ≤ 1024} only, so a coefficient like {@code 3866/3705} read as irrational
 * and the whole scheme as "not exact". Two of Perminov's published rational
 * serendipitous bases ({@code 2x4x6_m39}, {@code 2x5x7_m56}) were rejected that way on
 * every sync run — a {@code [FAIL]} line in a green job — which kept ⟨4,16,30⟩ and
 * ⟨8,14,20⟩ above his ranks. The second base also needed the scaled numerators to be
 * computed exactly: its common denominator (~2·10¹⁷) is beyond a double's 53 bits.</p>
 */
public class TestExactVerifierLargeDenominators {

	/** ⟨1,1,1⟩ as {@code r} products: {@code Σ u_k·v_k·w_k = 1}. */
	private static NonCubicBilinearAlgorithm scalar(double[] u, double[] v, double[] w) {
		return new NonCubicBilinearAlgorithm(1, 1, 1, new double[][] { u }, new double[][] { v }, new double[][] { w });
	}

	@Test
	public void a_denominator_above_1024_is_a_rational_not_an_irrational() {
		NonCubicBilinearAlgorithm exact = scalar(new double[] { 3866.0 / 3705 }, new double[] { 3705.0 / 3866 },
				new double[] { 1 });
		assertThat(Verifier.isExactNonCubic(exact)).isTrue();

		NonCubicBilinearAlgorithm off = scalar(new double[] { 3866.0 / 3705 }, new double[] { 3705.0 / 3866 + 1e-7 },
				new double[] { 1 });
		assertThat(Verifier.isExactNonCubic(off)).as("a perturbed coefficient is still rejected").isFalse();
	}

	@Test
	public void a_common_denominator_beyond_53_bits_is_scaled_exactly() {
		// Six terms (1/q)·q·(1/6) with six distinct primes q: the lcm of the denominators is
		// 6·∏q ≈ 7·10¹⁸ (63 bits). v·D in double arithmetic is off in its low digits.
		int[] primes = { 1009, 1013, 1019, 1021, 1031, 1033 };
		double[] u = new double[6], v = new double[6], w = new double[6];
		for (int k = 0; k < 6; k++) {
			u[k] = 1.0 / primes[k];
			v[k] = primes[k];
			w[k] = 1.0 / 6;
		}
		assertThat(Verifier.isExactNonCubic(scalar(u, v, w))).isTrue();

		w[5] = 1.0 / 7; // the six terms no longer sum to 1
		assertThat(Verifier.isExactNonCubic(scalar(u, v, w))).isFalse();
	}

	@Test
	public void an_irrational_coefficient_never_passes() {
		double s = Math.sqrt(2) / 2;
		// (√2/2)·(√2/2)·2 = 1 only in floating point. The wider recovery does find a
		// fraction for it (the convergent 470832/665857 is within 1e-12), and that is
		// harmless: 2·470832² = 665857² − 1, so the integer identity fails. A recovered
		// fraction can make an exact scheme fail; it cannot make an inexact one pass.
		assertThat(Verifier.isExactNonCubic(scalar(new double[] { s }, new double[] { s }, new double[] { 2 })))
				.isFalse();
	}
}
