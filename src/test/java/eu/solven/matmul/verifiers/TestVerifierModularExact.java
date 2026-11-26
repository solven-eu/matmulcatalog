package eu.solven.matmul.verifiers;

import static org.assertj.core.api.Assertions.assertThat;

import java.math.BigInteger;

import org.junit.jupiter.api.Test;

import eu.solven.matmul.NonCubicBilinearAlgorithm;

/**
 * The residue path of the exact verifier ({@link Verifier#modularDiscrepancy}) must decide
 * exactly what the BigInteger path decides. It exists because the BigInteger map does not
 * finish in useful time on a dense rational scheme of rank ~2000 (2026-11-24: the import of
 * He's ⟨14,16,16⟩:2110 was still accumulating after 30 minutes, and the 3-hourly sync would
 * have timed out on it every run) — and a faster verifier that accepted a wrong scheme
 * would be far worse than a slow one, hence the adversarial case below.
 */
public class TestVerifierModularExact {

	/** Naive ⟨3,4,5⟩ (rank 60) with a different rational scaling pushed through every product. */
	private static double[][][] scaledNaive() {
		int n = 3, m = 4, p = 5, r = n * m * p;
		double[][] u = new double[n * m][r], v = new double[m * p][r], w = new double[n * p][r];
		double[] a = { 3, -0.5, 7, 0.25, -5 }, b = { 2.0 / 3.0, -4, 1.0 / 6.0, 5 };
		int k = 0;
		for (int i = 0; i < n; i++) {
			for (int j = 0; j < m; j++) {
				for (int l = 0; l < p; l++) {
					double x = a[k % a.length], y = b[k % b.length];
					u[i * m + j][k] = x;
					v[j * p + l][k] = y;
					w[i * p + l][k] = 1.0 / (x * y);
					k++;
				}
			}
		}
		return new double[][][] { u, v, w };
	}

	private static NonCubicBilinearAlgorithm alg(double[][][] f) {
		return new NonCubicBilinearAlgorithm(3, 4, 5, f[0], f[1], f[2]);
	}

	private static String modular(NonCubicBilinearAlgorithm alg) {
		return Verifier.modularDiscrepancy(alg, Verifier.commonDenominator(alg));
	}

	@Test
	public void an_exact_rational_scheme_passes_on_both_paths() {
		NonCubicBilinearAlgorithm exact = alg(scaledNaive());
		assertThat(Verifier.commonDenominator(exact)).isGreaterThan(BigInteger.ONE);
		assertThat(Verifier.exactDiscrepancyNonCubic(exact, 0)).isNull();
		assertThat(modular(exact)).isNull();
	}

	@Test
	public void a_changed_coefficient_and_a_dropped_product_fail_on_both_paths() {
		double[][][] f = scaledNaive();
		f[1][7][7] += 1.0 / 3.0; // one V coefficient of product 7 off by a third
		assertThat(Verifier.exactDiscrepancyNonCubic(alg(f), 0)).isNotNull();
		assertThat(modular(alg(f))).isNotNull();

		double[][][] g = scaledNaive();
		for (double[] row : g[2]) row[11] = 0; // product 11 no longer reaches the output
		assertThat(Verifier.exactDiscrepancyNonCubic(alg(g), 0)).isNotNull();
		assertThat(modular(alg(g))).contains("matmul term");

		double[][][] h = scaledNaive();
		h[0][0][5] = 1; // a spurious monomial: product 5 also reads A[0][0]
		assertThat(Verifier.exactDiscrepancyNonCubic(alg(h), 0)).isNotNull();
		assertThat(modular(alg(h))).isNotNull();
	}

	/**
	 * A perturbation by a multiple of the first modulus (2³¹ − 1) is invisible modulo that
	 * modulus. The scheme must still be refused: the coefficient bound grows with the
	 * perturbation, so more moduli are taken until their lcm exceeds it.
	 */
	@Test
	public void a_perturbation_by_the_first_modulus_is_still_caught() {
		NonCubicBilinearAlgorithm exact = alg(scaledNaive());
		BigInteger d = Verifier.commonDenominator(exact);
		double[][][] f = scaledNaive();
		// scaled by D the coefficient moves by (2^31 - 1) * D: zero modulo the first modulus
		f[0][0][0] += (double) Integer.MAX_VALUE;
		NonCubicBilinearAlgorithm wrong = alg(f);
		assertThat(Verifier.commonDenominator(wrong)).isEqualTo(d);
		assertThat(Verifier.exactDiscrepancyNonCubic(wrong, 0)).isNotNull();
		assertThat(modular(wrong)).isNotNull().doesNotContain("(mod " + Integer.MAX_VALUE + ")");
	}
}
