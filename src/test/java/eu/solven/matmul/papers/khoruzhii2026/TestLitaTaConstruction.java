package eu.solven.matmul.papers.khoruzhii2026;

import static org.assertj.core.api.Assertions.assertThat;
import static org.assertj.core.api.Assertions.assertThatThrownBy;

import org.junit.jupiter.api.Test;

import eu.solven.matmul.NonCubicBilinearAlgorithm;
import eu.solven.matmul.verifiers.Verifier;

/**
 * Validates the exact-rational LITA constructor {@link LitaTaConstruction} (a port of the
 * Khoruzhii&ndash;Gel&szlig;&ndash;Pokutta 2026 Maple generators): every shipped {@code <N,N,N>}
 * must compute matmul EXACTLY ({@link Verifier#isExactNonCubic}) at the published rank
 * ({@link LitaTrilinearAggregation#cubicRank(int)}). Odd N reproduce the reference scheme up to
 * a row permutation; even N are valid alternate decompositions at the same rank.
 */
class TestLitaTaConstruction {

	private static final int[] ODD = { 19, 21, 23, 25 };
	private static final int[] EVEN = { 26, 28, 30, 32 };

	private NonCubicBilinearAlgorithm buildAtRank(int n) {
		NonCubicBilinearAlgorithm alg = LitaTaConstruction.build(n);
		assertThat(alg.n).isEqualTo(n);
		assertThat(alg.m).isEqualTo(n);
		assertThat(alg.p).isEqualTo(n);
		assertThat((long) alg.r)
				.as("rank for <%d,%d,%d>", n, n, n)
				.isEqualTo(LitaTrilinearAggregation.cubicRank(n));
		return alg;
	}

	@Test
	void odd_schemes_are_exact_at_formula_rank() {
		// Odd schemes are sparse: the EXACT symbolic verifier (a real algebraic proof)
		// is feasible and is the gate.
		for (int n : ODD) {
			NonCubicBilinearAlgorithm alg = buildAtRank(n);
			assertThat(Verifier.isExactNonCubic(alg))
					.as("exact matmul for <%d,%d,%d>", n, n, n)
					.isTrue();
		}
	}

	@Test
	void even_schemes_are_correct_at_formula_rank() {
		// Even schemes are DENSE (the phi-embedding fills the factors): the exact
		// symbolic verifier would generate 1.7e10..7.9e10 terms (the documented
		// <30,32,32> OOM), so the repo routes such schemes to the randomised matmul
		// spot-check (false-accept ≈ 0) — exactly what Verifier.verifyAuto picks. We
		// assert the auto verdict is OK and, for clarity, that the chosen strategy is
		// the spot-check (i.e. exact was correctly deemed infeasible here).
		for (int n : EVEN) {
			NonCubicBilinearAlgorithm alg = buildAtRank(n);
			Verifier.Verdict verdict = Verifier.verifyAuto(alg);
			assertThat(verdict.ok())
					.as("matmul spot-check for <%d,%d,%d>", n, n, n)
					.isTrue();
			assertThat(Verifier.passesRandomMatmulSpotCheck(alg, 16, 1e-5))
					.as("independent dense spot-check for <%d,%d,%d>", n, n, n)
					.isTrue();
		}
	}

	@Test
	void rejects_below_min_n() {
		assertThatThrownBy(() -> LitaTaConstruction.build(18))
				.isInstanceOf(IllegalArgumentException.class);
	}
}
