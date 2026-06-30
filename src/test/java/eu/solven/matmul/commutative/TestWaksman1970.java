package eu.solven.matmul.commutative;

import org.junit.jupiter.api.Tag;


import eu.solven.matmul.papers.waksman1970.Waksman1970;

import static org.assertj.core.api.Assertions.assertThat;

import org.junit.jupiter.api.Test;

import eu.solven.matmul.NonBilinearAlgorithm;
import eu.solven.matmul.verifiers.Verifier;

/**
 * Validates {@link Waksman1970} — the generic ⟨n,n,n⟩ commutative
 * matmul from Mezzarobba 2007 Fig.3 + closed-form — builds, verifies,
 * and matches the closed-form rank for every {@code n = 2..32}.
 *
 * <p>This is the canonical commutative baseline cited by DIS09 Table 4
 * via Mezzarobba 2007.</p>
 */
@Tag("slow")
public class TestWaksman1970 {

	@Test
	public void waksman_verifies_for_n_2_to_32() {
		// For n=32, algebraic residualNonBilinear is O(n^6·r) ≈ 10^13 ops
		// (impractical). Use the fast random-spot-check instead: O(samples·
		// r·n²). It is a probabilistic witness (Schwartz-Zippel: missing a
		// real implementation defect has measure-zero probability per sample
		// for fixed-coefficient schemes), NOT a proof. An exact-rational
		// cross-check (BigFraction over Q, independent re-derivation, n=2..32,
		// matches the closed-form rank) lives in
		// eu.solven.matmul.papers.waksman1970.VerifyWaksman1970.
		for (int n = 2; n <= 32; n++) {
			NonBilinearAlgorithm alg = Waksman1970.build(n);
			int expected = (n * n + 2 * n - 1) * (n / 2) + (n % 2) * n * n;

			assertThat(alg.n).isEqualTo(n);
			assertThat(alg.m).isEqualTo(n);
			assertThat(alg.p).isEqualTo(n);
			assertThat(alg.r).isEqualTo(expected);
			assertThat(Waksman1970.rank(n)).isEqualTo(expected);

			assertThat(Verifier.passesRandomMatmulSpotCheckNB(alg))
					.as("Waksman ⟨" + n + "," + n + "," + n + "⟩=" + expected + " spot-check")
					.isTrue();
		}
	}

	@Test
	public void waksman_n3_matches_fig3_exact_residual() {
		// At small n, also run the full algebraic O(n^6·r) check.
		assertThat(Waksman1970.rank(3)).isEqualTo(23);
		assertThat(Verifier.isExactNonBilinear(Waksman1970.build(3))).isTrue();
	}

	@Test
	public void waksman_n2_equals_strassen_count_over_commutative_ring() {
		assertThat(Waksman1970.rank(2)).isEqualTo(7);
		assertThat(Verifier.isExactNonBilinear(Waksman1970.build(2))).isTrue();
	}
}
