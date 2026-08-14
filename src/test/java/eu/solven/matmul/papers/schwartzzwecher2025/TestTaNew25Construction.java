package eu.solven.matmul.papers.schwartzzwecher2025;

import static org.assertj.core.api.Assertions.assertThat;
import static org.assertj.core.api.Assertions.assertThatThrownBy;

import org.junit.jupiter.api.Test;

import eu.solven.matmul.NonCubicBilinearAlgorithm;
import eu.solven.matmul.verifiers.Verifier;

/**
 * Regression guard for the Schwartz–Zwecher 2025 TA-New25 cubic constructor
 * ({@link TaNew25Construction}). The port was tensor-WRONG in two silent ways
 * during development — an asymmetric {@code B←L·B·Lᵀ} transform (must be the
 * uniform φ) and a transpose in family (b)'s C-form ({@code C*_{q,s}} vs
 * {@code C*_{s,q}}). Both reconstruct a valid RANK but the WRONG tensor, so a
 * rank-only check would have passed. Hence: assert the built scheme reconstructs
 * the matmul tensor EXACTLY (small even n0, arbiter of the whole construction).
 */
class TestTaNew25Construction {

	@Test
	void builds_exact_matmul_for_small_even_n0() {
		// n0=4 (γ=−2 dyadic), 6, 8 are small enough for the exact algebraic proof.
		for (int n0 : new int[] { 4, 6, 8 }) {
			NonCubicBilinearAlgorithm alg = TaNew25Construction.build(n0);
			assertThat(alg.n).isEqualTo(n0);
			assertThat(alg.m).isEqualTo(n0);
			assertThat(alg.p).isEqualTo(n0);
			assertThat((long) alg.r).as("rank ⟨%d³⟩ == tNew", n0).isEqualTo(TaNew25Construction.tNew(n0));
			assertThat(Verifier.isExactNonCubic(alg)).as("⟨%d³⟩ reconstructs matmul EXACTLY", n0).isTrue();
		}
	}

	@Test
	void rank_matches_closed_form_and_reproduces_held_imports() {
		// tNew = (4n³+45n²+122n+96)/12 — the held dense SZ imports at 28/30/32.
		assertThat(TaNew25Construction.tNew(28)).isEqualTo(10550L);
		assertThat(TaNew25Construction.tNew(30)).isEqualTo(12688L);
		assertThat(TaNew25Construction.tNew(32)).isEqualTo(15096L);
	}

	@Test
	void rejects_the_excluded_domain() {
		assertThatThrownBy(() -> TaNew25Construction.build(5)).isInstanceOf(IllegalArgumentException.class);
		assertThatThrownBy(() -> TaNew25Construction.build(16)).isInstanceOf(IllegalArgumentException.class);
	}
}
