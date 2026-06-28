package eu.solven.matmul;

import static org.assertj.core.api.Assertions.assertThat;

import org.junit.jupiter.api.Test;

/**
 * Guards the size/density-aware verification dispatcher ({@link Verifier#verifyAuto}):
 * it must route a small/sparse scheme to the EXACT symbolic proof, and a scheme whose
 * exact term-budget is exceeded to the memory-light random spot-check — while both still
 * return a correct verdict.
 */
public class TestVerifierAuto {

	@Test
	public void small_sparse_scheme_gets_exact_proof() {
		NonCubicBilinearAlgorithm naive = NonCubicBilinearAlgorithm.naive(2, 2, 2);
		Verifier.Verdict v = Verifier.verifyAuto(naive);
		assertThat(v.ok()).isTrue();
		assertThat(v.strategy()).isEqualTo(Verifier.VerifyStrategy.EXACT_SYMBOLIC);
		assertThat(v.isProof()).isTrue();
		// naive ⟨2,2,2⟩ has 8 rank-1 products, each a single (a,b,c) term ⇒ 8 generated terms.
		assertThat(v.estimatedTerms()).isEqualTo(8);
	}

	@Test
	public void over_budget_falls_back_to_spot_check() {
		// Force the fallback with a zero budget: even a tiny scheme then routes to the
		// memory-light spot-check, which must still confirm it.
		NonCubicBilinearAlgorithm naive = NonCubicBilinearAlgorithm.naive(2, 3, 4);
		Verifier.Verdict v = Verifier.verifyAuto(naive, 0L);
		assertThat(v.ok()).isTrue();
		assertThat(v.strategy()).isEqualTo(Verifier.VerifyStrategy.RANDOM_SPOT_CHECK);
		assertThat(v.isProof()).isFalse();
	}

	@Test
	public void estimate_short_circuits_at_cap() {
		NonCubicBilinearAlgorithm naive = NonCubicBilinearAlgorithm.naive(3, 3, 3);
		// 27 single-term products ⇒ estimate 27 when the cap allows it...
		assertThat(Verifier.estimateExactTerms(naive, 1_000)).isEqualTo(27);
		// ...and Long.MAX_VALUE once the running sum exceeds a tiny cap.
		assertThat(Verifier.estimateExactTerms(naive, 5)).isEqualTo(Long.MAX_VALUE);
	}
}
