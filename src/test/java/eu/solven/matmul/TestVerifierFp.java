package eu.solven.matmul;

import static org.assertj.core.api.Assertions.assertThat;

import java.io.File;

import org.junit.jupiter.api.Test;

import eu.solven.matmul.catalog.SchemeIO;

/**
 * Rational reduction mod a prime: a Q-exact scheme is F_p-valid exactly when
 * every coefficient's denominator is coprime to p ({@code num·den⁻¹ mod p}).
 * The two canonical directions:
 * <ul>
 *   <li>{@code 1/2} is representable in F₃ ({@code 2⁻¹ ≡ 2 mod 3}) but NOT F₂;</li>
 *   <li>{@code 1/3} is representable in F₂ ({@code 3 ≡ 1 mod 2}) but NOT F₃.</li>
 * </ul>
 */
public class TestVerifierFp {

	private static final String KNOWN = "src/main/resources/schemes/known/";

	@Test
	public void integer_scheme_reduces_mod_both_primes() throws Exception {
		NonCubicBilinearAlgorithm strassen = SchemeIO.readBilinear(
				new File(KNOWN, "section2/2x2x2-r7-strassen-db11bcc.json"));
		assertThat(Verifier.isExactNonCubicFp(strassen, 2)).as("Z ⇒ F₂").isTrue();
		assertThat(Verifier.isExactNonCubicFp(strassen, 3)).as("Z ⇒ F₃").isTrue();
	}

	@Test
	public void dyadic_scheme_is_F3_valid_but_not_F2_representable() throws Exception {
		// DPS 2025 ⟨4,4,4⟩=48: coefficients in {±1, ±1/2, ±1/4, ±1/8} — every
		// denominator a power of 2 → reduces mod 3, NOT mod 2.
		NonCubicBilinearAlgorithm dps = SchemeIO.readBilinear(
				new File(KNOWN, "section4/4x4x4-r48-dumas_pernet_sedoglavic_2025-929db4e.json"));
		assertThat(Verifier.isExactNonCubic(dps)).as("Q-exact (sanity)").isTrue();

		Verifier.FpReduction mod2 = Verifier.residualNonCubicFp(dps, 2);
		assertThat(mod2.representable()).as("1/2-style coefficients do not reduce mod 2").isFalse();

		Verifier.FpReduction mod3 = Verifier.residualNonCubicFp(dps, 3);
		assertThat(mod3.representable()).as("dyadic denominators are coprime to 3").isTrue();
		assertThat(mod3.exact()).as("Q⟨4,4,4⟩=48 verifies in F₃ via 1/2 ≡ 2").isTrue();
	}

	@Test
	public void third_scheme_is_F2_valid_but_not_F3_representable() throws Exception {
		// Synthetic 1/3 witness: per-term rescale of Strassen (U_k×3, V_k×1/3)
		// preserves exactness; coefficients gain denominator 3 → reduces mod 2
		// (1/3 ≡ 1), NOT mod 3.
		NonCubicBilinearAlgorithm strassen = SchemeIO.readBilinear(
				new File(KNOWN, "section2/2x2x2-r7-strassen-db11bcc.json"));
		double[][] u = strassen.denseU(), v = strassen.denseV(), w = strassen.denseW();
		for (int i = 0; i < u.length; i++) u[i][0] *= 3.0;
		for (int i = 0; i < v.length; i++) v[i][0] /= 3.0;
		NonCubicBilinearAlgorithm scaled = new NonCubicBilinearAlgorithm(2, 2, 2, u, v, w);
		assertThat(Verifier.isExactNonCubic(scaled)).as("rescale preserves exactness").isTrue();

		Verifier.FpReduction mod3 = Verifier.residualNonCubicFp(scaled, 3);
		assertThat(mod3.representable()).as("1/3 does not reduce mod 3").isFalse();

		Verifier.FpReduction mod2 = Verifier.residualNonCubicFp(scaled, 2);
		assertThat(mod2.representable()).as("denominator 3 is coprime to 2").isTrue();
		assertThat(mod2.exact()).as("verifies in F₂ via 1/3 ≡ 1").isTrue();
	}

	@Test
	public void fp_spotcheck_accepts_f2_only_scheme_that_char0_spotcheck_rejects() throws Exception {
		// AlphaTensor ⟨4,4,4⟩=47 computes matmul ONLY mod 2 (the char-0 optimum is 49).
		// The materialiser's verification must use the F_p spot-check for an F₂ sweep —
		// the char-0 (random-real) spot-check WRONGLY rejects this valid scheme.
		NonCubicBilinearAlgorithm at47 = SchemeIO.readBilinear(
				new File(KNOWN, "section4/4x4x4-r47-alphatensor_F2-258e5b7.json"));
		assertThat(Verifier.passesRandomMatmulSpotCheckFp(at47, 2))
				.as("AT ⟨4,4,4⟩=47 verifies over F₂ (fast spot-check)").isTrue();
		assertThat(Verifier.passesRandomMatmulSpotCheck(at47))
				.as("but the char-0 spot-check correctly rejects it (47 < 49 over char-0)").isFalse();
		// It is also not even representable / exact mod 3.
		assertThat(Verifier.passesRandomMatmulSpotCheckFp(at47, 3)).as("not an F₃ scheme").isFalse();
	}

	@Test
	public void fp_spotcheck_agrees_with_char0_on_an_integer_scheme() throws Exception {
		NonCubicBilinearAlgorithm strassen = SchemeIO.readBilinear(
				new File(KNOWN, "section2/2x2x2-r7-strassen-db11bcc.json"));
		assertThat(Verifier.passesRandomMatmulSpotCheck(strassen)).isTrue();
		assertThat(Verifier.passesRandomMatmulSpotCheckFp(strassen, 2)).isTrue();
		assertThat(Verifier.passesRandomMatmulSpotCheckFp(strassen, 3)).isTrue();
	}
}
