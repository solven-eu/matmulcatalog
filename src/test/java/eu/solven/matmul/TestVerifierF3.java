package eu.solven.matmul;

import static org.assertj.core.api.Assertions.assertThat;

import java.nio.file.Files;
import java.nio.file.Path;

import org.junit.jupiter.api.Test;

import eu.solven.matmul.catalog.SchemeIO;

/**
 * F₃ residual sanity. Any integer-coefficient bilinear scheme that
 * verifies over ℤ also verifies over F₃ (the trilinear identity holds
 * over ℤ, so reducing mod 3 preserves it). Strassen's ⟨2,2,2⟩=7 with
 * ±1 coefficients is the canonical witness.
 */
public class TestVerifierF3 {

	@Test
	public void integerSchemeAlsoVerifiesInF3() throws Exception {
		Path strassen = Path.of(
				"src/main/resources/schemes/known/section2/2x2x2-r7-strassen-db11bcc.json");
		NonCubicBilinearAlgorithm alg = SchemeIO.read(Files.readString(strassen));

		assertThat(Verifier.isExactNonCubic(alg))
				.as("ℤ verification (sanity)")
				.isTrue();
		assertThat(Verifier.isExactNonCubicF2(alg))
				.as("F₂ verification (sanity)")
				.isTrue();
		assertThat(Verifier.isExactNonCubicF3(alg))
				.as("F₃ verification — integer scheme should pass mod 3")
				.isTrue();
	}
}
