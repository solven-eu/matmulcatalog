package eu.solven.matmul.catalog;

import static org.assertj.core.api.Assertions.assertThat;

import java.io.File;

import org.junit.jupiter.api.Test;

import eu.solven.matmul.NonCubicBilinearAlgorithm;
import eu.solven.matmul.verifiers.Verifier;
import tools.jackson.databind.JsonNode;

/**
 * Guards the special-case Strassen–Winograd ⟨2,2,2⟩=7 addition counts (#185):
 * the flat bilinear count from {@code (U,V,W)} is 24, but the scheduled count
 * (when the CSE intermediates s1..s4 / t1..t4 are shared) is 15. The scheme
 * carries the scheduled figure as {@code scheduled_additions}; the flat figure
 * must agree with {@link Verifier#additionCount}. This distinction is the reason
 * raw addition counts are an unreliable cross-scheme comparison (see the
 * "Optimisation target" note in {@code paper/sections/wording.tex}).
 */
public class TestWinogradScheduledAdditions {

	private static final File WINOGRAD =
			eu.solven.matmul.catalog.SchemeResolver.byHint("src/main/resources/schemes/known/section2/winograd_1971-2x2x2_m7_a24.json");

	@Test
	public void winograd_222_scheduled_is_15_flat_is_24() throws Exception {
		JsonNode root = SchemeIO.parseJson(WINOGRAD);
		assertThat(root.has("scheduled_additions"))
				.as("Winograd ⟨2,2,2⟩ carries an explicit scheduled_additions field")
				.isTrue();
		assertThat(root.get("scheduled_additions").asInt())
				.as("Strassen-Winograd schedule = 15 (CSE-shared s1..s4,t1..t4)")
				.isEqualTo(15);

		NonCubicBilinearAlgorithm alg = SchemeIO.read(WINOGRAD);
		assertThat(alg.r).isEqualTo(7);
		assertThat(Verifier.additionCount(alg))
				.as("flat bilinear additionCount from U/V/W is 24 (> the 15 scheduled)")
				.isEqualTo(24);
	}
}
