package eu.solven.matmul;

import static org.assertj.core.api.Assertions.assertThat;

import java.io.File;
import java.util.Optional;

import org.junit.jupiter.api.Test;

import eu.solven.matmul.NonCubicBilinearAlgorithm;
import eu.solven.matmul.verifiers.Verifier;
import eu.solven.matmul.catalog.SchemeIO;

/**
 * Validates that {@link NonCubicBilinearAlgorithm#orientAs} produces a
 * correct algorithm for transpose-reachable axis orderings (not just
 * cyclic).
 */
public class TestOrientAsTranspose {

	@Test
	public void perminov_9_11_12_orient_as_9_12_11_passes_spot_check() throws Exception {
		NonCubicBilinearAlgorithm orig = SchemeIO.read(
				eu.solven.matmul.catalog.SchemeResolver.byHint("src/main/resources/schemes/known/section12/perminov_Q-9x11x12_m738_a39061.json"));
		System.out.printf("Original ⟨%d,%d,%d⟩=%d%n", orig.n, orig.m, orig.p, orig.r);
		assertThat(Verifier.passesRandomMatmulSpotCheck(orig))
				.as("⟨9,11,12⟩ original passes spot-check").isTrue();

		// Now orient as ⟨9,12,11⟩ — that's transpose∘cyclic² in S₃.
		Optional<NonCubicBilinearAlgorithm> oriented = orig.orientAs(9, 12, 11);
		assertThat(oriented).as("orientAs(9,12,11) returns a scheme").isPresent();
		NonCubicBilinearAlgorithm o = oriented.get();
		System.out.printf("Oriented ⟨%d,%d,%d⟩=%d%n", o.n, o.m, o.p, o.r);
		assertThat(o.n).isEqualTo(9);
		assertThat(o.m).isEqualTo(12);
		assertThat(o.p).isEqualTo(11);
		assertThat(o.r).isEqualTo(738);
		assertThat(Verifier.passesRandomMatmulSpotCheck(o))
				.as("⟨9,11,12⟩ oriented as ⟨9,12,11⟩ should also verify as matmul")
				.isTrue();
	}

	@Test
	public void strassen_222_transpose_passes_spot_check() throws Exception {
		// Cubic case: ⟨2,2,2⟩ transposed is still ⟨2,2,2⟩, sanity check.
		NonCubicBilinearAlgorithm s = SchemeIO.read(
				eu.solven.matmul.catalog.SchemeResolver.byHint("src/main/resources/schemes/known/section2/strassen-2x2x2_m7_a18.json"));
		assertThat(Verifier.passesRandomMatmulSpotCheck(s.transpose())).isTrue();
	}
}
