package eu.solven.matmul;

import static org.assertj.core.api.Assertions.assertThat;

import java.io.File;

import org.junit.jupiter.api.Test;

import eu.solven.matmul.NonCubicBilinearAlgorithm;
import eu.solven.matmul.verifiers.Verifier;
import eu.solven.matmul.catalog.SchemeIO;

/**
 * Cyclic-shift correctness: applying {@link NonCubicBilinearAlgorithm#cyclicShift}
 * to any algorithm yields a verifiable algorithm for the cyclically-shifted
 * format, and three applications return to the original.
 */
public class TestCyclicShift {

	@Test
	public void cyclic_shift_on_strassen_222_round_trips_and_verifies() throws Exception {
		NonCubicBilinearAlgorithm orig = SchemeIO.read(
				eu.solven.matmul.catalog.SchemeResolver.byHint("src/main/resources/schemes/known/section2/strassen-2x2x2_m7_a18.json"));
		assertThat(orig.n).isEqualTo(2);
		assertThat(Verifier.isExactNonCubic(orig)).isTrue();

		NonCubicBilinearAlgorithm shifted1 = orig.cyclicShift();
		assertThat(shifted1.n).isEqualTo(2);
		assertThat(shifted1.m).isEqualTo(2);
		assertThat(shifted1.p).isEqualTo(2);
		assertThat(Verifier.isExactNonCubic(shifted1)).isTrue();

		NonCubicBilinearAlgorithm shifted3 = shifted1.cyclicShift().cyclicShift();
		assertThat(Verifier.isExactNonCubic(shifted3)).isTrue();
	}

	@Test
	public void cyclic_shift_on_non_cubic_334_yields_343_and_433() throws Exception {
		// Find any ⟨3,3,4⟩ scheme.
		// Use the integer (Z) variant; skip F2 since it doesn't verify under real arithmetic.
		File f = eu.solven.matmul.catalog.SchemeResolver.byHint("src/main/resources/schemes/known/section4/alphatensor_Z-3x3x4_m29_a148.json");
		NonCubicBilinearAlgorithm orig = SchemeIO.read(f);
		assertThat(orig.n).isEqualTo(3);
		assertThat(orig.m).isEqualTo(3);
		assertThat(orig.p).isEqualTo(4);
		assertThat(Verifier.isExactNonCubic(orig)).isTrue();

		NonCubicBilinearAlgorithm shift1 = orig.cyclicShift();
		assertThat(shift1.n).isEqualTo(3);
		assertThat(shift1.m).isEqualTo(4);
		assertThat(shift1.p).isEqualTo(3);
		assertThat(shift1.r).isEqualTo(orig.r);
		assertThat(Verifier.isExactNonCubic(shift1)).isTrue();

		NonCubicBilinearAlgorithm shift2 = shift1.cyclicShift();
		assertThat(shift2.n).isEqualTo(4);
		assertThat(shift2.m).isEqualTo(3);
		assertThat(shift2.p).isEqualTo(3);
		assertThat(shift2.r).isEqualTo(orig.r);
		assertThat(Verifier.isExactNonCubic(shift2)).isTrue();

		// Third shift returns to original format ⟨3,3,4⟩.
		NonCubicBilinearAlgorithm shift3 = shift2.cyclicShift();
		assertThat(shift3.n).isEqualTo(3);
		assertThat(shift3.m).isEqualTo(3);
		assertThat(shift3.p).isEqualTo(4);
		assertThat(Verifier.isExactNonCubic(shift3)).isTrue();
	}

	@Test
	public void orient_as_reaches_all_cyclic_orientations() throws Exception {
		File f = eu.solven.matmul.catalog.SchemeResolver.byHint("src/main/resources/schemes/known/section4/alphatensor_Z-3x3x4_m29_a148.json");
		NonCubicBilinearAlgorithm orig = SchemeIO.read(f);

		assertThat(orig.orientAs(3, 3, 4)).isPresent();
		assertThat(orig.orientAs(3, 4, 3)).isPresent();
		assertThat(orig.orientAs(4, 3, 3)).isPresent();

		// Non-cyclic orientations (would need transpose): not implemented yet.
		// ⟨3,3,4⟩ doesn't have a non-cyclic distinct orientation since two dims match,
		// so all 3 orderings ARE cyclic. Use ⟨2,3,4⟩ for an actually-non-cyclic test.
	}
}
