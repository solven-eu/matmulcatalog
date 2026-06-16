package eu.solven.matmul;

import static org.assertj.core.api.Assertions.assertThat;

import java.io.File;

import org.junit.jupiter.api.Test;

import eu.solven.matmul.papers.laderman1976.Laderman23;
import eu.solven.matmul.NonCubicBilinearAlgorithm;
import eu.solven.matmul.Verifier;
import eu.solven.matmul.catalog.SchemeIO;

/**
 * Sanity test for the fast randomised verifier
 * {@link Verifier#passesRandomMatmulSpotCheck}: must accept known-good
 * schemes (Strassen, Laderman) and reject a deliberately-corrupted scheme.
 */
public class TestRandomSpotCheck {

	@Test
	public void strassen_passes_spot_check() throws Exception {
		NonCubicBilinearAlgorithm s = SchemeIO.read(
				eu.solven.matmul.catalog.SchemeResolver.byHint("src/main/resources/schemes/known/section2/strassen-2x2x2_m7_a18.json"));
		assertThat(Verifier.passesRandomMatmulSpotCheck(s)).isTrue();
	}

	@Test
	public void laderman_passes_spot_check() {
		NonCubicBilinearAlgorithm l = NonCubicBilinearAlgorithm.fromCubic(Laderman23.get());
		assertThat(Verifier.passesRandomMatmulSpotCheck(l)).isTrue();
	}

	@Test
	public void corrupted_strassen_fails_spot_check() throws Exception {
		NonCubicBilinearAlgorithm s = SchemeIO.read(
				eu.solven.matmul.catalog.SchemeResolver.byHint("src/main/resources/schemes/known/section2/strassen-2x2x2_m7_a18.json"));
		// Flip a single sign in W to break the algorithm.
		double[][] srcU = s.denseU();
		double[][] srcV = s.denseV();
		double[][] srcW = s.denseW();
		double[][] badW = new double[srcW.length][];
		for (int i = 0; i < srcW.length; i++) badW[i] = srcW[i].clone();
		badW[0][0] = -badW[0][0] + 1; // arbitrary corruption
		NonCubicBilinearAlgorithm bad = new NonCubicBilinearAlgorithm(s.n, s.m, s.p, srcU, srcV, badW);
		assertThat(Verifier.passesRandomMatmulSpotCheck(bad)).isFalse();
	}
}
