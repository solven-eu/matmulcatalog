package eu.solven.matmul.f2.sat;

import static org.assertj.core.api.Assertions.assertThat;

import java.util.Optional;

import org.junit.jupiter.api.Test;

import eu.solven.matmul.f2.sat.SatMatmulPipeline;

/**
 * Tests the non-cubic ⟨n, m, p⟩ Z/2 search pipeline. Validates the SAT path
 * at the known upper bound {@code R(⟨2,2,3⟩) ≤ 11} (Hopcroft–Kerr 1971,
 * tight over characteristic-0 fields; Z/2 may permit lower, which is one
 * reason to keep the search pipeline general).
 */
public class TestSatNonCubic {

	@Test
	public void denseZ2_223_r11_isSat() {
		int[][][] target = SatMatmulPipeline.z2DenseMatmulTensor(2, 2, 3);
		Optional<double[][][]> result =
				SatMatmulPipeline.findZ2NonCubicDecomposition(2, 2, 3, 11, target);
		assertThat(result).isPresent();
		assertThat(SatMatmulPipeline.verifyZ2NonCubic(result.get(), target)).isTrue();
	}

	/** Verifies the target tensor shape: dim(U) = n·m, dim(V) = m·p, dim(W) = n·p. */
	@Test
	public void targetTensorShape() {
		int[][][] target223 = SatMatmulPipeline.z2DenseMatmulTensor(2, 2, 3);
		assertThat(target223.length).isEqualTo(4);          // dimU = 2·2
		assertThat(target223[0].length).isEqualTo(6);       // dimV = 2·3
		assertThat(target223[0][0].length).isEqualTo(6);    // dimW = 2·3

		int[][][] target233 = SatMatmulPipeline.z2DenseMatmulTensor(2, 3, 3);
		assertThat(target233.length).isEqualTo(6);          // dimU = 2·3
		assertThat(target233[0].length).isEqualTo(9);       // dimV = 3·3
		assertThat(target233[0][0].length).isEqualTo(6);    // dimW = 2·3
	}
}
