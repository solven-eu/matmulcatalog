package eu.solven.matmul.f2.sat;

import static org.assertj.core.api.Assertions.assertThat;

import java.util.Optional;

import org.junit.jupiter.api.Test;

import eu.solven.matmul.BilinearAlgorithm;
import eu.solven.matmul.f2.sat.SatMatmulPipeline;

public class TestSatMatmulPipeline {

	// ---- dense ⟨2,2,2⟩ over Z/2 -------------------------------------------

	@Test
	public void denseZ2_2x2_r7_isSat() {
		int[][][] target = SatMatmulPipeline.z2DenseMatmulTensor(2);
		Optional<BilinearAlgorithm> result = SatMatmulPipeline.findZ2Decomposition(2, 7, target);
		assertThat(result).isPresent();
		assertThat(SatMatmulPipeline.verifyZ2(result.get(), target)).isTrue();
	}

	@Test
	public void denseZ2_2x2_r6_isUnsat() {
		// R(⟨2,2,2⟩) = 7 over any ring (Hopcroft–Kerr 1971); Z/2 included.
		int[][][] target = SatMatmulPipeline.z2DenseMatmulTensor(2);
		Optional<BilinearAlgorithm> result = SatMatmulPipeline.findZ2Decomposition(2, 6, target);
		assertThat(result).isEmpty();
	}

	// ---- upper-triangular ⟨2,2,2⟩ over Z/2 ---------------------------------

	@Test
	public void triangularZ2_2x2_r4_isSat() {
		// A and B upper-triangular: C = AB needs exactly 4 multiplications.
		int[][][] target = SatMatmulPipeline.z2UpperTriangularMatmulTensor(2);
		Optional<BilinearAlgorithm> result = SatMatmulPipeline.findZ2Decomposition(2, 4, target);
		assertThat(result).isPresent();
		assertThat(SatMatmulPipeline.verifyZ2(result.get(), target)).isTrue();
	}

	@Test
	public void triangularZ2_2x2_r3_isUnsat() {
		// 4 mults is tight for upper-triangular ⟨2,2,2⟩ — r=3 should be UNSAT.
		int[][][] target = SatMatmulPipeline.z2UpperTriangularMatmulTensor(2);
		Optional<BilinearAlgorithm> result = SatMatmulPipeline.findZ2Decomposition(2, 3, target);
		assertThat(result).isEmpty();
	}
}
