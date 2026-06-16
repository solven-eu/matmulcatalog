package eu.solven.matmul.f2.sat;

import static org.assertj.core.api.Assertions.assertThat;

import java.util.Optional;

import org.junit.jupiter.api.Test;

import eu.solven.matmul.BilinearAlgorithm;
import eu.solven.matmul.f2.sat.SatMatmulPipeline;

/**
 * Phase 1.5 sanity tests: restricted ⟨3,3,3⟩ over Z/2 at n=3 scale, with
 * inputs/outputs constrained to a small position subset. Validates that the
 * SAT pipeline works correctly at the n=3 index space without paying the
 * full ⟨3,3,3⟩ search cost.
 *
 * <p><b>Diagonal + top-right (0, 2) case</b>: positions S = {(0,0), (0,2),
 * (1,1), (2,2)}. The restricted matmul has structure:
 * <pre>
 *   C[0,0] = A[0,0]·B[0,0]
 *   C[0,2] = A[0,0]·B[0,2] + A[0,2]·B[2,2]
 *   C[1,1] = A[1,1]·B[1,1]
 *   C[2,2] = A[2,2]·B[2,2]
 * </pre>
 * Provable rank by direct-sum decomposition: the {0, 2} sub-block is
 * triangular ⟨2,2,2⟩ with rank 4; the isolated cell C[1,1] adds 1.
 * So R = 5 exactly.</p>
 */
public class TestSatRestricted3x3 {

	@Test
	public void diagonalPlusTopRight_r5_isSat() {
		boolean[] allowed = SatMatmulPipeline.diagonalPlusOne(3, 0, 2);
		int[][][] target = SatMatmulPipeline.z2RestrictedMatmulTensor(3, allowed);
		Optional<BilinearAlgorithm> result =
				SatMatmulPipeline.findZ2RestrictedDecomposition(3, 5, target, allowed);
		assertThat(result).isPresent();
		assertThat(SatMatmulPipeline.verifyZ2(result.get(), target)).isTrue();
	}

	@Test
	public void diagonalPlusTopRight_r4_isUnsat() {
		// R = 5 by direct-sum: triangular ⟨2,2,2⟩ sub-block (rank 4) +
		// isolated C[1,1]=A[1,1]·B[1,1] (rank 1). So r=4 must be UNSAT.
		boolean[] allowed = SatMatmulPipeline.diagonalPlusOne(3, 0, 2);
		int[][][] target = SatMatmulPipeline.z2RestrictedMatmulTensor(3, allowed);
		Optional<BilinearAlgorithm> result =
				SatMatmulPipeline.findZ2RestrictedDecomposition(3, 4, target, allowed);
		assertThat(result).isEmpty();
	}
}
