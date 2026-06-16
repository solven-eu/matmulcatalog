package eu.solven.matmul;

import static org.assertj.core.api.Assertions.assertThat;

import org.junit.jupiter.api.Test;

import eu.solven.matmul.NonCubicBilinearAlgorithm;
import eu.solven.matmul.Verifier;

/**
 * Empirically confirms the S₃ orbit equivalence for the user's
 * concrete test case: ⟨1,5,1⟩ (dot product) ↔ ⟨1,1,5⟩ (scalar·vector)
 * ↔ ⟨5,1,1⟩ (vector·scalar). All three have rank 5; a tensor
 * decomposition of one re-orients to a valid decomposition of any
 * other in the orbit, AND the re-oriented algorithm runs correctly
 * on inputs of its new format.
 */
public class TestSmallOrientationOrbit {

	/** Trivial ⟨1,5,1⟩ rank-5 algorithm: c[0,0] = Σ A[0,j]·B[j,0]. */
	private static NonCubicBilinearAlgorithm buildDotProduct() {
		int n = 1, m = 5, p = 1, r = 5;
		double[][] U = new double[n * m][r]; // [5][5]
		double[][] V = new double[m * p][r]; // [5][5]
		double[][] W = new double[n * p][r]; // [1][5]
		for (int j = 0; j < 5; j++) {
			U[j][j] = 1;  // A[0,j] used by product k=j
			V[j][j] = 1;  // B[j,0] used by product k=j
			W[0][j] = 1;  // C[0,0] = Σ products
		}
		return new NonCubicBilinearAlgorithm(n, m, p, U, V, W);
	}

	@Test
	public void orbit_of_1_5_1() {
		NonCubicBilinearAlgorithm dot151 = buildDotProduct();
		assertThat(Verifier.passesRandomMatmulSpotCheck(dot151))
				.as("⟨1,5,1⟩ original (dot product) verifies").isTrue();

		// Cyclic shift: ⟨1,5,1⟩ → ⟨5,1,1⟩
		NonCubicBilinearAlgorithm a511 = dot151.cyclicShift();
		assertThat(a511.n).isEqualTo(5);
		assertThat(a511.m).isEqualTo(1);
		assertThat(a511.p).isEqualTo(1);
		assertThat(a511.r).isEqualTo(5);
		assertThat(Verifier.passesRandomMatmulSpotCheck(a511))
				.as("⟨5,1,1⟩ re-oriented from ⟨1,5,1⟩ verifies as 5×1·1×1=5×1 column-vector·scalar")
				.isTrue();

		// Second cyclic shift: ⟨5,1,1⟩ → ⟨1,1,5⟩
		NonCubicBilinearAlgorithm a115 = a511.cyclicShift();
		assertThat(a115.n).isEqualTo(1);
		assertThat(a115.m).isEqualTo(1);
		assertThat(a115.p).isEqualTo(5);
		assertThat(a115.r).isEqualTo(5);
		assertThat(Verifier.passesRandomMatmulSpotCheck(a115))
				.as("⟨1,1,5⟩ re-oriented from ⟨1,5,1⟩ verifies as 1×1·1×5=1×5 scalar·row-vector")
				.isTrue();

		// Third cyclic = back to original
		NonCubicBilinearAlgorithm back = a115.cyclicShift();
		assertThat(back.n).isEqualTo(1);
		assertThat(back.m).isEqualTo(5);
		assertThat(back.p).isEqualTo(1);

		// orientAs covers the full S₃ — should find ⟨1,1,5⟩ from ⟨1,5,1⟩.
		var oriented = dot151.orientAs(1, 1, 5);
		assertThat(oriented).isPresent();
		assertThat(Verifier.passesRandomMatmulSpotCheck(oriented.get())).isTrue();
	}
}
