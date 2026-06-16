package eu.solven.matmul;

import org.junit.jupiter.api.Test;

import eu.solven.matmul.NonCubicBilinearAlgorithm;

/**
 * Dump the U/V/W of {@code ⟨1,5,1⟩} and its cyclic shifts so the
 * "what does the re-oriented algorithm actually compute" question is
 * answered explicitly.
 */
public class TestSmallOrientationOrbitDump {

	@Test
	public void dump_151_orbit() {
		NonCubicBilinearAlgorithm dot = buildDot();
		System.out.println("\n=== Original ⟨1,5,1⟩ (dot product) ===");
		System.out.println("M_k = A[0,k] · B[k,0]  for k=0..4");
		dump(dot);

		System.out.println("\n=== Cyclic shift → ⟨5,1,1⟩ ===");
		dump(dot.cyclicShift());

		System.out.println("\n=== Cyclic²  → ⟨1,1,5⟩ ===");
		dump(dot.cyclicShift().cyclicShift());
	}

	private static NonCubicBilinearAlgorithm buildDot() {
		int n = 1, m = 5, p = 1, r = 5;
		double[][] U = new double[n * m][r];
		double[][] V = new double[m * p][r];
		double[][] W = new double[n * p][r];
		for (int j = 0; j < 5; j++) {
			U[j][j] = 1;
			V[j][j] = 1;
			W[0][j] = 1;
		}
		return new NonCubicBilinearAlgorithm(n, m, p, U, V, W);
	}

	private static void dump(NonCubicBilinearAlgorithm a) {
		double[][] srcU = a.denseU();
		double[][] srcV = a.denseV();
		double[][] srcW = a.denseW();
		System.out.printf("  n=%d m=%d p=%d r=%d%n", a.n, a.m, a.p, a.r);
		for (int k = 0; k < a.r; k++) {
			StringBuilder alpha = new StringBuilder();
			for (int i = 0; i < a.n; i++) {
				for (int j = 0; j < a.m; j++) {
					double v = srcU[i * a.m + j][k];
					if (v != 0) {
						alpha.append((v > 0 ? "+" : "")).append((int) v)
								.append("·A[").append(i).append(",").append(j).append("] ");
					}
				}
			}
			StringBuilder beta = new StringBuilder();
			for (int j = 0; j < a.m; j++) {
				for (int l = 0; l < a.p; l++) {
					double v = srcV[j * a.p + l][k];
					if (v != 0) {
						beta.append((v > 0 ? "+" : "")).append((int) v)
								.append("·B[").append(j).append(",").append(l).append("] ");
					}
				}
			}
			StringBuilder cs = new StringBuilder();
			for (int i = 0; i < a.n; i++) {
				for (int l = 0; l < a.p; l++) {
					double v = srcW[i * a.p + l][k];
					if (v != 0) {
						cs.append("C[").append(i).append(",").append(l).append("]")
								.append((v > 0 ? "+=" : "-=")).append("M ");
					}
				}
			}
			System.out.printf("  M%d = (%s) · (%s)   used in: %s%n", k, alpha, beta, cs);
		}
	}
}
