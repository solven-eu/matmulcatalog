package eu.solven.matmul.research;

import eu.solven.matmul.NonBilinearAlgorithm;
import eu.solven.matmul.Verifier;
import eu.solven.matmul.papers.makarov1986.Makarov22;

/**
 * Search for a single sign-flip in one γ_k's definition (one entry of
 * Ua/Ub/Va/Vb) that brings the residual to 0. Complements
 * {@link MakarovSearch} which flips signs in the W output combinations.
 */
public final class MakarovGammaSignSearch {
	private MakarovGammaSignSearch() {}

	public static void main(String[] args) {
		NonBilinearAlgorithm base = Makarov22.buildDefault();
		int rank = base.r;
		int dimA = base.n * base.m, dimB = base.m * base.p;

		int hits = 0;
		System.out.println("Single-sign-flip in γ_k definitions (Ua / Ub / Va / Vb):");

		// Try flipping each non-zero entry in each of the 4 factor matrices
		for (int sideLabel = 0; sideLabel < 4; sideLabel++) {
			double[][] target;
			int dim;
			String label;
			switch (sideLabel) {
				case 0 -> { target = base.Ua; dim = dimA; label = "Ua"; }
				case 1 -> { target = base.Ub; dim = dimB; label = "Ub"; }
				case 2 -> { target = base.Va; dim = dimA; label = "Va"; }
				case 3 -> { target = base.Vb; dim = dimB; label = "Vb"; }
				default -> throw new IllegalStateException();
			}
			for (int row = 0; row < dim; row++) {
				for (int k = 0; k < rank; k++) {
					double orig = target[row][k];
					if (orig == 0) continue;
					target[row][k] = -orig;
					double r = Verifier.residualNonBilinear(base);
					if (r < 1e-9) {
						System.out.printf("  HIT: flip sign of %s[row=%d][γ%d] (was %+g) → residual %.4g%n",
								label, row, k + 1, orig, r);
						hits++;
					}
					target[row][k] = orig;  // restore
				}
			}
		}
		System.out.printf("Total γ-definition single-flip hits: %d%n", hits);
	}
}
