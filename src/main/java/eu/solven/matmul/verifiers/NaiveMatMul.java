package eu.solven.matmul.verifiers;
import eu.solven.matmul.NonCubicBilinearAlgorithm;
import eu.solven.matmul.BilinearAlgorithm;

/**
 * The textbook n³-multiplication algorithm:
 *
 *   M_{(i,l,j)} = A[i][l] · B[l][j]
 *   C[i][j]     = Σ_l M_{(i,l,j)}
 *
 * Rank = n³. For n=2: 8 mults. For n=3: 27 mults. Used as a baseline and as a
 * sanity check for {@link Verifier}.
 */
public class NaiveMatMul {

	/**
	 * The textbook {@code n·m·p}-product naive scheme for the rectangular shape
	 * {@code ⟨n,m,p⟩}: {@code M_{(i,l,j)} = A[i][l]·B[l][j]}, {@code C[i][j]=Σ_l}.
	 * Rank {@code n·m·p}. For any shape with a unit axis this is also the OPTIMAL
	 * bilinear rank (no Strassen saving when a dim is 1), so it is the canonical
	 * constructor for degenerate Kronecker factors (e.g. {@code ⟨1,1,k⟩} in
	 * {@code ⟨1,1,3⟩⊗⟨3,3,6⟩=⟨3,3,18⟩}).
	 */
	public static NonCubicBilinearAlgorithm ofNonCubic(int n, int m, int p) {
		int r = n * m * p;
		double[][] U = new double[n * m][r];
		double[][] V = new double[m * p][r];
		double[][] W = new double[n * p][r];
		int k = 0;
		for (int i = 0; i < n; i++) {
			for (int l = 0; l < m; l++) {
				for (int j = 0; j < p; j++) {
					U[i * m + l][k] = 1.0;
					V[l * p + j][k] = 1.0;
					W[i * p + j][k] = 1.0;
					k++;
				}
			}
		}
		return new NonCubicBilinearAlgorithm(n, m, p, U, V, W);
	}

	public static BilinearAlgorithm of(int n) {
		int n2 = n * n;
		int r = n * n * n;
		double[][] U = new double[n2][r];
		double[][] V = new double[n2][r];
		double[][] W = new double[n2][r];

		int k = 0;
		for (int i = 0; i < n; i++) {
			for (int l = 0; l < n; l++) {
				for (int j = 0; j < n; j++) {
					U[i * n + l][k] = 1.0;
					V[l * n + j][k] = 1.0;
					W[i * n + j][k] = 1.0;
					k++;
				}
			}
		}
		return new BilinearAlgorithm(n, U, V, W);
	}
}
