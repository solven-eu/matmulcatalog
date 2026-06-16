package eu.solven.matmul;

/**
 * A bilinear algorithm for n×n matrix multiplication, expressed as a rank-r tensor
 * decomposition. For k = 0..r-1:
 *
 *   M_k = (Σ_{a,b} U[a*n+b][k] · A[a][b]) · (Σ_{c,d} V[c*n+d][k] · B[c][d])
 *
 *   C[i][j] = Σ_k W[i*n+j][k] · M_k
 *
 * Indexing convention: row-major flatten of n×n matrices to length-n² vectors.
 * Index (a,b) → a*n + b.
 */
public class BilinearAlgorithm {

	public final int n;
	public final int r;
	public final double[][] U; // shape [n²][r]
	public final double[][] V; // shape [n²][r]
	public final double[][] W; // shape [n²][r]

	public BilinearAlgorithm(int n, double[][] U, double[][] V, double[][] W) {
		int n2 = n * n;
		if (U.length != n2 || V.length != n2 || W.length != n2) {
			throw new IllegalArgumentException("U/V/W must have n² = " + n2 + " rows");
		}
		int rank = U[0].length;
		for (int a = 0; a < n2; a++) {
			if (U[a].length != rank || V[a].length != rank || W[a].length != rank) {
				throw new IllegalArgumentException("U/V/W must all have rank columns");
			}
		}
		this.n = n;
		this.r = rank;
		this.U = U;
		this.V = V;
		this.W = W;
	}
}
