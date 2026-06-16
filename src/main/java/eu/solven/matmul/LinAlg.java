package eu.solven.matmul;

/**
 * Minimal dense-matrix helpers for ALS: matrix–matrix multiply and Gauss–Jordan
 * inverse with partial pivoting. Adequate because all matrices that appear in the
 * ALS Gram solves here are tiny (rank-r square, r ≤ ~30).
 */
public final class LinAlg {

	private LinAlg() {}

	/** Returns A·B for A: m×k, B: k×n. */
	public static double[][] mul(double[][] A, double[][] B) {
		int m = A.length;
		int k = A[0].length;
		int n = B[0].length;
		double[][] C = new double[m][n];
		for (int i = 0; i < m; i++) {
			double[] Ai = A[i];
			double[] Ci = C[i];
			for (int p = 0; p < k; p++) {
				double aip = Ai[p];
				if (aip == 0.0) continue;
				double[] Bp = B[p];
				for (int j = 0; j < n; j++) {
					Ci[j] += aip * Bp[j];
				}
			}
		}
		return C;
	}

	/**
	 * Returns the inverse of a square matrix via Gauss–Jordan elimination with
	 * partial pivoting. Caller is responsible for adding any ridge regularization
	 * to keep the matrix invertible — this routine throws on singular input.
	 */
	public static double[][] invert(double[][] A) {
		int n = A.length;
		double[][] aug = new double[n][2 * n];
		for (int i = 0; i < n; i++) {
			System.arraycopy(A[i], 0, aug[i], 0, n);
			aug[i][n + i] = 1.0;
		}
		for (int col = 0; col < n; col++) {
			int pivot = col;
			for (int row = col + 1; row < n; row++) {
				if (Math.abs(aug[row][col]) > Math.abs(aug[pivot][col])) {
					pivot = row;
				}
			}
			if (Math.abs(aug[pivot][col]) < 1e-14) {
				throw new ArithmeticException("singular matrix at column " + col);
			}
			double[] tmp = aug[col]; aug[col] = aug[pivot]; aug[pivot] = tmp;
			double diag = aug[col][col];
			double[] pivRow = aug[col];
			for (int j = 0; j < 2 * n; j++) pivRow[j] /= diag;
			for (int row = 0; row < n; row++) {
				if (row == col) continue;
				double f = aug[row][col];
				if (f == 0.0) continue;
				double[] target = aug[row];
				for (int j = 0; j < 2 * n; j++) target[j] -= f * pivRow[j];
			}
		}
		double[][] inv = new double[n][n];
		for (int i = 0; i < n; i++) System.arraycopy(aug[i], n, inv[i], 0, n);
		return inv;
	}
}
