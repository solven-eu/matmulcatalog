package eu.solven.matmul.papers.strassen1969;

import eu.solven.matmul.BilinearAlgorithm;

/**
 * Strassen's 1969 algorithm: 7 multiplications for 2×2 matrix multiplication.
 *
 *   M1 = (A11 + A22)(B11 + B22)
 *   M2 = (A21 + A22) · B11
 *   M3 = A11 · (B12 − B22)
 *   M4 = A22 · (B21 − B11)
 *   M5 = (A11 + A12) · B22
 *   M6 = (A21 − A11)(B11 + B12)
 *   M7 = (A12 − A22)(B21 + B22)
 *
 *   C11 = M1 + M4 − M5 + M7
 *   C12 = M3 + M5
 *   C21 = M2 + M4
 *   C22 = M1 − M2 + M3 + M6
 *
 * Indexing: row-major, so Aij is at U[(i-1)*2 + (j-1)].
 */
public class Strassen7 {

	private static final int A11 = 0, A12 = 1, A21 = 2, A22 = 3;
	private static final int B11 = 0, B12 = 1, B21 = 2, B22 = 3;
	private static final int C11 = 0, C12 = 1, C21 = 2, C22 = 3;

	public static BilinearAlgorithm get() {
		double[][] U = new double[4][7];
		double[][] V = new double[4][7];
		double[][] W = new double[4][7];

		// M1 = (A11 + A22)(B11 + B22)
		U[A11][0] = 1; U[A22][0] = 1;
		V[B11][0] = 1; V[B22][0] = 1;
		// M2 = (A21 + A22) · B11
		U[A21][1] = 1; U[A22][1] = 1;
		V[B11][1] = 1;
		// M3 = A11 · (B12 − B22)
		U[A11][2] = 1;
		V[B12][2] = 1; V[B22][2] = -1;
		// M4 = A22 · (B21 − B11)
		U[A22][3] = 1;
		V[B21][3] = 1; V[B11][3] = -1;
		// M5 = (A11 + A12) · B22
		U[A11][4] = 1; U[A12][4] = 1;
		V[B22][4] = 1;
		// M6 = (A21 − A11)(B11 + B12)
		U[A21][5] = 1; U[A11][5] = -1;
		V[B11][5] = 1; V[B12][5] = 1;
		// M7 = (A12 − A22)(B21 + B22)
		U[A12][6] = 1; U[A22][6] = -1;
		V[B21][6] = 1; V[B22][6] = 1;

		// C11 = M1 + M4 − M5 + M7
		W[C11][0] = 1; W[C11][3] = 1; W[C11][4] = -1; W[C11][6] = 1;
		// C12 = M3 + M5
		W[C12][2] = 1; W[C12][4] = 1;
		// C21 = M2 + M4
		W[C21][1] = 1; W[C21][3] = 1;
		// C22 = M1 − M2 + M3 + M6
		W[C22][0] = 1; W[C22][1] = -1; W[C22][2] = 1; W[C22][5] = 1;

		return new BilinearAlgorithm(2, U, V, W);
	}
}
