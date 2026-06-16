package eu.solven.matmul.papers.laderman1976;

import eu.solven.matmul.BilinearAlgorithm;

/**
 * Laderman's 1976 algorithm: 23 multiplications for 3×3 matrix multiplication.
 *
 * Source: J. Laderman, "A noncommutative algorithm for multiplying (3×3) matrices using
 * 23 multiplications", Bull. AMS 82(1):126–128, 1976. Reproduced in the form used by
 * Heun/Smirnov (Heun, "A New General-Purpose Method to Multiply 3×3 Matrices Using
 * Only 23 Multiplications", arXiv:1108.2830, §2.4) — signs absorbed into a few of the
 * products and into c11/c31, equivalent to Laderman's original up to ± rearrangement.
 *
 * Indexing: row-major, a_{ij} at U[(i-1)*3 + (j-1)]. Same for B and C.
 */
public class Laderman23 {

	private static final int A11 = 0, A12 = 1, A13 = 2;
	private static final int A21 = 3, A22 = 4, A23 = 5;
	private static final int A31 = 6, A32 = 7, A33 = 8;

	private static final int B11 = 0, B12 = 1, B13 = 2;
	private static final int B21 = 3, B22 = 4, B23 = 5;
	private static final int B31 = 6, B32 = 7, B33 = 8;

	private static final int C11 = 0, C12 = 1, C13 = 2;
	private static final int C21 = 3, C22 = 4, C23 = 5;
	private static final int C31 = 6, C32 = 7, C33 = 8;

	public static BilinearAlgorithm get() {
		double[][] U = new double[9][23];
		double[][] V = new double[9][23];
		double[][] W = new double[9][23];

		// P1 = (a11 - a12 - a13 + a21 - a22 - a32 - a33) · (-b22)
		U[A11][0] = 1; U[A12][0] = -1; U[A13][0] = -1;
		U[A21][0] = 1; U[A22][0] = -1; U[A32][0] = -1; U[A33][0] = -1;
		V[B22][0] = -1;
		// P2 = (a11 + a21) · (b12 + b22)
		U[A11][1] = 1; U[A21][1] = 1;
		V[B12][1] = 1; V[B22][1] = 1;
		// P3 = a22 · (b11 - b12 + b21 - b22 - b23 + b31 - b33)
		U[A22][2] = 1;
		V[B11][2] = 1; V[B12][2] = -1; V[B21][2] = 1; V[B22][2] = -1;
		V[B23][2] = -1; V[B31][2] = 1; V[B33][2] = -1;
		// P4 = (-a11 - a21 + a22) · (-b11 + b12 + b22)
		U[A11][3] = -1; U[A21][3] = -1; U[A22][3] = 1;
		V[B11][3] = -1; V[B12][3] = 1; V[B22][3] = 1;
		// P5 = (-a21 + a22) · (-b11 + b12)
		U[A21][4] = -1; U[A22][4] = 1;
		V[B11][4] = -1; V[B12][4] = 1;
		// P6 = a11 · (-b11)
		U[A11][5] = 1;
		V[B11][5] = -1;
		// P7 = (a11 + a31 + a32) · (b11 - b13 + b23)
		U[A11][6] = 1; U[A31][6] = 1; U[A32][6] = 1;
		V[B11][6] = 1; V[B13][6] = -1; V[B23][6] = 1;
		// P8 = (a11 + a31) · (-b13 + b23)
		U[A11][7] = 1; U[A31][7] = 1;
		V[B13][7] = -1; V[B23][7] = 1;
		// P9 = (a31 + a32) · (b11 - b13)
		U[A31][8] = 1; U[A32][8] = 1;
		V[B11][8] = 1; V[B13][8] = -1;
		// P10 = (a11 + a12 - a13 - a22 + a23 + a31 + a32) · b23
		U[A11][9] = 1; U[A12][9] = 1; U[A13][9] = -1; U[A22][9] = -1;
		U[A23][9] = 1; U[A31][9] = 1; U[A32][9] = 1;
		V[B23][9] = 1;
		// P11 = a32 · (-b11 + b13 + b21 - b22 - b23 - b31 + b32)
		U[A32][10] = 1;
		V[B11][10] = -1; V[B13][10] = 1; V[B21][10] = 1; V[B22][10] = -1;
		V[B23][10] = -1; V[B31][10] = -1; V[B32][10] = 1;
		// P12 = (a13 + a32 + a33) · (b22 + b31 - b32)
		U[A13][11] = 1; U[A32][11] = 1; U[A33][11] = 1;
		V[B22][11] = 1; V[B31][11] = 1; V[B32][11] = -1;
		// P13 = (a13 + a33) · (-b22 + b32)
		U[A13][12] = 1; U[A33][12] = 1;
		V[B22][12] = -1; V[B32][12] = 1;
		// P14 = a13 · b31
		U[A13][13] = 1;
		V[B31][13] = 1;
		// P15 = (-a32 - a33) · (-b31 + b32)
		U[A32][14] = -1; U[A33][14] = -1;
		V[B31][14] = -1; V[B32][14] = 1;
		// P16 = (a13 + a22 - a23) · (b23 - b31 + b33)
		U[A13][15] = 1; U[A22][15] = 1; U[A23][15] = -1;
		V[B23][15] = 1; V[B31][15] = -1; V[B33][15] = 1;
		// P17 = (-a13 + a23) · (b23 + b33)
		U[A13][16] = -1; U[A23][16] = 1;
		V[B23][16] = 1; V[B33][16] = 1;
		// P18 = (a22 - a23) · (b31 - b33)
		U[A22][17] = 1; U[A23][17] = -1;
		V[B31][17] = 1; V[B33][17] = -1;
		// P19 = a12 · b21
		U[A12][18] = 1;
		V[B21][18] = 1;
		// P20 = a23 · b32
		U[A23][19] = 1;
		V[B32][19] = 1;
		// P21 = a21 · b13
		U[A21][20] = 1;
		V[B13][20] = 1;
		// P22 = a31 · b12
		U[A31][21] = 1;
		V[B12][21] = 1;
		// P23 = a33 · b33
		U[A33][22] = 1;
		V[B33][22] = 1;

		// c11 = -P6 + P14 + P19
		W[C11][5] = -1; W[C11][13] = 1; W[C11][18] = 1;
		// c12 = P1 - P4 + P5 - P6 - P12 + P14 + P15
		W[C12][0] = 1; W[C12][3] = -1; W[C12][4] = 1; W[C12][5] = -1;
		W[C12][11] = -1; W[C12][13] = 1; W[C12][14] = 1;
		// c13 = -P6 - P7 + P9 + P10 + P14 + P16 + P18
		W[C13][5] = -1; W[C13][6] = -1; W[C13][8] = 1; W[C13][9] = 1;
		W[C13][13] = 1; W[C13][15] = 1; W[C13][17] = 1;
		// c21 = P2 + P3 + P4 + P6 + P14 + P16 + P17
		W[C21][1] = 1; W[C21][2] = 1; W[C21][3] = 1; W[C21][5] = 1;
		W[C21][13] = 1; W[C21][15] = 1; W[C21][16] = 1;
		// c22 = P2 + P4 - P5 + P6 + P20
		W[C22][1] = 1; W[C22][3] = 1; W[C22][4] = -1; W[C22][5] = 1; W[C22][19] = 1;
		// c23 = P14 + P16 + P17 + P18 + P21
		W[C23][13] = 1; W[C23][15] = 1; W[C23][16] = 1; W[C23][17] = 1; W[C23][20] = 1;
		// c31 = P6 + P7 - P8 + P11 + P12 + P13 - P14
		W[C31][5] = 1; W[C31][6] = 1; W[C31][7] = -1; W[C31][10] = 1;
		W[C31][11] = 1; W[C31][12] = 1; W[C31][13] = -1;
		// c32 = P12 + P13 - P14 - P15 + P22
		W[C32][11] = 1; W[C32][12] = 1; W[C32][13] = -1; W[C32][14] = -1; W[C32][21] = 1;
		// c33 = P6 + P7 - P8 - P9 + P23
		W[C33][5] = 1; W[C33][6] = 1; W[C33][7] = -1; W[C33][8] = -1; W[C33][22] = 1;

		return new BilinearAlgorithm(3, U, V, W);
	}
}
