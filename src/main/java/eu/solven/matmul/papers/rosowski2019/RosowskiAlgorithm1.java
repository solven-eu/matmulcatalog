package eu.solven.matmul.papers.rosowski2019;

import eu.solven.matmul.NonBilinearAlgorithm;

/**
 * Rosowski 2019 Algorithm 1 — explicit ⟨n,3,3⟩ = 6n + 3 NON-BILINEAR
 * algorithm over a commutative ring. Generalises {@link Rosowski21}
 * (which is the n=3 cubic case): the 3 B-only products
 * {@code b12·b21, b13·b31, b23·b32} are computed once and reused
 * across all {@code n} rows of A; each row adds 6 products. See
 * {@code references/rosowski-algorithms.md} for the paper's exposition.
 *
 * <p>For each row {@code i ∈ [1..n]}:</p>
 * <pre>
 *   q1(i) = (a_{i,2} + b12)(a_{i,1} + b21)
 *   q2(i) = (a_{i,3} + b13)(a_{i,1} + b31)
 *   q3(i) = (a_{i,3} + b23)(a_{i,2} + b32)
 *   q4(i) = a_{i,1} · (b11 − b12 − b13 − a_{i,2} − a_{i,3})
 *   q5(i) = a_{i,2} · (b22 − b21 − b23 − a_{i,1} − a_{i,3})
 *   q6(i) = a_{i,3} · (b33 − b31 − b32 − a_{i,1} − a_{i,2})
 * </pre>
 *
 * <p>Plus 3 shared B-only products {@code s1 = b12·b21}, {@code s2 = b13·b31},
 * {@code s3 = b23·b32}. Outputs:</p>
 * <pre>
 *   C[i,1] = q4(i) + q1(i) + q2(i) − s1 − s2
 *   C[i,2] = q5(i) + q1(i) + q3(i) − s1 − s3
 *   C[i,3] = q6(i) + q2(i) + q3(i) − s2 − s3
 * </pre>
 *
 * <p>Total rank: {@code 6n + 3}. For n = 3 this is 21, matching Rosowski's
 * Corollary 1.</p>
 */
public final class RosowskiAlgorithm1 {

	private RosowskiAlgorithm1() {}

	/**
	 * Builds the explicit non-bilinear scheme for {@code ⟨n, 3, 3⟩}.
	 *
	 * @param n number of rows of A (and rows of C); must be ≥ 1
	 * @return a non-bilinear algorithm with {@code r = 6n + 3} products
	 */
	public static NonBilinearAlgorithm build(int n) {
		if (n < 1) throw new IllegalArgumentException("n must be ≥ 1, got " + n);
		final int m = 3, p = 3;
		final int r = 6 * n + 3;

		double[][] Ua = new double[n * m][r];
		double[][] Ub = new double[m * p][r];
		double[][] Va = new double[n * m][r];
		double[][] Vb = new double[m * p][r];
		double[][] W  = new double[n * p][r];

		for (int iRow = 1; iRow <= n; iRow++) {
			int base = (iRow - 1) * 6;  // 6 products per row, 0-indexed offset
			int q1 = base, q2 = base + 1, q3 = base + 2;
			int q4 = base + 3, q5 = base + 4, q6 = base + 5;

			// q1 = (a_{i,2} + b12)(a_{i,1} + b21)
			Ua[a(iRow, 2)][q1] = 1; Ub[b(1, 2)][q1] = 1;
			Va[a(iRow, 1)][q1] = 1; Vb[b(2, 1)][q1] = 1;

			// q2 = (a_{i,3} + b13)(a_{i,1} + b31)
			Ua[a(iRow, 3)][q2] = 1; Ub[b(1, 3)][q2] = 1;
			Va[a(iRow, 1)][q2] = 1; Vb[b(3, 1)][q2] = 1;

			// q3 = (a_{i,3} + b23)(a_{i,2} + b32)
			Ua[a(iRow, 3)][q3] = 1; Ub[b(2, 3)][q3] = 1;
			Va[a(iRow, 2)][q3] = 1; Vb[b(3, 2)][q3] = 1;

			// q4 = a_{i,1} · (b11 − b12 − b13 − a_{i,2} − a_{i,3})
			Ua[a(iRow, 1)][q4] = 1;
			Vb[b(1, 1)][q4] = 1; Vb[b(1, 2)][q4] = -1; Vb[b(1, 3)][q4] = -1;
			Va[a(iRow, 2)][q4] = -1; Va[a(iRow, 3)][q4] = -1;

			// q5 = a_{i,2} · (b22 − b21 − b23 − a_{i,1} − a_{i,3})
			Ua[a(iRow, 2)][q5] = 1;
			Vb[b(2, 2)][q5] = 1; Vb[b(2, 1)][q5] = -1; Vb[b(2, 3)][q5] = -1;
			Va[a(iRow, 1)][q5] = -1; Va[a(iRow, 3)][q5] = -1;

			// q6 = a_{i,3} · (b33 − b31 − b32 − a_{i,1} − a_{i,2})
			Ua[a(iRow, 3)][q6] = 1;
			Vb[b(3, 3)][q6] = 1; Vb[b(3, 1)][q6] = -1; Vb[b(3, 2)][q6] = -1;
			Va[a(iRow, 1)][q6] = -1; Va[a(iRow, 2)][q6] = -1;

			// Outputs C[i, l] for l = 1..3
			int s1 = 6 * n, s2 = 6 * n + 1, s3 = 6 * n + 2;
			// C[i,1] = q4 + q1 + q2 − s1 − s2
			W[cRow(iRow, 1)][q4] += 1; W[cRow(iRow, 1)][q1] += 1; W[cRow(iRow, 1)][q2] += 1;
			W[cRow(iRow, 1)][s1] -= 1; W[cRow(iRow, 1)][s2] -= 1;
			// C[i,2] = q5 + q1 + q3 − s1 − s3
			W[cRow(iRow, 2)][q5] += 1; W[cRow(iRow, 2)][q1] += 1; W[cRow(iRow, 2)][q3] += 1;
			W[cRow(iRow, 2)][s1] -= 1; W[cRow(iRow, 2)][s3] -= 1;
			// C[i,3] = q6 + q2 + q3 − s2 − s3
			W[cRow(iRow, 3)][q6] += 1; W[cRow(iRow, 3)][q2] += 1; W[cRow(iRow, 3)][q3] += 1;
			W[cRow(iRow, 3)][s2] -= 1; W[cRow(iRow, 3)][s3] -= 1;
		}

		// Shared B-only products s1, s2, s3 — last 3 indices.
		int s1 = 6 * n, s2 = 6 * n + 1, s3 = 6 * n + 2;
		Ub[b(1, 2)][s1] = 1; Vb[b(2, 1)][s1] = 1;
		Ub[b(1, 3)][s2] = 1; Vb[b(3, 1)][s2] = 1;
		Ub[b(2, 3)][s3] = 1; Vb[b(3, 2)][s3] = 1;

		return new NonBilinearAlgorithm(n, m, p, Ua, Ub, Va, Vb, W);
	}

	// (1-indexed paper convention) → 0-indexed flat array index helpers.
	/** Row index in {@code Ua/Va} for {@code A[i, j]} (1-indexed). */
	private static int a(int i, int j) { return (i - 1) * 3 + (j - 1); }
	/** Row index in {@code Ub/Vb} for {@code B[i, j]} (1-indexed). */
	private static int b(int i, int j) { return (i - 1) * 3 + (j - 1); }
	/** Row index in {@code W} for {@code C[i, l]} (1-indexed; uses p = 3). */
	private static int cRow(int i, int l) { return (i - 1) * 3 + (l - 1); }
}
