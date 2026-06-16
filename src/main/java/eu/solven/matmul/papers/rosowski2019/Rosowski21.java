package eu.solven.matmul.papers.rosowski2019;

import eu.solven.matmul.NonBilinearAlgorithm;

/**
 * Rosowski 2019 Corollary 1 — explicit ⟨3,3,3⟩ = 21 NON-BILINEAR
 * algorithm over a commutative ring. Reference:
 * {@code references/rosowski-algorithms.md}.
 *
 * <p>Per-product structure (paper notation, 1-indexed):</p>
 * <pre>
 * p1  := (a12+b12)(a11+b21)        p7  := (a22+b12)(a21+b21)        p13 := (a32+b12)(a31+b21)
 * p2  := (a13+b13)(a11+b31)        p8  := (a23+b13)(a21+b31)        p14 := (a33+b13)(a31+b31)
 * p3  := (a13+b23)(a12+b32)        p9  := (a23+b23)(a22+b32)        p15 := (a33+b23)(a32+b32)
 * p4  := a11 · (b11−b12−b13−a12−a13)  p10 := a21 · (b11−b12−b13−a22−a23)  p16 := a31 · (b11−b12−b13−a32−a33)
 * p5  := a12 · (b22−b21−b23−a11−a13)  p11 := a22 · (b22−b21−b23−a21−a23)  p17 := a32 · (b22−b21−b23−a31−a33)
 * p6  := a13 · (b33−b31−b32−a11−a12)  p12 := a23 · (b33−b31−b32−a21−a22)  p18 := a33 · (b33−b31−b32−a31−a32)
 *
 * p19 := b12 · b21
 * p20 := b13 · b31
 * p21 := b23 · b32
 * </pre>
 *
 * <p>Output (i = 1, 2, 3):</p>
 * <pre>
 * (AB)_{i,1} = p_{6i−2} + p_{6i−5} + p_{6i−4} − p19 − p20
 * (AB)_{i,2} = p_{6i−1} + p_{6i−5} + p_{6i−3} − p19 − p21
 * (AB)_{i,3} = p_{6i}   + p_{6i−4} + p_{6i−3} − p20 − p21
 * </pre>
 */
public final class Rosowski21 {

	private Rosowski21() {}

	/**
	 * Returns the explicit non-bilinear algorithm. 0-indexed throughout
	 * ({@code a_{ij}} in the paper → {@code A[i-1, j-1]} → flat index
	 * {@code (i-1)*3 + (j-1)}).
	 */
	public static NonBilinearAlgorithm build() {
		int n = 3, m = 3, p = 3, r = 21;
		double[][] Ua = new double[n * m][r];
		double[][] Ub = new double[m * p][r];
		double[][] Va = new double[n * m][r];
		double[][] Vb = new double[m * p][r];
		double[][] W  = new double[n * p][r];

		// Helpers: a_{ij}/b_{ij} flatten (1-indexed → 0-indexed).
		// a/b: flat = (i-1)*3 + (j-1)
		// Cell coefficient setters:
		//   uA(k, i, j, c) — A[i,j] coefficient in α_k
		//   uB(k, i, j, c) — B[i,j] coefficient in α_k
		//   vA(k, i, j, c) — A[i,j] coefficient in β_k
		//   vB(k, i, j, c) — B[i,j] coefficient in β_k

		// Block 1: i=1 (paper) — products p1..p6 in our k=0..5
		// p1 = (a12+b12)(a11+b21)
		Ua[idx(1,2)][0] = 1; Ub[idx(1,2)][0] = 1;
		Va[idx(1,1)][0] = 1; Vb[idx(2,1)][0] = 1;
		// p2 = (a13+b13)(a11+b31)
		Ua[idx(1,3)][1] = 1; Ub[idx(1,3)][1] = 1;
		Va[idx(1,1)][1] = 1; Vb[idx(3,1)][1] = 1;
		// p3 = (a13+b23)(a12+b32)
		Ua[idx(1,3)][2] = 1; Ub[idx(2,3)][2] = 1;
		Va[idx(1,2)][2] = 1; Vb[idx(3,2)][2] = 1;
		// p4 = a11 · (b11 − b12 − b13 − a12 − a13)
		Ua[idx(1,1)][3] = 1;
		Vb[idx(1,1)][3] = 1; Vb[idx(1,2)][3] = -1; Vb[idx(1,3)][3] = -1;
		Va[idx(1,2)][3] = -1; Va[idx(1,3)][3] = -1;
		// p5 = a12 · (b22 − b21 − b23 − a11 − a13)
		Ua[idx(1,2)][4] = 1;
		Vb[idx(2,2)][4] = 1; Vb[idx(2,1)][4] = -1; Vb[idx(2,3)][4] = -1;
		Va[idx(1,1)][4] = -1; Va[idx(1,3)][4] = -1;
		// p6 = a13 · (b33 − b31 − b32 − a11 − a12)
		Ua[idx(1,3)][5] = 1;
		Vb[idx(3,3)][5] = 1; Vb[idx(3,1)][5] = -1; Vb[idx(3,2)][5] = -1;
		Va[idx(1,1)][5] = -1; Va[idx(1,2)][5] = -1;

		// Block 2: i=2 — products p7..p12 in our k=6..11
		Ua[idx(2,2)][6] = 1; Ub[idx(1,2)][6] = 1;
		Va[idx(2,1)][6] = 1; Vb[idx(2,1)][6] = 1;
		Ua[idx(2,3)][7] = 1; Ub[idx(1,3)][7] = 1;
		Va[idx(2,1)][7] = 1; Vb[idx(3,1)][7] = 1;
		Ua[idx(2,3)][8] = 1; Ub[idx(2,3)][8] = 1;
		Va[idx(2,2)][8] = 1; Vb[idx(3,2)][8] = 1;
		Ua[idx(2,1)][9] = 1;
		Vb[idx(1,1)][9] = 1; Vb[idx(1,2)][9] = -1; Vb[idx(1,3)][9] = -1;
		Va[idx(2,2)][9] = -1; Va[idx(2,3)][9] = -1;
		Ua[idx(2,2)][10] = 1;
		Vb[idx(2,2)][10] = 1; Vb[idx(2,1)][10] = -1; Vb[idx(2,3)][10] = -1;
		Va[idx(2,1)][10] = -1; Va[idx(2,3)][10] = -1;
		Ua[idx(2,3)][11] = 1;
		Vb[idx(3,3)][11] = 1; Vb[idx(3,1)][11] = -1; Vb[idx(3,2)][11] = -1;
		Va[idx(2,1)][11] = -1; Va[idx(2,2)][11] = -1;

		// Block 3: i=3 — products p13..p18 in our k=12..17
		Ua[idx(3,2)][12] = 1; Ub[idx(1,2)][12] = 1;
		Va[idx(3,1)][12] = 1; Vb[idx(2,1)][12] = 1;
		Ua[idx(3,3)][13] = 1; Ub[idx(1,3)][13] = 1;
		Va[idx(3,1)][13] = 1; Vb[idx(3,1)][13] = 1;
		Ua[idx(3,3)][14] = 1; Ub[idx(2,3)][14] = 1;
		Va[idx(3,2)][14] = 1; Vb[idx(3,2)][14] = 1;
		Ua[idx(3,1)][15] = 1;
		Vb[idx(1,1)][15] = 1; Vb[idx(1,2)][15] = -1; Vb[idx(1,3)][15] = -1;
		Va[idx(3,2)][15] = -1; Va[idx(3,3)][15] = -1;
		Ua[idx(3,2)][16] = 1;
		Vb[idx(2,2)][16] = 1; Vb[idx(2,1)][16] = -1; Vb[idx(2,3)][16] = -1;
		Va[idx(3,1)][16] = -1; Va[idx(3,3)][16] = -1;
		Ua[idx(3,3)][17] = 1;
		Vb[idx(3,3)][17] = 1; Vb[idx(3,1)][17] = -1; Vb[idx(3,2)][17] = -1;
		Va[idx(3,1)][17] = -1; Va[idx(3,2)][17] = -1;

		// Shared B-only products p19..p21 — k=18..20
		Ub[idx(1,2)][18] = 1; Vb[idx(2,1)][18] = 1;  // p19 = b12 * b21
		Ub[idx(1,3)][19] = 1; Vb[idx(3,1)][19] = 1;  // p20 = b13 * b31
		Ub[idx(2,3)][20] = 1; Vb[idx(3,2)][20] = 1;  // p21 = b23 * b32

		// W matrix — for each output C[i,l] (W row (i-1)*3 + (l-1)), set coefficients.
		// C[1,1] = p4+p1+p2−p19−p20
		setW(W, 1, 1, new int[]{4, 1, 2}, new int[]{19, 20});
		// C[1,2] = p5+p1+p3−p19−p21
		setW(W, 1, 2, new int[]{5, 1, 3}, new int[]{19, 21});
		// C[1,3] = p6+p2+p3−p20−p21
		setW(W, 1, 3, new int[]{6, 2, 3}, new int[]{20, 21});
		// C[2,1] = p10+p7+p8−p19−p20
		setW(W, 2, 1, new int[]{10, 7, 8}, new int[]{19, 20});
		// C[2,2] = p11+p7+p9−p19−p21
		setW(W, 2, 2, new int[]{11, 7, 9}, new int[]{19, 21});
		// C[2,3] = p12+p8+p9−p20−p21
		setW(W, 2, 3, new int[]{12, 8, 9}, new int[]{20, 21});
		// C[3,1] = p16+p13+p14−p19−p20
		setW(W, 3, 1, new int[]{16, 13, 14}, new int[]{19, 20});
		// C[3,2] = p17+p13+p15−p19−p21
		setW(W, 3, 2, new int[]{17, 13, 15}, new int[]{19, 21});
		// C[3,3] = p18+p14+p15−p20−p21
		setW(W, 3, 3, new int[]{18, 14, 15}, new int[]{20, 21});

		return new NonBilinearAlgorithm(n, m, p, Ua, Ub, Va, Vb, W);
	}

	/** Paper {@code a_{ij}/b_{ij}} (1-indexed) → flat row index (0-indexed). */
	private static int idx(int i, int j) { return (i - 1) * 3 + (j - 1); }

	/** Set output row for {@code C_{i,l}} (1-indexed paper, 0-indexed array). */
	private static void setW(double[][] W, int i, int l, int[] plus, int[] minus) {
		int row = (i - 1) * 3 + (l - 1);
		for (int pk : plus)  W[row][pk - 1] += 1;
		for (int mk : minus) W[row][mk - 1] -= 1;
	}
}
