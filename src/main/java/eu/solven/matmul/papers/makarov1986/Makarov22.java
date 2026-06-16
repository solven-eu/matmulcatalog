package eu.solven.matmul.papers.makarov1986;

import eu.solven.matmul.NonBilinearAlgorithm;

/**
 * Makarov's explicit ⟨3,3,3⟩ = 22 NON-BILINEAR commutative scheme.
 * One better than Laderman 1976's bilinear ⟨3,3,3⟩=23 — achievable
 * because Makarov's factors mix entries of A and B inside each
 * rank-1 product (requires scalar commutativity; does NOT lift to
 * recursive matmul over non-commutative rings).
 *
 * <p>Sourced from Makarov 1986, "An algorithm for multiplication of
 * 3×3 matrices", <em>USSR Comput. Math. Math. Phys.</em>
 * 26(2):293–294 (DOI:10.1016/0041-5553(86)90203-X). The 22 products
 * γ_1..γ_22 and the output combinations c_{i,j} below transcribe
 * Makarov 1986's M_1..M_22 and r_1..r_9 respectively, using
 * {@code a_{ij}=A[i,j]} and {@code b_{ij}=B[i,j]} in place of
 * Makarov's a/b/c/k mixed naming.</p>
 *
 * <p>Note: Islam 2009 §3.3.1 transcribes Makarov's algorithm with a
 * single-index typo in γ18 (the second factor reads {@code b32}
 * where the Russian original has {@code k_6 = b_{2,3} = b23}).
 * Verified by total-residual = 0 sympy / Verifier check after
 * applying the index correction.</p>
 */
public final class Makarov22 {

	private Makarov22() {}

	private static int idx(int i, int j) { return (i - 1) * 3 + (j - 1); }

	/**
	 * Builds the algorithm with the W output combinations supplied by
	 * {@code outputs}, one per (i,j) ∈ {1,2,3}². Each {@code outputs[i-1][j-1]}
	 * is a {@code WSpec}: list of signed γ indices that sum to c_{i,j}.
	 *
	 * <p>Pass {@link #defaultOutputs()} to get Islam 2009's exact transcribed
	 * formula; pass a hand-edited variant to test typo hypotheses.</p>
	 */
	public static NonBilinearAlgorithm build(WSpec[][] outputs) {
		final int n = 3, m = 3, p = 3, r = 22;
		double[][] Ua = new double[n * m][r];
		double[][] Ub = new double[m * p][r];
		double[][] Va = new double[n * m][r];
		double[][] Vb = new double[m * p][r];
		double[][] W  = new double[n * p][r];

		// γ1 = (b13 + a13 − a23)(a11 + b31 − b32 + b33)
		Ua[idx(1,3)][0] = 1;  Ua[idx(2,3)][0] = -1;
		Ub[idx(1,3)][0] = 1;
		Va[idx(1,1)][0] = 1;
		Vb[idx(3,1)][0] = 1;  Vb[idx(3,2)][0] = -1;  Vb[idx(3,3)][0] = 1;
		// γ2 = (b12 + a12 + a22)(a21 − b21 + b22 − b23)
		Ua[idx(1,2)][1] = 1;  Ua[idx(2,2)][1] = 1;
		Ub[idx(1,2)][1] = 1;
		Va[idx(2,1)][1] = 1;
		Vb[idx(2,1)][1] = -1; Vb[idx(2,2)][1] = 1;   Vb[idx(2,3)][1] = -1;
		// γ3 = (b12 + a12 + a32)(a31 − b21 + b22 − b23)
		Ua[idx(1,2)][2] = 1;  Ua[idx(3,2)][2] = 1;
		Ub[idx(1,2)][2] = 1;
		Va[idx(3,1)][2] = 1;
		Vb[idx(2,1)][2] = -1; Vb[idx(2,2)][2] = 1;   Vb[idx(2,3)][2] = -1;
		// γ4 = (b13 − a23 − a33)(a31 − b31 + b32 − b33)
		Ua[idx(2,3)][3] = -1; Ua[idx(3,3)][3] = -1;
		Ub[idx(1,3)][3] = 1;
		Va[idx(3,1)][3] = 1;
		Vb[idx(3,1)][3] = -1; Vb[idx(3,2)][3] = 1;   Vb[idx(3,3)][3] = -1;
		// γ5 = (b11 − a13 + a23)·a11
		Ua[idx(1,3)][4] = -1; Ua[idx(2,3)][4] = 1;
		Ub[idx(1,1)][4] = 1;
		Va[idx(1,1)][4] = 1;
		// γ6 = (b11 + a12 + a22)·a21
		Ua[idx(1,2)][5] = 1;  Ua[idx(2,2)][5] = 1;
		Ub[idx(1,1)][5] = 1;
		Va[idx(2,1)][5] = 1;
		// γ7 = (b11 + a12 + a32 + a23 + a33)·a31
		Ua[idx(1,2)][6] = 1;  Ua[idx(3,2)][6] = 1;  Ua[idx(2,3)][6] = 1;  Ua[idx(3,3)][6] = 1;
		Ub[idx(1,1)][6] = 1;
		Va[idx(3,1)][6] = 1;
		// γ8 = b12·(a11 + b21 − b22 + b23)
		Ub[idx(1,2)][7] = 1;
		Va[idx(1,1)][7] = 1;
		Vb[idx(2,1)][7] = 1;  Vb[idx(2,2)][7] = -1; Vb[idx(2,3)][7] = 1;
		// γ9 = b13·(a21 + b31 − b32 + b33)
		Ub[idx(1,3)][8] = 1;
		Va[idx(2,1)][8] = 1;
		Vb[idx(3,1)][8] = 1;  Vb[idx(3,2)][8] = -1; Vb[idx(3,3)][8] = 1;
		// γ10 = a12·b21
		Ua[idx(1,2)][9] = 1;
		Vb[idx(2,1)][9] = 1;
		// γ11 = a23·b31
		Ua[idx(2,3)][10] = 1;
		Vb[idx(3,1)][10] = 1;
		// γ12 = (a13 − a23)(a11 + b31)
		Ua[idx(1,3)][11] = 1;  Ua[idx(2,3)][11] = -1;
		Va[idx(1,1)][11] = 1;
		Vb[idx(3,1)][11] = 1;
		// γ13 = (a12 + a22)(b21 − a21)
		Ua[idx(1,2)][12] = 1;  Ua[idx(2,2)][12] = 1;
		Va[idx(2,1)][12] = -1;
		Vb[idx(2,1)][12] = 1;
		// γ14 = (a12 + b12)(b21 − b22 + b23)
		Ua[idx(1,2)][13] = 1;
		Ub[idx(1,2)][13] = 1;
		Vb[idx(2,1)][13] = 1;  Vb[idx(2,2)][13] = -1; Vb[idx(2,3)][13] = 1;
		// γ15 = a22·b23
		Ua[idx(2,2)][14] = 1;
		Vb[idx(2,3)][14] = 1;
		// γ16 = (b13 − a23)(b31 − b32 + b33)
		Ua[idx(2,3)][15] = -1;
		Ub[idx(1,3)][15] = 1;
		Vb[idx(3,1)][15] = 1;  Vb[idx(3,2)][15] = -1; Vb[idx(3,3)][15] = 1;
		// γ17 = a23·b32
		Ua[idx(2,3)][16] = 1;
		Vb[idx(3,2)][16] = 1;
		// γ18 = (a32 − a23 − a33)·b23
		// (Islam 2009 §3.3.1 transcribes the second factor as b32 — that is
		// a typo; cross-checked against Makarov 1986 Russian original where
		// M18 = (b_3 − c_2 − c_3)·k_6 and k_6 maps to b_{2,3} in our
		// notation, not k_8 = b_{3,2}.)
		Ua[idx(3,2)][17] = 1;  Ua[idx(2,3)][17] = -1;  Ua[idx(3,3)][17] = -1;
		Vb[idx(2,3)][17] = 1;
		// γ19 = (a13 + a33 − a12 − a32)·b32
		Ua[idx(1,3)][18] = 1;  Ua[idx(3,3)][18] = 1;  Ua[idx(1,2)][18] = -1; Ua[idx(3,2)][18] = -1;
		Vb[idx(3,2)][18] = 1;
		// γ20 = (a12 + a32)(b21 − a31 + b23 + b32)
		Ua[idx(1,2)][19] = 1;  Ua[idx(3,2)][19] = 1;
		Va[idx(3,1)][19] = -1;
		Vb[idx(2,1)][19] = 1;  Vb[idx(2,3)][19] = 1;  Vb[idx(3,2)][19] = 1;
		// γ21 = (a23 + a33)(a31 + b23 − b31 + b32)
		Ua[idx(2,3)][20] = 1;  Ua[idx(3,3)][20] = 1;
		Va[idx(3,1)][20] = 1;
		Vb[idx(2,3)][20] = 1;  Vb[idx(3,1)][20] = -1; Vb[idx(3,2)][20] = 1;
		// γ22 = (a23 + a33 − a12 − a32)(b23 + b32)
		Ua[idx(2,3)][21] = 1;  Ua[idx(3,3)][21] = 1;  Ua[idx(1,2)][21] = -1; Ua[idx(3,2)][21] = -1;
		Vb[idx(2,3)][21] = 1;  Vb[idx(3,2)][21] = 1;

		// W from the supplied output specs.
		for (int i = 1; i <= 3; i++) {
			for (int j = 1; j <= 3; j++) {
				int row = idx(i, j);
				for (WSpec.Term t : outputs[i - 1][j - 1].terms) {
					W[row][t.gammaIdx - 1] += t.sign;
				}
			}
		}

		return new NonBilinearAlgorithm(n, m, p, Ua, Ub, Va, Vb, W);
	}

	/** Build with Islam 2009's exact output formulas (residual ≈ 4.9). */
	public static NonBilinearAlgorithm buildDefault() {
		return build(defaultOutputs());
	}

	/** Islam 2009 §3.3.1 transcription. */
	public static WSpec[][] defaultOutputs() {
		return new WSpec[][] {
			{
				WSpec.of(+5, +10, +11, +12),                                  // c11
				WSpec.of(+8, +10, -14, +17, -18, +19, -22),                   // c12
				WSpec.of(+1, -11, -12, -16, +17, -18, +19, -22),              // c13
			},
			{
				WSpec.of(+6, -10, +11, +13),                                  // c21
				WSpec.of(+2, -10, +13, +14, +15, +17),                        // c22
				WSpec.of(+9, -11, +15, -16, +17),                             // c23
			},
			{
				WSpec.of(+7, -10, -11, +20, -21, +22),                        // c31
				WSpec.of(+3, -10, +14, -17, +18, +20, +22),                   // c32
				WSpec.of(+4, +11, +16, -17, +18, +21),                        // c33
			},
		};
	}

	/** A signed term in a c_{i,j} output combination. */
	public record WSpec(java.util.List<Term> terms) {
		public record Term(int sign, int gammaIdx) {}
		public static WSpec of(int... signedGammas) {
			java.util.List<Term> ts = new java.util.ArrayList<>();
			for (int g : signedGammas) ts.add(new Term(g > 0 ? 1 : -1, Math.abs(g)));
			return new WSpec(ts);
		}
		/** Return a new WSpec with the sign of γ_{flipIdx} negated. */
		public WSpec flipSign(int flipIdx) {
			java.util.List<Term> nts = new java.util.ArrayList<>();
			boolean found = false;
			for (Term t : terms) {
				if (t.gammaIdx == flipIdx) {
					nts.add(new Term(-t.sign, t.gammaIdx));
					found = true;
				} else nts.add(t);
			}
			if (!found) return null;
			return new WSpec(nts);
		}
	}
}
