package eu.solven.matmul.papers.dis2009;

import java.util.ArrayList;
import java.util.List;
import java.util.function.BiFunction;
import java.util.function.IntFunction;

import eu.solven.matmul.NonCubicBilinearAlgorithm;

/**
 * Concrete factor-matrix construction for {@link PanTrilinearAggregation#build}.
 *
 * <p>Direct port of Islam's Magma {@code TA.mgm} {@code ProductEven} /
 * {@code ProductOdd} functions, verified against
 * {@code references/islam2009/sympy/04c_lemma4_from_magma.py} and
 * {@code 05_lemma4_odd_from_magma.py}: the sympy versions check both
 * symbolic correctness ({@code C_TA == A·B}) and multiplication count
 * exactly matches {@link PanTrilinearAggregation#cubicBound(int)} for
 * {@code n ∈ {4, 6}} (even) and {@code nn ∈ {3, 5}} (odd).</p>
 *
 * <p>The algorithm builds the U/V/W factor matrices column by column:
 * each invocation of {@link #mul(LinForm, LinForm, LinForm)} corresponds
 * to one rank-1 atom and appends one column. Linear forms are
 * represented as dense {@code double[]} sparse-vector views; the only
 * non-integer coefficients in the construction come from the padding
 * matrix R (denominators of {@code n+1} where {@code n} is the input
 * block size, ≤ 32 here), well within {@code double} precision.</p>
 */
final class PanTrilinearAggregationBuilder {

	private final int nn;     // input matrix side length
	private final int n;      // block size (nn/2 for even, (nn+1)/2 for odd)
	private final int m;      // padded block size = n + 1; full padded matrix is 2m × 2m
	private final int nInput; // = nn²  (input size for A and B)
	private final int mPadC;  // = m·m  (padded C-variable count)

	private final List<double[]> uCols = new ArrayList<>();
	private final List<double[]> vCols = new ArrayList<>();
	/** Padded c-coefficients (m·m space). Remapped to nn² output space at end. */
	private final List<double[]> wColsPadded = new ArrayList<>();

	/** Mask of which padded c-variables are zeroed (algorithm never touches them). */
	private final boolean[] cZeroed;

	private PanTrilinearAggregationBuilder(int nn, int n) {
		this.nn = nn;
		this.n = n;
		this.m = n + 1;
		this.nInput = nn * nn;
		this.mPadC = (2 * m) * (2 * m);
		this.cZeroed = new boolean[mPadC];
	}

	// ──────────────────────────────────────────────────────────────────
	//  Linear-form algebra (LinForm = double[] sparse-vector view).
	// ──────────────────────────────────────────────────────────────────

	/** Dense linear-form vector of fixed length {@code len}. */
	static record LinForm(double[] coeffs) {
		static LinForm zero(int len) { return new LinForm(new double[len]); }

		/** Returns a NEW LinForm = this + other. */
		LinForm add(LinForm o) {
			double[] r = coeffs.clone();
			for (int i = 0; i < r.length; i++) r[i] += o.coeffs[i];
			return new LinForm(r);
		}

		LinForm sub(LinForm o) {
			double[] r = coeffs.clone();
			for (int i = 0; i < r.length; i++) r[i] -= o.coeffs[i];
			return new LinForm(r);
		}

		LinForm neg() {
			double[] r = new double[coeffs.length];
			for (int i = 0; i < r.length; i++) r[i] = -coeffs[i];
			return new LinForm(r);
		}

		LinForm scale(double s) {
			double[] r = new double[coeffs.length];
			for (int i = 0; i < r.length; i++) r[i] = coeffs[i] * s;
			return new LinForm(r);
		}

		boolean isZero() {
			for (double c : coeffs) if (Math.abs(c) > 1e-12) return false;
			return true;
		}
	}

	/** Matrix of linear forms (each entry is a LinForm of the same length). */
	static final class FormMatrix {
		final int rows, cols;
		final LinForm[][] entries;

		FormMatrix(int rows, int cols, int formLen) {
			this.rows = rows;
			this.cols = cols;
			this.entries = new LinForm[rows][cols];
			for (int i = 0; i < rows; i++)
				for (int j = 0; j < cols; j++)
					entries[i][j] = LinForm.zero(formLen);
		}

		LinForm at(int r, int c) { return entries[r][c]; }

		FormMatrix add(FormMatrix o) {
			FormMatrix out = new FormMatrix(rows, cols, entries[0][0].coeffs.length);
			for (int i = 0; i < rows; i++)
				for (int j = 0; j < cols; j++)
					out.entries[i][j] = entries[i][j].add(o.entries[i][j]);
			return out;
		}

		FormMatrix sub(FormMatrix o) {
			FormMatrix out = new FormMatrix(rows, cols, entries[0][0].coeffs.length);
			for (int i = 0; i < rows; i++)
				for (int j = 0; j < cols; j++)
					out.entries[i][j] = entries[i][j].sub(o.entries[i][j]);
			return out;
		}

		FormMatrix neg() {
			FormMatrix out = new FormMatrix(rows, cols, entries[0][0].coeffs.length);
			for (int i = 0; i < rows; i++)
				for (int j = 0; j < cols; j++)
					out.entries[i][j] = entries[i][j].neg();
			return out;
		}

		/** Returns L · this where L is a numeric {@code (rows_L × this.rows)} matrix. */
		FormMatrix leftMul(double[][] L) {
			int rl = L.length;
			FormMatrix out = new FormMatrix(rl, cols, entries[0][0].coeffs.length);
			for (int i = 0; i < rl; i++)
				for (int j = 0; j < cols; j++)
					for (int k = 0; k < rows; k++) {
						double l = L[i][k];
						if (l == 0) continue;
						out.entries[i][j] = out.entries[i][j].add(entries[k][j].scale(l));
					}
			return out;
		}

		/** Returns this · R where R is a numeric {@code (this.cols × cols_R)} matrix. */
		FormMatrix rightMul(double[][] R) {
			int cR = R[0].length;
			FormMatrix out = new FormMatrix(rows, cR, entries[0][0].coeffs.length);
			for (int i = 0; i < rows; i++)
				for (int j = 0; j < cR; j++)
					for (int k = 0; k < cols; k++) {
						double r = R[k][j];
						if (r == 0) continue;
						out.entries[i][j] = out.entries[i][j].add(entries[i][k].scale(r));
					}
			return out;
		}

		/** Build from a per-entry linear-form supplier. */
		static FormMatrix build(int rows, int cols, int formLen,
				BiFunction<Integer, Integer, LinForm> supplier) {
			FormMatrix out = new FormMatrix(rows, cols, formLen);
			for (int i = 0; i < rows; i++)
				for (int j = 0; j < cols; j++)
					out.entries[i][j] = supplier.apply(i, j);
			return out;
		}

		/** Swap rows {@code r1} and {@code r2} (1-based). */
		FormMatrix swapRows(int r1, int r2) {
			FormMatrix out = new FormMatrix(rows, cols, entries[0][0].coeffs.length);
			for (int i = 0; i < rows; i++)
				for (int j = 0; j < cols; j++)
					out.entries[i][j] = entries[i][j];
			LinForm[] tmp = out.entries[r1 - 1];
			out.entries[r1 - 1] = out.entries[r2 - 1];
			out.entries[r2 - 1] = tmp;
			return out;
		}

		FormMatrix swapCols(int c1, int c2) {
			FormMatrix out = new FormMatrix(rows, cols, entries[0][0].coeffs.length);
			for (int i = 0; i < rows; i++)
				for (int j = 0; j < cols; j++)
					out.entries[i][j] = entries[i][j];
			for (int i = 0; i < rows; i++) {
				LinForm tmp = out.entries[i][c1 - 1];
				out.entries[i][c1 - 1] = out.entries[i][c2 - 1];
				out.entries[i][c2 - 1] = tmp;
			}
			return out;
		}
	}

	// ──────────────────────────────────────────────────────────────────
	//  Atom emission (= one column of U, V, W).
	// ──────────────────────────────────────────────────────────────────

	/**
	 * Emit one rank-1 atom: {@code (a) × (b) × (c)} where a, b, c are
	 * linear forms in A entries, B entries, and padded-C entries
	 * respectively. Zero atoms are dropped (matches the Magma {@code Mul}
	 * counter behaviour).
	 */
	private void mul(LinForm a, LinForm b, LinForm c) {
		if (a.isZero() || b.isZero() || c.isZero()) return;
		uCols.add(a.coeffs());
		vCols.add(b.coeffs());
		wColsPadded.add(c.coeffs());
	}

	// ──────────────────────────────────────────────────────────────────
	//  Padding matrices L, R  (from Magma TA.mgm lines 46-58).
	// ──────────────────────────────────────────────────────────────────

	private static double[][] buildL(int sz) {
		// L = [I_sz ; -1^T_sz] of size (sz+1) × sz.
		double[][] L = new double[sz + 1][sz];
		for (int i = 0; i < sz; i++) L[i][i] = 1;
		for (int j = 0; j < sz; j++) L[sz][j] = -1;
		return L;
	}

	private static double[][] buildR(int sz) {
		// R = [I - (1/(sz+1)) u u^T  | -(1/(sz+1)) u]  size sz × (sz+1).
		double inv = 1.0 / (sz + 1);
		double[][] R = new double[sz][sz + 1];
		for (int i = 0; i < sz; i++)
			for (int j = 0; j < sz; j++)
				R[i][j] = (i == j ? 1 : 0) - inv;
		for (int i = 0; i < sz; i++) R[i][sz] = -inv;
		return R;
	}

	// ──────────────────────────────────────────────────────────────────
	//  Set up A^{i,j}, B^{i,j} as FormMatrix in the nn²-input space.
	// ──────────────────────────────────────────────────────────────────

	/**
	 * Build the (rows, cols) sub-block of input A (or B) as a FormMatrix
	 * whose {@code (r,c)} entry is the unit basis vector for input
	 * position {@code (rowOffset+r, colOffset+c)} of the nn×nn matrix.
	 */
	private FormMatrix inputBlock(int rows, int cols, int rowOffset, int colOffset) {
		return FormMatrix.build(rows, cols, nInput, (i, j) -> {
			double[] coeffs = new double[nInput];
			int row = rowOffset + i;
			int col = colOffset + j;
			coeffs[row * nn + col] = 1.0;
			return new LinForm(coeffs);
		});
	}

	/** Pad a c-block into (n+1) × (n+1), inserting at (1,1) with trailing zeros. */
	private FormMatrix paddedCBlock(int blockRow, int blockCol, int origRows, int origCols) {
		// The padded m × m C-variable space has blocks at row offsets {0, m} and col offsets {0, m}.
		// "blockRow" / "blockCol" are 0 or 1 (which of the four C-blocks).
		int rOff = blockRow * m;
		int cOff = blockCol * m;
		FormMatrix M = new FormMatrix(m, m, mPadC);
		for (int i = 0; i < origRows; i++)
			for (int j = 0; j < origCols; j++) {
				int idx = (rOff + i) * (2 * m) + (cOff + j);
				double[] coeffs = new double[mPadC];
				coeffs[idx] = 1.0;
				M.entries[i][j] = new LinForm(coeffs);
			}
		return M;
	}

	private void markCZeroed(int paddedIdx) {
		cZeroed[paddedIdx] = true;
	}

	// ──────────────────────────────────────────────────────────────────
	//  Operators s_0, s_1, s_2, u_1, u_2, u_2', u_3, u_4
	//  (faithful ports of TA.mgm lines 116-192).
	// ──────────────────────────────────────────────────────────────────

	private void s0(FormMatrix A, FormMatrix B, FormMatrix C) {
		for (int i = 1; i < m; i++) {   // 1 ≤ i < m  (0-indexed: 0 ≤ i-1 < m-1)
			mul(A.at(i - 1, i - 1), B.at(i - 1, i - 1), C.at(i - 1, i - 1));
		}
	}

	private void s1(FormMatrix A, FormMatrix B, FormMatrix C) {
		for (int i = 1; i <= m; i++) for (int j = 1; j <= m; j++) for (int k = 1; k <= m; k++) {
			if (!propS1Improved(i, j, k)) continue;
			LinForm f1 = A.at(i - 1, j - 1).add(A.at(j - 1, k - 1)).add(A.at(k - 1, i - 1));
			LinForm f2 = B.at(j - 1, k - 1).add(B.at(k - 1, i - 1)).add(B.at(i - 1, j - 1));
			LinForm f3 = C.at(k - 1, i - 1).add(C.at(i - 1, j - 1)).add(C.at(j - 1, k - 1));
			mul(f1, f2, f3);
		}
	}

	private void s2(FormMatrix A, FormMatrix B, FormMatrix C,
			FormMatrix U, FormMatrix V, FormMatrix W,
			FormMatrix X, FormMatrix Y, FormMatrix Z) {
		for (int i = 1; i <= m; i++) for (int j = 1; j <= m; j++) for (int k = 1; k <= m; k++) {
			if (!propS2Improved(i, j, k)) continue;
			LinForm f1 = A.at(i - 1, j - 1).add(U.at(j - 1, k - 1)).add(X.at(k - 1, i - 1));
			LinForm f2 = B.at(j - 1, k - 1).add(V.at(k - 1, i - 1)).add(Y.at(i - 1, j - 1));
			LinForm f3 = C.at(k - 1, i - 1).add(W.at(i - 1, j - 1)).add(Z.at(j - 1, k - 1));
			mul(f1, f2, f3);
		}
	}

	private void u1(FormMatrix A, FormMatrix B, FormMatrix C, FormMatrix W, FormMatrix Z) {
		for (int j = 1; j < m; j++) for (int i = 1; i < m; i++) {
			if (i == j) continue;
			LinForm outerA = A.at(i - 1, j - 1);
			LinForm outerB = B.at(i - 1, j - 1);
			LinForm inner = W.at(i - 1, j - 1).scale(m);
			for (int k = 1; k <= m; k++) {
				inner = inner.add(C.at(k - 1, i - 1)).add(Z.at(j - 1, k - 1));
			}
			// u-family operators are SUBTRACTED from the polynomial in the
			// sympy reference (`product -= u_op(...)`), so we negate the
			// inner c-form to flip the sign of this atom's W contribution.
			mul(outerA, outerB, inner.neg());
		}
	}

	private void u2(FormMatrix A, FormMatrix B, FormMatrix D, FormMatrix W, FormMatrix Z) {
		for (int i = 1; i < m; i++) {
			LinForm outerA = A.at(i - 1, i - 1);
			LinForm outerB = B.at(i - 1, i - 1);
			LinForm inner = LinForm.zero(mPadC);
			for (int k = 1; k <= m; k++) inner = inner.add(D.at(k - 1, i - 1));
			inner = inner.add(W.at(i - 1, i - 1).scale(m));
			for (int k = 1; k <= m; k++) inner = inner.add(Z.at(i - 1, k - 1));
			// u-family operators are SUBTRACTED from the polynomial in the
			// sympy reference (`product -= u_op(...)`), so we negate the
			// inner c-form to flip the sign of this atom's W contribution.
			mul(outerA, outerB, inner.neg());
		}
	}

	private void u2Prime(FormMatrix A, FormMatrix B, FormMatrix C) {
		for (int i = 1; i < m; i++) {
			LinForm outerA = A.at(i - 1, i - 1);
			LinForm outerB = B.at(i - 1, i - 1);
			LinForm inner = C.at(i - 1, i - 1).scale(m - 9);
			for (int k = 1; k <= m; k++) {
				inner = inner.add(C.at(k - 1, i - 1)).add(C.at(i - 1, k - 1));
			}
			// u-family operators are SUBTRACTED from the polynomial in the
			// sympy reference (`product -= u_op(...)`), so we negate the
			// inner c-form to flip the sign of this atom's W contribution.
			mul(outerA, outerB, inner.neg());
		}
	}

	private void u3(FormMatrix A, FormMatrix Y, FormMatrix C) {
		for (int i = 1; i < m; i++) {
			LinForm outerA = A.at(i - 1, m - 1);
			LinForm outerB = Y.at(i - 1, m - 1);
			LinForm inner = LinForm.zero(mPadC);
			for (int k = 1; k <= m; k++) inner = inner.add(C.at(k - 1, i - 1));
			// u-family operators are SUBTRACTED from the polynomial in the
			// sympy reference (`product -= u_op(...)`), so we negate the
			// inner c-form to flip the sign of this atom's W contribution.
			mul(outerA, outerB, inner.neg());
		}
	}

	private void u4(FormMatrix A, FormMatrix Y, FormMatrix Z) {
		for (int j = 1; j < m; j++) {
			LinForm outerA = A.at(m - 1, j - 1);
			LinForm outerB = Y.at(m - 1, j - 1);
			LinForm inner = LinForm.zero(mPadC);
			for (int k = 1; k <= m; k++) inner = inner.add(Z.at(j - 1, k - 1));
			// u-family operators are SUBTRACTED from the polynomial in the
			// sympy reference (`product -= u_op(...)`), so we negate the
			// inner c-form to flip the sign of this atom's W contribution.
			mul(outerA, outerB, inner.neg());
		}
	}

	private boolean propS1(int i, int j, int k) {
		boolean b1 = (0 <= i) && (i <= j) && (j < k);
		boolean b2 = (0 <= k) && (k < j) && (j <= i);
		return b1 || b2;
	}

	private boolean propSizeM(int i, int j, int k) {
		int eqM = 0;
		if (i == m) eqM++;
		if (j == m) eqM++;
		if (k == m) eqM++;
		return eqM < 2;
	}

	private boolean propS1Improved(int i, int j, int k) {
		return propS1(i, j, k) && propSizeM(i, j, k);
	}

	private boolean propS2Improved(int i, int j, int k) {
		return propSizeM(i, j, k) && !(i == j && j == k);
	}

	// ──────────────────────────────────────────────────────────────────
	//  Output assembly: convert padded-W to nn²-output-W.
	// ──────────────────────────────────────────────────────────────────

	/**
	 * For each output position {@code (i, l) ∈ [0,nn)²}, return the
	 * (padded-c row, padded-c col) that holds its coefficient per the
	 * Magma's {@code extr + Transpose + HorizontalJoin + VerticalJoin
	 * + Submatrix} assembly chain.
	 */
	private int[] outputToPaddedCIndex() {
		// outputMap[outIdx] = paddedIdx where outIdx = i·nn + l.
		int[] outputMap = new int[nn * nn];
		int padFull = 2 * m;
		for (int i = 0; i < nn; i++) {
			for (int l = 0; l < nn; l++) {
				int blockRow = (i < n) ? 0 : 1;
				int blockCol = (l < n) ? 0 : 1;
				int localRow = (i < n) ? i : (i - n);
				int localCol = (l < n) ? l : (l - n);
				// extr[blockRow][blockCol] = out_C[blockRow*m : ...].T, so
				// after transpose, extr[i_local, l_local] = out_C[block, l_local, block, i_local].
				// extr_blockR_blockC[localR, localC] = out_C[blockR*m + localC, blockC*m + localR]
				// (the .T in sympy swaps the two indices). Output blockRow → padded ROW
				// offset; output blockCol → padded COL offset; locals swap inside.
				int paddedRow = blockRow * m + localCol;
				int paddedCol = blockCol * m + localRow;
				outputMap[i * nn + l] = paddedRow * padFull + paddedCol;
			}
		}
		return outputMap;
	}

	private double[] remapPaddedToOutput(double[] wPadded, int[] outputMap) {
		double[] w = new double[nn * nn];
		for (int o = 0; o < w.length; o++) {
			int padIdx = outputMap[o];
			w[o] = wPadded[padIdx];
		}
		return w;
	}

	// ──────────────────────────────────────────────────────────────────
	//  Even-case driver (port of Magma TA.mgm ProductEven, lines 200-296).
	// ──────────────────────────────────────────────────────────────────

	static NonCubicBilinearAlgorithm buildEven(int nn) {
		PanTrilinearAggregationBuilder b = new PanTrilinearAggregationBuilder(nn, nn / 2);
		b.runEven();
		return b.finish();
	}

	static NonCubicBilinearAlgorithm buildOdd(int nn) {
		PanTrilinearAggregationBuilder b = new PanTrilinearAggregationBuilder(nn, (nn + 1) / 2);
		b.runOdd();
		return b.finish();
	}

	private void runEven() {
		// Zero out the padded c-variables the Magma drops (lines 211-216).
		int padFull = 2 * m;
		// Note: Magma uses 1-based indices; n+1 in Magma = n in Java; m in Magma = padFull - 1 in Java.
		int colNplus1 = n;             // Magma `n+1` (1-based) → 0-based index n
		int colM     = padFull - 1;    // Magma `m`   (1-based) → 0-based padFull - 1
		int rowNplus1 = n;
		int rowM     = padFull - 1;
		for (int i = 0; i < nn; i++) {
			markCZeroed(i * padFull + colNplus1);
			markCZeroed(i * padFull + colM);
			markCZeroed(rowNplus1 * padFull + i);
			markCZeroed(rowM * padFull + i);
		}

		// Four (nn/2 × nn/2) sub-blocks of A and B as FormMatrix.
		FormMatrix a11 = inputBlock(n, n, 0,        0);
		FormMatrix a12 = inputBlock(n, n, 0,        n);
		FormMatrix a21 = inputBlock(n, n, n,        0);
		FormMatrix a22 = inputBlock(n, n, n,        n);
		// B uses the same nInput space but separately — we have to model A and B as
		// the SAME input space (we'll separate them at U/V emission). Reuse helper.
		// For Pan TA, A and B input spaces are DIFFERENT (V has B's coefficients).
		// We'll use two separate FormMatrix trees: a* for A-input, b* for B-input.
		// → done below by ALSO calling inputBlock for the B side; the LinForm
		// length is the same (nInput = nn²), and uCols / vCols are stored separately.

		FormMatrix b11 = inputBlock(n, n, 0,        0);
		FormMatrix b12 = inputBlock(n, n, 0,        n);
		FormMatrix b21 = inputBlock(n, n, n,        0);
		FormMatrix b22 = inputBlock(n, n, n,        n);

		// C blocks (padded, all 4 are (n+1)² placeholders into the m·m padded space).
		FormMatrix C11 = paddedCBlock(0, 0, m, m);
		FormMatrix C12 = paddedCBlock(0, 1, m, m);
		FormMatrix C21 = paddedCBlock(1, 0, m, m);
		FormMatrix C22 = paddedCBlock(1, 1, m, m);

		// Zero-sum padding: enlarge each block via L · X · R or L · X · L^T.
		double[][] Ln = buildL(n);
		double[][] Rn = buildR(n);
		double[][] LnT = transpose(Ln);
		// IMPORTANT: A uses L · A · R; B uses L · B · L^T.
		FormMatrix A11 = a11.leftMul(Ln).rightMul(Rn);
		FormMatrix A12 = a12.leftMul(Ln).rightMul(Rn);
		FormMatrix A21 = a21.leftMul(Ln).rightMul(Rn);
		FormMatrix A22 = a22.leftMul(Ln).rightMul(Rn);
		FormMatrix B11 = b11.leftMul(Ln).rightMul(LnT);
		FormMatrix B12 = b12.leftMul(Ln).rightMul(LnT);
		FormMatrix B21 = b21.leftMul(Ln).rightMul(LnT);
		FormMatrix B22 = b22.leftMul(Ln).rightMul(LnT);

		// 4 × s_0
		s0(A12.sub(A11).add(A21), B21.add(B12).add(B11), C11.sub(C12).add(C21));
		s0(A12.add(A21).sub(A22), B22.add(B12).add(B21), C12.add(C22).sub(C21));
		// 2 × s_1
		s1(A11, B11, C11);
		s1(A22, B22, C22);
		// 2 × s_2
		s2(A12, B21, C11,  A11.neg(), B12, C12.neg(),  A21, B11, C21);
		s2(A12, B22, C12,  A21, B12, C22,  A22.neg(), B21, C21.neg());
		// 8 × u_1
		u1(A11,          B11,  C11,  C11,         C11);
		u1(A11.neg(),    B21,  C12.neg(), C21,    C11);
		u1(A12,          B21,  C12,  C22,         C21.neg());
		u1(A22.neg(),    B12,  C21.neg(), C12,    C22);
		u1(A12,          B11,  C11,  C12.neg(),   C21);
		u1(A21,          B12,  C21,  C11,         C12.neg());
		u1(A21,          B22,  C22,  C21.neg(),   C12);
		u1(A22,          B22,  C22,  C22,         C22);
		// 2 × u_2' + 6 × u_2
		u2Prime(A11, B11, C11);
		u2(A11, B21, C12, C21.neg(), C11.neg());
		u2(A12, B11, C11, C12.neg(), C21);
		u2(A21, B12, C21, C11, C12.neg());
		u2(A22, B12, C21, C12.neg(), C22.neg());
		u2(A21, B22, C22, C21.neg(), C12);
		u2(A12, B21, C12, C22, C21.neg());
		u2Prime(A22, B22, C22);
		// 4 × u_3
		u3(A11.add(A12), B11, C11);
		u3(A11.add(A12), B21, C12);
		u3(A21.add(A22), B12, C21);
		u3(A21.add(A22), B22, C22);
		// 4 × u_4
		u4(A11, B11.sub(B21), C11);
		u4(A12, B11.sub(B21), C21);
		u4(A21, B22.sub(B12), C12);
		u4(A22, B22.sub(B12), C22);
	}

	// ──────────────────────────────────────────────────────────────────
	//  Odd-case driver (port of Magma TA.mgm ProductOdd, lines 304-440).
	// ──────────────────────────────────────────────────────────────────

	private void runOdd() {
		int padFull = 2 * m;
		// Odd case uses rectangular original blocks: a11 (n×n), a12 (n×(nn-n)), a21 ((nn-n)×n), a22 ((nn-n)×(nn-n)).
		int n2 = nn - n;
		FormMatrix a11_raw = inputBlock(n,  n,  0, 0);
		FormMatrix a12_raw = inputBlock(n,  n2, 0, n);
		FormMatrix a21_raw = inputBlock(n2, n,  n, 0);
		FormMatrix a22_raw = inputBlock(n2, n2, n, n);
		FormMatrix b11_raw = inputBlock(n,  n,  0, 0);
		FormMatrix b12_raw = inputBlock(n,  n2, 0, n);
		FormMatrix b21_raw = inputBlock(n2, n,  n, 0);
		FormMatrix b22_raw = inputBlock(n2, n2, n, n);

		// Enlarge each. EnlargeLeft for A uses L_rows · M · R_cols, where rows/cols
		// are the matrix's own dims (NOT necessarily n).
		FormMatrix a11_en = enlargeLeftForA(a11_raw);
		FormMatrix a12_en = enlargeLeftForA(a12_raw);
		FormMatrix a21_en = enlargeLeftForA(a21_raw);
		FormMatrix a22_en = enlargeLeftForA(a22_raw);
		FormMatrix b11_en = enlargeRightForB(b11_raw);
		FormMatrix b12_en = enlargeRightForB(b12_raw);
		FormMatrix b21_en = enlargeRightForB(b21_raw);
		FormMatrix b22_en = enlargeRightForB(b22_raw);

		// Pad each enlarged block to (m × m) = (n+1)² with zeros, then swap to put
		// the zero rows/cols at the second-last position (per Magma odd-case logic).
		FormMatrix A11 = padToSquare(a11_en, m);
		FormMatrix A12 = padToSquare(a12_en, m).swapCols(n, n + 1);
		FormMatrix A21 = padToSquare(a21_en, m).swapRows(n, n + 1);
		FormMatrix A22 = padToSquare(a22_en, m).swapRows(n, n + 1).swapCols(n, n + 1);
		FormMatrix B11 = padToSquare(b11_en, m);
		FormMatrix B12 = padToSquare(b12_en, m).swapCols(n, n + 1);
		FormMatrix B21 = padToSquare(b21_en, m).swapRows(n, n + 1);
		FormMatrix B22 = padToSquare(b22_en, m).swapRows(n, n + 1).swapCols(n, n + 1);

		// C blocks: rectangular originals padded to (m × m) with zeros.
		FormMatrix C11 = paddedCBlock(0, 0, n + 1, n + 1);
		FormMatrix C12 = paddedCBlock(0, 1, n + 1, n2 + 1);
		FormMatrix C21 = paddedCBlock(1, 0, n2 + 1, n + 1);
		FormMatrix C22 = paddedCBlock(1, 1, n2 + 1, n2 + 1);

		// Same operator sequence as even, with the s_2 cyclic shift the Magma uses.
		s0(A12.sub(A11).add(A21), B21.add(B12).add(B11), C11.sub(C12).add(C21));
		s0(A12.add(A21).sub(A22), B22.add(B12).add(B21), C12.add(C22).sub(C21));
		s1(A11, B11, C11);
		s1(A22, B22, C22);
		// s_2 cyclic shift (Magma odd-case lines 401-402).
		s2(A11.neg(), B12, C12.neg(),  A21, B11, C21,  A12, B21, C11);
		s2(A22.neg(), B21, C21.neg(),  A12, B22, C12,  A21, B12, C22);

		// In Magma odd-case u_2 block precedes u_1 (lines 404-411 vs 413-420).
		u2Prime(A11, B11, C11);
		u2(A11.neg(),    B21,  C12.neg(), C21,    C11);
		u2(A12,          B11,  C11,  C12.neg(),   C21);
		u2(A21,          B12,  C21,  C11,         C12.neg());
		u2(A22.neg(),    B12,  C21.neg(), C12,    C22);
		u2(A21,          B22,  C22,  C21.neg(),   C12);
		u2(A12,          B21,  C12,  C22,         C21.neg());
		u2Prime(A22, B22, C22);

		u1(A11,          B11,  C11,  C11,         C11);
		u1(A11.neg(),    B21,  C12.neg(), C21,    C11);
		u1(A12,          B21,  C12,  C22,         C21.neg());
		u1(A22.neg(),    B12,  C21.neg(), C12,    C22);
		u1(A12,          B11,  C11,  C12.neg(),   C21);
		u1(A21,          B12,  C21,  C11,         C12.neg());
		u1(A21,          B22,  C22,  C21.neg(),   C12);
		u1(A22,          B22,  C22,  C22,         C22);

		u3(A11.add(A12), B11, C11);
		u3(A11.add(A12), B21, C12);
		u3(A21.add(A22), B12, C21);
		u3(A21.add(A22), B22, C22);

		u4(A11, B11.sub(B21), C11);
		u4(A12, B11.sub(B21), C21);
		u4(A21, B22.sub(B12), C12);
		u4(A22, B22.sub(B12), C22);
	}

	private FormMatrix enlargeLeftForA(FormMatrix M) {
		// L_rows · M · R_cols.
		double[][] L = buildL(M.rows);
		double[][] R = buildR(M.cols);
		return M.leftMul(L).rightMul(R);
	}

	private FormMatrix enlargeRightForB(FormMatrix M) {
		// L_rows · M · L_cols^T.
		double[][] L1 = buildL(M.rows);
		double[][] L2T = transpose(buildL(M.cols));
		return M.leftMul(L1).rightMul(L2T);
	}

	private FormMatrix padToSquare(FormMatrix M, int targetSize) {
		int formLen = M.entries[0][0].coeffs.length;
		FormMatrix out = new FormMatrix(targetSize, targetSize, formLen);
		for (int i = 0; i < M.rows && i < targetSize; i++)
			for (int j = 0; j < M.cols && j < targetSize; j++)
				out.entries[i][j] = M.entries[i][j];
		return out;
	}

	private static double[][] transpose(double[][] M) {
		double[][] T = new double[M[0].length][M.length];
		for (int i = 0; i < M.length; i++)
			for (int j = 0; j < M[0].length; j++)
				T[j][i] = M[i][j];
		return T;
	}

	// ──────────────────────────────────────────────────────────────────
	//  Finalisation: post-process W and assemble NonCubicBilinearAlgorithm.
	// ──────────────────────────────────────────────────────────────────

	private NonCubicBilinearAlgorithm finish() {
		int r = uCols.size();
		int[] outputMap = outputToPaddedCIndex();

		double[][] U = new double[nInput][r];
		double[][] V = new double[nInput][r];
		double[][] W = new double[nn * nn][r];

		for (int k = 0; k < r; k++) {
			double[] uCol = uCols.get(k);
			double[] vCol = vCols.get(k);
			double[] wRemapped = remapPaddedToOutput(wColsPadded.get(k), outputMap);
			for (int i = 0; i < nInput; i++) U[i][k] = uCol[i];
			for (int i = 0; i < nInput; i++) V[i][k] = vCol[i];
			for (int i = 0; i < nn * nn; i++) W[i][k] = wRemapped[i];
		}

		return new NonCubicBilinearAlgorithm(nn, nn, nn, U, V, W);
	}
}
