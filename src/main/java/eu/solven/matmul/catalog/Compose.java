package eu.solven.matmul.catalog;

import java.util.List;

import eu.solven.matmul.BilinearAlgorithm;
import eu.solven.matmul.ComplexNonCubicBilinearAlgorithm;
import eu.solven.matmul.FactorMatrix;
import eu.solven.matmul.NonCubicBilinearAlgorithm;
import eu.solven.matmul.SparseFactorMatrix;
import eu.solven.matmul.papers.strassen1969.Strassen7;

/**
 * Constructive recursive composition of bilinear matmul algorithms via the
 * Kronecker (tensor) product. Combines algorithms for cubic formats
 * {@code ⟨n_1, n_1, n_1⟩} and {@code ⟨n_2, n_2, n_2⟩} into one for
 * {@code ⟨n_1·n_2, n_1·n_2, n_1·n_2⟩} at rank {@code r_1·r_2}.
 *
 * <p>Index convention: the outer algorithm sees the matrix as an
 * {@code n_2 × n_2} grid of blocks (each block an {@code n_1 × n_1} matrix),
 * and the inner algorithm handles each block multiplication. Global flatten
 * stays row-major: {@code a = (I_block · n_1 + i_block) · n + (J_block · n_1 + j_block)}
 * where {@code n = n_1 · n_2}.</p>
 *
 * <p><b>Correctness</b>: the Kronecker construction is mathematically exact
 * by the bilinear algebra theorem
 * {@code R(⟨n,m,p⟩ ⊗ ⟨n',m',p'⟩) ≤ R(⟨n,m,p⟩) · R(⟨n',m',p'⟩)} (in fact equal
 * when both decompositions are realised). Tests verify exactness at
 * {@code ⟨4,4,4⟩} and {@code ⟨8,8,8⟩}; larger compositions
 * ({@code ⟨16,16,16⟩}, {@code ⟨32,32,32⟩}) are correct by construction —
 * the verifier's `O(n⁶ · r)` cost makes it impractical to re-verify there.</p>
 *
 * <p><b>Scope</b>: cubic-to-cubic only for now. Non-cubic composition (e.g.
 * mixing {@code ⟨2,2,3⟩} with {@code ⟨3,3,2⟩} to get {@code ⟨6,6,6⟩}) would
 * need a non-cubic {@code BilinearAlgorithm} class first — see TODO at the
 * bottom of {@code SMALL_MATMUL_CATALOG.md}.</p>
 */
public final class Compose {

	private Compose() {}

	/**
	 * Cubic Kronecker composition. Convenience wrapper over {@link #kroneckerGeneral}.
	 */
	public static BilinearAlgorithm kronecker(BilinearAlgorithm outer, BilinearAlgorithm inner) {
		NonCubicBilinearAlgorithm composed = kroneckerGeneral(
				NonCubicBilinearAlgorithm.fromCubic(outer),
				NonCubicBilinearAlgorithm.fromCubic(inner));
		return composed.asCubic();
	}

	/**
	 * General Kronecker composition. Combines
	 * {@code outer : ⟨n_2, m_2, p_2⟩} (handling the top-level block structure)
	 * with {@code inner : ⟨n_1, m_1, p_1⟩} (handling each block's multiplication)
	 * to produce {@code ⟨n_1·n_2, m_1·m_2, p_1·p_2⟩} at rank {@code r_outer · r_inner}.
	 *
	 * <p>Indexing: the global row-major flatten of any factor matrix decomposes
	 * uniquely as {@code (I_block · dim_inner + i_block) · (cols_outer · cols_inner)
	 * + (J_block · cols_inner + j_block)}, where the four sub-indices range over
	 * the outer/inner row and column extents of the appropriate slot.</p>
	 */
	public static NonCubicBilinearAlgorithm kroneckerGeneral(
			NonCubicBilinearAlgorithm outer, NonCubicBilinearAlgorithm inner) {
		int n = outer.n * inner.n;
		int m = outer.m * inner.m;
		int p = outer.p * inner.p;
		int rOuter = outer.r, rInner = inner.r;

		// Build each factor sparse, directly from the operands' sparse columns —
		// no dense double[n*m][rTotal] intermediate (~413 MB at ⟨32,32,32⟩). The
		// "cols" argument is the factor's inner-dim: U is (n,m), V is (m,p), W is (n,p).
		FactorMatrix U = kronFactor(outer.u(), inner.u(), outer.m, inner.m, rOuter, rInner);
		FactorMatrix V = kronFactor(outer.v(), inner.v(), outer.p, inner.p, rOuter, rInner);
		FactorMatrix W = kronFactor(outer.w(), inner.w(), outer.p, inner.p, rOuter, rInner);

		return NonCubicBilinearAlgorithm.fromFactors(n, m, p, U, V, W);
	}

	/** Snapshot of a factor's non-zeros grouped by column (product), each column
	 *  in ascending row order. */
	private record ColSnapshot(int[][] rows, double[][] vals) {}

	private static ColSnapshot snapshot(FactorMatrix f) {
		int cols = f.cols();
		int[] cnt = new int[cols];
		f.forEachNonZero((row, col, val) -> cnt[col]++);
		int[][] rows = new int[cols][];
		double[][] vals = new double[cols][];
		for (int c = 0; c < cols; c++) {
			rows[c] = new int[cnt[c]];
			vals[c] = new double[cnt[c]];
		}
		int[] w = new int[cols];
		// forEachNonZero is column-major for the sparse backing, so rows land sorted.
		f.forEachNonZero((row, col, val) -> {
			int i = w[col]++;
			rows[col][i] = row;
			vals[col][i] = val;
		});
		return new ColSnapshot(rows, vals);
	}

	/**
	 * One factor of {@code outer ⊗ inner}, built sparse. {@code colsOuter} /
	 * {@code colsInner} are the factor's inner-dim on each operand (e.g. {@code m}
	 * for U); the result row {@code a = i·cols + j} with {@code i = Ib·rowsInner+ib},
	 * {@code j = Jb·colsInner+jb} is the standard Kronecker row mapping.
	 */
	private static FactorMatrix kronFactor(FactorMatrix oF, FactorMatrix iF,
			int colsOuter, int colsInner, int rOuter, int rInner) {
		ColSnapshot O = snapshot(oF);
		ColSnapshot I = snapshot(iF);
		int rowsInner = iF.rows() / colsInner;
		int cols = colsOuter * colsInner;          // result factor inner-dim
		int resultRows = oF.rows() * iF.rows();
		int[][] colRows = new int[rOuter * rInner][];
		double[][] colVals = new double[rOuter * rInner][];
		for (int ko = 0; ko < rOuter; ko++) {
			int[] oR = O.rows()[ko];
			double[] oV = O.vals()[ko];
			for (int ki = 0; ki < rInner; ki++) {
				int[] iR = I.rows()[ki];
				double[] iV = I.vals()[ki];
				int len = oR.length * iR.length;
				int[] cr = new int[len];
				double[] cv = new double[len];
				int w = 0;
				for (int oi = 0; oi < oR.length; oi++) {
					int aOuter = oR[oi];
					int Ib = aOuter / colsOuter, Jb = aOuter % colsOuter;
					double uo = oV[oi];
					for (int ii = 0; ii < iR.length; ii++) {
						int aInner = iR[ii];
						int ib = aInner / colsInner, jb = aInner % colsInner;
						cr[w] = (Ib * rowsInner + ib) * cols + (Jb * colsInner + jb);
						cv[w] = uo * iV[ii];
						w++;
					}
				}
				colRows[ko * rInner + ki] = cr;
				colVals[ko * rInner + ki] = cv;
			}
		}
		return SparseFactorMatrix.fromColumns(resultRows, colRows, colVals);
	}

	/**
	 * Projection (Perminov draft Def 2.8 / meta-flip-graph {@code Project};
	 * = padding+DCE run in reverse). Restrict a scheme for {@code ⟨n,m,p⟩} to the
	 * kept indices on each axis, yielding an exact scheme for
	 * {@code ⟨|keepN|,|keepM|,|keepP|⟩}, and dead-code-eliminate any product whose
	 * restricted U, V, OR W is all-zero.
	 *
	 * <p>Sound because the matmul tensor restricted to index subsets
	 * {@code Kn×Km, Km×Kp, Kn×Kp} is exactly the {@code ⟨|Kn|,|Km|,|Kp|⟩} matmul
	 * tensor; a product contributing nothing to the restricted tensor (any factor
	 * zero) is removed. Which indices are dropped — and on which axis (incl. the
	 * contracted axis {@code m}) — changes how many products DCE, hence the rank;
	 * the caller enumerates choices (see {@link ProjectionSearch}). The result is
	 * still re-checked by {@link eu.solven.matmul.Verifier#isExactNonCubic}.</p>
	 *
	 * @param keepN sorted kept row indices in {@code [0,n)} (A-rows / C-rows)
	 * @param keepM sorted kept inner indices in {@code [0,m)} (A-cols / B-rows)
	 * @param keepP sorted kept col indices in {@code [0,p)} (B-cols / C-cols)
	 */
	public static NonCubicBilinearAlgorithm project(
			NonCubicBilinearAlgorithm a, int[] keepN, int[] keepM, int[] keepP) {
		int n2 = keepN.length, m2 = keepM.length, p2 = keepP.length;
		double[][] srcU = a.denseU(), srcV = a.denseV(), srcW = a.denseW();
		java.util.List<double[]> uc = new java.util.ArrayList<>();
		java.util.List<double[]> vc = new java.util.ArrayList<>();
		java.util.List<double[]> wc = new java.util.ArrayList<>();
		for (int l = 0; l < a.r; l++) {
			double[] u = new double[n2 * m2];
			boolean unz = false;
			for (int i = 0; i < n2; i++)
				for (int j = 0; j < m2; j++) {
					double v = srcU[keepN[i] * a.m + keepM[j]][l];
					u[i * m2 + j] = v;
					if (v != 0) unz = true;
				}
			if (!unz) continue;
			double[] vv = new double[m2 * p2];
			boolean vnz = false;
			for (int j = 0; j < m2; j++)
				for (int k = 0; k < p2; k++) {
					double v = srcV[keepM[j] * a.p + keepP[k]][l];
					vv[j * p2 + k] = v;
					if (v != 0) vnz = true;
				}
			if (!vnz) continue;
			double[] ww = new double[n2 * p2];
			boolean wnz = false;
			for (int i = 0; i < n2; i++)
				for (int k = 0; k < p2; k++) {
					double v = srcW[keepN[i] * a.p + keepP[k]][l];
					ww[i * p2 + k] = v;
					if (v != 0) wnz = true;
				}
			if (!wnz) continue;
			uc.add(u);
			vc.add(vv);
			wc.add(ww);
		}
		int r = uc.size();
		double[][] U = new double[n2 * m2][r], V = new double[m2 * p2][r], W = new double[n2 * p2][r];
		for (int c = 0; c < r; c++) {
			for (int row = 0; row < n2 * m2; row++) U[row][c] = uc.get(c)[row];
			for (int row = 0; row < m2 * p2; row++) V[row][c] = vc.get(c)[row];
			for (int row = 0; row < n2 * p2; row++) W[row][c] = wc.get(c)[row];
		}
		return new NonCubicBilinearAlgorithm(n2, m2, p2, U, V, W);
	}

	/** Convenience: keep-set = all indices except {@code drop}. */
	public static int[] keepExcept(int dim, int... drop) {
		java.util.Set<Integer> d = new java.util.HashSet<>();
		for (int x : drop) d.add(x);
		int[] keep = new int[dim - d.size()];
		int c = 0;
		for (int i = 0; i < dim; i++) if (!d.contains(i)) keep[c++] = i;
		return keep;
	}

	/**
	 * Block-column concatenation: combine algorithms for {@code ⟨n,m,p1⟩}
	 * and {@code ⟨n,m,p2⟩} into one for {@code ⟨n,m,p1+p2⟩} by stacking
	 * the right operands. The B inputs become {@code [B1|B2]} (m×(p1+p2)),
	 * the C outputs become {@code [C1|C2]} (n×(p1+p2)). Rank is
	 * {@code r1+r2} — strictly additive. Useful when no single-shot
	 * scheme is known for the combined shape but parts are catalogued.
	 *
	 * <p>FMM ⟨2,10,16⟩=248 is decomposed by FMM itself as
	 * {@code ⟨2,10,6⟩+⟨2,10,10⟩}; this method is the bilinear formal
	 * version of that block split.</p>
	 *
	 * @throws IllegalArgumentException if {@code left.n != right.n} or
	 *     {@code left.m != right.m}
	 */
	public static NonCubicBilinearAlgorithm concatRight(
			NonCubicBilinearAlgorithm left, NonCubicBilinearAlgorithm right) {
		if (left.n != right.n || left.m != right.m) {
			throw new IllegalArgumentException(
					"concatRight: ⟨" + left.n + "," + left.m + "," + left.p
							+ "⟩ and ⟨" + right.n + "," + right.m + "," + right.p
							+ "⟩ — n and m must match");
		}
		int n = left.n, m = left.m, p = left.p + right.p;
		int r1 = left.r, r2 = right.r, r = r1 + r2;

		double[][] leftU = left.denseU(), leftV = left.denseV(), leftW = left.denseW();
		double[][] rightU = right.denseU(), rightV = right.denseV(), rightW = right.denseW();

		double[][] U = new double[n * m][r];
		double[][] V = new double[m * p][r];
		double[][] W = new double[n * p][r];

		// U is over A (n×m) — same shape for both halves. Left products
		// (k ∈ [0, r1)) use left.U[ab][k]; right (k ∈ [r1, r)) use right.U[ab][k-r1].
		for (int ab = 0; ab < n * m; ab++) {
			System.arraycopy(leftU[ab], 0, U[ab], 0, r1);
			System.arraycopy(rightU[ab], 0, U[ab], r1, r2);
		}

		// V is over B (m×p). Left half acts on B's first p1 columns; right
		// half on the last p2 columns. Flat index for B[j, l]: j*p + l.
		// Left product k contributes left.V[j*left.p + l'][k] to V[j*p + l'][k]
		// for l' ∈ [0, left.p); right product k' contributes
		// right.V[j*right.p + l''][k'] to V[j*p + (left.p + l'')][r1 + k'].
		int p1 = left.p, p2 = right.p;
		for (int j = 0; j < m; j++) {
			for (int l = 0; l < p1; l++) {
				int srcRow = j * p1 + l;
				int dstRow = j * p + l;
				System.arraycopy(leftV[srcRow], 0, V[dstRow], 0, r1);
			}
			for (int l = 0; l < p2; l++) {
				int srcRow = j * p2 + l;
				int dstRow = j * p + (p1 + l);
				System.arraycopy(rightV[srcRow], 0, V[dstRow], r1, r2);
			}
		}

		// W is over C (n×p). Same column-stacking as V.
		for (int i = 0; i < n; i++) {
			for (int l = 0; l < p1; l++) {
				int srcRow = i * p1 + l;
				int dstRow = i * p + l;
				System.arraycopy(leftW[srcRow], 0, W[dstRow], 0, r1);
			}
			for (int l = 0; l < p2; l++) {
				int srcRow = i * p2 + l;
				int dstRow = i * p + (p1 + l);
				System.arraycopy(rightW[srcRow], 0, W[dstRow], r1, r2);
			}
		}

		return new NonCubicBilinearAlgorithm(n, m, p, U, V, W);
	}

	/**
	 * Block-row concatenation: combine algorithms for {@code ⟨n1,m,p⟩}
	 * and {@code ⟨n2,m,p⟩} into one for {@code ⟨n1+n2,m,p⟩}. Same
	 * additive-rank idea, splitting A's rows instead of C's columns.
	 *
	 * @throws IllegalArgumentException if {@code top.m != bottom.m} or
	 *     {@code top.p != bottom.p}
	 */
	public static NonCubicBilinearAlgorithm concatBelow(
			NonCubicBilinearAlgorithm top, NonCubicBilinearAlgorithm bottom) {
		if (top.m != bottom.m || top.p != bottom.p) {
			throw new IllegalArgumentException(
					"concatBelow: m and p must match");
		}
		int m = top.m, p = top.p, n = top.n + bottom.n;
		int r1 = top.r, r2 = bottom.r, r = r1 + r2;
		int n1 = top.n, n2 = bottom.n;

		double[][] topU = top.denseU(), topV = top.denseV(), topW = top.denseW();
		double[][] bottomU = bottom.denseU(), bottomV = bottom.denseV(), bottomW = bottom.denseW();

		double[][] U = new double[n * m][r];
		double[][] V = new double[m * p][r];
		double[][] W = new double[n * p][r];

		// U: top half acts on A's first n1 rows; bottom on last n2.
		for (int i = 0; i < n1; i++) {
			for (int j = 0; j < m; j++) {
				int srcRow = i * m + j;
				int dstRow = i * m + j;
				System.arraycopy(topU[srcRow], 0, U[dstRow], 0, r1);
			}
		}
		for (int i = 0; i < n2; i++) {
			for (int j = 0; j < m; j++) {
				int srcRow = i * m + j;
				int dstRow = (n1 + i) * m + j;
				System.arraycopy(bottomU[srcRow], 0, U[dstRow], r1, r2);
			}
		}

		// V: top and bottom both use full B (m×p) for their respective halves.
		for (int ab = 0; ab < m * p; ab++) {
			System.arraycopy(topV[ab], 0, V[ab], 0, r1);
			System.arraycopy(bottomV[ab], 0, V[ab], r1, r2);
		}

		// W: top contributes to rows [0, n1); bottom to rows [n1, n).
		for (int i = 0; i < n1; i++) {
			for (int l = 0; l < p; l++) {
				int srcRow = i * p + l;
				int dstRow = i * p + l;
				System.arraycopy(topW[srcRow], 0, W[dstRow], 0, r1);
			}
		}
		for (int i = 0; i < n2; i++) {
			for (int l = 0; l < p; l++) {
				int srcRow = i * p + l;
				int dstRow = (n1 + i) * p + l;
				System.arraycopy(bottomW[srcRow], 0, W[dstRow], r1, r2);
			}
		}

		return new NonCubicBilinearAlgorithm(n, m, p, U, V, W);
	}

	/**
	 * Inner-dimension (contraction) sum: combine algorithms for
	 * {@code ⟨n,m1,p⟩} and {@code ⟨n,m2,p⟩} into one for
	 * {@code ⟨n,m1+m2,p⟩}. Unlike {@link #concatRight} / {@link #concatBelow}
	 * this does <em>not</em> tile the output into disjoint regions — the
	 * shared inner (contraction) dimension is split ({@code A = [A1|A2]},
	 * {@code B = [B1;B2]}) and the two sub-products are <strong>accumulated</strong>:
	 * {@code C = A1·B1 + A2·B2}. Both halves write the full {@code n×p}
	 * output; the bilinear reconstruction {@code C = Σ_k W[:,k]·prod_k}
	 * sums them automatically. Rank is {@code r1+r2} — strictly additive —
	 * by the direct-sum decomposition of the matmul tensor along its
	 * <em>middle</em> mode ({@code ⟨n,m1+m2,p⟩ = ⟨n,m1,p⟩ ⊕_m ⟨n,m2,p⟩}).
	 *
	 * <p>This is the third sibling of {@code concatRight} (p-axis tile) and
	 * {@code concatBelow} (n-axis tile): the m-axis (contraction) variant.
	 * Note neither operand is reused (each half owns its own A-columns and
	 * B-rows) and C is summed rather than concatenated — so it is a
	 * <em>sum</em>, not a geometric concat. (The older
	 * {@code ConcatSplitSearch} docstring claimed the m-axis "needs a real
	 * matmul algorithm"; that conflated tiling with additive composition —
	 * the m-split is a valid rank-{@code r1+r2} construction, just summed.)</p>
	 *
	 * @throws IllegalArgumentException if {@code left.n != right.n} or
	 *     {@code left.p != right.p}
	 */
	public static NonCubicBilinearAlgorithm concatInner(
			NonCubicBilinearAlgorithm left, NonCubicBilinearAlgorithm right) {
		if (left.n != right.n || left.p != right.p) {
			throw new IllegalArgumentException(
					"concatInner: ⟨" + left.n + "," + left.m + "," + left.p
							+ "⟩ and ⟨" + right.n + "," + right.m + "," + right.p
							+ "⟩ — n and p must match");
		}
		int n = left.n, p = left.p, m = left.m + right.m;
		int m1 = left.m, m2 = right.m;
		int r1 = left.r, r2 = right.r, r = r1 + r2;

		double[][] leftU = left.denseU(), leftV = left.denseV(), leftW = left.denseW();
		double[][] rightU = right.denseU(), rightV = right.denseV(), rightW = right.denseW();

		double[][] U = new double[n * m][r];
		double[][] V = new double[m * p][r];
		double[][] W = new double[n * p][r];

		// U over A (n×m): A's columns split. Left products own columns
		// [0, m1); right products own [m1, m). Flat index for A[i, j]: i*m + j.
		for (int i = 0; i < n; i++) {
			for (int j = 0; j < m1; j++) {
				System.arraycopy(leftU[i * m1 + j], 0, U[i * m + j], 0, r1);
			}
			for (int j = 0; j < m2; j++) {
				System.arraycopy(rightU[i * m2 + j], 0, U[i * m + (m1 + j)], r1, r2);
			}
		}

		// V over B (m×p): B's rows split the same way. Flat index for B[j, l]: j*p + l.
		for (int j = 0; j < m1; j++) {
			for (int l = 0; l < p; l++) {
				System.arraycopy(leftV[j * p + l], 0, V[j * p + l], 0, r1);
			}
		}
		for (int j = 0; j < m2; j++) {
			for (int l = 0; l < p; l++) {
				System.arraycopy(rightV[j * p + l], 0, V[(m1 + j) * p + l], r1, r2);
			}
		}

		// W over C (n×p): SHARED — both halves write every output entry; the
		// reconstruction's Σ_k naturally accumulates A1·B1 + A2·B2.
		for (int c = 0; c < n * p; c++) {
			System.arraycopy(leftW[c], 0, W[c], 0, r1);
			System.arraycopy(rightW[c], 0, W[c], r1, r2);
		}

		return new NonCubicBilinearAlgorithm(n, m, p, U, V, W);
	}

	/**
	 * Left-fold composition over a chain. {@code chain.get(0)} is the outermost
	 * level; the result has matrix dimension {@code ∏ n_i} and rank {@code ∏ r_i}.
	 */
	public static BilinearAlgorithm chain(List<BilinearAlgorithm> chain) {
		if (chain.isEmpty()) {
			throw new IllegalArgumentException("chain must be non-empty");
		}
		BilinearAlgorithm acc = chain.get(0);
		for (int i = 1; i < chain.size(); i++) {
			acc = kronecker(acc, chain.get(i));
		}
		return acc;
	}

	/**
	 * Complex Kronecker composition. Mirrors {@link #kroneckerGeneral} for
	 * {@link ComplexNonCubicBilinearAlgorithm}: each composed factor entry is
	 * the <i>complex</i> product
	 * {@code (a + bi)(c + di) = (ac − bd) + (ad + bc)i} of the outer and
	 * inner entries.
	 *
	 * <p>Use case: {@code AlphaEvolve ⟨4,4,4⟩=48 ⊗ AlphaEvolve ⟨4,4,4⟩=48}
	 * yields {@code ⟨16,16,16⟩=2304} over {@code C} — matching fmm-lille's
	 * best-known.</p>
	 */
	public static ComplexNonCubicBilinearAlgorithm kroneckerComplex(
			ComplexNonCubicBilinearAlgorithm outer,
			ComplexNonCubicBilinearAlgorithm inner) {
		int n = outer.n * inner.n;
		int m = outer.m * inner.m;
		int p = outer.p * inner.p;
		int rOuter = outer.r, rInner = inner.r, rTotal = rOuter * rInner;

		double[][] uRe = new double[n * m][rTotal], uIm = new double[n * m][rTotal];
		double[][] vRe = new double[m * p][rTotal], vIm = new double[m * p][rTotal];
		double[][] wRe = new double[n * p][rTotal], wIm = new double[n * p][rTotal];

		fillKroneckerComplex(uRe, uIm, outer.uRe, outer.uIm, inner.uRe, inner.uIm,
				outer.n, outer.m, inner.n, inner.m, rOuter, rInner);
		fillKroneckerComplex(vRe, vIm, outer.vRe, outer.vIm, inner.vRe, inner.vIm,
				outer.m, outer.p, inner.m, inner.p, rOuter, rInner);
		fillKroneckerComplex(wRe, wIm, outer.wRe, outer.wIm, inner.wRe, inner.wIm,
				outer.n, outer.p, inner.n, inner.p, rOuter, rInner);

		return new ComplexNonCubicBilinearAlgorithm(n, m, p, uRe, uIm, vRe, vIm, wRe, wIm);
	}

	private static void fillKroneckerComplex(double[][] dstRe, double[][] dstIm,
			double[][] outerRe, double[][] outerIm,
			double[][] innerRe, double[][] innerIm,
			int rowsOuter, int colsOuter, int rowsInner, int colsInner,
			int rOuter, int rInner) {
		int cols = colsOuter * colsInner;
		for (int Ib = 0; Ib < rowsOuter; Ib++) {
			for (int Jb = 0; Jb < colsOuter; Jb++) {
				int aOuter = Ib * colsOuter + Jb;
				for (int ib = 0; ib < rowsInner; ib++) {
					for (int jb = 0; jb < colsInner; jb++) {
						int aInner = ib * colsInner + jb;
						int i = Ib * rowsInner + ib;
						int j = Jb * colsInner + jb;
						int a = i * cols + j;
						for (int ko = 0; ko < rOuter; ko++) {
							double or = outerRe[aOuter][ko], oi = outerIm[aOuter][ko];
							if (or == 0.0 && oi == 0.0) continue;
							for (int ki = 0; ki < rInner; ki++) {
								double ir = innerRe[aInner][ki], ii = innerIm[aInner][ki];
								dstRe[a][ko * rInner + ki] = or * ir - oi * ii;
								dstIm[a][ko * rInner + ki] = or * ii + oi * ir;
							}
						}
					}
				}
			}
		}
	}

	/**
	 * Block-split constructor for cubic {@code ⟨n,n,n⟩} with axis split
	 * {@code n = u + v}. Builds the composed algorithm constructively from
	 * 8 sub-products (the 2×2×2 block tensor with mixed shapes
	 * {@code ⟨u/v, u/v, u/v⟩}) drawn from our catalog via {@code lookup}.
	 *
	 * <p>Resulting rank: {@code R(⟨u,u,u⟩) + 3·R(⟨u,u,v⟩) + 3·R(⟨u,v,v⟩) +
	 * R(⟨v,v,v⟩)} — the **naïve** block decomposition. Sedoglavic 2017's
	 * algebraic identity drops the trailing {@code ⟨v,v,v⟩} term via
	 * cross-block sharing (separate construction; see ROADMAP).</p>
	 *
	 * <p>For {@code ⟨7,7,7⟩} with {@code u=4, v=3}: yields
	 * {@code 49 + 3·29 + 3·38 + 23 = 273} multiplications (vs Sedoglavic's
	 * 250 with the algebraic trick).</p>
	 *
	 * @throws IllegalStateException if any required sub-format is missing
	 *                               from {@code lookup}.
	 */
	public static NonCubicBilinearAlgorithm blockSplitCubic(int n, int u, int v,
			eu.solven.matmul.catalog.Recombination.AlgorithmLookup lookup) {
		if (u + v != n || u < 1 || v < 1) {
			throw new IllegalArgumentException("require u + v = n, both ≥ 1; got u=" + u + ", v=" + v + ", n=" + n);
		}
		int[] axisSize = { u, v };
		int[] axisOffset = { 0, u };

		// Look up all 8 sub-products; track total rank.
		NonCubicBilinearAlgorithm[] subs = new NonCubicBilinearAlgorithm[8];
		int totalRank = 0;
		int subIdx = 0;
		for (int ai = 0; ai < 2; ai++) {
			for (int mi = 0; mi < 2; mi++) {
				for (int bi = 0; bi < 2; bi++) {
					int aa = axisSize[ai], mid = axisSize[mi], bb = axisSize[bi];
					NonCubicBilinearAlgorithm sub = lookup.find(aa, mid, bb)
							.orElseThrow(() -> new IllegalStateException(
									"blockSplitCubic: missing ⟨" + aa + "," + mid + "," + bb + "⟩"));
					subs[subIdx++] = sub;
					totalRank += sub.r;
				}
			}
		}

		double[][] U = new double[n * n][totalRank];
		double[][] V = new double[n * n][totalRank];
		double[][] W = new double[n * n][totalRank];

		int kOff = 0;
		subIdx = 0;
		for (int ai = 0; ai < 2; ai++) {
			int aaOff = axisOffset[ai];
			int aa = axisSize[ai];
			for (int mi = 0; mi < 2; mi++) {
				int midOff = axisOffset[mi];
				int mid = axisSize[mi];
				for (int bi = 0; bi < 2; bi++) {
					int bbOff = axisOffset[bi];
					int bb = axisSize[bi];
					NonCubicBilinearAlgorithm sub = subs[subIdx];
					int subR = sub.r;
					double[][] subU = sub.denseU(), subV = sub.denseV(), subW = sub.denseW();

					// Embed sub.U (shape (aa, mid)) into the global (aaOff..aaOff+aa) × (midOff..midOff+mid) sub-block.
					for (int i = 0; i < aa; i++) {
						for (int j = 0; j < mid; j++) {
							int aSubFlat = i * mid + j;
							int aGlobalFlat = (aaOff + i) * n + (midOff + j);
							double[] dstRow = U[aGlobalFlat];
							double[] srcRow = subU[aSubFlat];
							for (int k = 0; k < subR; k++) dstRow[kOff + k] = srcRow[k];
						}
					}
					// Embed sub.V (shape (mid, bb))
					for (int i = 0; i < mid; i++) {
						for (int j = 0; j < bb; j++) {
							int bSubFlat = i * bb + j;
							int bGlobalFlat = (midOff + i) * n + (bbOff + j);
							double[] dstRow = V[bGlobalFlat];
							double[] srcRow = subV[bSubFlat];
							for (int k = 0; k < subR; k++) dstRow[kOff + k] = srcRow[k];
						}
					}
					// Embed sub.W (shape (aa, bb))
					for (int i = 0; i < aa; i++) {
						for (int j = 0; j < bb; j++) {
							int cSubFlat = i * bb + j;
							int cGlobalFlat = (aaOff + i) * n + (bbOff + j);
							double[] dstRow = W[cGlobalFlat];
							double[] srcRow = subW[cSubFlat];
							for (int k = 0; k < subR; k++) dstRow[kOff + k] = srcRow[k];
						}
					}

					kOff += subR;
					subIdx++;
				}
			}
		}

		return new NonCubicBilinearAlgorithm(n, n, n, U, V, W);
	}

	/** Non-cubic-aware chain: left-fold across general {@code ⟨n,m,p⟩} factors. */
	public static NonCubicBilinearAlgorithm chainGeneral(List<NonCubicBilinearAlgorithm> chain) {
		if (chain.isEmpty()) {
			throw new IllegalArgumentException("chain must be non-empty");
		}
		NonCubicBilinearAlgorithm acc = chain.get(0);
		for (int i = 1; i < chain.size(); i++) {
			acc = kroneckerGeneral(acc, chain.get(i));
		}
		return acc;
	}

	/**
	 * Build {@code Strassen^k} — `k` levels of Strassen ⟨2,2,2⟩ → algorithm for
	 * {@code ⟨2^k, 2^k, 2^k⟩} at rank {@code 7^k}. Examples:
	 * {@code k=1} → ⟨2,2,2⟩ at 7;
	 * {@code k=5} → ⟨32,32,32⟩ at 16,807.
	 */
	public static BilinearAlgorithm strassenPower(int k) {
		if (k < 1) throw new IllegalArgumentException("k must be ≥ 1, got " + k);
		BilinearAlgorithm acc = Strassen7.get();
		for (int i = 1; i < k; i++) {
			acc = kronecker(acc, Strassen7.get());
		}
		return acc;
	}
}
