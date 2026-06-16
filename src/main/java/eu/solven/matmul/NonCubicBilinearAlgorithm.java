package eu.solven.matmul;

/**
 * A bilinear algorithm for non-cubic matrix multiplication {@code ⟨n, m, p⟩}:
 * compute {@code C = A·B} where {@code A} is {@code n×m}, {@code B} is
 * {@code m×p}, {@code C} is {@code n×p}.
 *
 * <p>For {@code k = 0..r-1}:</p>
 * <pre>
 *   M_k = (Σ_{i,j} U[i·m + j][k] · A[i][j]) · (Σ_{j,l} V[j·p + l][k] · B[j][l])
 *   C[i][l] = Σ_k W[i·p + l][k] · M_k
 * </pre>
 *
 * <p>Row-major flatten throughout. When {@code n == m == p}, this is
 * equivalent to {@link BilinearAlgorithm}; use {@link #asCubic} / {@link #fromCubic}
 * to convert.</p>
 */
public class NonCubicBilinearAlgorithm {

	public final int n;
	public final int m;
	public final int p;
	public final int r;

	// Factor matrices, stored COLUMN-MAJOR SPARSE (task #7). A factor is ~3–5 %
	// dense, so this is ~15–25× smaller than the historical dense double[][] —
	// the memory fix for caching many large schemes (e.g. PanTA cubes). Consumers
	// read via u()/v()/w() (sparse-efficient: dotColumn / axpyColumn /
	// forEachInColumn) or, where a dense array is genuinely needed, denseU()/…
	// which materialise on demand.
	private final FactorMatrix uMat;
	private final FactorMatrix vMat;
	private final FactorMatrix wMat;

	/** Factor matrix {@code U} (coefficients of A-entries per product). */
	public FactorMatrix u() { return uMat; }

	/** Factor matrix {@code V} (coefficients of B-entries per product). */
	public FactorMatrix v() { return vMat; }

	/** Factor matrix {@code W} (output combination of products per C-entry). */
	public FactorMatrix w() { return wMat; }

	// Dense materialisation hatch (task #6): for builders / sibling-class interop
	// (BilinearAlgorithm, asCubic) that genuinely need a double[][]. Once the
	// backing is sparse (task #7) these allocate on demand — callers should grab
	// one once and reuse, not call per element. Hot paths must use u()/v()/w().
	/** Dense {@code [n·m][r]} view of {@code U} (materialised on demand). */
	public double[][] denseU() { return ((SparseFactorMatrix) uMat).toDense(); }
	/** Dense {@code [m·p][r]} view of {@code V} (materialised on demand). */
	public double[][] denseV() { return ((SparseFactorMatrix) vMat).toDense(); }
	/** Dense {@code [n·p][r]} view of {@code W} (materialised on demand). */
	public double[][] denseW() { return ((SparseFactorMatrix) wMat).toDense(); }

	/** The ⟨n,m,p⟩ shape of this algorithm (orientation-significant). */
	public Shape shape() {
		return Shape.of(n, m, p);
	}

	public NonCubicBilinearAlgorithm(int n, int m, int p, double[][] U, double[][] V, double[][] W) {
		int dimU = n * m, dimV = m * p, dimW = n * p;
		if (U.length != dimU) {
			throw new IllegalArgumentException("U rows must be n·m = " + dimU + ", got " + U.length);
		}
		if (V.length != dimV) {
			throw new IllegalArgumentException("V rows must be m·p = " + dimV + ", got " + V.length);
		}
		if (W.length != dimW) {
			throw new IllegalArgumentException("W rows must be n·p = " + dimW + ", got " + W.length);
		}
		int rank = U[0].length;
		for (int i = 0; i < dimU; i++) {
			if (U[i].length != rank) throw new IllegalArgumentException("U row " + i + " has wrong rank");
		}
		for (int i = 0; i < dimV; i++) {
			if (V[i].length != rank) throw new IllegalArgumentException("V row " + i + " has wrong rank");
		}
		for (int i = 0; i < dimW; i++) {
			if (W[i].length != rank) throw new IllegalArgumentException("W row " + i + " has wrong rank");
		}
		this.n = n;
		this.m = m;
		this.p = p;
		this.r = rank;
		this.uMat = SparseFactorMatrix.fromDense(U);
		this.vMat = SparseFactorMatrix.fromDense(V);
		this.wMat = SparseFactorMatrix.fromDense(W);
	}

	/**
	 * Build directly from (already-sparse) factor matrices, with no dense
	 * round-trip. Used by the orientation operators ({@link #cyclicShift},
	 * {@link #transpose}), which only re-index rows of an existing scheme.
	 */
	private NonCubicBilinearAlgorithm(int n, int m, int p,
			FactorMatrix u, FactorMatrix v, FactorMatrix w) {
		if (u.rows() != n * m || v.rows() != m * p || w.rows() != n * p) {
			throw new IllegalArgumentException("factor row counts do not match ⟨"
					+ n + "," + m + "," + p + "⟩");
		}
		if (v.cols() != u.cols() || w.cols() != u.cols()) {
			throw new IllegalArgumentException("U/V/W must share the rank (column count)");
		}
		this.n = n;
		this.m = m;
		this.p = p;
		this.r = u.cols();
		this.uMat = u;
		this.vMat = v;
		this.wMat = w;
	}

	/**
	 * Build from already-sparse factor matrices, with no dense round-trip — for
	 * composition operators (e.g. the Kronecker product) that produce sparse
	 * factors directly. Factors must be immutable and sized for ⟨n,m,p⟩.
	 */
	public static NonCubicBilinearAlgorithm fromFactors(int n, int m, int p,
			FactorMatrix u, FactorMatrix v, FactorMatrix w) {
		return new NonCubicBilinearAlgorithm(n, m, p, u, v, w);
	}

	/**
	 * The elementary "naive" scheme for ⟨n,m,p⟩: one product per (i,j,l) triple,
	 * rank n·m·p, every factor column a single {@code 1}. Built directly as sparse
	 * factor matrices — NO dense {@code double[][]} round-trip (a width-1 block like
	 * ⟨1,32,32⟩ would otherwise allocate an ~8 MB mostly-zero matrix). This is the
	 * single source of truth for trivial/width-1 sub-blocks, shared by recombination
	 * construct and lineage replay so the two agree bit-for-bit.
	 */
	public static NonCubicBilinearAlgorithm naive(int n, int m, int p) {
		int r = n * m * p;
		int[][] uCols = new int[r][], vCols = new int[r][], wCols = new int[r][];
		double[][] vals = new double[r][];
		double[] one = { 1.0 };
		int k = 0;
		for (int i = 0; i < n; i++)
			for (int j = 0; j < m; j++)
				for (int l = 0; l < p; l++) {
					uCols[k] = new int[] { i * m + j };
					vCols[k] = new int[] { j * p + l };
					wCols[k] = new int[] { i * p + l };
					vals[k] = one;   // immutable, safe to share across columns
					k++;
				}
		return new NonCubicBilinearAlgorithm(n, m, p,
				SparseFactorMatrix.fromColumns(n * m, uCols, vals),
				SparseFactorMatrix.fromColumns(m * p, vCols, vals),
				SparseFactorMatrix.fromColumns(n * p, wCols, vals));
	}

	public boolean isCubic() {
		return n == m && m == p;
	}

	public int dimU() { return n * m; }
	public int dimV() { return m * p; }
	public int dimW() { return n * p; }

	/** Lift a cubic algorithm into the non-cubic representation. */
	public static NonCubicBilinearAlgorithm fromCubic(BilinearAlgorithm alg) {
		return new NonCubicBilinearAlgorithm(alg.n, alg.n, alg.n, alg.U, alg.V, alg.W);
	}

	/**
	 * Project to a cubic algorithm. Throws if {@code !isCubic()}.
	 */
	public BilinearAlgorithm asCubic() {
		if (!isCubic()) {
			throw new IllegalStateException(
					"not cubic: ⟨" + n + "," + m + "," + p + "⟩");
		}
		return new BilinearAlgorithm(n, denseU(), denseV(), denseW());
	}

	/**
	 * Returns the algorithm for the cyclically-shifted format
	 * {@code ⟨m, p, n⟩} derived from this one.
	 *
	 * <p>Matmul tensor identity: {@code T(A, B, C) = T(B, C, A)} under
	 * cyclic re-labelling of the trilinear slots (since
	 * {@code Σ A_{ij} B_{jk} C_{ik}} is invariant under cyclic relabelling
	 * of the index roles). For factor matrices that means
	 * {@code U_new = V}, {@code V_new = W} with an index-transpose, and
	 * {@code W_new = U} with an index-transpose.</p>
	 *
	 * <p>Three applications return to the original format
	 * {@code ⟨n,m,p⟩}. Used by block-decomposition constructors to
	 * obtain mixed-shape sub-algorithms (e.g. derive {@code ⟨4,4,3⟩},
	 * {@code ⟨4,3,4⟩}, {@code ⟨3,4,4⟩} from a single catalog entry).</p>
	 */
	public NonCubicBilinearAlgorithm cyclicShift() {
		int newN = m, newM = p, newP = n;
		SparseFactorMatrix su = (SparseFactorMatrix) uMat;
		SparseFactorMatrix sv = (SparseFactorMatrix) vMat;
		SparseFactorMatrix sw = (SparseFactorMatrix) wMat;

		// U' = V — both over the (m, p) flatten, same indexing; reuse (immutable).
		FactorMatrix newU = sv;
		// V' from W: new row k·n+i, old row i·p+k (i ∈ [n], k ∈ [p]).
		FactorMatrix newV = sw.mapRows(newM * newP, old -> (old % p) * n + (old / p));
		// W' from U: new row j·n+i, old row i·m+j (i ∈ [n], j ∈ [m]).
		FactorMatrix newW = su.mapRows(newN * newP, old -> (old % m) * n + (old / m));

		return new NonCubicBilinearAlgorithm(newN, newM, newP, newU, newV, newW);
	}

	/**
	 * Returns the algorithm for the transposed format {@code ⟨p, m, n⟩}.
	 *
	 * <p>Matmul tensor identity: {@code (AB)^T = B^T A^T}. If {@code A} is
	 * {@code n×m} and {@code B} is {@code m×p}, then {@code C^T = B^T·A^T}
	 * is a {@code ⟨p, m, n⟩} matmul where {@code A_new = B^T} (p×m) and
	 * {@code B_new = A^T} (m×n). The factor matrices map as:</p>
	 *
	 * <ul>
	 *   <li>{@code U_new[k·m + j] = V[j·p + k]} (swap U and V with index transpose)</li>
	 *   <li>{@code V_new[j·n + i] = U[i·m + j]}</li>
	 *   <li>{@code W_new[k·n + i] = W[i·p + k]} (W index transposed)</li>
	 * </ul>
	 *
	 * <p>Combined with {@link #cyclicShift} (the cyclic Z/3 in S₃) this
	 * transpose-involution generates the full S₃ action on tensor slots —
	 * 6 orientations total for cubic bases. Different orientations of the
	 * same algorithm have different U/V/W column-support multisets per
	 * axis, which can produce different sub-shape distributions when used
	 * as outer bases under min-reduction.</p>
	 */
	public NonCubicBilinearAlgorithm transpose() {
		int newN = p, newM = m, newP = n;
		SparseFactorMatrix su = (SparseFactorMatrix) uMat;
		SparseFactorMatrix sv = (SparseFactorMatrix) vMat;
		SparseFactorMatrix sw = (SparseFactorMatrix) wMat;

		// U' from V: new row k·m+j, old row j·p+k (j ∈ [m], k ∈ [p]).
		FactorMatrix newU = sv.mapRows(newN * newM, old -> (old % p) * m + (old / p));
		// V' from U: new row j·n+i, old row i·m+j (i ∈ [n], j ∈ [m]).
		FactorMatrix newV = su.mapRows(newM * newP, old -> (old % m) * n + (old / m));
		// W' from W: new row k·n+i, old row i·p+k (i ∈ [n], k ∈ [p]).
		FactorMatrix newW = sw.mapRows(newN * newP, old -> (old % p) * n + (old / p));

		return new NonCubicBilinearAlgorithm(newN, newM, newP, newU, newV, newW);
	}

	/**
	 * Returns the algorithm for the requested cyclic orientation of this
	 * algorithm's format. If {@code (newN, newM, newP)} is reachable via
	 * 0–2 cyclic shifts, the corresponding algorithm is returned;
	 * otherwise empty (non-cyclic permutations like swap-axis are not yet
	 * implemented).
	 */
	public java.util.Optional<NonCubicBilinearAlgorithm> orientAs(int newN, int newM, int newP) {
		// Try all 6 elements of the S₃ orbit: 3 cyclic shifts of this, plus
		// 3 cyclic shifts of this.transpose(). Together these cover the full
		// symmetric group action on matmul tensor slots (⟨9,11,12⟩ ↔ ⟨9,12,11⟩
		// requires the transpose+cyclic² composition, not just cyclic).
		NonCubicBilinearAlgorithm cur = this;
		for (int i = 0; i < 3; i++) {
			if (cur.n == newN && cur.m == newM && cur.p == newP) {
				return java.util.Optional.of(cur);
			}
			cur = cur.cyclicShift();
		}
		NonCubicBilinearAlgorithm t = this.transpose();
		for (int i = 0; i < 3; i++) {
			if (t.n == newN && t.m == newM && t.p == newP) {
				return java.util.Optional.of(t);
			}
			t = t.cyclicShift();
		}
		return java.util.Optional.empty();
	}
}
