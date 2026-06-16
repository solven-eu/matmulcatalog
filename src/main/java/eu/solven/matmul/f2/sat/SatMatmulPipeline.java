package eu.solven.matmul.f2.sat;

import java.io.IOException;
import java.util.List;
import java.util.Optional;

import eu.solven.matmul.BilinearAlgorithm;
import eu.solven.matmul.Verifier;

/**
 * End-to-end pipeline: build the Z/2 CNF encoding, hand to SAT4J, decode any
 * satisfying assignment back into a {@link BilinearAlgorithm}.
 *
 * Provides factory helpers for the common target tensors (dense and triangular
 * matmul over Z/2). For arbitrary tensors callers can supply their own.
 */
public final class SatMatmulPipeline {

	private SatMatmulPipeline() {}

	/**
	 * Search for a rank-{@code r} decomposition of {@code target} over Z/2.
	 * Uses the default encoder (column lex-ordering on) and no BreakID
	 * preprocessing.
	 *
	 * @return the decomposition (entries 0.0 or 1.0) if SAT, empty if UNSAT.
	 */
	public static Optional<BilinearAlgorithm> findZ2Decomposition(int n, int r, int[][][] target) {
		return findZ2Decomposition(n, r, target, false);
	}

	/**
	 * Search for a rank-{@code r} decomposition with optional BreakID
	 * preprocessing.
	 *
	 * <p><b>Symmetry-breaking modes</b> (selected by {@code useBreakId} and
	 * BreakID availability):</p>
	 * <ul>
	 *   <li>{@code useBreakId == false}: use the encoder's hand-coded column
	 *       lex-ordering (always available).</li>
	 *   <li>{@code useBreakId == true} and BreakID is on PATH: <b>disable</b>
	 *       the hand-coded lex-ordering and let BreakID generate the symmetry-
	 *       breaking clauses automatically (it should find a superset of what
	 *       we hand-coded, including sign flips and the transpose action).</li>
	 *   <li>{@code useBreakId == true} but BreakID NOT on PATH: fall back to
	 *       hand-coded lex-ordering — without some form of symmetry breaking,
	 *       SAT4J cannot solve hard ranks.</li>
	 * </ul>
	 */
	public static Optional<BilinearAlgorithm> findZ2Decomposition(
			int n, int r, int[][][] target, boolean useBreakId) {
		boolean breakIdAvailable = useBreakId && BreakIdBridge.isAvailable();

		// If we'll use BreakID, skip the hand-coded lex (BreakID does its own).
		// Otherwise, hand-coded lex is essential — without it SAT4J hangs.
		Z2CnfEncoder encoder = new Z2CnfEncoder(n, r, target, !breakIdAvailable);
		int varCount = encoder.getVarCount();
		List<int[]> clauses = encoder.getClauses();

		if (breakIdAvailable) {
			try {
				Cnf.ReadResult augmented = BreakIdBridge.preprocess(varCount, clauses);
				varCount = augmented.varCount;
				clauses = augmented.clauses;
			} catch (IOException e) {
				// BreakID failed mid-run; the CNF has no hand-coded lex so
				// SAT4J would hang. Re-encode with hand-coded lex as fallback.
				encoder = new Z2CnfEncoder(n, r, target, true);
				varCount = encoder.getVarCount();
				clauses = encoder.getClauses();
			}
		}

		Optional<boolean[]> model = solveWithBestAvailable(varCount, clauses);
		return model.map(encoder::decode);
	}

	/** Prefer kissat over SAT4J when available — kissat is orders of magnitude faster on large CNFs. */
	private static Optional<boolean[]> solveWithBestAvailable(int varCount, List<int[]> clauses) {
		if (KissatSolver.isAvailable()) {
			return KissatSolver.solve(varCount, clauses);
		}
		return Sat4jSolver.solve(varCount, clauses);
	}

	/**
	 * Dense {@code n×n} matmul tensor over Z/2: same shape as
	 * {@link Verifier#intMatmulTensor(int)} but treated as a Z/2 target.
	 */
	public static int[][][] z2DenseMatmulTensor(int n) {
		return Verifier.intMatmulTensor(n);
	}

	/**
	 * Restricted matmul tensor over Z/2: the dense tensor with every entry
	 * involving an "off-restriction" position zeroed out. The
	 * {@code allowedPositions} mask (length n²) selects which flatten indices
	 * are non-zero in A, B, and C. Used for sub-problems like diagonal-only,
	 * triangular, or "diagonal + a single off-diagonal entry."
	 *
	 * <p>Returns both the restricted target tensor and a corresponding
	 * "force-zero" mask suitable for {@link Z2CnfEncoder}'s 5-arg constructor:
	 * positions <b>not</b> in {@code allowedPositions} should have their
	 * U/V/W rows forced to 0 in the encoder.</p>
	 */
	public static int[][][] z2RestrictedMatmulTensor(int n, boolean[] allowedPositions) {
		int n2 = n * n;
		int[][][] T = z2DenseMatmulTensor(n);
		for (int a = 0; a < n2; a++) {
			for (int b = 0; b < n2; b++) {
				for (int c = 0; c < n2; c++) {
					if (!allowedPositions[a] || !allowedPositions[b] || !allowedPositions[c]) {
						T[a][b][c] = 0;
					}
				}
			}
		}
		return T;
	}

	/** Helper: positions on the diagonal plus an arbitrary single extra. */
	public static boolean[] diagonalPlusOne(int n, int extraRow, int extraCol) {
		boolean[] mask = new boolean[n * n];
		for (int i = 0; i < n; i++) mask[i * n + i] = true;
		mask[extraRow * n + extraCol] = true;
		return mask;
	}

	/** Search the restricted-positions case. */
	public static java.util.Optional<BilinearAlgorithm> findZ2RestrictedDecomposition(
			int n, int r, int[][][] target, boolean[] allowedPositions) {
		boolean[] forceZero = new boolean[allowedPositions.length];
		for (int i = 0; i < forceZero.length; i++) forceZero[i] = !allowedPositions[i];
		Z2CnfEncoder encoder = new Z2CnfEncoder(n, r, target, true, forceZero);
		Optional<boolean[]> model = solveWithBestAvailable(encoder.getVarCount(), encoder.getClauses());
		return model.map(encoder::decode);
	}

	/**
	 * Upper-triangular matmul tensor over Z/2: the dense tensor with every
	 * entry involving a strictly-lower-triangular flatten index zeroed out.
	 * The non-zero entries correspond exactly to the contributions
	 * {@code A[i, l] · B[l, j]} → {@code C[i, j]} where each of (i,l), (l,j),
	 * (i,j) is on or above the diagonal.
	 *
	 * For {@code n=2}: lower-triangular flatten indices are just {@code (1, 0) = 2}.
	 * For {@code n=3}: indices {@code (1,0)=3, (2,0)=6, (2,1)=7}.
	 */
	public static int[][][] z2UpperTriangularMatmulTensor(int n) {
		int n2 = n * n;
		int[][][] T = z2DenseMatmulTensor(n);
		boolean[] isLowerTri = new boolean[n2];
		for (int i = 0; i < n; i++) {
			for (int j = 0; j < i; j++) {
				isLowerTri[i * n + j] = true;
			}
		}
		for (int a = 0; a < n2; a++) {
			for (int b = 0; b < n2; b++) {
				for (int c = 0; c < n2; c++) {
					if (isLowerTri[a] || isLowerTri[b] || isLowerTri[c]) {
						T[a][b][c] = 0;
					}
				}
			}
		}
		return T;
	}

	// ───────────────────────────────────────────────────────────────────────────
	// Non-cubic ⟨n, m, p⟩ helpers
	// ───────────────────────────────────────────────────────────────────────────

	/**
	 * Dense non-cubic matmul tensor over Z/2 for {@code ⟨n, m, p⟩}: same shape as
	 * {@link Verifier#intMatmulTensor(int, int, int)} but treated as a Z/2 target.
	 */
	public static int[][][] z2DenseMatmulTensor(int n, int m, int p) {
		return Verifier.intMatmulTensor(n, m, p);
	}

	/**
	 * Search for a rank-{@code r} decomposition of a non-cubic {@code ⟨n, m, p⟩}
	 * matmul tensor over Z/2. Returns the raw {@code {U, V, W}} factors
	 * (shapes {@code [n·m][r], [m·p][r], [n·p][r]}) if SAT, empty if UNSAT.
	 */
	public static Optional<double[][][]> findZ2NonCubicDecomposition(
			int n, int m, int p, int r, int[][][] target) {
		int dimU = n * m, dimV = m * p, dimW = n * p;
		Z2CnfEncoder encoder = new Z2CnfEncoder(dimU, dimV, dimW, r, target, true, null);
		Optional<boolean[]> model = solveWithBestAvailable(encoder.getVarCount(), encoder.getClauses());
		return model.map(encoder::decodeRaw);
	}

	/**
	 * Verify a non-cubic decomposition against a Z/2 target: returns true iff
	 * {@code ⊕_k U[a, k] · V[b, k] · W[c, k] = target[a, b, c]} (mod 2) for all
	 * {@code (a, b, c) ∈ [dimU] × [dimV] × [dimW]}.
	 */
	public static boolean verifyZ2NonCubic(double[][][] uvw, int[][][] target) {
		double[][] U = uvw[0], V = uvw[1], W = uvw[2];
		int dimU = U.length, dimV = V.length, dimW = W.length;
		int r = U[0].length;
		for (int a = 0; a < dimU; a++) {
			for (int b = 0; b < dimV; b++) {
				for (int c = 0; c < dimW; c++) {
					int sum = 0;
					for (int k = 0; k < r; k++) {
						int u = (int) Math.round(U[a][k]);
						int v = (int) Math.round(V[b][k]);
						int w = (int) Math.round(W[c][k]);
						sum ^= (u & v & w);
					}
					if (sum != target[a][b][c]) return false;
				}
			}
		}
		return true;
	}

	/**
	 * Verify a decomposition against a Z/2 target: returns true iff
	 * {@code ⊕_k U[a, k] · V[b, k] · W[c, k] = target[a, b, c]} (mod 2) for all (a, b, c).
	 */
	public static boolean verifyZ2(BilinearAlgorithm alg, int[][][] target) {
		int n2 = alg.n * alg.n;
		for (int a = 0; a < n2; a++) {
			for (int b = 0; b < n2; b++) {
				for (int c = 0; c < n2; c++) {
					int sum = 0;
					for (int k = 0; k < alg.r; k++) {
						int u = (int) Math.round(alg.U[a][k]);
						int v = (int) Math.round(alg.V[b][k]);
						int w = (int) Math.round(alg.W[c][k]);
						sum ^= (u & v & w);
					}
					if (sum != target[a][b][c]) return false;
				}
			}
		}
		return true;
	}
}
