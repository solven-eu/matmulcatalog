package eu.solven.matmul.catalog;

import java.util.ArrayList;
import java.util.List;

import org.apache.commons.math3.fitting.leastsquares.LeastSquaresBuilder;
import org.apache.commons.math3.fitting.leastsquares.LeastSquaresOptimizer;
import org.apache.commons.math3.fitting.leastsquares.LevenbergMarquardtOptimizer;
import org.apache.commons.math3.fitting.leastsquares.MultivariateJacobianFunction;
import org.apache.commons.math3.linear.Array2DRowRealMatrix;
import org.apache.commons.math3.linear.ArrayRealVector;
import org.apache.commons.math3.linear.RealMatrix;
import org.apache.commons.math3.linear.RealVector;
import org.apache.commons.math3.util.Pair;

import eu.solven.matmul.NonCubicBilinearAlgorithm;

/**
 * Continuous-isotropy equivalence for rank-7 {@code ⟨2,2,2⟩} schemes (#168).
 *
 * <p>The matmul structure tensor for our index convention (Û_k is n×m on
 * {@code (A[i][j])}, V̂_k is m×p on {@code (B[j][l])}, Ŵ_k is n×p on
 * {@code (C[i][l])}, satisfying {@code Σ_k Û V̂ Ŵ = δ_{ii'}δ_{j'j''}δ_{ll'}}) is
 * invariant under {@code GL₂(K)³ = GL(n)×GL(m)×GL(p)} acting as:</p>
 * <pre>
 *   Û_k ↦ P · Û_k · Qᵀ
 *   V̂_k ↦ Q⁻ᵀ · V̂_k · Sᵀ
 *   Ŵ_k ↦ P⁻ᵀ · Ŵ_k · S⁻¹
 * </pre>
 * <p>derived by requiring each shared index space to contract to a δ: the
 * n-space couples Û-left ({@code P}) with Ŵ-left ({@code P⁻ᵀ}), the m-space
 * couples Û-right ({@code Qᵀ}) with V̂-left ({@code Q⁻ᵀ}), the p-space couples
 * V̂-right ({@code Sᵀ}) with Ŵ-right ({@code Sᵀ})… {@code Σ_i (P e_i)⊗(P⁻ᵀ e_i)
 * = Σ_i e_i⊗e_i}, etc. On top sits the discrete {@code S₇} relabelling of the
 * seven products (and, for the cubic shape, an {@code S₃} axis permutation —
 * not needed to relate two schemes already at {@code ⟨2,2,2⟩}, which this MVP
 * omits).</p>
 *
 * <p>De Groote 1978 / Burichenko 2014: over {@code ℂ} every rank-7
 * decomposition lies in a SINGLE {@code GL₂³} orbit, so Strassen and
 * Strassen–Winograd are the same scheme up to continuous coordinate change.
 * This class verifies that constructively.</p>
 *
 * <h3>Algorithm</h3>
 * <ol>
 *   <li><b>Invariant pruning.</b> For each product column the triple of 2×2
 *       matrix-ranks {@code (rank Û_k, rank V̂_k, rank Ŵ_k)} is a
 *       {@code GL₂³}-invariant. A valid relabelling {@code π} must map each
 *       column of {@code a} to a column of {@code b} of equal signature — this
 *       collapses the {@code 7! = 5040} permutations to the product of the
 *       per-signature group factorials (a few dozen for Strassen).</li>
 *   <li><b>Levenberg–Marquardt.</b> For each surviving {@code π} we solve the
 *       12-parameter least-squares problem {@code ‖action(a; P,Q,R)∘π − b‖²}
 *       (84 residuals) starting from {@code P=Q=R=I} plus a small
 *       deterministic perturbation. Residual {@code < TOL} ⇒ equivalent.</li>
 * </ol>
 *
 * <p><b>Scope note.</b> Collapsing equivalent entries is for cataloguing /
 * reporting only — the discrete search pool must KEEP all representatives,
 * because the discrete axis-flip / sign-flip enumeration does not span the
 * continuous orbit (see #108, #110, and the task #168 search-pool note).</p>
 */
public final class Isotropy222 {

	private Isotropy222() {}

	/** Residual threshold below which two schemes are declared equivalent. */
	public static final double TOL = 1e-9;

	private static final int R = 7;

	/** A scheme's seven products as 2×2 matrices, {@code m[k][slot][row][col]}
	 *  with {@code slot} 0=U, 1=V, 2=W. */
	public record Scheme222(double[][][][] m) {
		public static Scheme222 of(NonCubicBilinearAlgorithm alg) {
			if (alg.n != 2 || alg.m != 2 || alg.p != 2 || alg.r != 7) {
				throw new IllegalArgumentException(
						"Isotropy222 expects ⟨2,2,2⟩=7, got ⟨" + alg.n + "," + alg.m + "," + alg.p
								+ "⟩=" + alg.r);
			}
			double[][] srcU = alg.denseU();
			double[][] srcV = alg.denseV();
			double[][] srcW = alg.denseW();
			double[][][][] m = new double[R][3][2][2];
			for (int k = 0; k < R; k++) {
				for (int i = 0; i < 2; i++) {
					for (int j = 0; j < 2; j++) {
						m[k][0][i][j] = srcU[i * 2 + j][k];
						m[k][1][i][j] = srcV[i * 2 + j][k];
						m[k][2][i][j] = srcW[i * 2 + j][k];
					}
				}
			}
			return new Scheme222(m);
		}
	}

	/** A change-of-basis witnessing {@code a ≅ b}: {@code P,Q,R ∈ GL₂}, product
	 *  relabelling {@code perm} (so that {@code action(a[perm[k]]) ≈ b[k]}), and
	 *  the achieved least-squares residual. */
	public record Transform(double[][] P, double[][] Q, double[][] Rm, int[] perm, double residual) {}

	/** Try to prove {@code a ≅ b} under {@code GL₂³ ⋊ S₇}; empty if none found. */
	public static java.util.Optional<Transform> findEquivalence(Scheme222 a, Scheme222 b) {
		int[][] sigA = signatures(a);
		int[][] sigB = signatures(b);
		for (int[] perm : candidatePermutations(sigA, sigB)) {
			Transform t = solve(a, b, perm);
			if (t != null && t.residual() < TOL) {
				return java.util.Optional.of(t);
			}
		}
		return java.util.Optional.empty();
	}

	/** Convenience overload on raw algorithms. */
	public static java.util.Optional<Transform> findEquivalence(
			NonCubicBilinearAlgorithm a, NonCubicBilinearAlgorithm b) {
		return findEquivalence(Scheme222.of(a), Scheme222.of(b));
	}

	// ── invariant signatures + permutation pruning ──────────────────────────

	/** Per-column {@code (rank Û_k, rank V̂_k, rank Ŵ_k)} — a GL₂³ invariant. */
	private static int[][] signatures(Scheme222 s) {
		int[][] sig = new int[R][3];
		for (int k = 0; k < R; k++) {
			for (int slot = 0; slot < 3; slot++) {
				sig[k][slot] = rank2x2(s.m()[k][slot]);
			}
		}
		return sig;
	}

	private static int rank2x2(double[][] x) {
		double a = x[0][0], b = x[0][1], c = x[1][0], d = x[1][1];
		double fro = Math.abs(a) + Math.abs(b) + Math.abs(c) + Math.abs(d);
		if (fro < 1e-12) return 0;
		double det = a * d - b * c;
		return Math.abs(det) > 1e-9 * (1 + fro) ? 2 : 1;
	}

	/** All {@code π} mapping each b-column k to an a-column {@code π[k]} of equal
	 *  signature. Enumerated as the cartesian product of per-signature
	 *  bijections — tiny for real schemes. */
	private static List<int[]> candidatePermutations(int[][] sigA, int[][] sigB) {
		// Group a-columns by signature key.
		java.util.Map<Integer, List<Integer>> byKeyA = new java.util.HashMap<>();
		for (int k = 0; k < R; k++) byKeyA.computeIfAbsent(key(sigA[k]), x -> new ArrayList<>()).add(k);
		java.util.Map<Integer, List<Integer>> byKeyB = new java.util.HashMap<>();
		for (int k = 0; k < R; k++) byKeyB.computeIfAbsent(key(sigB[k]), x -> new ArrayList<>()).add(k);

		// Multisets must match, else no permutation can exist.
		if (!byKeyA.keySet().equals(byKeyB.keySet())) return List.of();
		for (Integer key : byKeyA.keySet()) {
			if (byKeyA.get(key).size() != byKeyB.getOrDefault(key, List.of()).size()) return List.of();
		}

		// Build perm[bCol] = aCol by independently permuting within each group.
		List<int[]> out = new ArrayList<>();
		List<Integer> keys = new ArrayList<>(byKeyA.keySet());
		buildPerms(keys, 0, byKeyA, byKeyB, new int[R], out);
		return out;
	}

	private static void buildPerms(List<Integer> keys, int ki,
			java.util.Map<Integer, List<Integer>> byKeyA,
			java.util.Map<Integer, List<Integer>> byKeyB,
			int[] perm, List<int[]> out) {
		if (ki == keys.size()) { out.add(perm.clone()); return; }
		int key = keys.get(ki);
		List<Integer> aCols = byKeyA.get(key);
		List<Integer> bCols = byKeyB.get(key);
		for (int[] sub : permutationsOf(aCols.size())) {
			for (int t = 0; t < bCols.size(); t++) perm[bCols.get(t)] = aCols.get(sub[t]);
			buildPerms(keys, ki + 1, byKeyA, byKeyB, perm, out);
		}
	}

	private static int key(int[] sig) { return (sig[0] * 3 + sig[1]) * 3 + sig[2]; }

	private static List<int[]> permutationsOf(int n) {
		List<int[]> res = new ArrayList<>();
		permute(new int[n], new boolean[n], 0, res);
		return res;
	}

	private static void permute(int[] cur, boolean[] used, int d, List<int[]> res) {
		if (d == cur.length) { res.add(cur.clone()); return; }
		for (int v = 0; v < cur.length; v++) {
			if (used[v]) continue;
			used[v] = true; cur[d] = v;
			permute(cur, used, d + 1, res);
			used[v] = false;
		}
	}

	// ── Levenberg–Marquardt solve for one fixed permutation ─────────────────

	private static Transform solve(Scheme222 a, Scheme222 b, int[] perm) {
		// Parameters: P(4), Q(4), R(4) row-major.
		MultivariateJacobianFunction model = params -> {
			double[] x = params.toArray();
			double[] res = residualAt(a, b, perm, x);
			// Numeric Jacobian (84×12) by central differences — robust and the
			// problem is tiny, so analytic derivatives aren't worth the risk.
			double[][] jac = new double[res.length][12];
			double h = 1e-7;
			for (int c = 0; c < 12; c++) {
				double[] xp = x.clone(); xp[c] += h;
				double[] xm = x.clone(); xm[c] -= h;
				double[] rp = residualAt(a, b, perm, xp);
				double[] rm = residualAt(a, b, perm, xm);
				for (int rIdx = 0; rIdx < res.length; rIdx++) jac[rIdx][c] = (rp[rIdx] - rm[rIdx]) / (2 * h);
			}
			return new Pair<RealVector, RealMatrix>(new ArrayRealVector(res), new Array2DRowRealMatrix(jac));
		};

		double bestRes = Double.MAX_VALUE;
		double[] bestX = null;
		// A few deterministic starts: identity, then small fixed perturbations.
		double[][] starts = identityStarts(perm.length);
		for (double[] start : starts) {
			try {
				LeastSquaresOptimizer.Optimum opt = new LevenbergMarquardtOptimizer().optimize(
						new LeastSquaresBuilder()
								.start(start)
								.model(model)
								.target(new double[3 * 4 * R])
								.lazyEvaluation(false)
								.maxEvaluations(2000)
								.maxIterations(2000)
								.build());
				double rms = opt.getResiduals().getNorm();
				if (rms < bestRes) {
					bestRes = rms;
					bestX = opt.getPoint().toArray();
					if (bestRes < TOL) break;
				}
			} catch (RuntimeException e) {
				// Singular step / non-convergence for this start — try the next.
			}
		}
		if (bestX == null) return null;
		return new Transform(mat(bestX, 0), mat(bestX, 4), mat(bestX, 8), perm.clone(), bestRes);
	}

	private static double[] residualAt(Scheme222 a, Scheme222 b, int[] perm, double[] x) {
		double[] res = new double[3 * 4 * R];
		int idx = 0;
		double[][] P = mat(x, 0), Q = mat(x, 4), S = mat(x, 8);
		// Matmul-tensor symmetry for our (U:n×m, V:m×p, W:n×p) convention. An
		// element of a⊗b transforms as X ↦ G·X·Hᵀ. The n-space couples U (P) with
		// W (P⁻ᵀ); m couples U (Q) with V (Q⁻ᵀ); p couples V (S) with W (S⁻ᵀ):
		//   Û ↦ P·Û·Qᵀ ;  V̂ ↦ Q⁻ᵀ·V̂·Sᵀ ;  Ŵ ↦ P⁻ᵀ·Ŵ·(S⁻ᵀ)ᵀ = P⁻ᵀ·Ŵ·S⁻¹.
		double[][] Qt = transpose(Q), St = transpose(S);
		double[][] Pit = transpose(inv(P)), Qit = transpose(inv(Q)), Si = inv(S);
		for (int k = 0; k < R; k++) {
			int ak = perm[k];
			idx = diff(res, idx, mul(mul(P, a.m()[ak][0]), Qt), b.m()[k][0]);
			idx = diff(res, idx, mul(mul(Qit, a.m()[ak][1]), St), b.m()[k][1]);
			idx = diff(res, idx, mul(mul(Pit, a.m()[ak][2]), Si), b.m()[k][2]);
		}
		return res;
	}

	/** Deterministic LM starts (no {@code Math.random()} — keeps runs
	 *  reproducible). A nonconvex 12-parameter problem needs several basins:
	 *  identity, sign-flips, a 90° rotation, and small perturbations of each. */
	private static double[][] identityStarts(int n) {
		// Per-2×2 seed bases to combine across the three blocks (P,Q,S).
		double[][] bases = {
				{ 1, 0, 0, 1 },    // I
				{ -1, 0, 0, 1 },   // diag(-1, 1)
				{ 1, 0, 0, -1 },   // diag(1, -1)
				{ 0, 1, -1, 0 },   // 90° rotation
				{ 1, 1, 0, 1 },    // shear
		};
		List<double[]> out = new ArrayList<>();
		// All identical-base triples (P=Q=S=base) — cheap and covers the common
		// "same coordinate change on all three axes" case.
		for (double[] base : bases) {
			double[] x = new double[12];
			for (int blk = 0; blk < 3; blk++) System.arraycopy(base, 0, x, blk * 4, 4);
			out.add(x);
		}
		// A few mixed triples to break the symmetry of the above.
		int[][] mixes = { { 0, 1, 2 }, { 3, 0, 3 }, { 1, 3, 0 }, { 2, 2, 3 } };
		for (int[] mix : mixes) {
			double[] x = new double[12];
			for (int blk = 0; blk < 3; blk++) System.arraycopy(bases[mix[blk]], 0, x, blk * 4, 4);
			out.add(x);
		}
		// Small deterministic perturbations of the identity start.
		for (int s = 1; s <= 3; s++) {
			double[] x = out.get(0).clone();
			for (int c = 0; c < 12; c++) x[c] += 0.05 * Math.sin(1.0 + c + 3.0 * s);
			out.add(x);
		}
		return out.toArray(new double[0][]);
	}

	// ── 2×2 linear algebra helpers ──────────────────────────────────────────

	private static double[][] mat(double[] x, int off) {
		return new double[][] { { x[off], x[off + 1] }, { x[off + 2], x[off + 3] } };
	}

	private static double[][] mul(double[][] a, double[][] b) {
		return new double[][] {
				{ a[0][0] * b[0][0] + a[0][1] * b[1][0], a[0][0] * b[0][1] + a[0][1] * b[1][1] },
				{ a[1][0] * b[0][0] + a[1][1] * b[1][0], a[1][0] * b[0][1] + a[1][1] * b[1][1] } };
	}

	private static double[][] transpose(double[][] a) {
		return new double[][] { { a[0][0], a[1][0] }, { a[0][1], a[1][1] } };
	}

	private static double[][] inv(double[][] a) {
		double det = a[0][0] * a[1][1] - a[0][1] * a[1][0];
		if (Math.abs(det) < 1e-12) det = Math.copySign(1e-12, det == 0 ? 1 : det);
		double id = 1.0 / det;
		return new double[][] { { a[1][1] * id, -a[0][1] * id }, { -a[1][0] * id, a[0][0] * id } };
	}

	private static int diff(double[] res, int idx, double[][] got, double[][] want) {
		res[idx++] = got[0][0] - want[0][0];
		res[idx++] = got[0][1] - want[0][1];
		res[idx++] = got[1][0] - want[1][0];
		res[idx++] = got[1][1] - want[1][1];
		return idx;
	}
}
