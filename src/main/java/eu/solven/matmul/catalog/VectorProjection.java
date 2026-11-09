package eu.solven.matmul.catalog;

import java.util.ArrayList;
import java.util.BitSet;
import java.util.LinkedHashMap;
import java.util.List;
import java.util.Locale;
import java.util.Map;
import java.util.Optional;

import eu.solven.matmul.NonCubicBilinearAlgorithm;
import eu.solven.matmul.verifiers.Verifier;

import lombok.extern.slf4j.Slf4j;

/**
 * Projection <em>along a direction</em> — FMM-Lille's {@code ProjOfAlongMode} — expressed
 * as a change of basis followed by the ordinary coordinate projection
 * ({@link ProjectionSearch}).
 *
 * <p>{@link ProjectionSearch} drops a <em>coordinate</em> of one axis and dead-code-eliminates
 * the products whose support on that axis, in either carrying matrix, is that single
 * coordinate. That is basis-dependent: Smirnov's ⟨3,3,6⟩=40 has coordinate margins
 * {@code [0,0,0]} in every representative we hold (and in FMM's own artifact), yet FMM
 * prices {@code ⟨9,18,18⟩ = ⟨3,3,6⟩ ⊗ ⟨3,6,3⟩} down to ⟨9,17,18⟩:1591 and ⟨9,17,17⟩:1574
 * (2026-11, issue #23). The products it kills are those whose factor on the mode is
 * <b>rank-1 with a common direction</b> {@code d} — after a change of basis that sends
 * {@code d} to a coordinate they die under a plain coordinate drop. The matmul tensor is
 * invariant under {@code GL_n × GL_m × GL_p} ({@code A → H·A·G, B → G⁻¹·B·K, C → H·C·K}), so
 * the rebased scheme is an exact scheme for the same format, GL-equivalent to the original
 * — the same rank, a different (margin-rich) <em>representative</em>: the projection
 * sibling of the bud-representative caveat (fmm-gap skill §4c).</p>
 *
 * <p>The mechanics, for the m-axis ({@code j}): a product dies when coordinate {@code j} is
 * dropped iff {@code u_l} (n×m) is supported on column {@code j} only, or {@code v_l} (m×p) on
 * row {@code j} only. Under {@code A → A·G, B → G⁻¹·B} the scheme becomes
 * {@code u'_l = u_l·G⁻ᵀ, v'_l = Gᵀ·v_l} and {@code u'_l} is single-column-{@code j} iff
 * {@code u_l = a·dᵀ} with {@code d ∝ G·e_j} (column {@code j} of {@code G}: the
 * <b>primal</b> side), while {@code v'_l} is single-row-{@code j} iff {@code v_l = d·bᵀ} with
 * {@code d ∝ G⁻ᵀ·e_j} (row {@code j} of {@code G⁻¹}: the <b>dual</b> side). So one basis can
 * serve a primal direction {@code d_p} and a dual direction {@code d_q} with
 * {@code d_q·d_p = 1}: {@code G = [basis of ker d_qᵀ | d_p]}. The n-axis (primal W-rows, dual
 * U-rows, {@code H = G⁻¹}) and the p-axis (primal V-columns, dual W-columns, {@code K = G})
 * follow the same pattern.</p>
 *
 * <p>Every result is re-verified ({@link Verifier#isExactNonCubic}); a transform that does
 * not verify is dropped, never returned. Kills are counted on distinct products. Tier of
 * the choice: <b>bound</b> — the most populous (primal, dual) pair per axis, greedily per
 * axis; no claim that it is the best basis.</p>
 */
@Slf4j
public final class VectorProjection {
	private VectorProjection() {}

	/** Tolerance for "zero" / "proportional" on the double-valued factor matrices. */
	static final double EPS = 1e-9;

	/** One rebased axis: the GL-equivalent scheme, the axis, how many products the drop of
	 *  coordinate {@code dim−1} of that axis now kills, and the basis used. */
	public record Rebased(NonCubicBilinearAlgorithm scheme, int axis, int kills, double[][] basis) {}

	/** A common rank-1 direction on one side of an axis and the products carrying it. */
	record Group(double[] dir, BitSet terms) {}

	/**
	 * Greedy rebase of every axis in turn (n, then m, then p — the three group actions
	 * commute, each acts on its own index space). Returns the input itself when no axis
	 * offers a direction, so callers can test identity by reference.
	 */
	public static NonCubicBilinearAlgorithm rebaseAll(NonCubicBilinearAlgorithm a) {
		NonCubicBilinearAlgorithm cur = a;
		for (int axis = 0; axis < 3; axis++) {
			Optional<Rebased> r = rebase(cur, axis);
			if (r.isPresent()) {
				cur = r.get().scheme();
			}
		}
		return cur;
	}

	/**
	 * Rebase one axis ({@code 0 = n, 1 = m, 2 = p}) so that the most populous rank-1
	 * direction pair becomes coordinate {@code dim−1}. Empty when no product is rank-1 along
	 * that axis on either side, or when the transform does not verify.
	 */
	public static Optional<Rebased> rebase(NonCubicBilinearAlgorithm a, int axis) {
		int dim = switch (axis) {
			case 0 -> a.n;
			case 1 -> a.m;
			default -> a.p;
		};
		double[][] U = a.denseU(), V = a.denseV(), W = a.denseW();
		Map<String, Group> primal = new LinkedHashMap<>(), dual = new LinkedHashMap<>();
		for (int l = 0; l < a.r; l++) {
			double[][] u = term(U, l, a.n, a.m), v = term(V, l, a.m, a.p), w = term(W, l, a.n, a.p);
			switch (axis) {
				case 0 -> { // n: primal = W rows, dual = U rows
					collect(primal, rowDirection(w), l);
					collect(dual, rowDirection(u), l);
				}
				case 1 -> { // m: primal = U columns, dual = V rows
					collect(primal, colDirection(u), l);
					collect(dual, rowDirection(v), l);
				}
				default -> { // p: primal = V columns, dual = W columns
					collect(primal, colDirection(v), l);
					collect(dual, colDirection(w), l);
				}
			}
		}
		if (primal.isEmpty() && dual.isEmpty()) {
			return Optional.empty();
		}
		// The most populous compatible (primal, dual) pair; a lone side pairs with a unit vector.
		Group bestP = null, bestQ = null;
		int best = 0;
		for (Group p : primal.values()) {
			int solo = p.terms().cardinality();
			if (solo > best) {
				best = solo;
				bestP = p;
				bestQ = null;
			}
			for (Group q : dual.values()) {
				if (Math.abs(dot(p.dir(), q.dir())) <= EPS) {
					continue;
				}
				BitSet both = (BitSet) p.terms().clone();
				both.or(q.terms());
				if (both.cardinality() > best) {
					best = both.cardinality();
					bestP = p;
					bestQ = q;
				}
			}
		}
		for (Group q : dual.values()) {
			int solo = q.terms().cardinality();
			if (solo > best) {
				best = solo;
				bestP = null;
				bestQ = q;
			}
		}
		if (best == 0) {
			return Optional.empty();
		}
		double[] dp = bestP == null ? null : bestP.dir().clone();
		double[] dq = bestQ == null ? null : bestQ.dir().clone();
		if (dp == null) {
			int i = argmaxAbs(dq);
			dp = new double[dim];
			dp[i] = 1.0 / dq[i];
		} else if (dq == null) {
			int i = argmaxAbs(dp);
			dq = new double[dim];
			dq[i] = 1.0 / dp[i];
		} else {
			double s = dot(dq, dp);
			for (int i = 0; i < dim; i++) {
				dp[i] /= s;
			}
		}
		double[][] G = basisWith(dp, dq);
		NonCubicBilinearAlgorithm out;
		try {
			out = changeBasis(a, axis, G);
		} catch (RuntimeException e) {
			log.warn("rebase axis {} of ⟨{},{},{}⟩=r{} failed: {}", axis, a.n, a.m, a.p, a.r, e.toString());
			return Optional.empty();
		}
		if (!Verifier.isExactNonCubic(out)) {
			log.warn("rebase axis {} of ⟨{},{},{}⟩=r{} does not verify — dropped", axis, a.n, a.m, a.p, a.r);
			return Optional.empty();
		}
		int margin = ProjectionSearch.axisMargins(out)[axis];
		if (margin < best) {
			log.warn("rebase axis {} of ⟨{},{},{}⟩: expected ≥{} kills, coordinate margin is {} — dropped", axis, a.n,
					a.m, a.p, best, margin);
			return Optional.empty();
		}
		return Optional.of(new Rebased(out, axis, best, G));
	}

	/**
	 * The {@code GL} action on one axis with matrix {@code G} (columns = new basis), i.e. the
	 * scheme for the same format in the coordinates where column {@code j} of {@code G} is
	 * {@code e_j}: m-axis {@code u' = u·G⁻ᵀ, v' = Gᵀ·v}; n-axis {@code u' = Gᵀ·u, w' = G⁻¹·w};
	 * p-axis {@code v' = v·G⁻ᵀ, w' = w·G}. Exact for any invertible {@code G}; callers verify.
	 */
	public static NonCubicBilinearAlgorithm changeBasis(NonCubicBilinearAlgorithm a, int axis, double[][] G) {
		double[][] Ginv = invert(G);
		double[][] Gt = transpose(G), GinvT = transpose(Ginv);
		double[][] U = a.denseU(), V = a.denseV(), W = a.denseW();
		double[][] U2 = new double[U.length][a.r], V2 = new double[V.length][a.r], W2 = new double[W.length][a.r];
		for (int l = 0; l < a.r; l++) {
			double[][] u = term(U, l, a.n, a.m), v = term(V, l, a.m, a.p), w = term(W, l, a.n, a.p);
			switch (axis) {
				case 0 -> {
					u = mul(Gt, u);
					w = mul(Ginv, w);
				}
				case 1 -> {
					u = mul(u, GinvT);
					v = mul(Gt, v);
				}
				default -> {
					v = mul(v, GinvT);
					w = mul(w, G);
				}
			}
			store(U2, l, u);
			store(V2, l, v);
			store(W2, l, w);
		}
		return new NonCubicBilinearAlgorithm(a.n, a.m, a.p, clean(U2), clean(V2), clean(W2));
	}

	// ---- rank-1 directions -------------------------------------------------------------

	/** {@code X = d·bᵀ} → the normalised column direction {@code d} (rows), else null. */
	static double[] rowDirection(double[][] X) {
		int rows = X.length, cols = X[0].length;
		int r0 = -1, c0 = -1;
		outer:
		for (int r = 0; r < rows; r++) {
			for (int c = 0; c < cols; c++) {
				if (Math.abs(X[r][c]) > EPS) {
					r0 = r;
					c0 = c;
					break outer;
				}
			}
		}
		if (r0 < 0) {
			return null;
		}
		double[] d = new double[rows];
		for (int r = 0; r < rows; r++) {
			d[r] = X[r][c0] / X[r0][c0];
		}
		for (int r = 0; r < rows; r++) {
			for (int c = 0; c < cols; c++) {
				if (Math.abs(X[r][c] - d[r] * X[r0][c]) > EPS) {
					return null;
				}
			}
		}
		return d;
	}

	/** {@code X = a·dᵀ} → the normalised row direction {@code d} (columns), else null. */
	static double[] colDirection(double[][] X) {
		double[] d = rowDirection(transpose(X));
		return d;
	}

	private static void collect(Map<String, Group> into, double[] dir, int term) {
		if (dir == null) {
			return;
		}
		StringBuilder key = new StringBuilder();
		for (double x : dir) {
			key.append(String.format(Locale.ROOT, "%.9f", Math.abs(x) < EPS ? 0.0 : x)).append(',');
		}
		into.computeIfAbsent(key.toString(), k -> new Group(dir, new BitSet())).terms().set(term);
	}

	/** {@code G = [basis of ker(dqᵀ) | dp]} with {@code dq·dp = 1}: column {@code dim−1} is the
	 *  primal direction, row {@code dim−1} of {@code G⁻¹} the dual one. */
	static double[][] basisWith(double[] dp, double[] dq) {
		int dim = dp.length;
		double[][] G = new double[dim][dim];
		int i0 = argmaxAbs(dq);
		int col = 0;
		for (int k = 0; k < dim; k++) {
			if (k == i0) {
				continue;
			}
			G[k][col] += 1.0;
			G[i0][col] -= dq[k] / dq[i0];
			col++;
		}
		for (int i = 0; i < dim; i++) {
			G[i][dim - 1] = dp[i];
		}
		return G;
	}

	// ---- small dense linear algebra ----------------------------------------------------

	static double[][] term(double[][] F, int l, int rows, int cols) {
		double[][] t = new double[rows][cols];
		for (int i = 0; i < rows; i++) {
			for (int j = 0; j < cols; j++) {
				t[i][j] = F[i * cols + j][l];
			}
		}
		return t;
	}

	private static void store(double[][] F, int l, double[][] t) {
		int cols = t[0].length;
		for (int i = 0; i < t.length; i++) {
			for (int j = 0; j < cols; j++) {
				F[i * cols + j][l] = t[i][j];
			}
		}
	}

	static double[][] mul(double[][] A, double[][] B) {
		int n = A.length, k = B.length, m = B[0].length;
		if (A[0].length != k) {
			throw new IllegalArgumentException("mul: " + n + "x" + A[0].length + " · " + k + "x" + m);
		}
		double[][] C = new double[n][m];
		for (int i = 0; i < n; i++) {
			for (int t = 0; t < k; t++) {
				double a = A[i][t];
				if (a == 0) {
					continue;
				}
				for (int j = 0; j < m; j++) {
					C[i][j] += a * B[t][j];
				}
			}
		}
		return C;
	}

	static double[][] transpose(double[][] A) {
		double[][] T = new double[A[0].length][A.length];
		for (int i = 0; i < A.length; i++) {
			for (int j = 0; j < A[0].length; j++) {
				T[j][i] = A[i][j];
			}
		}
		return T;
	}

	/** Gauss–Jordan with partial pivoting; throws on a singular matrix. */
	static double[][] invert(double[][] G) {
		int n = G.length;
		double[][] a = new double[n][2 * n];
		for (int i = 0; i < n; i++) {
			System.arraycopy(G[i], 0, a[i], 0, n);
			a[i][n + i] = 1.0;
		}
		for (int c = 0; c < n; c++) {
			int piv = c;
			for (int r = c + 1; r < n; r++) {
				if (Math.abs(a[r][c]) > Math.abs(a[piv][c])) {
					piv = r;
				}
			}
			if (Math.abs(a[piv][c]) <= EPS) {
				throw new IllegalArgumentException("singular basis matrix");
			}
			double[] tmp = a[c];
			a[c] = a[piv];
			a[piv] = tmp;
			double d = a[c][c];
			for (int j = 0; j < 2 * n; j++) {
				a[c][j] /= d;
			}
			for (int r = 0; r < n; r++) {
				if (r == c) {
					continue;
				}
				double f = a[r][c];
				if (f == 0) {
					continue;
				}
				for (int j = 0; j < 2 * n; j++) {
					a[r][j] -= f * a[c][j];
				}
			}
		}
		double[][] inv = new double[n][n];
		for (int i = 0; i < n; i++) {
			System.arraycopy(a[i], n, inv[i], 0, n);
		}
		return inv;
	}

	/** Snap floating noise to exact zero / nearby small rationals (the scheme coefficients are
	 *  small rationals; a product of two such basis matrices stays one). */
	private static double[][] clean(double[][] F) {
		for (double[] row : F) {
			for (int j = 0; j < row.length; j++) {
				double x = row[j];
				if (Math.abs(x) < EPS) {
					row[j] = 0.0;
					continue;
				}
				double r = Math.rint(x * 4096.0) / 4096.0; // denominators up to 4096
				if (Math.abs(r - x) < 1e-9) {
					row[j] = r;
				}
			}
		}
		return F;
	}

	static double dot(double[] a, double[] b) {
		double s = 0;
		for (int i = 0; i < a.length; i++) {
			s += a[i] * b[i];
		}
		return s;
	}

	static int argmaxAbs(double[] v) {
		int best = 0;
		for (int i = 1; i < v.length; i++) {
			if (Math.abs(v[i]) > Math.abs(v[best])) {
				best = i;
			}
		}
		return best;
	}

	/** Per-axis direction census — how many products are rank-1 along each axis on the primal and
	 *  dual side, and the size of the largest common-direction group: the diagnostic behind a
	 *  "why does FMM project this and we cannot" question. */
	public static String census(NonCubicBilinearAlgorithm a) {
		StringBuilder sb = new StringBuilder();
		String[] names = { "n", "m", "p" };
		for (int axis = 0; axis < 3; axis++) {
			Optional<Rebased> r = rebase(a, axis);
			sb.append(names[axis]).append(": best direction kills ")
					.append(r.map(Rebased::kills).orElse(0)).append("; ");
		}
		return sb.toString();
	}

	/** Convenience: the per-axis kills a full rebase would expose, without applying it. */
	public static int[] directionMargins(NonCubicBilinearAlgorithm a) {
		int[] out = new int[3];
		for (int axis = 0; axis < 3; axis++) {
			out[axis] = rebase(a, axis).map(Rebased::kills).orElse(0);
		}
		return out;
	}

	/** All distinct rebased variants worth offering as projection parents: the original
	 *  (identity), each single-axis rebase, and the greedy all-axes one — de-duplicated by
	 *  content hash. Small schemes only (factors), so the cost is negligible. */
	public static List<NonCubicBilinearAlgorithm> variants(NonCubicBilinearAlgorithm a) {
		Map<String, NonCubicBilinearAlgorithm> out = new LinkedHashMap<>();
		out.put(SchemeIO.contentHash(a), a);
		for (int axis = 0; axis < 3; axis++) {
			rebase(a, axis).ifPresent(r -> out.putIfAbsent(SchemeIO.contentHash(r.scheme()), r.scheme()));
		}
		NonCubicBilinearAlgorithm all = rebaseAll(a);
		out.putIfAbsent(SchemeIO.contentHash(all), all);
		return new ArrayList<>(out.values());
	}
}
