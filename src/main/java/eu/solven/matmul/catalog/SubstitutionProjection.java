package eu.solven.matmul.catalog;

import java.util.ArrayList;
import java.util.List;

import eu.solven.matmul.NonCubicBilinearAlgorithm;

/**
 * Projection by <em>substitution</em> — FMM-Lille's {@code ProjOfAlongMode}
 * ({@code https://fmm.univ-lille.fr/projalongmode.html}, Lemma 1): from a decomposition of
 * {@code ⟨n,m,p+1⟩} of rank {@code r} one constructs {@code ⟨n,m,p⟩} of rank {@code ≤ r − m};
 * likewise {@code ⟨n,m+1,p⟩ → ⟨n,m,p⟩} at {@code r − n} and {@code ⟨n+1,m,p⟩ → ⟨n,m,p⟩} at
 * {@code r − p}. It is the constructive face of the substitution method.
 *
 * <p>Coordinate projection ({@link Compose#project}) sets the dropped slice's variables to
 * zero and removes only the products whose factor lived on that slice alone — for a Kronecker
 * product of generic factors that is usually none ({@code ⟨9,18,18⟩ = ⟨3,3,6⟩ ⊗ ⟨3,6,3⟩}: 0
 * on every axis, in every representative we hold and in FMM's own artifact). Substitution
 * sets them to <em>linear forms in the kept variables</em> instead. For the p-case: write
 * {@code B̃ = [B' | λ(B')]} with {@code λ : R^{m×p} → R^m} linear; since the scheme computes
 * {@code A·B̃} for every {@code B̃}, its first {@code p} output columns are {@code A·B'} for
 * every {@code B'}. Each product's B-form becomes {@code ṽ_l = v_l|_{B'} + Σ_j c_l[j]·λ_j}
 * where {@code c_l[j]} is the coefficient of the dropped variable {@code b_{j,p+1}} in
 * {@code v_l}. Pick {@code m} products whose vectors {@code c_l} are linearly independent and
 * solve the {@code m×m} system so that exactly those {@code ṽ_l} vanish: {@code λ = −C⁻¹·V_rest}.
 * The scheme on the remaining products, with the dropped output column removed, is an exact
 * {@code ⟨n,m,p⟩} decomposition. Products that vanish incidentally are removed too.</p>
 *
 * <p>The m- and n-cases are the same lemma on a cyclic orientation of the tensor
 * ({@link NonCubicBilinearAlgorithm#orientByPerm}), so one implementation serves all three.
 * The result is <b>exact by construction</b> and is re-verified by every writer; its rank is a
 * <b>bound</b> (which {@code m} products are killed is the deterministic greedy choice below —
 * the first {@code m} independent ones — so a lineage node {@code SubstProject(child, axis,
 * index)} replays bit-for-bit). Dropping the axes in a different order can give a different
 * rank: {@code ⟨9,18,18⟩:1600 → ⟨9,17,17⟩} is {@code 1600−9−17 = 1574} (FMM's order) or
 * {@code 1600−18−9 = 1573}.</p>
 *
 * <p>2026-11-09: the 61 FMM index moves of issue #23 were all this lemma applied to products
 * we already held; the catalog's projection machinery, being coordinate-only, could not see
 * them (research/DISCOVERIES_PENDING_ANALYSIS.md).</p>
 */
public final class SubstitutionProjection {
	private SubstitutionProjection() {}

	static final double EPS = 1e-9;

	/** Perms (target-axis ← source-axis) that send source axis {@code a} into the p slot. */
	private static final int[][] TO_P = {
			{ 1, 2, 0 }, // n → p   (target n ← source m, target m ← source p, target p ← source n)
			{ 2, 0, 1 }, // m → p   (target n ← source p, target m ← source n, target p ← source m)
			{ 0, 1, 2 }, // p → p
	};

	/** How many products the lemma removes for a drop on {@code axis}: p-drop {@code m}, m-drop {@code n}, n-drop {@code p}. */
	public static int guaranteedKills(NonCubicBilinearAlgorithm a, int axis) {
		return switch (axis) {
			case 0 -> a.p;
			case 1 -> a.n;
			default -> a.m;
		};
	}

	/**
	 * Drop index {@code index} of {@code axis} ({@code 0 = n, 1 = m, 2 = p}) by substitution.
	 * Returns an exact scheme for the format with that axis one shorter, of rank
	 * {@code ≤ r − guaranteedKills}. Callers re-verify ({@code Verifier.isExactNonCubic}).
	 */
	public static NonCubicBilinearAlgorithm drop(NonCubicBilinearAlgorithm a, int axis, int index) {
		if (axis < 0 || axis > 2) {
			throw new IllegalArgumentException("axis must be 0 (n), 1 (m) or 2 (p): " + axis);
		}
		int dim = switch (axis) {
			case 0 -> a.n;
			case 1 -> a.m;
			default -> a.p;
		};
		if (index < 0 || index >= dim) {
			throw new IllegalArgumentException("index " + index + " out of [0," + dim + ") on axis " + axis);
		}
		if (dim < 2) {
			throw new IllegalArgumentException("cannot drop the only index of axis " + axis);
		}
		if (axis == 2) {
			return dropP(a, index);
		}
		int[] perm = TO_P[axis];
		NonCubicBilinearAlgorithm oriented = a.orientByPerm(perm);
		NonCubicBilinearAlgorithm dropped = dropP(oriented, index);
		return dropped.orientByPerm(inverse(perm));
	}

	/** The index on {@code axis} whose drop leaves the fewest products (ties → lowest index). */
	public static int bestIndex(NonCubicBilinearAlgorithm a, int axis) {
		int dim = switch (axis) {
			case 0 -> a.n;
			case 1 -> a.m;
			default -> a.p;
		};
		int best = 0, bestRank = Integer.MAX_VALUE;
		for (int i = 0; i < dim; i++) {
			int r = drop(a, axis, i).r;
			if (r < bestRank) {
				bestRank = r;
				best = i;
			}
		}
		return best;
	}

	private static int[] inverse(int[] perm) {
		int[] inv = new int[perm.length];
		for (int t = 0; t < perm.length; t++) {
			inv[perm[t]] = t;
		}
		return inv;
	}

	/** The p-case proper: drop column {@code k} of B and of C. */
	static NonCubicBilinearAlgorithm dropP(NonCubicBilinearAlgorithm a, int k) {
		int n = a.n, m = a.m, p = a.p, r = a.r, p2 = p - 1;
		double[][] U = a.denseU(), V = a.denseV(), W = a.denseW();
		// c_l[j] = coefficient of b_{j,k} in v_l ;  vRest_l = v_l without column k (m × p2, row-major)
		double[][] c = new double[r][m];
		double[][] vRest = new double[r][m * p2];
		for (int l = 0; l < r; l++) {
			for (int j = 0; j < m; j++) {
				c[l][j] = V[j * p + k][l];
				int t = 0;
				for (int col = 0; col < p; col++) {
					if (col == k) {
						continue;
					}
					vRest[l][j * p2 + t++] = V[j * p + col][l];
				}
			}
		}
		// Greedy: the first m products whose c_l are linearly independent (deterministic → replayable).
		List<Integer> chosen = new ArrayList<>();
		List<double[]> echelon = new ArrayList<>();
		for (int l = 0; l < r && chosen.size() < m; l++) {
			double[] v = c[l].clone();
			for (double[] e : echelon) {
				int piv = pivot(e);
				double f = v[piv] / e[piv];
				if (f != 0) {
					for (int j = 0; j < m; j++) {
						v[j] -= f * e[j];
					}
				}
			}
			if (norm(v) > EPS) {
				echelon.add(v);
				chosen.add(l);
			}
		}
		int q = chosen.size(); // = m unless the decomposition is degenerate on this slice
		// Solve C·Λ = −V_rest over the chosen rows: Λ (q × m·p2) with Λ_j the form for b_{j,k}.
		// When q < m the unchosen directions get λ = 0 (plain coordinate drop on them).
		double[][] C = new double[q][m];
		double[][] R = new double[q][m * p2];
		for (int i = 0; i < q; i++) {
			C[i] = c[chosen.get(i)];
			R[i] = vRest[chosen.get(i)];
		}
		double[][] lambda = solveLambda(C, R, m, p2); // m × (m·p2), rows for unchosen directions zero
		// ṽ_l = vRest_l + Σ_j c_l[j] λ_j ; drop products whose ṽ vanishes (the chosen ones by construction).
		List<Integer> alive = new ArrayList<>();
		double[][] vNew = new double[r][];
		for (int l = 0; l < r; l++) {
			double[] v = vRest[l].clone();
			for (int j = 0; j < m; j++) {
				double cj = c[l][j];
				if (cj == 0) {
					continue;
				}
				double[] lj = lambda[j];
				for (int t = 0; t < v.length; t++) {
					v[t] += cj * lj[t];
				}
			}
			if (norm(v) > EPS) {
				snap(v);
				vNew[l] = v;
				alive.add(l);
			}
		}
		int r2 = alive.size();
		double[][] U2 = new double[n * m][r2], V2 = new double[m * p2][r2], W2 = new double[n * p2][r2];
		for (int col = 0; col < r2; col++) {
			int l = alive.get(col);
			for (int row = 0; row < n * m; row++) {
				U2[row][col] = U[row][l];
			}
			for (int row = 0; row < m * p2; row++) {
				V2[row][col] = vNew[l][row];
			}
			for (int i = 0; i < n; i++) {
				int t = 0;
				for (int colP = 0; colP < p; colP++) {
					if (colP == k) {
						continue;
					}
					W2[i * p2 + t++][col] = W[i * p + colP][l];
				}
			}
		}
		return new NonCubicBilinearAlgorithm(n, m, p2, U2, V2, W2);
	}

	/** {@code Λ = −C⁺·R} on the chosen directions: Gauss–Jordan on {@code [C | −R]} (C is q×m of
	 *  rank q); the solution is expressed in the pivot columns, other λ rows stay zero. */
	private static double[][] solveLambda(double[][] C, double[][] R, int m, int p2) {
		int q = C.length, w = m * p2;
		double[][] aug = new double[q][m + w];
		for (int i = 0; i < q; i++) {
			System.arraycopy(C[i], 0, aug[i], 0, m);
			for (int t = 0; t < w; t++) {
				aug[i][m + t] = -R[i][t];
			}
		}
		int[] pivotCol = new int[q];
		int row = 0;
		for (int col = 0; col < m && row < q; col++) {
			int piv = row;
			for (int i = row + 1; i < q; i++) {
				if (Math.abs(aug[i][col]) > Math.abs(aug[piv][col])) {
					piv = i;
				}
			}
			if (Math.abs(aug[piv][col]) <= EPS) {
				continue;
			}
			double[] tmp = aug[row];
			aug[row] = aug[piv];
			aug[piv] = tmp;
			double d = aug[row][col];
			for (int t = 0; t < m + w; t++) {
				aug[row][t] /= d;
			}
			for (int i = 0; i < q; i++) {
				if (i == row) {
					continue;
				}
				double f = aug[i][col];
				if (f == 0) {
					continue;
				}
				for (int t = 0; t < m + w; t++) {
					aug[i][t] -= f * aug[row][t];
				}
			}
			pivotCol[row] = col;
			row++;
		}
		double[][] lambda = new double[m][w];
		for (int i = 0; i < row; i++) {
			System.arraycopy(aug[i], m, lambda[pivotCol[i]], 0, w);
		}
		return lambda;
	}

	private static int pivot(double[] e) {
		for (int j = 0; j < e.length; j++) {
			if (Math.abs(e[j]) > EPS) {
				return j;
			}
		}
		throw new IllegalStateException("zero echelon row");
	}

	private static double norm(double[] v) {
		double s = 0;
		for (double x : v) {
			s = Math.max(s, Math.abs(x));
		}
		return s;
	}

	/** Snap floating noise: exact zeros, and small rationals (denominators up to 4096). */
	private static void snap(double[] v) {
		for (int i = 0; i < v.length; i++) {
			double x = v[i];
			if (Math.abs(x) < EPS) {
				v[i] = 0.0;
				continue;
			}
			double rr = Math.rint(x * 4096.0) / 4096.0;
			if (Math.abs(rr - x) < 1e-9) {
				v[i] = rr;
			}
		}
	}
}
