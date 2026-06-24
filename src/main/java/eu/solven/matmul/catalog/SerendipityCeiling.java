package eu.solven.matmul.catalog;

import eu.solven.matmul.recombination.Recombination;

/**
 * CERTIFIED ceiling on the serendipitous saving achievable by single-type
 * bud structure — the "how far are we from optimal" bracket (user 2026-06-11).
 * Optimality tier: <strong>optimal-within-scope</strong> — the scope is
 * (single bud type, current catalog enlarged-block prices, partition
 * combinatorics + the spanning/divisibility constraints below). The ceiling is
 * an upper bound on savings that may well be UNREACHABLE (it ignores whether
 * the class slice-matrices are realizable at their class sizes); the catalog
 * value is a lower bound; the optimum lies in between.
 *
 * <p>Constraints used (both are theorems, not heuristics):</p>
 * <ul>
 *   <li><b>Spanning</b>: the distinct directions of one factor's columns must
 *       span that factor's flattening row space, which for the matmul tensor
 *       has FULL rank — so #classes ≥ n·m (U), m·p (V), n·p (W).</li>
 *   <li><b>Contracted-axis divisibility</b>: if the direction count equals the
 *       flattening rank exactly, the directions are a basis and every class
 *       matrix is a slice combination — whose rank is a positive multiple of
 *       the contracted dimension (a λ-combination of W-slices is m disjoint
 *       copies of λ). Every class then needs ≥ m terms, total ≥ m·n·p; when
 *       that exceeds the scheme rank, the minimal class count rises by one.</li>
 * </ul>
 *
 * <p>σ prices come from the catalog ({@code σ(k) = k·R(inner) −
 * R(k-enlarged)}), so the ceiling moves when the catalog improves — it is a
 * bound on what serendipitous products can deliver TODAY, not an absolute.</p>
 */
public final class SerendipityCeiling {

	private SerendipityCeiling() {}

	/** Per-axis result: ceiling on savings and the matching floor on the
	 *  predicted product rank ({@code rank·R(inner) − maxSavings}). */
	public record AxisCeiling(long maxSavings, long productRankFloor, int minClasses) {}

	/**
	 * Ceiling for buds of {@code type} on a rank-{@code r} base ⟨n,m,p⟩
	 * composed with inner ⟨n2,m2,p2⟩, priced by {@code lookup}.
	 */
	public static AxisCeiling forAxis(SerendipitousBudProduct.BudType type, int n, int m, int p,
			int r, FieldAwareLookup lookup, int n2, int m2, int p2) {
		long inner = lookup.findRank(n2, m2, p2);
		if (inner >= Recombination.SotaResolver.UNKNOWN_RANK) {
			return new AxisCeiling(0, Long.MAX_VALUE / 4, 0);
		}
		int minClasses = minClasses(type, n, m, p, r);
		long[] sigma = sigmaTable(type, r, lookup, n2, m2, p2);
		long best = maxSavings(r, minClasses, sigma);
		return new AxisCeiling(best, (long) r * inner - best, minClasses);
	}

	/** Spanning bound (flattening rank) + contracted-axis divisibility bump. */
	public static int minClasses(SerendipitousBudProduct.BudType type, int n, int m, int p, int r) {
		int flat = switch (type) { case U -> n * m; case V -> m * p; case W -> n * p; };
		int contracted = switch (type) { case U -> p; case V -> n; case W -> m; };
		return flat + ((long) contracted * flat > r ? 1 : 0);
	}

	/** σ(k) = k·R(inner) − R(k-enlarged along the axis), floored at 0; index k. */
	public static long[] sigmaTable(SerendipitousBudProduct.BudType type, int r,
			FieldAwareLookup lookup, int n2, int m2, int p2) {
		long inner = lookup.findRank(n2, m2, p2);
		long[] sigma = new long[r + 1];
		if (inner >= Recombination.SotaResolver.UNKNOWN_RANK) {
			return sigma;
		}
		for (int k = 2; k <= r; k++) {
			long enlarged = switch (type) {
				case U -> lookup.findRank(n2, m2, k * p2);
				case V -> lookup.findRank(k * n2, m2, p2);
				case W -> lookup.findRank(n2, k * m2, p2);
			};
			// No catalog price for the enlarged block → the bud is unbuildable
			// today → contributes 0 within this scope.
			sigma[k] = enlarged >= Recombination.SotaResolver.UNKNOWN_RANK ? 0
					: Math.max(0, k * inner - enlarged);
		}
		return sigma;
	}

	/**
	 * Max {@code Σ σ(k_c)} over partitions of {@code r} into at least
	 * {@code minClasses} parts. Small DP: {@code f[c][s]} = best savings with
	 * {@code c} classes (capped at {@code minClasses}) summing to {@code s}.
	 */
	static long maxSavings(int r, int minClasses, long[] sigma) {
		int cap = Math.min(minClasses, r);
		long[][] f = new long[cap + 1][r + 1];
		for (long[] row : f) {
			java.util.Arrays.fill(row, Long.MIN_VALUE / 4);
		}
		f[0][0] = 0;
		for (int c = 0; c < cap + 1; c++) {
			for (int s = 0; s <= r; s++) {
				if (f[c][s] < 0) {
					continue;
				}
				for (int k = 1; s + k <= r; k++) {
					int c2 = Math.min(c + 1, cap);
					long v = f[c][s] + (k < sigma.length ? sigma[k] : 0);
					if (v > f[c2][s + k]) {
						f[c2][s + k] = v;
					}
				}
			}
		}
		return Math.max(0, f[cap][r]);
	}
}
