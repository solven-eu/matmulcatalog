package eu.solven.matmul.search;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.Comparator;
import java.util.LinkedHashMap;
import java.util.List;
import java.util.Map;

import eu.solven.matmul.NonCubicBilinearAlgorithm;
import eu.solven.matmul.catalog.Recombination.SotaResolver;

/**
 * Fast analytical scoring of axis-flip mask variants for an outer scheme at a
 * fixed allocation, without rebuilding factor matrices.
 *
 * <p>Equivalence exploited: applying axis-flip mask {@code (sA, sB, sC)} to a
 * scheme and evaluating at allocation {@code (allocA, allocB, allocC)} is
 * exactly equivalent to evaluating the canonical scheme at allocation
 * {@code (reverse-if-sA(allocA), reverse-if-sB(allocB), reverse-if-sC(allocC))}.
 * The block-support structure is fixed; only the allocation labels change.
 *
 * <p>So instead of calling
 * {@link eu.solven.matmul.catalog.Recombination#recombineWithAllocation} eight
 * times (one per mask) and re-scanning U/V/W each time, we extract per-product
 * block supports ONCE from the canonical scheme then re-score the 8 mask-
 * permuted allocations in O(r) each.
 *
 * <p><b>Scope</b>: any outer scheme, but typical use is rank-7 ⟨2,2,2⟩ where
 * the axis-flip orbit dominates. Larger axes still work — reversing a length-3
 * allocation permutes blocks {@code [0,1,2] → [2,1,0]} which corresponds to
 * the {@code J} (anti-diagonal) axis-flip in {@link
 * eu.solven.matmul.SymmetryTransforms#applyAxisFlip}.
 *
 * <p><b>What this module does NOT model</b>: zero-peel (γ5 reduction). This
 * scorer treats {@code alloc} as effective sizes. Callers that depend on peel
 * should request top-K candidates (K ≥ 2-3) here and run brute-force
 * peel-aware materialisation on each to pick the actual winner. The
 * analytical scorer is a fast pruner; peel-aware {@code recombineWithAllocation}
 * remains the source of truth.
 */
public final class AnalyticalMaskSearch {

	private AnalyticalMaskSearch() {}

	/**
	 * Per-product row/column block supports of a scheme's U/V/W factors.
	 * Extracted once per canonical scheme, reused across many mask × allocation
	 * scorings.
	 */
	public static final class SchemeSupports {
		public final int n, m, p, r;
		/** {@code uRowSupport[k]} = sorted distinct A-row indices {@code i} with U[i·m+j][k] ≠ 0 for some j. */
		public final int[][] uRowSupport;
		/** {@code uColSupport[k]} = sorted distinct A-col indices {@code j} with U[i·m+j][k] ≠ 0 for some i. */
		public final int[][] uColSupport;
		public final int[][] vRowSupport;
		public final int[][] vColSupport;
		public final int[][] wRowSupport;
		public final int[][] wColSupport;

		private SchemeSupports(int n, int m, int p, int r,
				int[][] uR, int[][] uC, int[][] vR, int[][] vC, int[][] wR, int[][] wC) {
			this.n = n; this.m = m; this.p = p; this.r = r;
			this.uRowSupport = uR; this.uColSupport = uC;
			this.vRowSupport = vR; this.vColSupport = vC;
			this.wRowSupport = wR; this.wColSupport = wC;
		}

		public static SchemeSupports extract(NonCubicBilinearAlgorithm alg) {
			int n = alg.n, m = alg.m, p = alg.p, r = alg.r;
			int[][] uR = new int[r][];
			int[][] uC = new int[r][];
			int[][] vR = new int[r][];
			int[][] vC = new int[r][];
			int[][] wR = new int[r][];
			int[][] wC = new int[r][];
			for (int k = 0; k < r; k++) {
				uR[k] = supportRows(alg.denseU(), k, n, m);
				uC[k] = supportCols(alg.denseU(), k, n, m);
				vR[k] = supportRows(alg.denseV(), k, m, p);
				vC[k] = supportCols(alg.denseV(), k, m, p);
				wR[k] = supportRows(alg.denseW(), k, n, p);
				wC[k] = supportCols(alg.denseW(), k, n, p);
			}
			return new SchemeSupports(n, m, p, r, uR, uC, vR, vC, wR, wC);
		}

		private static int[] supportRows(double[][] factor, int rank, int rows, int cols) {
			boolean[] hit = new boolean[rows];
			int count = 0;
			for (int i = 0; i < rows; i++) {
				for (int j = 0; j < cols; j++) {
					if (factor[i * cols + j][rank] != 0.0) {
						if (!hit[i]) { hit[i] = true; count++; }
						break;
					}
				}
			}
			int[] out = new int[count];
			int idx = 0;
			for (int i = 0; i < rows; i++) if (hit[i]) out[idx++] = i;
			return out;
		}

		private static int[] supportCols(double[][] factor, int rank, int rows, int cols) {
			boolean[] hit = new boolean[cols];
			int count = 0;
			for (int j = 0; j < cols; j++) {
				for (int i = 0; i < rows; i++) {
					if (factor[i * cols + j][rank] != 0.0) {
						if (!hit[j]) { hit[j] = true; count++; }
						break;
					}
				}
			}
			int[] out = new int[count];
			int idx = 0;
			for (int j = 0; j < cols; j++) if (hit[j]) out[idx++] = j;
			return out;
		}
	}

	/**
	 * Compute the per-product sub-shape array for {@code supports} at the
	 * given allocation. Matches the semantics of
	 * {@link eu.solven.matmul.catalog.Recombination#recombineWithAllocation}
	 * with {@code peel = null}: each sub-dim is the {@code min} of the two
	 * relevant {@code max-over-support} views (U vs W for A, U vs V for B,
	 * V vs W for C).
	 *
	 * @return {@code shapes[k] = [subA, subB, subC]} for product {@code k}
	 */
	public static int[][] shapesAt(SchemeSupports supports, int[] allocA, int[] allocB, int[] allocC) {
		if (allocA.length != supports.n || allocB.length != supports.m || allocC.length != supports.p) {
			throw new IllegalArgumentException(String.format(
					"alloc lengths (%d,%d,%d) must equal scheme dims (%d,%d,%d)",
					allocA.length, allocB.length, allocC.length, supports.n, supports.m, supports.p));
		}
		int[][] shapes = new int[supports.r][3];
		for (int k = 0; k < supports.r; k++) {
			int subA_U = maxIndexed(allocA, supports.uRowSupport[k]);
			int subB_U = maxIndexed(allocB, supports.uColSupport[k]);
			int subB_V = maxIndexed(allocB, supports.vRowSupport[k]);
			int subC_V = maxIndexed(allocC, supports.vColSupport[k]);
			int subA_W = maxIndexed(allocA, supports.wRowSupport[k]);
			int subC_W = maxIndexed(allocC, supports.wColSupport[k]);
			shapes[k][0] = Math.min(subA_U, subA_W);
			shapes[k][1] = Math.min(subB_U, subB_V);
			shapes[k][2] = Math.min(subC_V, subC_W);
		}
		return shapes;
	}

	private static int maxIndexed(int[] alloc, int[] indices) {
		int max = 0;
		for (int i : indices) if (alloc[i] > max) max = alloc[i];
		return max;
	}

	/** Sum {@code sota.getRank(shape)} over a shape multiset. */
	public static long costOf(int[][] shapes, SotaResolver sota) {
		long total = 0;
		for (int[] s : shapes) total += sota.getRank(s[0], s[1], s[2]);
		return total;
	}

	/** A single mask evaluation result. */
	public static final class MaskCandidate {
		public final int mask;
		public final long cost;
		public final int[][] shapes;

		public MaskCandidate(int mask, long cost, int[][] shapes) {
			this.mask = mask;
			this.cost = cost;
			this.shapes = shapes;
		}

		@Override
		public String toString() {
			return String.format("mask=%d cost=%d shapes=%s",
					mask, cost, Arrays.deepToString(shapes));
		}
	}

	/**
	 * Score all 8 axis-flip mask variants of {@code canonical} at the given
	 * allocation. Returns the top-{@code k} cheapest, after de-duplication by
	 * shape multiset (two masks that produce the same multiset are merged,
	 * keeping the smaller mask label for traceability).
	 *
	 * <p>Equivalence: mask bit 0 = reverse A axis, bit 1 = reverse B, bit 2 =
	 * reverse C. Reversing the allocation is equivalent to applying the
	 * matching axis-flip transform to the canonical scheme — see class doc.
	 */
	public static List<MaskCandidate> topKMasks(NonCubicBilinearAlgorithm canonical,
			int[] allocA, int[] allocB, int[] allocC, SotaResolver sota, int k) {
		return topKMasks(SchemeSupports.extract(canonical), allocA, allocB, allocC, sota, k);
	}

	/**
	 * Variant that accepts pre-extracted supports — avoids re-scanning the
	 * canonical scheme on every call. Use this when iterating many
	 * allocations of the same base; extract once per base, reuse per
	 * allocation. Cuts per-allocation cost from O(r·dim²) extract +
	 * 8·O(r·dim) score to just 8·O(r·dim) score.
	 */
	public static List<MaskCandidate> topKMasks(SchemeSupports supports,
			int[] allocA, int[] allocB, int[] allocC, SotaResolver sota, int k) {
		if (k <= 0) {
			throw new IllegalArgumentException("k must be > 0, got " + k);
		}
		Map<Long, MaskCandidate> byMultiset = new LinkedHashMap<>();
		for (int mask = 0; mask < 8; mask++) {
			int[] aA = ((mask & 1) != 0) ? reverse(allocA) : allocA;
			int[] aB = ((mask & 2) != 0) ? reverse(allocB) : allocB;
			int[] aC = ((mask & 4) != 0) ? reverse(allocC) : allocC;
			int[][] shapes = shapesAt(supports, aA, aB, aC);
			long cost = costOf(shapes, sota);
			MaskCandidate cand = new MaskCandidate(mask, cost, shapes);
			long key = multisetHash(shapes);
			MaskCandidate prev = byMultiset.get(key);
			if (prev == null || cand.cost < prev.cost
					|| (cand.cost == prev.cost && cand.mask < prev.mask)) {
				byMultiset.put(key, cand);
			}
		}
		List<MaskCandidate> deduped = new ArrayList<>(byMultiset.values());
		deduped.sort(Comparator.<MaskCandidate>comparingLong(c -> c.cost)
				.thenComparingInt(c -> c.mask));
		if (deduped.size() <= k) return deduped;
		return new ArrayList<>(deduped.subList(0, k));
	}

	private static int[] reverse(int[] a) {
		int n = a.length;
		int[] r = new int[n];
		for (int i = 0; i < n; i++) r[i] = a[n - 1 - i];
		return r;
	}

	/**
	 * Order-invariant 64-bit hash of the multiset of sub-shape triples.
	 * Uses addition of well-mixed per-element hashes — commutative and
	 * associative (independent of order) AND multiplicity-aware (two
	 * copies count twice, unlike XOR which would cancel them). O(r) with
	 * no allocation and no sort, vs the previous O(r log r) sort-then-FNV.
	 *
	 * <p>The per-element mix is the Murmur3-64 finalizer (high-quality
	 * avalanche); collisions across distinct multisets are statistically
	 * negligible at search sizes &lt; 1M entries.
	 */
	private static long multisetHash(int[][] shapes) {
		long h = 0;
		int n = shapes.length;
		for (int i = 0; i < n; i++) {
			int[] x = shapes[i];
			long packed = ((long) x[0] << 32) | ((long) x[1] << 16) | (long) x[2];
			h += murmur3Finalizer(packed);
		}
		return h;
	}

	private static long murmur3Finalizer(long x) {
		x ^= x >>> 33;
		x *= 0xff51afd7ed558ccdL;
		x ^= x >>> 33;
		x *= 0xc4ceb9fe1a85ec53L;
		x ^= x >>> 33;
		return x;
	}

	/**
	 * Canonical string key for a shape multiset (sort by lex order on tuples).
	 * Two mask variants that produce the same multiset get the same key, so the
	 * scorer collapses them rather than reporting redundant candidates.
	 *
	 * @deprecated kept for tests; use {@link #multisetHash} in hot paths.
	 */
	@Deprecated
	static String multisetKey(int[][] shapes) {
		String[] s = new String[shapes.length];
		for (int i = 0; i < shapes.length; i++) {
			int[] x = shapes[i];
			s[i] = x[0] + "," + x[1] + "," + x[2];
		}
		Arrays.sort(s);
		return String.join("|", s);
	}
}
