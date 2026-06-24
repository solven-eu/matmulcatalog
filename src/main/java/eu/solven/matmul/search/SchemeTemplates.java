package eu.solven.matmul.search;

import java.util.IdentityHashMap;
import java.util.Map;
import java.util.concurrent.atomic.AtomicLong;

import eu.solven.matmul.NonCubicBilinearAlgorithm;
import eu.solven.matmul.recombination.Recombination;

/**
 * Cache of per-scheme "symbolic templates" — the per-product block-support
 * bitmasks for U, V, W factors. Once extracted, the per-call cost of
 * computing sub-product shapes drops from O(r · dim²) (re-scanning all
 * factor entries) to O(r · dim) (just looking up cached bitmasks).
 *
 * <p>This is the symbolic-template layer per task #125: cache scope is
 * per (scheme identity). Allocation and peel are NOT part of the cache
 * key — they're query-time parameters. The templates capture the
 * intrinsic structure of the scheme; concrete (a, b, c) values
 * instantiate the shapes at query time.
 *
 * <p>Subsumed work from {@link Recombination#processAdditions}: rather
 * than scanning the factor matrix each call, we scan ONCE during
 * template construction and store 6 bitmasks per product (one per
 * U-row, U-col, V-row, V-col, W-row, W-col). Each bitmask encodes
 * "which block-row/col indices have a nonzero entry in this product's
 * U/V/W column".
 */
public final class SchemeTemplates {

	private static final Map<NonCubicBilinearAlgorithm, SchemeTemplates> CACHE = new IdentityHashMap<>();
	private static final AtomicLong cacheHits = new AtomicLong();
	private static final AtomicLong cacheMisses = new AtomicLong();

	/** Lookup-or-build. Identity-keyed so the same scheme reference always hits. */
	public static SchemeTemplates forScheme(NonCubicBilinearAlgorithm scheme) {
		synchronized (CACHE) {
			SchemeTemplates t = CACHE.get(scheme);
			if (t != null) { cacheHits.incrementAndGet(); return t; }
			cacheMisses.incrementAndGet();
			t = new SchemeTemplates(scheme);
			CACHE.put(scheme, t);
			return t;
		}
	}

	public static String stats() {
		long h = cacheHits.get(), m = cacheMisses.get();
		long total = h + m;
		double rate = total == 0 ? 0.0 : 100.0 * h / total;
		return String.format("SchemeTemplates: %d hits, %d misses (%.1f%% hit rate, %d schemes cached)",
				h, m, rate, m);
	}

	// ── Per-scheme extracted bitmasks ──

	private final int n, m, p, r;
	/** Per-product U-row block-support bitmask. Bit i = 1 iff U[i*m+j][k] ≠ 0 for some j. */
	private final int[] uRowMask;
	private final int[] uColMask;
	private final int[] vRowMask;
	private final int[] vColMask;
	private final int[] wRowMask;
	private final int[] wColMask;

	private SchemeTemplates(NonCubicBilinearAlgorithm scheme) {
		this.n = scheme.n;
		this.m = scheme.m;
		this.p = scheme.p;
		this.r = scheme.r;
		this.uRowMask = new int[r];
		this.uColMask = new int[r];
		this.vRowMask = new int[r];
		this.vColMask = new int[r];
		this.wRowMask = new int[r];
		this.wColMask = new int[r];
		extractAxisMasks(scheme.denseU(), n, m, uRowMask, uColMask);
		extractAxisMasks(scheme.denseV(), m, p, vRowMask, vColMask);
		extractAxisMasks(scheme.denseW(), n, p, wRowMask, wColMask);
	}

	/** For each product k, scan factor[i*cols+j][k] and OR i into rowMask[k], j into colMask[k] for nonzeros. */
	private static void extractAxisMasks(double[][] factor, int rows, int cols,
			int[] rowMask, int[] colMask) {
		int r = rowMask.length;
		for (int k = 0; k < r; k++) {
			int rm = 0, cm = 0;
			for (int i = 0; i < rows; i++) {
				for (int j = 0; j < cols; j++) {
					if (factor[i * cols + j][k] != 0.0) {
						rm |= (1 << i);
						cm |= (1 << j);
					}
				}
			}
			rowMask[k] = rm;
			colMask[k] = cm;
		}
	}

	private static int maxOverMask(int[] alloc, int mask) {
		int max = 0;
		int idx = 0;
		while (mask != 0) {
			if ((mask & 1) != 0 && alloc[idx] > max) max = alloc[idx];
			mask >>>= 1;
			idx++;
		}
		return max;
	}

	/**
	 * Instantiate the templates at concrete (alloc, peel) values, returning total rank
	 * = Σ sota.getRank(sub-shape) for each of {@code r} sub-products.
	 *
	 * <p>Effective per-axis allocations are {@code effA[i] = allocA[i] - (peelA == null ? 0 : peelA[i])},
	 * mirroring {@link Recombination#recombineWithAllocation}.
	 */
	public long totalRank(int[] allocA, int[] allocB, int[] allocC,
			int[] peelA, int[] peelB, int[] peelC,
			Recombination.SotaResolver sota) {
		int[] effA = applyPeel(allocA, peelA);
		int[] effB = applyPeel(allocB, peelB);
		int[] effC = applyPeel(allocC, peelC);
		long total = 0;
		for (int k = 0; k < r; k++) {
			int uRow = maxOverMask(effA, uRowMask[k]);
			int uCol = maxOverMask(effB, uColMask[k]);
			int vRow = maxOverMask(effB, vRowMask[k]);
			int vCol = maxOverMask(effC, vColMask[k]);
			int wRow = maxOverMask(effA, wRowMask[k]);
			int wCol = maxOverMask(effC, wColMask[k]);
			int subA = Math.min(uRow, wRow);
			int subB = Math.min(uCol, vRow);
			int subC = Math.min(vCol, wCol);
			total += sota.getRank(subA, subB, subC);
		}
		return total;
	}

	private static int[] applyPeel(int[] alloc, int[] peel) {
		if (peel == null) return alloc;
		int[] eff = new int[alloc.length];
		for (int i = 0; i < alloc.length; i++) eff[i] = alloc[i] - peel[i];
		return eff;
	}
}
