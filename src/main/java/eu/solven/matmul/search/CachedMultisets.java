package eu.solven.matmul.search;

import java.util.Arrays;
import java.util.IdentityHashMap;
import java.util.Map;
import java.util.concurrent.atomic.AtomicLong;

import eu.solven.matmul.NonCubicBilinearAlgorithm;
import eu.solven.matmul.catalog.Recombination;

/**
 * Caches per-(scheme, allocation, peel) shape multisets for use by
 * recombination-driven search loops (FrontierClosure, BlockSplitSearch,
 * SchemeSweep).
 *
 * <p><b>Why this exists</b>: {@link Recombination#recombineWithAllocation}
 * re-computes the per-product effective sub-shapes on every call (~O(r ·
 * dim²) per call for support extraction + max/min reductions). When the
 * same (scheme, allocation, peel) tuple is queried many times — typical
 * in a frontier sweep where catalog SOTAs evolve but the scheme/allocation
 * inventory does not — caching the resulting shape multiset gives a clean
 * ~r-x speedup. The catalog-SOTA-dependent cost is recomputed per query
 * (cheap lookups), but the structural shape multiset is reused.
 *
 * <p><b>What's cacheable</b>: the {@code int[][]} {@code smallMatrixSizes}
 * from {@link Recombination.Result} — purely structural, depends only on
 * (scheme block-supports, allocation, peel). The total rank is recomputed
 * by summing {@code sota.getRank(shape)} over the cached multiset.
 *
 * <p><b>Cache key</b>: identity-based on the scheme reference (cheap), plus
 * value-based on the allocation and peel arrays. Same scheme object +
 * same allocation arrays + same peel arrays → cache hit.
 *
 * <p><b>Next step beyond this</b> (per user's idea): cache symbolic
 * multiset templates per (scheme, allocation-pattern-type, peel-pattern)
 * such that one cache entry covers all concrete (a, b, c, …) values at
 * the same structural pattern. This class is the simpler concrete-key
 * version; the symbolic-template layer can wrap it later.
 */
public final class CachedMultisets {

	/** Cap per-scheme entries to avoid OOM during long sweeps. Sized so that
	 *  ~5 schemes × ~200k entries × ~100 bytes/entry = ~100 MB heap. */
	private static final int MAX_ENTRIES_PER_SCHEME = 200_000;

	/** Per-scheme cache. The inner map is keyed by (alloc + peel) value-equality. */
	private final Map<NonCubicBilinearAlgorithm, Map<Key, int[][]>> bySchemeIdentity =
			new IdentityHashMap<>();
	private final AtomicLong hits = new AtomicLong();
	private final AtomicLong misses = new AtomicLong();
	private final AtomicLong evicted = new AtomicLong();

	/** Per-scheme cache of the extracted block-supports — reused across
	 *  all (allocation, peel) calls for that scheme. */
	private final Map<NonCubicBilinearAlgorithm, AnalyticalMaskSearch.SchemeSupports> supportsCache =
			new IdentityHashMap<>();

	/**
	 * Look up (or compute and cache) the shape multiset for
	 * {@code Recombination.recombineWithAllocation(scheme, ..., allocA,
	 * allocB, allocC, peelA, peelB, peelC)}. Returns a flat
	 * {@code int[r][3]} mirroring {@link Recombination.Result#smallMatrixSizes}.
	 *
	 * <p>The {@code sota} parameter is only used on a cache miss to drive
	 * the {@link Recombination#recombineWithAllocation} call. Cache hits
	 * skip the sota entirely. Callers do their own per-query summing.
	 */
	public int[][] getMultiset(NonCubicBilinearAlgorithm scheme,
			int[] allocA, int[] allocB, int[] allocC,
			int[] peelA, int[] peelB, int[] peelC,
			Recombination.SotaResolver sota) {
		Map<Key, int[][]> perAlloc;
		synchronized (this) {
			perAlloc = bySchemeIdentity.computeIfAbsent(scheme, s -> new java.util.HashMap<>());
		}
		Key k = new Key(allocA, allocB, allocC, peelA, peelB, peelC);
		int[][] cached;
		synchronized (perAlloc) {
			cached = perAlloc.get(k);
		}
		if (cached != null) {
			hits.incrementAndGet();
			return cached;
		}
		misses.incrementAndGet();
		Recombination.Result r = (peelA == null && peelB == null && peelC == null)
				? Recombination.recombineWithAllocation(scheme, sota, allocA, allocB, allocC)
				: Recombination.recombineWithAllocation(scheme, sota, allocA, allocB, allocC,
						peelA, peelB, peelC);
		synchronized (perAlloc) {
			if (perAlloc.size() < MAX_ENTRIES_PER_SCHEME) {
				perAlloc.put(k, r.smallMatrixSizes);
			} else {
				evicted.incrementAndGet();
			}
		}
		return r.smallMatrixSizes;
	}

	/** Compute total rank from a cached multiset using the given sota. */
	public static long costOf(int[][] multiset, Recombination.SotaResolver sota) {
		long total = 0;
		for (int[] s : multiset) {
			total += sota.getRank(s[0], s[1], s[2]);
		}
		return total;
	}

	public long hits() { return hits.get(); }
	public long misses() { return misses.get(); }

	public String stats() {
		long h = hits.get(), m = misses.get(), e = evicted.get();
		long total = h + m;
		double rate = total == 0 ? 0.0 : 100.0 * h / total;
		return String.format("CachedMultisets: %d hits, %d misses, %d not-cached-after-cap (%.1f%% hit rate)",
				h, m, e, rate);
	}

	/** Value-based key on allocation + peel arrays. */
	private static final class Key {
		final int[] aA, aB, aC, pA, pB, pC;
		final int hash;
		Key(int[] aA, int[] aB, int[] aC, int[] pA, int[] pB, int[] pC) {
			this.aA = aA.clone();
			this.aB = aB.clone();
			this.aC = aC.clone();
			this.pA = pA == null ? null : pA.clone();
			this.pB = pB == null ? null : pB.clone();
			this.pC = pC == null ? null : pC.clone();
			int h = Arrays.hashCode(this.aA);
			h = 31 * h + Arrays.hashCode(this.aB);
			h = 31 * h + Arrays.hashCode(this.aC);
			h = 31 * h + Arrays.hashCode(this.pA);
			h = 31 * h + Arrays.hashCode(this.pB);
			h = 31 * h + Arrays.hashCode(this.pC);
			this.hash = h;
		}
		@Override public int hashCode() { return hash; }
		@Override public boolean equals(Object o) {
			if (!(o instanceof Key k)) return false;
			return Arrays.equals(aA, k.aA) && Arrays.equals(aB, k.aB) && Arrays.equals(aC, k.aC)
					&& Arrays.equals(pA, k.pA) && Arrays.equals(pB, k.pB) && Arrays.equals(pC, k.pC);
		}
	}
}
