package eu.solven.matmul.search;

import eu.solven.matmul.recombination.BlockSplitSearch;

import eu.solven.matmul.NonCubicBilinearAlgorithm;

import eu.solven.matmul.recombination.Recombination;

import java.util.HashMap;
import java.util.List;
import java.util.Map;
import java.util.Optional;

/**
 * SOTA resolver that performs the DIS09-style <strong>recursive
 * closure</strong>: for each requested shape {@code ⟨a,b,c⟩}, returns
 * the min over (i) direct catalog lookup, (ii) any closed-form bound
 * via {@link CitedBound}, (iii) the best <em>composition</em>
 * found by {@link BlockSplitSearch#findBestStrategy} — where the
 * sub-problem ranks consulted by the composition come from this very
 * resolver, recursively.
 *
 * <p>This is the "fill the {@code T[a,b,c]} table bottom-up" pattern
 * from DIS09 §3.4: each cell is the minimum across all recipes that
 * read smaller cells, with memoisation preventing recomputation.</p>
 *
 * <p>Termination: the search shrinks shape monotonically (composition
 * sub-shapes are strictly smaller in at least one axis), so the
 * recursion bottoms out at small shapes that hit the catalog
 * directly. Memoised in {@link #cache}; bound by recursion depth on
 * the maximum axis (effectively cubic worst-case in {@code MAX_DIM}).</p>
 */
public final class RecursiveClosureSota implements Recombination.SotaResolver {

	private final Recombination.AlgorithmLookup lookup;
	private final List<BlockSplitSearch.NamedBase> pool;
	private final boolean balancedOnly;
	private final boolean useFormulaBounds;
	private final Map<Long, Integer> cache = new HashMap<>();
	// Guard against re-entrant cycles (shouldn't happen — shape strictly decreases —
	// but defensive).
	private final java.util.Set<Long> inProgress = new java.util.HashSet<>();

	public RecursiveClosureSota(Recombination.AlgorithmLookup lookup,
			List<BlockSplitSearch.NamedBase> pool,
			boolean balancedOnly,
			boolean useFormulaBounds) {
		this.lookup = lookup;
		this.pool = pool;
		this.balancedOnly = balancedOnly;
		this.useFormulaBounds = useFormulaBounds;
	}

	@Override
	public int getRank(int a, int b, int c) {
		if (a <= 0 || b <= 0 || c <= 0) return 0;
		long key = key(a, b, c);
		Integer hit = cache.get(key);
		if (hit != null) return hit;
		if (inProgress.contains(key)) {
			// Shouldn't happen with monotonically-shrinking shape, but if it does,
			// bail out to direct lookup to break the cycle.
			return directRank(a, b, c);
		}
		inProgress.add(key);
		try {
			int best = computeBest(a, b, c);
			cache.put(key, best);
			return best;
		} finally {
			inProgress.remove(key);
		}
	}

	private int computeBest(int a, int b, int c) {
		// (1) Direct catalog rank — never invalid; floor of a*b*c (naive).
		int best = directRank(a, b, c);
		if (best <= a + b + c) {
			// Already at a trivial lower bound; no point searching further.
			return best;
		}

		// (2) Closed-form bounds (Pan TA, HK, Waksman, …) if enabled.
		if (useFormulaBounds) {
			int formulaBest = new CitedBound(lookup, false).getRank(a, b, c);
			if (formulaBest > 0 && formulaBest < best) best = formulaBest;
		}

		// (3) Recursive composition: try both recombination over pool and concat
		//     splits, but consult THIS resolver for the sub-shape ranks.
		Optional<BlockSplitSearch.NonCubicStrategy> picked =
				BlockSplitSearch.findBestStrategy(a, b, c, pool, this, balancedOnly);
		if (picked.isPresent() && picked.get().rank() < best) {
			best = (int) Math.min(picked.get().rank(), Recombination.SotaResolver.UNKNOWN_RANK);
		}
		return best;
	}

	private int directRank(int a, int b, int c) {
		// findRank (NOT find): find() skips lineage-only STUBS (maxDim>16), so the entire
		// 17–32 band reads as NAIVE a*b*c here — making sota blind to every stub's real
		// rank and systematically undervaluing concat/composition decompositions that use
		// them (e.g. ⟨25,26,32⟩ never valued the 6+26 concat 2402+8552=10954). findRank sees
		// stubs and already excludes corrupted ones; the write-guard keeps stubs replayable.
		int r = lookup.findRank(a, b, c);
		int naive = a * b * c;
		return (r > 0 && r < naive) ? r : naive;
	}

	private static long key(int a, int b, int c) {
		// Canonical sorted key; rank is invariant under cyclic + transpose symmetry.
		int[] s = { a, b, c };
		java.util.Arrays.sort(s);
		return ((long) s[0] << 32) | ((long) s[1] << 16) | (long) s[2];
	}

	/** Visible for debugging: how many shapes have been resolved into cache. */
	public int cacheSize() { return cache.size(); }
}
