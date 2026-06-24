package eu.solven.matmul.search;

import eu.solven.matmul.recombination.Recombination;

import eu.solven.matmul.papers.hopcroftkerr1971.HopcroftKerrBound;

import eu.solven.matmul.papers.waksman1970.WaksmanBound;

import eu.solven.matmul.papers.dis2009.PanTrilinearAggregation;

import java.util.Optional;

import eu.solven.matmul.NonCubicBilinearAlgorithm;

/**
 * SOTA resolver that augments a catalog-backed lookup with
 * <strong>closed-form rank bounds</strong> from the literature:
 * Pan trilinear-aggregation (NC cubic), Hopcroft-Kerr {@code ⟨a,2,c⟩}
 * family, Waksman commutative cubic. For each query ⟨a,b,c⟩ returns
 * the {@code min} over:
 *
 * <ol>
 *   <li>the catalog-best (via the wrapped lookup),</li>
 *   <li>Pan TA cubic formula (if {@code a=b=c}),</li>
 *   <li>Hopcroft-Kerr formula (if any axis is 2),</li>
 *   <li>(optional) Waksman commutative — only when the caller flags
 *       the target as commutative-acceptable.</li>
 * </ol>
 *
 * <p>This lets {@link Recombination#findBestMultiBaseSplit}-style
 * searches discover better outer allocations even when the catalog
 * doesn't carry a constructive scheme for the bound. The
 * MATERIALIZATION path is separate ({@link Recombination#constructFromResult})
 * and still requires a real scheme via the wrapped lookup — so a
 * formula bound that's lower than the catalog's best scheme will
 * IMPROVE the rank-only search prediction but will NOT change the
 * materialised output until the corresponding scheme is added.</p>
 *
 * <p>Field discipline: Pan TA and Hopcroft-Kerr are non-commutative
 * and apply to R/Q/Z/C/F₂. Waksman is gated on
 * {@code includeCommutative=true}.</p>
 */
public final class CitedBound implements Recombination.SotaResolver {

	private final Recombination.AlgorithmLookup lookup;
	private final boolean includeCommutative;

	public CitedBound(Recombination.AlgorithmLookup lookup, boolean includeCommutative) {
		this.lookup = lookup;
		this.includeCommutative = includeCommutative;
	}

	/** Default: NC bounds only (Pan TA + Hopcroft-Kerr); no commutative fallback. */
	public CitedBound(Recombination.AlgorithmLookup lookup) {
		this(lookup, false);
	}

	@Override
	public int getRank(int a, int b, int c) {
		if (a == 0 || b == 0 || c == 0) return 0;

		// (1) Catalog-best rank. Use findRank (NOT find): findRank canonicalises the
		// orientation (sorts the axes — rank is S₃-invariant, so ⟨6,3,3⟩ resolves to
		// the catalog's ⟨3,3,6⟩=40) AND sees stub schemes, whereas find() requires the
		// exact orientation and skips stubs. Using find() here silently returned the
		// naive a·b·c for any non-sorted factor, which made the Kronecker pass miss
		// e.g. ⟨12,12,15⟩ = ⟨2,4,5⟩⊗⟨6,3,3⟩ (it scored ⟨6,3,3⟩ as 54 not 40). (#201)
		int catRank = lookup.findRank(a, b, c);
		int best = (catRank > 0 && catRank < Recombination.SotaResolver.UNKNOWN_RANK) ? catRank : a * b * c;

		// (2) Pan TA (non-commutative, cubic square only).
		if (a == b && b == c) {
			long taBound = PanTrilinearAggregation.cubicBound(a);
			if (taBound < best && taBound > 0) best = (int) taBound;
		}

		// (3) Hopcroft-Kerr ⟨a,2,c⟩ family (NC, applies when one axis is 2).
		long hkBound = HopcroftKerrBound.forShape(a, b, c);
		if (hkBound > 0 && hkBound < best) best = (int) hkBound;

		// (4) Waksman commutative cubic (only if caller opted in).
		if (includeCommutative) {
			long wak = WaksmanBound.forShape(a, b, c);
			if (wak > 0 && wak < best) best = (int) wak;
		}

		return best;
	}
}
