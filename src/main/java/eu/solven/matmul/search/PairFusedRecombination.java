package eu.solven.matmul.search;

import eu.solven.matmul.catalog.Recombination;

import eu.solven.matmul.papers.pan1978.PanPairProduct;
import eu.solven.matmul.papers.dis2009.PanTrilinearAggregation;

import java.util.Optional;

import eu.solven.matmul.NonCubicBilinearAlgorithm;

/**
 * Predicts the rank of a Pan-pair-fused recombination over a Strassen-like
 * base whose outer products all share the same cubic sub-shape
 * {@code ⟨k,k,k⟩} (i.e. the standard balanced {@code [k,k]³} allocation
 * over an outer Strassen-style base).
 *
 * <p>Pair-fusion model: of the {@code N} outer products, pair {@code ⌊N/2⌋}
 * of them (each pair costing {@link PanPairProduct#rank} =
 * {@code k³+3k²}) and leave 1 solo (cost {@code R(⟨k,k,k⟩)}) if N is odd,
 * else 0 solos. Strassen {@code N=7} ⇒ 3 pairs + 1 solo — what
 * {@link eu.solven.matmul.research.MaterializeViaPanPair} hard-codes.</p>
 *
 * <p>Limited to:</p>
 * <ul>
 *   <li>Outer base must produce uniform cubic sub-shapes (e.g. Strassen
 *       with {@code [k,k]³} allocation).</li>
 *   <li>Target must be cubic {@code ⟨n,n,n⟩} with {@code n = 2k} (Strassen)
 *       or {@code n = 3k} (Laderman) etc.</li>
 *   <li>Pair-fusion is profitable iff
 *       {@code 2·R(⟨k,k,k⟩) > k³+3k²} (the user's "MaterializeViaPanPair
 *       profitability table").</li>
 * </ul>
 *
 * <p>For non-uniform sub-shapes (irregular allocations, non-cubic
 * targets) the cyclic-pair compatibility {@code ⟨a,b,c⟩+⟨b,c,a⟩} needs
 * dedicated handling — this MVP only covers the all-cubic case.</p>
 *
 * <p>{@link #predict} weighs <em>only</em> SOLO vs PAIR_FUSED. To compare
 * <em>all</em> trilinear-aggregation options — including the full
 * single-product TA ({@link PanTrilinearAggregation}) — use
 * {@link #chooseBest}, which returns the cheapest {@link TaStrategy} and
 * the per-option cost breakdown. See {@link TaStrategy} for the
 * pair-vs-full-TA crossover (full TA wins only for even {@code k ≥ 18}).</p>
 */
public final class PairFusedRecombination {

	private PairFusedRecombination() {}

	public record Prediction(
			int targetN, int outerProducts, int subK,
			long pairCost, long soloCost, int pairs, int solos, long totalRank) {
		public String breakdown() {
			return String.format(
					"pair-fused[%d pairs × ⟨%d,%d,%d⟩=%d + %d solo × ⟨%d,%d,%d⟩=%d] = %d",
					pairs, subK, subK, subK, pairCost, solos, subK, subK, subK, soloCost, totalRank);
		}
	}

	/**
	 * Predict the rank of pair-fused recombination for cubic target
	 * {@code ⟨n,n,n⟩ = ⟨k,k,k⟩^N} where the outer base has {@code N}
	 * products. Returns empty when not profitable (pair-cost ≥ 2·solo)
	 * or when N is too small (≤ 1).
	 */
	public static Optional<Prediction> predict(int targetN, int subK, int outerProducts,
			Recombination.SotaResolver sota) {
		if (outerProducts < 2 || subK < 1) return Optional.empty();
		if (subK * (long) outerProducts == 0) return Optional.empty();
		long soloRank = sota.getRank(subK, subK, subK);
		if (soloRank <= 0) return Optional.empty();
		long pairRank = PanPairProduct.rank(subK, subK, subK);
		// Pairing each pair: 1 pair = 2 products at total rank pairRank;
		// not pairing: 2 products at total 2*soloRank. Profitable iff
		// pairRank < 2*soloRank.
		if (pairRank >= 2L * soloRank) return Optional.empty();
		int pairs = outerProducts / 2;
		int solos = outerProducts % 2;
		long total = pairs * pairRank + solos * soloRank;
		return Optional.of(new Prediction(targetN, outerProducts, subK,
				pairRank, soloRank, pairs, solos, total));
	}

	/**
	 * Convenience: predict pair-fusion for a balanced {@code Base ⊗ [k,k]³}
	 * Strassen-style recombination. Calls {@link #predict} with the base's
	 * own product count.
	 */
	public static Optional<Prediction> predictBalancedCubic(int targetN,
			NonCubicBilinearAlgorithm outerBase, Recombination.SotaResolver sota) {
		if (outerBase.n != outerBase.m || outerBase.m != outerBase.p) return Optional.empty();
		int baseDim = outerBase.n;
		if (targetN % baseDim != 0) return Optional.empty();
		int subK = targetN / baseDim;
		return predict(targetN, subK, outerBase.r, sota);
	}

	// ------------------------------------------------------------------
	// Multiple trilinear-aggregation (TA) options
	// ------------------------------------------------------------------

	/**
	 * The trilinear-aggregation strategies the catalog can apply to a group
	 * of {@code N} same-shape cubic leaves {@code ⟨k,k,k⟩}.
	 *
	 * <ul>
	 *   <li>{@link #SOLO} — no aggregation; each leaf costs its catalog rank
	 *       {@code R(⟨k,k,k⟩)}.</li>
	 *   <li>{@link #PAIR_FUSED} — Pan's <em>two-product</em> cross-fusion
	 *       ({@link PanPairProduct}): a cyclic pair {@code ⟨k,k,k⟩+⟨k,k,k⟩}
	 *       costs {@code k³+3k²} (per product {@code ≈ k³/2}). This is the
	 *       only TA that fuses <em>separate</em> products.</li>
	 *   <li>{@link #FULL_TA} — Pan/Islam <em>single-product</em> TA
	 *       ({@link PanTrilinearAggregation#bestPanTaBound}): each leaf is
	 *       recomputed internally aggregating its scalar products in triples
	 *       (per product {@code ≈ k³/3}).</li>
	 * </ul>
	 *
	 * <p><strong>There is no "fuse three separate products" primitive.</strong>
	 * The {@code /3} of {@code FULL_TA} comes from a single product's own
	 * 3-fold cyclic index symmetry, not from fusing three matmuls. So the two
	 * genuine choices when several same-shape leaves arrive are PAIR_FUSED
	 * (shallow, cross-product, low {@code O(k²)} overhead) versus FULL_TA
	 * (deep, per-product, {@code ≈4×} the {@code O(k²)} overhead but a better
	 * leading coefficient).</p>
	 *
	 * <p><strong>Crossover.</strong> {@code FULL_TA}'s {@code ≈k³/3} only
	 * overcomes its heavier corrections at large leaves: it beats
	 * {@code PAIR_FUSED} only for <em>even</em> {@code k ≥ 18} (odd {@code k}
	 * carries Islam's heavier {@code 15k²} term, so pairs stay ahead in any
	 * practical range). For small leaves — e.g. the {@code ⟨7,7,7⟩} leaves of
	 * a {@code ⟨14,14,14⟩} or padded-{@code ⟨17,17,17⟩} recombination —
	 * {@code PAIR_FUSED} (245/product) dominates both {@code SOLO} (249) and
	 * {@code FULL_TA} (390). This is the algebraic reason the catalog
	 * pair-fuses rather than full-TAs its leaves.</p>
	 */
	public enum TaStrategy {
		SOLO, PAIR_FUSED, FULL_TA
	}

	/**
	 * The cost of each TA option for {@code outerProducts} leaves of cubic
	 * shape {@code ⟨subK,subK,subK⟩}, plus which one wins. All totals are
	 * <em>upper-bound cost models</em> (not proven-optimal); {@code FULL_TA}
	 * and {@code PAIR_FUSED} are materialisable via
	 * {@link PanTrilinearAggregation#build} / {@link PanPairProduct#build}
	 * respectively, while {@code SOLO} defers to the catalog rank.
	 */
	public record TaChoice(
			int subK, int outerProducts,
			long soloTotal, long pairTotal, long fullTaTotal,
			TaStrategy best, long bestTotal) {
		public String breakdown() {
			return String.format(
					"TA options for %d×⟨%d,%d,%d⟩: solo=%d, pair-fused=%d, full-TA=%d → best=%s (%d)",
					outerProducts, subK, subK, subK,
					soloTotal, pairTotal, fullTaTotal, best, bestTotal);
		}
	}

	/**
	 * Compare all available trilinear-aggregation options ({@link TaStrategy})
	 * for a group of {@code outerProducts} cubic leaves {@code ⟨subK,subK,subK⟩}
	 * and return the cheapest. Unlike {@link #predict}, which only weighs
	 * SOLO against PAIR_FUSED, this also weighs the full single-product TA.
	 *
	 * <p>The per-leaf solo cost is the catalog rank {@code R(⟨k,k,k⟩)}.
	 * Pairing is only applied when a pair is cheaper than two solos
	 * ({@code k³+3k² < 2·R}); otherwise PAIR_FUSED falls back to all-solo.
	 * FULL_TA replaces <em>every</em> leaf by {@link PanTrilinearAggregation#bestPanTaBound}.</p>
	 *
	 * @return empty when inputs are degenerate or the catalog has no rank for
	 *         {@code ⟨subK,subK,subK⟩}
	 */
	public static Optional<TaChoice> chooseBest(int subK, int outerProducts,
			Recombination.SotaResolver sota) {
		if (outerProducts < 1 || subK < 1) return Optional.empty();
		long soloRank = sota.getRank(subK, subK, subK);
		if (soloRank <= 0) return Optional.empty();

		long soloTotal = (long) outerProducts * soloRank;

		long pairRank = PanPairProduct.rank(subK, subK, subK);
		long pairTotal;
		if (pairRank < 2L * soloRank) {
			int pairs = outerProducts / 2;
			int solos = outerProducts % 2;
			pairTotal = pairs * pairRank + solos * soloRank;
		} else {
			// pairing never helps here — equivalent to all-solo
			pairTotal = soloTotal;
		}

		long fullTaBound = PanTrilinearAggregation.bestPanTaBound(subK);
		long fullTaTotal = fullTaBound > 0
				? (long) outerProducts * fullTaBound
				: Long.MAX_VALUE;

		TaStrategy best = TaStrategy.SOLO;
		long bestTotal = soloTotal;
		if (pairTotal < bestTotal) {
			best = TaStrategy.PAIR_FUSED;
			bestTotal = pairTotal;
		}
		if (fullTaTotal < bestTotal) {
			best = TaStrategy.FULL_TA;
			bestTotal = fullTaTotal;
		}
		return Optional.of(new TaChoice(subK, outerProducts,
				soloTotal, pairTotal, fullTaTotal, best, bestTotal));
	}
}
