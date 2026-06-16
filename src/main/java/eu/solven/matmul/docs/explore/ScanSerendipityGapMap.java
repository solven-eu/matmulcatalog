package eu.solven.matmul.docs.explore;

import eu.solven.matmul.catalog.Recombination;

import java.util.ArrayList;
import java.util.Comparator;
import java.util.List;

import eu.solven.matmul.algebra.Field;
import eu.solven.matmul.catalog.FieldAwareLookup;
import lombok.extern.slf4j.Slf4j;

/**
 * Serendipity GAP MAP over the catalog (default {@code Field.Q}, shapes up to
 * dim 32). For every known shape {@code I = ⟨n,m,p⟩} and every bud type/size
 * {@code k}, prices the single-bud fusion saving
 *
 * <pre>
 *   σ_ax(k) = k·R(I) − R(k-enlarged along ax)
 * </pre>
 *
 * (U enlarges p, V enlarges n, W enlarges m — the
 * {@code SerendipitousBudProduct.costOf} convention), purely from
 * {@link FieldAwareLookup#findRank} — NO scheme generation. This is the cheap
 * "where is the catalog rank sub-additive, hence where could a bud win big"
 * potential layer (companion to {@code references/BUD_STRUCTURE_THEORY.md} and
 * {@code SERENDIPITY_RANK_TRADEOFF.md}).
 *
 * <p><strong>Potential, not realizability.</strong> A big σ means a base with
 * such a bud <em>would</em> save that much — it does NOT say a low-rank base
 * carrying that bud exists (the ⟨6,8,9⟩/⟨3,2,3⟩ rigidity wall). Feed the top
 * rows to {@code ScanSerendipityFrontier} / {@code StructuredWAls} for the
 * (expensive) realizability decision.</p>
 *
 * <p><strong>Catalog-bounded.</strong> σ is only computed where BOTH R(I) and
 * R(k-enlarged) are in the catalog. Enlarging the large axis of a big shape
 * usually leaves the catalog (e.g. ⟨32,32,64⟩) → σ is "unpriceable", NOT zero.
 * So absolute σ is visible mainly on the smaller axes of big shapes; the report
 * counts how many candidates were dropped as unpriceable.</p>
 *
 * <p>Args (all optional, positional): {@code maxDim=32 kMax=4 topN=40}.</p>
 */
@Slf4j
public class ScanSerendipityGapMap {

	enum BudAxis { U, V, W }

	record Gap(int n, int m, int p, int rI, BudAxis axis, int k, int en, int em, int ep,
			int enlargedRank, long sigma) {
		/** Fractional saving σ / (k·R(I)) — strips the trivial size-scaling. */
		double fraction() {
			return (double) sigma / ((long) k * rI);
		}

		/** Saving per extra product spent (the σ-density s of the tradeoff doc). */
		double density() {
			return (double) sigma / (k - 1);
		}
	}

	public static void main(String[] args) {
		int maxDim = args.length > 0 ? Integer.parseInt(args[0]) : 32;
		int kMax = args.length > 1 ? Integer.parseInt(args[1]) : 4;
		int topN = args.length > 2 ? Integer.parseInt(args[2]) : 40;
		Field field = Field.Q;

		FieldAwareLookup lookup = new FieldAwareLookup(field);
		final int UNKNOWN = Recombination.SotaResolver.UNKNOWN_RANK;

		List<Gap> gaps = new ArrayList<>();
		long shapesKnown = 0;
		long attempts = 0;

		for (int n = 1; n <= maxDim; n++) {
			for (int m = n; m <= maxDim; m++) {
				for (int p = m; p <= maxDim; p++) {
					int rI = lookup.findRank(n, m, p);
					if (rI >= UNKNOWN) {
						continue;
					}
					shapesKnown++;
					for (int k = 2; k <= kMax; k++) {
						// U enlarges p, V enlarges n, W enlarges m.
						attempts += 3;
						tryGap(lookup, gaps, BudAxis.U, n, m, p, rI, k, n, m, k * p, UNKNOWN);
						tryGap(lookup, gaps, BudAxis.V, n, m, p, rI, k, k * n, m, p, UNKNOWN);
						tryGap(lookup, gaps, BudAxis.W, n, m, p, rI, k, n, k * m, p, UNKNOWN);
					}
				}
			}
		}

		log.info("=== Serendipity gap map ({}), shapes≤{}, k≤{} ===", field, maxDim, kMax);
		log.info("known shapes: {} | priced σ-candidates: {} | unpriceable (enlarged off-catalog): {}",
				shapesKnown, gaps.size(), attempts - gaps.size());

		log.info("");
		log.info("--- TOP {} by ABSOLUTE σ (raw rank saved by one bud of size k) ---", topN);
		log.info("    (the user-expected ranking: bigger shapes dominate by scale)");
		printTop(gaps, topN, Comparator.comparingLong(Gap::sigma).reversed());

		log.info("");
		log.info("--- TOP {} by FRACTIONAL σ  (σ / k·R(I): genuine sub-additivity, scale-free) ---", topN);
		printTop(gaps, topN, Comparator.comparingDouble(Gap::fraction).reversed());
	}

	private static void tryGap(FieldAwareLookup lookup, List<Gap> gaps, BudAxis axis,
			int n, int m, int p, int rI, int k, int en, int em, int ep, int unknown) {
		int enlarged = lookup.findRank(en, em, ep);
		if (enlarged >= unknown) {
			return;
		}
		long sigma = (long) k * rI - enlarged;
		gaps.add(new Gap(n, m, p, rI, axis, k, en, em, ep, enlarged, sigma));
	}

	private static void printTop(List<Gap> gaps, int topN, Comparator<Gap> order) {
		gaps.stream().sorted(order).limit(topN).forEach(g -> log.info(
				"  σ={}  ({}%)  base ⟨{},{},{}⟩ R={}  {}-bud k={} → ⟨{},{},{}⟩ R={}  density={}",
				g.sigma(), Math.round(g.fraction() * 100), g.n(), g.m(), g.p(), g.rI(),
				g.axis(), g.k(), g.en(), g.em(), g.ep(), g.enlargedRank(),
				Math.round(g.density() * 10) / 10.0));
	}
}
