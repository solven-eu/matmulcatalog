package eu.solven.matmul.docs.explore;

import eu.solven.matmul.NonCubicBilinearAlgorithm;
import eu.solven.matmul.algebra.Field;
import eu.solven.matmul.catalog.FieldAwareLookup;
import eu.solven.matmul.catalog.SerendipitousBudProduct;
import eu.solven.matmul.search.flip.FlipGraphWalk;
import eu.solven.matmul.search.flip.FlipObjective;
import eu.solven.matmul.search.flip.FlipObjectives;
import eu.solven.matmul.search.flip.FlipScheme;
import lombok.extern.slf4j.Slf4j;

/**
 * Phase-3 decision probe: does maximizing the budScore PROXY via weighted flip
 * walks actually improve the PREDICTED serendipitous-product cost
 * ({@link SerendipitousBudProduct#serendipitousCost}, the real currency)?
 *
 * <p>For the ⟨3,3,3⟩ base: walk under several (wRank, wBud) tradeoffs, then
 * compare {@code serendipitousCost(base, ·, 3,3,3)} (predicts ⟨9,9,9⟩) and
 * {@code serendipitousCost(base, ·, 2,2,2)} (predicts ⟨6,6,6⟩) against the
 * seed's predictions and the catalog SOTA at the product shapes.</p>
 *
 * <ul>
 *   <li>If walked bases predict BELOW catalog → harvesting wins is already
 *       possible; phase 3 should prioritize a direct-cost objective + a
 *       multi-seed sweep to mine these systematically.</li>
 *   <li>If budScore rises but predicted cost doesn't beat catalog → the proxy
 *       is insufficient; phase 3 MUST start with the oracle-backed
 *       direct-cost objective.</li>
 * </ul>
 */
@Slf4j
public class ProbeFlipBudHarvest {

	public static void main(String[] args) {
		FieldAwareLookup zLookup = new FieldAwareLookup(Field.Z);
		FieldAwareLookup qLookup = new FieldAwareLookup(Field.Q);
		NonCubicBilinearAlgorithm seedAlg = zLookup.find(3, 3, 3).orElseThrow();
		FlipScheme seed = FlipScheme.of(seedAlg);
		report(qLookup, "SEED", seedAlg, seed);
		log.info("catalog SOTA at product shapes: ⟨9,9,9⟩={} ⟨6,6,6⟩={}",
				qLookup.findRank(9, 9, 9), qLookup.findRank(6, 6, 6));

		long[][] weightConfigs = { { 5, 1 }, { 10, 1 }, { 3, 1 } };
		for (long[] wc : weightConfigs) {
			FlipObjective obj = FlipObjectives.weighted(wc[0], wc[1], 0);
			FlipGraphWalk.Result best = null;
			for (long rng = 42; rng < 45; rng++) {
				FlipGraphWalk.Config cfg = new FlipGraphWalk.Config(
						250_000, rng, 1, 0.0, 4_000, 6, 0, 0.05, false);
				FlipGraphWalk.Result r = FlipGraphWalk.walk(seed, obj, cfg);
				if (best == null || r.bestCost() < best.bestCost()) {
					best = r;
				}
			}
			FlipScheme s = best.best();
			report(qLookup, "walk w=(" + wc[0] + "," + wc[1] + ")", s.toAlgorithm(), s);
		}

		// Probe v2 (the phase-3a question): walk the DIRECT predicted-product
		// cost. Fewer steps — each evaluation runs the full greedy decomposition
		// + oracle lookups — but the gradient now points at the real currency.
		for (int[] inner : new int[][] { { 3, 3, 3 }, { 2, 2, 2 } }) {
			FlipObjective direct = FlipObjectives.serendipitous(
					qLookup, inner[0], inner[1], inner[2]);
			FlipGraphWalk.Result best = null;
			for (long rng = 42; rng < 45; rng++) {
				FlipGraphWalk.Config cfg = new FlipGraphWalk.Config(
						40_000, rng, 1, 0.0, 3_000, 6, 10_000, 0.05, false);
				FlipGraphWalk.Result r = FlipGraphWalk.walk(seed, direct, cfg);
				if (best == null || r.bestCost() < best.bestCost()) {
					best = r;
				}
			}
			FlipScheme s = best.best();
			report(qLookup, "DIRECT inner=⟨" + inner[0] + "," + inner[1] + "," + inner[2] + "⟩",
					s.toAlgorithm(), s);
		}
	}

	private static void report(FieldAwareLookup lookup, String label,
			NonCubicBilinearAlgorithm alg, FlipScheme s) {
		long c999 = SerendipitousBudProduct.serendipitousCost(alg, lookup, 3, 3, 3);
		long c666 = SerendipitousBudProduct.serendipitousCost(alg, lookup, 2, 2, 2);
		log.info("{}: rank={} budScore={} margin={} → predicted ⟨9,9,9⟩≤{} ⟨6,6,6⟩≤{}",
				label, s.rank(), FlipObjectives.budScore(s), FlipObjectives.projectionMargin(s),
				c999 >= Long.MAX_VALUE / 8 ? "unknown" : c999,
				c666 >= Long.MAX_VALUE / 8 ? "unknown" : c666);
	}
}
