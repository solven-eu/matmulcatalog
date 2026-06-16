package eu.solven.matmul.docs.explore;

import eu.solven.matmul.catalog.Recombination;

import java.util.EnumSet;

import eu.solven.matmul.NonCubicBilinearAlgorithm;
import eu.solven.matmul.Verifier;
import eu.solven.matmul.algebra.Field;
import eu.solven.matmul.catalog.FieldAwareLookup;
import eu.solven.matmul.catalog.SerendipitousBudProduct;
import eu.solven.matmul.search.flip.FlipGraphWalk;
import eu.solven.matmul.search.flip.FlipObjective;
import eu.solven.matmul.search.flip.FlipObjectives;
import eu.solven.matmul.search.flip.FlipScheme;
import lombok.extern.slf4j.Slf4j;

/**
 * DEMONSTRATION probe (user 2026-06-11): improve the smallest catalog shape
 * whose SOTA relies on a serendipitous product, by metaflip-walking its BASE
 * under the direct product-cost objective.
 *
 * <p>Target: {@code ⟨6,8,9⟩ = 296 = Serendipitous(2x4x3@82c9bd ⊗ ⟨3,2,3⟩)}
 * (smallest serendipitous-reliant SOTA, {@code derived/section9}). Walk
 * ⟨2,4,3⟩-class bases minimizing {@code serendipitousCost(·, ⟨3,2,3⟩)}; if a
 * base predicts {@code < 296}, BUILD the product from the cost-matching
 * decomposition, verify it, and report. A verified build below 296 is a new
 * ⟨6,8,9⟩ upper bound — the end-to-end validation of the bud-base approach
 * (same playbook as Kauers–Moosbauer–Wood's "divide less, conquer more").</p>
 */
@Slf4j
public class ProbeSerendipitousBaseWalk {

	public static void main(String[] args) {
		FieldAwareLookup zLookup = new FieldAwareLookup(Field.Z);
		FieldAwareLookup qLookup = new FieldAwareLookup(Field.Q);
		int n2 = 3;
		int m2 = 2;
		int p2 = 3;
		long current = qLookup.findRank(6, 8, 9);
		log.info("current ⟨6,8,9⟩ catalog rank = {} (record base: 2x4x3@82c9bd ⊗ ⟨3,2,3⟩)",
				current);

		NonCubicBilinearAlgorithm seedAlg = zLookup.find(2, 4, 3).orElseThrow();
		FlipScheme seed = FlipScheme.of(seedAlg);
		long seedCost = SerendipitousBudProduct.serendipitousCost(seedAlg, qLookup, n2, m2, p2);
		log.info("seed ⟨2,4,3⟩: rank={} budScore={} → predicted ⟨6,8,9⟩ ≤ {}",
				seed.rank(), FlipObjectives.budScore(seed), seedCost);

		// serendipitousCost is a PLATEAU objective: it only moves when a complete
		// priced bud forms. The σ-table for THIS factorization
		// (ProbeSerendipitySigmaTable) says only the W axis pays — σ_U=σ_V=0 ∀k,
		// σ_W = {2:1, 3:5, 4:5} — so the plateau tie-break must be the PRICED
		// W-axis class potential (Σ σ_W(k) over W direction-classes, overlap
		// ignored), not the axis-blind budScore (which drifted to worthless
		// U/V classes in the previous attempt).
		long[] sigmaW = new long[8];
		for (int k = 2; k < 8; k++) {
			long enlarged = qLookup.findRank(n2, k * m2, p2);
			sigmaW[k] = enlarged >= Recombination.SotaResolver.UNKNOWN_RANK ? 0
					: Math.max(0, (long) k * qLookup.findRank(n2, m2, p2) - enlarged);
		}
		FlipObjective direct = FlipObjectives.serendipitous(qLookup, n2, m2, p2);
		FlipObjective obj = new FlipObjective() {
			@Override
			public long cost(FlipScheme s) {
				long wPotential = 0;
				int[] wSizes = SerendipitousBudProduct
						.independentClassSizes(s.toAlgorithm())[2];
				for (int size : wSizes) {
					if (size >= 2) {
						wPotential += sigmaW[Math.min(size, 7)];
					}
				}
				return direct.cost(s) * 1_000 - wPotential;
			}

			@Override
			public String describe() {
				return direct.describe() + " then W-axis σ-potential↑";
			}
		};
		FlipScheme bestBase = seed;
		long bestCost = obj.cost(seed);
		for (long rng = 42; rng < 48; rng++) {
			int cap = rng < 45 ? 1 : 2;  // half ternary, half |coef|≤2 for freedom
			FlipGraphWalk.Config cfg = new FlipGraphWalk.Config(
					150_000, rng, cap, 0.01, 8_000, 2, 75_000, 0.05, false);
			FlipGraphWalk.Result r = FlipGraphWalk.walk(seed, obj, cfg);
			log.info("walk rng={} cap={}: predicted ⟨6,8,9⟩ ≤ {} (base rank={} budScore={})",
					rng, cap, Math.ceilDiv(r.bestCost(), 1_000), r.best().rank(),
					FlipObjectives.budScore(r.best()));
			if (r.bestCost() < bestCost) {
				bestCost = r.bestCost();
				bestBase = r.best();
			}
		}
		bestCost = Math.ceilDiv(bestCost, 1_000);  // cost = direct·1000 − potential
		log.info("best walked base: rank={} budScore={} predicted ⟨6,8,9⟩ ≤ {} (current {})",
				bestBase.rank(), FlipObjectives.budScore(bestBase), bestCost, current);
		if (bestCost >= current) {
			log.info("no predicted improvement below {} — demonstration not achieved with"
					+ " this seed/budget", current);
			return;
		}

		// Build from the ordering whose cost matches the prediction, then verify.
		NonCubicBilinearAlgorithm baseAlg = bestBase.toAlgorithm();
		SerendipitousBudProduct.BudDecomposition bestDec = null;
		long bestDecCost = Long.MAX_VALUE / 4;
		for (SerendipitousBudProduct.BudType[] order : SerendipitousBudProduct.ALL_ORDERINGS) {
			SerendipitousBudProduct.BudDecomposition dec =
					SerendipitousBudProduct.findBuds(baseAlg, order);
			long c = SerendipitousBudProduct.costOf(dec, qLookup, n2, m2, p2);
			if (c < bestDecCost) {
				bestDecCost = c;
				bestDec = dec;
			}
		}
		NonCubicBilinearAlgorithm product = SerendipitousBudProduct.productFromDecomposition(
				baseAlg, bestDec, qLookup, n2, m2, p2,
				EnumSet.allOf(SerendipitousBudProduct.BudType.class));
		boolean exact = Verifier.isExactNonCubic(product);
		log.info("BUILT ⟨{},{},{}⟩ rank={} exact={} — {} vs current catalog {}",
				product.n, product.m, product.p, product.r, exact,
				product.r < current ? "IMPROVEMENT (bound)" : "no improvement",
				current);
	}
}
