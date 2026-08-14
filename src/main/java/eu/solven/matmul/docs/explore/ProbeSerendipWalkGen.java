package eu.solven.matmul.docs.explore;

import java.util.EnumSet;

import eu.solven.matmul.NonCubicBilinearAlgorithm;
import eu.solven.matmul.algebra.Field;
import eu.solven.matmul.catalog.FieldAwareLookup;
import eu.solven.matmul.catalog.SchemeIO;
import eu.solven.matmul.catalog.SerendipitousBudProduct;
import eu.solven.matmul.recombination.Recombination;
import eu.solven.matmul.search.flip.FlipGraphWalk;
import eu.solven.matmul.search.flip.FlipObjective;
import eu.solven.matmul.search.flip.FlipObjectives;
import eu.solven.matmul.search.flip.FlipScheme;
import eu.solven.matmul.verifiers.Verifier;
import lombok.extern.slf4j.Slf4j;

/**
 * Generalised serendipitous-base walk (the parameterised sibling of
 * {@link ProbeSerendipitousBaseWalk}). Metaflip-walks a BASE minimising the
 * direct serendipitous product-cost against an INNER, with an all-axes
 * σ-potential tie-break to break the plateau; if a base predicts below the
 * target's catalog SOTA, builds and verifies.
 *
 * <p>Args: {@code tN tM tP  bN bM bP  iN iM iP}. Example (⟨6,9,12⟩ via
 * ⟨2,3,4⟩ ⊗ ⟨3,3,3⟩): {@code 6 9 12  2 3 4  3 3 3}. The ⟨3,3,3⟩-class inner
 * pays σ≈6 per k=2 bud (vs σ=1 for ⟨3,2,3⟩), so a single extra bud is a −6.</p>
 */
@Slf4j
public class ProbeSerendipWalkGen {

	public static void main(String[] args) throws Exception {
		int tN = Integer.parseInt(args[0]), tM = Integer.parseInt(args[1]), tP = Integer.parseInt(args[2]);
		int bN = Integer.parseInt(args[3]), bM = Integer.parseInt(args[4]), bP = Integer.parseInt(args[5]);
		int n2 = Integer.parseInt(args[6]), m2 = Integer.parseInt(args[7]), p2 = Integer.parseInt(args[8]);
		FieldAwareLookup zLookup = new FieldAwareLookup(Field.Z);
		FieldAwareLookup qLookup = new FieldAwareLookup(Field.Q);
		long current = qLookup.findRank(tN, tM, tP);
		NonCubicBilinearAlgorithm seedAlg = args.length > 9
				? SchemeIO.read(new java.io.File(args[9]))
				: zLookup.find(bN, bM, bP).orElseThrow();
		FlipScheme seed = FlipScheme.of(seedAlg);
		long seedCost = SerendipitousBudProduct.serendipitousCost(seedAlg, qLookup, n2, m2, p2);
		log.info("target ⟨{},{},{}⟩ current={} | seed base ⟨{},{},{}⟩ rank={} → predicted ≤ {}",
				tN, tM, tP, current, bN, bM, bP, seed.rank(), seedCost);

		// Per-axis σ table. U-bud enlarges p2 (index 0), V-bud enlarges n2 (index 1),
		// W-bud enlarges m2 (index 2) — matching SerendipitousBudProduct.independentClassSizes.
		long RI = qLookup.findRank(n2, m2, p2);
		long[][] sigma = new long[3][8];
		for (int k = 2; k < 8; k++) {
			sigma[0][k] = sig(qLookup, RI, k, n2, m2, k * p2);
			sigma[1][k] = sig(qLookup, RI, k, k * n2, m2, p2);
			sigma[2][k] = sig(qLookup, RI, k, n2, k * m2, p2);
		}
		log.info("σ/bud (k=2) : U={} V={} W={}", sigma[0][2], sigma[1][2], sigma[2][2]);

		FlipObjective direct = FlipObjectives.serendipitous(qLookup, n2, m2, p2);
		FlipObjective obj = new FlipObjective() {
			@Override
			public long cost(FlipScheme s) {
				long potential = 0;
				int[][] sizes = SerendipitousBudProduct.independentClassSizes(s.toAlgorithm());
				for (int ax = 0; ax < 3; ax++) {
					for (int size : sizes[ax]) {
						if (size >= 2) {
							potential += sigma[ax][Math.min(size, 7)];
						}
					}
				}
				return direct.cost(s) * 1_000 - potential;
			}

			@Override
			public String describe() {
				return direct.describe() + " then all-axes σ-potential↑";
			}
		};

		FlipScheme bestBase = seed;
		long bestCost = obj.cost(seed);
		for (long rng = 42; rng < 50; rng++) {
			int cap = rng < 46 ? 1 : 2;
			FlipGraphWalk.Config cfg = new FlipGraphWalk.Config(
					200_000, rng, cap, 0.01, 8_000, 2, 100_000, 0.05, false);
			FlipGraphWalk.Result r = FlipGraphWalk.walk(seed, obj, cfg);
			long pred = Math.ceilDiv(r.bestCost(), 1_000);
			log.info("walk rng={} cap={}: predicted ⟨{},{},{}⟩ ≤ {} (base rank={})",
					rng, cap, tN, tM, tP, pred, r.best().rank());
			if (r.bestCost() < bestCost) {
				bestCost = r.bestCost();
				bestBase = r.best();
			}
		}
		long predicted = Math.ceilDiv(bestCost, 1_000);
		log.info("best walked base rank={} predicted ⟨{},{},{}⟩ ≤ {} (current {})",
				bestBase.rank(), tN, tM, tP, predicted, current);
		if (predicted >= current) {
			log.info("NO predicted improvement below {} — not achieved with this seed/budget", current);
			return;
		}
		// Build + verify.
		NonCubicBilinearAlgorithm baseAlg = bestBase.toAlgorithm();
		SerendipitousBudProduct.BudDecomposition bestDec = null;
		long bestDecCost = Long.MAX_VALUE / 4;
		for (SerendipitousBudProduct.BudType[] order : SerendipitousBudProduct.ALL_ORDERINGS) {
			SerendipitousBudProduct.BudDecomposition dec = SerendipitousBudProduct.findBuds(baseAlg, order);
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
		log.info("BUILT ⟨{},{},{}⟩ rank={} exact={} — {} vs current {}",
				product.n, product.m, product.p, product.r, exact,
				product.r < current ? "*** IMPROVEMENT ***" : "no improvement", current);
	}

	private static long sig(FieldAwareLookup q, long RI, int k, int n, int m, int p) {
		long enl = q.findRank(n, m, p);
		return enl >= Recombination.SotaResolver.UNKNOWN_RANK ? 0 : Math.max(0, (long) k * RI - enl);
	}
}
