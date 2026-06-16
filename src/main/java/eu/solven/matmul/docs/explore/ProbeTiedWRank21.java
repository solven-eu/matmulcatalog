package eu.solven.matmul.docs.explore;

import eu.solven.matmul.NonCubicBilinearAlgorithm;
import eu.solven.matmul.Verifier;
import eu.solven.matmul.catalog.FieldAwareLookup;
import eu.solven.matmul.algebra.Field;
import eu.solven.matmul.catalog.SerendipitousBudProduct;
import eu.solven.matmul.search.als.StructuredWAls;
import lombok.extern.slf4j.Slf4j;

/**
 * Rank-TRADEABLE follow-up to {@link ProbeTiedWDecision}: at base rank 21
 * (⟨2,4,3⟩, one above optimal), predicted ⟨6,8,9⟩ cost is 21·15 − σ_W, so any
 * W-profile with σ ≥ 20 beats the current 296. σ_W = {2:1, 3:5, 4:5} — triples
 * are the efficient unit. Ladder from the rank-21 ceiling (7 triples → 280,
 * exactly at the 7-class lower bound) down to the break-even (4 triples → 295).
 */
@Slf4j
public class ProbeTiedWRank21 {

	public static void main(String[] args) {
		run("7 triples (→280, ceiling)", profile(21, 3, 3, 3, 3, 3, 3, 3));
		run("6 triples (→285)", profile(21, 3, 3, 3, 3, 3, 3));
		run("5 triples (→290)", profile(21, 3, 3, 3, 3, 3));
		run("4 triples + 2 pairs (→293)", profile(21, 3, 3, 3, 3, 2, 2));
		run("4 triples (→295)", profile(21, 3, 3, 3, 3));
	}

	private static void run(String label, int[] classOf) {
		double best = Double.MAX_VALUE;
		double bestMaxAbs = 0;
		for (long seed = 0; seed < 300; seed++) {
			StructuredWAls.Result r = StructuredWAls.solve(2, 4, 3, classOf, seed, 4_000,
					null, null, null);
			if (r.residual() < best) {
				best = r.residual();
				bestMaxAbs = r.maxAbs();
			}
			if (r.solved()) {
				NonCubicBilinearAlgorithm alg = StructuredWAls.expand(2, 4, 3, classOf, r);
				boolean exact = Verifier.passesRandomMatmulSpotCheck(alg);
				long cost = SerendipitousBudProduct.serendipitousCost(alg,
						new FieldAwareLookup(Field.Q), 3, 2, 3);
				log.info("{}: SOLVED over ℝ — residual={} maxAbs={} spotCheck={} → predicted"
						+ " ⟨6,8,9⟩ ≤ {} (rationalization pending)",
						label, r.residual(), r.maxAbs(), exact, cost);
				return;
			}
		}
		log.info("{}: NO convergence in 300 restarts (best residual {} at maxAbs {}) — "
				+ "evidence of infeasibility, NOT a proof", label, best, bestMaxAbs);
	}

	private static int[] profile(int rank, int... sizes) {
		int[] out = new int[rank];
		int l = 0;
		int cls = 0;
		for (int s : sizes) {
			for (int k = 0; k < s; k++) {
				out[l++] = cls;
			}
			cls++;
		}
		while (l < rank) {
			out[l++] = cls++;
		}
		return out;
	}
}
