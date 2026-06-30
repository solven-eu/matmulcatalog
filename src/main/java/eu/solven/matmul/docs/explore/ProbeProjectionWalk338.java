package eu.solven.matmul.docs.explore;

import java.util.LinkedHashMap;
import java.util.Map;

import eu.solven.matmul.NonCubicBilinearAlgorithm;
import eu.solven.matmul.verifiers.Verifier;
import eu.solven.matmul.algebra.Field;
import eu.solven.matmul.catalog.Compose;
import eu.solven.matmul.catalog.FieldAwareLookup;
import eu.solven.matmul.catalog.ProjectionSearch;
import eu.solven.matmul.search.flip.FlipGraphWalk;
import eu.solven.matmul.search.flip.FlipObjectives;
import eu.solven.matmul.search.flip.FlipScheme;
import lombok.extern.slf4j.Slf4j;

/**
 * Margin-walk demo on a REAL projection-achieved catalog value (user
 * 2026-06-12: "pick a shape with a good projection and improve it by searching
 * a better projecting base"). Target: ⟨3,3,7⟩ = 49 (room 28 over the
 * flattening LB), currently MATCHED by projecting the rank-55 ⟨3,3,8⟩ closure
 * scheme (rational, unwalkable). Any walkable ⟨3,3,8⟩ base with projected rank
 * ≤ 48 improves the catalog.
 *
 * <p>Seeds (exact-int, two distinct basins):</p>
 * <ul>
 *   <li><b>block-diagonal</b>: best integer ⟨3,3,7⟩ ⊕ naive ⟨3,3,1⟩ via
 *       {@code concatRight} → rank 58, projects to 49 by construction (the
 *       "free margin" seed — starts AT the catalog tie, needs −1);</li>
 *   <li><b>Perminov ZT r=56</b> (projects to 51, needs −3 but 2 lower rank).</li>
 * </ul>
 *
 * <p>Objective: {@code FlipObjectives.projectedTo(3,3,7)} — the exact
 * projected rank (rank/margin exchange rate is 1). Success: cost ≤ 48 →
 * materialize the projection ({@code ProjectionSearch.bestFor}), exact-verify,
 * report for catalog entry.</p>
 */
@Slf4j
public class ProbeProjectionWalk338 {

	public static void main(String[] args) {
		FieldAwareLookup q = new FieldAwareLookup(Field.Q);
		NonCubicBilinearAlgorithm t337 = q.find(3, 3, 7).orElseThrow();
		NonCubicBilinearAlgorithm blockDiag = Compose.concatRight(t337,
				NonCubicBilinearAlgorithm.naive(3, 3, 1));
		log.info("block-diagonal seed: r={} proj→⟨3,3,7⟩={} (catalog 49)",
				blockDiag.r, ProjectionSearch.projectedRank(blockDiag, 3, 3, 7, 1));

		Map<String, NonCubicBilinearAlgorithm> seeds = new LinkedHashMap<>();
		seeds.put("blockDiag_r" + blockDiag.r, blockDiag);
		q.findByHash(3, 3, 8, "5109c4b").map(FieldAwareLookup.WithSource::alg)
				.ifPresent(a -> seeds.put("perminov_ZT_r56", a));

		var objective = FlipObjectives.projectedTo(3, 3, 7, 1);
		long bestOverall = Long.MAX_VALUE;
		for (Map.Entry<String, NonCubicBilinearAlgorithm> e : seeds.entrySet()) {
			for (long rng = 1; rng <= 6; rng++) {
				FlipScheme seed = FlipScheme.of(e.getValue());
				FlipGraphWalk.Config cfg = new FlipGraphWalk.Config(400_000, rng, 2, 0.02,
						5_000, 2, 100_000, 0.02, false);
				FlipGraphWalk.Result res = FlipGraphWalk.walk(seed, objective, cfg);
				bestOverall = Math.min(bestOverall, res.bestCost());
				log.info("{} rng={}: bestProj={} (seed {}) rank={} steps={} restarts={}",
						e.getKey(), rng, res.bestCost(), res.seedCost(),
						res.best().rank(), res.steps(), res.restarts());
				if (res.bestCost() <= 48) {
					NonCubicBilinearAlgorithm base = res.best().toAlgorithm();
					boolean exactBase = res.best().isExactIntTensor();
					var hit = ProjectionSearch.bestFor(3, 3, 7,
							java.util.List.of(base), 49, 1);
					log.info("WIN candidate: base r={} exactInt={} → projection {}",
							base.r, exactBase, hit.map(h -> h.rank() + " (verified spot-check "
									+ Verifier.passesRandomMatmulSpotCheck(h.scheme()) + ")")
									.orElse("FAILED TO MATERIALIZE"));
					return;
				}
			}
		}
		log.info("no improvement: best projected rank {} vs catalog 49 over {} seeds × 6 walks "
				+ "(evidence about THESE basins only)", bestOverall, seeds.size());
	}
}
