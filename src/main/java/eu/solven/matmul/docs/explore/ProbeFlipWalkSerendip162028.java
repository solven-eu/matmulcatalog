package eu.solven.matmul.docs.explore;

import java.io.File;

import eu.solven.matmul.NonCubicBilinearAlgorithm;
import eu.solven.matmul.catalog.FieldAwareLookup;
import eu.solven.matmul.catalog.SchemeIO;
import eu.solven.matmul.catalog.SerendipitousBudProduct;
import eu.solven.matmul.search.flip.FlipGraphWalk;
import eu.solven.matmul.search.flip.FlipObjective;
import eu.solven.matmul.search.flip.FlipObjectives;
import eu.solven.matmul.search.flip.FlipScheme;
import lombok.extern.slf4j.Slf4j;

/**
 * fmm-gap 2026-07-08, ⟨16,20,28⟩: FMM's 4944 = (⟨4,5,7⟩:104 − 17) ⊗ ⟨4,4,4⟩:48
 * + 8·⟨4,4,8⟩:96 implies a rank-104 representative absorbing 17 terms into 8
 * ⟨4,4,8⟩ blocks (one block eats a 3-term class — span compression). Neither
 * our catalog 104s (best greedy σ → 4989) nor FMM's own published base
 * ({2=7,3=1} U-classes → 4989) shows it. Walk the flip graph from all seeds
 * under the DIRECT serendipitous objective (inner ⟨4,4,4⟩), hunting a
 * representative predicting < 4986 (our incumbent). Success target: ≤ 4944.
 */
@Slf4j
public final class ProbeFlipWalkSerendip162028 {
	private ProbeFlipWalkSerendip162028() {}

	public static void main(String[] args) throws Exception {
		int steps = args.length > 0 ? Integer.parseInt(args[0]) : 30_000;
		int seeds = args.length > 1 ? Integer.parseInt(args[1]) : 4;
		FieldAwareLookup q = new FieldAwareLookup("Q");
		String[][] bases = {
				{ "ours-505458b", "src/main/resources/schemes/known/section7/4x5x7-r104-perminov_c924_ZT-505458b.json" },
				{ "ours-655248c", "src/main/resources/schemes/known/section7/4x5x7-r104-perminov_cr400_fv163_cn931_ZT_reduced-655248c.json" },
				{ "fmm-a1410", "src/main/resources/schemes/bud-bases/section7/fmm-lille_4x5x7_r104_a1410.json" } };
		FlipObjective obj = FlipObjectives.serendipitous(q, 4, 4, 4);
		for (String[] b : bases) {
			NonCubicBilinearAlgorithm alg = SchemeIO.read(new File(b[1]));
			long seedCost = SerendipitousBudProduct.serendipitousCost(alg, q, 4, 4, 4);
			log.info("{} seed: rank={} predicted ⟨16,20,28⟩ ≤ {}", b[0], alg.r, seedCost);
			FlipScheme seed = FlipScheme.of(alg);
			FlipGraphWalk.Result best = null;
			for (long rng = 42; rng < 42 + seeds; rng++) {
				FlipGraphWalk.Config cfg = new FlipGraphWalk.Config(
						steps, rng, 1, 0.0, Math.max(3_000, steps / 10), 6, 10_000, 0.05, false);
				FlipGraphWalk.Result r = FlipGraphWalk.walk(seed, obj, cfg);
				log.info("  [{} rng={}] walked bestCost={}", b[0], rng, r.bestCost());
				if (best == null || r.bestCost() < best.bestCost()) best = r;
			}
			NonCubicBilinearAlgorithm walked = best.best().toAlgorithm();
			long walkedCost = SerendipitousBudProduct.serendipitousCost(walked, q, 4, 4, 4);
			log.info("{} BEST after walk: rank={} predicted ⟨16,20,28⟩ ≤ {} (seed was {})",
					b[0], walked.r, walkedCost, seedCost);
			if (walkedCost < 4986) {
				File out = new File("target/flip-walk-4x5x7-" + b[0] + "-cost" + walkedCost + ".json");
				SchemeIO.write(walked, out);
				log.info("  WIN candidate written to {} — verify + productViaBudsBest + persist", out);
			}
		}
	}
}
