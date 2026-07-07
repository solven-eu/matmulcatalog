package eu.solven.matmul.docs.explore;

import java.io.File;

import eu.solven.matmul.NonCubicBilinearAlgorithm;
import eu.solven.matmul.SymmetryTransforms;
import eu.solven.matmul.catalog.FieldAwareLookup;
import eu.solven.matmul.catalog.SchemeIO;
import eu.solven.matmul.catalog.SerendipitousBudProduct;
import eu.solven.matmul.search.flip.FlipGraphWalk;
import eu.solven.matmul.search.flip.FlipObjective;
import eu.solven.matmul.search.flip.FlipObjectives;
import eu.solven.matmul.search.flip.FlipScheme;
import lombok.extern.slf4j.Slf4j;

/**
 * fmm-react 2026-07-06, ⟨20,24,25⟩ follow-up: FMM's 6466 = (⟨5,12,5⟩:204 − 42)
 * ⊗ ⟨4,2,5⟩ + patches implies a rank-204 representative with bud profile
 * 16×W2 + 2×W3 + 1×V4 (σ-savings 62). Neither our catalog 204 nor FMM's own
 * published 204 carries it (both: 18×W2 → 54). This walks the flip graph from
 * both seeds under the DIRECT serendipitous-cost objective (inner ⟨4,2,5⟩),
 * hunting a representative predicting < 6474. Success target: ≤ 6466.
 */
@Slf4j
public final class ProbeFlipWalkSerendip202425 {
	private ProbeFlipWalkSerendip202425() {}

	public static void main(String[] args) throws Exception {
		int steps = args.length > 0 ? Integer.parseInt(args[0]) : 30_000;
		int seeds = args.length > 1 ? Integer.parseInt(args[1]) : 4;
		FieldAwareLookup q = new FieldAwareLookup("Q");
		String[][] bases = {
				{ "ours-61a6cb7", "src/main/resources/schemes/known/section12/5x5x12-r204-perminov_ZT-61a6cb7.json" },
				{ "fmm-a2326", "src/main/resources/schemes/bud-bases/section12/fmm-lille_5x5x12_r204_a2326.json" } };
		FlipObjective obj = FlipObjectives.serendipitous(q, 4, 2, 5);
		for (String[] b : bases) {
			NonCubicBilinearAlgorithm raw = SchemeIO.read(new File(b[1]));
			// Orient ⟨5,5,12⟩ → ⟨5,12,5⟩ so the ⊗⟨4,2,5⟩ product hits ⟨20,24,25⟩.
			NonCubicBilinearAlgorithm alg = SymmetryTransforms.permuteAxes(raw, "ABC->ACB");
			long seedCost = SerendipitousBudProduct.serendipitousCost(alg, q, 4, 2, 5);
			log.info("{} seed: rank={} predicted ⟨20,24,25⟩ ≤ {}", b[0], alg.r, seedCost);
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
			long walkedCost = SerendipitousBudProduct.serendipitousCost(walked, q, 4, 2, 5);
			log.info("{} BEST after walk: rank={} predicted ⟨20,24,25⟩ ≤ {} (seed was {})",
					b[0], walked.r, walkedCost, seedCost);
			if (walkedCost < 6474) {
				File out = new File("target/flip-walk-5x12x5-" + b[0] + "-cost" + walkedCost + ".json");
				SchemeIO.write(walked, out);
				log.info("  WIN candidate written to {} — verify + productViaBudsBest + persist", out);
			}
		}
	}
}
