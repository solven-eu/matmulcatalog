package eu.solven.matmul.docs.explore;

import java.util.Arrays;

import eu.solven.matmul.NonCubicBilinearAlgorithm;
import eu.solven.matmul.Verifier;
import eu.solven.matmul.algebra.Field;
import eu.solven.matmul.catalog.FieldAwareLookup;
import eu.solven.matmul.catalog.SerendipitousBudProduct;
import eu.solven.matmul.search.als.StructuredWAls;
import lombok.extern.slf4j.Slf4j;

/**
 * WARM local repair for the reversed ⟨6,8,9⟩ orientation (2026-06-12): a
 * V-TRIPLE in a rank-15 ⟨3,2,3⟩ base gives ⟨6,8,9⟩ ≤ 294 (σ_V(3)=6 against
 * inner ⟨2,4,3⟩; a quad gives ≤ 293). The catalog base is bud-free, and cold
 * ALS has no power at tight rank — but warm repair (keep the exact base
 * factors, tie a V-class, descend locally) retains power near the incumbent.
 *
 * <p>Tying V of ⟨3,2,3⟩ = tying W of the doubly-rotated ⟨3,3,2⟩ scheme
 * (cyclicShift²), so we warm-start the W-tied solver on the rotated factors
 * and map any solution back via {@code expandTied(V, …)} semantics (here:
 * cyclicShift once from ⟨3,3,2⟩).</p>
 */
@Slf4j
public class ProbeTiedVRepair323 {

	public static void main(String[] args) {
		FieldAwareLookup q = new FieldAwareLookup(Field.Q);
		NonCubicBilinearAlgorithm base = q.find(3, 2, 3).orElseThrow();
		// ⟨3,2,3⟩ → ⟨2,3,3⟩ → ⟨3,3,2⟩: W of the rotated scheme is V of the base.
		NonCubicBilinearAlgorithm rot = base.cyclicShift().cyclicShift();
		int r = rot.r;
		double[][] wu = rot.denseU();
		double[][] wv = rot.denseV();
		double[][] ww = rot.denseW();

		double bestRes = Double.MAX_VALUE;
		int[] bestTriple = null;
		int tried = 0;
		long t0 = System.currentTimeMillis();
		for (int x = 0; x < r; x++) {
			for (int y = x + 1; y < r; y++) {
				for (int z = y + 1; z < r; z++) {
					int[] classOf = classMap(r, x, y, z);
					double[][] warmW = tiedW(ww, classOf);
					StructuredWAls.Result res = StructuredWAls.solve(3, 3, 2, classOf, 11,
							1_500, wu, wv, warmW);
					tried++;
					if (res.residual() < bestRes) {
						bestRes = res.residual();
						bestTriple = new int[] { x, y, z };
					}
					if (res.solved()) {
						NonCubicBilinearAlgorithm alg = StructuredWAls
								.expand(3, 3, 2, classOf, res).cyclicShift();
						boolean exact = Verifier.passesRandomMatmulSpotCheck(alg);
						long cost = SerendipitousBudProduct.serendipitousCost(alg, q, 2, 4, 3);
						log.info("V-triple {} SOLVED over ℝ — residual={} maxAbs={} "
								+ "spotCheck={} → serendipitousCost ⟨6,8,9⟩ ≤ {} "
								+ "(rationalization pending)",
								Arrays.toString(new int[] { x, y, z }), res.residual(),
								res.maxAbs(), exact, cost);
						return;
					}
					if (tried % 50 == 0) {
						log.info("[progress] {}/455 triples, best residual {} at {} ({}ms)",
								tried, bestRes, Arrays.toString(bestTriple),
								System.currentTimeMillis() - t0);
					}
				}
			}
		}
		log.info("V-triple warm repair: NO convergence over {} triples (best residual {} "
				+ "at {}) — local evidence around THIS base only, not a proof",
				tried, bestRes, Arrays.toString(bestTriple));
	}

	private static int[] classMap(int r, int x, int y, int z) {
		int[] out = new int[r];
		out[x] = 0;
		out[y] = 0;
		out[z] = 0;
		int cls = 1;
		for (int l = 0; l < r; l++) {
			if (l != x && l != y && l != z) {
				out[l] = cls++;
			}
		}
		return out;
	}

	private static double[][] tiedW(double[][] w, int[] classOf) {
		int d = Arrays.stream(classOf).max().orElseThrow() + 1;
		double[][] out = new double[w.length][d];
		for (int l = classOf.length - 1; l >= 0; l--) {
			for (int c = 0; c < w.length; c++) {
				out[c][classOf[l]] = w[c][l];
			}
		}
		return out;
	}
}
