package eu.solven.matmul.docs.explore;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;

import eu.solven.matmul.NonCubicBilinearAlgorithm;
import eu.solven.matmul.Verifier;
import eu.solven.matmul.algebra.Field;
import eu.solven.matmul.catalog.FieldAwareLookup;
import eu.solven.matmul.catalog.SerendipitousBudProduct;
import eu.solven.matmul.search.als.StructuredWAls;
import lombok.extern.slf4j.Slf4j;

/**
 * CONSTRUCTIVE decision runs for the ⟨2,4,3⟩ W-structure question
 * (user 2026-06-11: "constructive — faster to get a solution, or announce
 * there is none"). Each target is a W-class SIZE PROFILE at rank 20; tied-W
 * ALS either constructs a real solution over ℝ (→ predicted ⟨6,8,9⟩ drop) or
 * accumulates evidence of infeasibility.
 *
 * <p>Profiles and their ⟨6,8,9⟩ payoffs (σ_W = {2:1, 3:5, 4:5}):
 * (2⁴,1¹²)=296 (the record, sanity check) · (2⁵,1¹⁰)=295 ·
 * (3,2³,1¹¹)=292 · (3,2⁴,1⁹)=291.</p>
 */
@Slf4j
public class ProbeTiedWDecision {

	public static void main(String[] args) {
		FieldAwareLookup z = new FieldAwareLookup(Field.Z);
		NonCubicBilinearAlgorithm known = z.find(2, 4, 3).orElseThrow();

		// 0. Sanity: the KNOWN structure, warm-started — must solve instantly.
		int[] knownClasses = classMapFromScheme(known);
		double[][] warmW = tiedWFromScheme(known, knownClasses);
		StructuredWAls.Result sanity = StructuredWAls.solve(2, 4, 3, knownClasses, 1, 50,
				known.denseU(), known.denseV(), warmW);
		log.info("sanity (known 4-pair structure, warm): residual={} solved={}",
				sanity.residual(), sanity.solved());

		// 1..3. De-novo profiles, random restarts. Big iteration budget is
		// affordable: the solver stall-exits swamped restarts early.
		runProfile("5 pairs (→295)", profile(20, 2, 2, 2, 2, 2), 300, 4_000);
		runProfile("triple + 3 pairs (→292)", profile(20, 3, 2, 2, 2), 300, 4_000);
		runProfile("triple + 4 pairs (→291)", profile(20, 3, 2, 2, 2, 2), 300, 4_000);

		// 4. Warm local repair: keep U,V and the 4 known pairs; tie one extra
		// pair among the 12 singletons (all 66 choices).
		List<Integer> singles = new ArrayList<>();
		int maxCls = Arrays.stream(knownClasses).max().orElseThrow();
		int[] count = new int[maxCls + 1];
		for (int c : knownClasses) {
			count[c]++;
		}
		for (int l = 0; l < knownClasses.length; l++) {
			if (count[knownClasses[l]] == 1) {
				singles.add(l);
			}
		}
		double best = Double.MAX_VALUE;
		for (int x = 0; x < singles.size(); x++) {
			for (int y = x + 1; y < singles.size(); y++) {
				int[] cls = knownClasses.clone();
				cls[singles.get(y)] = cls[singles.get(x)];
				int[] compact = compact(cls);
				StructuredWAls.Result r = StructuredWAls.solve(2, 4, 3, compact, 7, 2_000,
						known.denseU(), known.denseV(), null);
				best = Math.min(best, r.residual());
				if (r.solved()) {
					report("warm 5th-pair (" + singles.get(x) + "," + singles.get(y) + ")",
							compact, r);
					return;
				}
			}
		}
		log.info("warm 5th-pair repair: NO convergence over 66 pairings (best residual {})",
				best);
	}

	private static void runProfile(String label, int[] classOf, int restarts, int iters) {
		double best = Double.MAX_VALUE;
		double bestMaxAbs = 0;
		for (long seed = 0; seed < restarts; seed++) {
			StructuredWAls.Result r = StructuredWAls.solve(2, 4, 3, classOf, seed, iters,
					null, null, null);
			if (r.residual() < best) {
				best = r.residual();
				bestMaxAbs = r.maxAbs();
			}
			if (r.solved()) {
				report(label, classOf, r);
				return;
			}
		}
		log.info("{}: NO convergence in {} restarts (best residual {} at maxAbs {}) — "
				+ "evidence of infeasibility, NOT a proof", label, restarts, best, bestMaxAbs);
	}

	private static void report(String label, int[] classOf, StructuredWAls.Result r) {
		NonCubicBilinearAlgorithm alg = StructuredWAls.expand(2, 4, 3, classOf, r);
		boolean exact = Verifier.passesRandomMatmulSpotCheck(alg);
		FieldAwareLookup q = new FieldAwareLookup(Field.Q);
		long cost = SerendipitousBudProduct.serendipitousCost(alg, q, 3, 2, 3);
		log.info("{}: SOLVED over ℝ — residual={} maxAbs={} spotCheck={} → predicted"
				+ " ⟨6,8,9⟩ ≤ {} (rationalization pending)",
				label, r.residual(), r.maxAbs(), exact, cost);
	}

	/** Size profile → class map: explicit sizes first, then singletons to rank. */
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

	/** Class map of an existing scheme: group products by W direction. */
	private static int[] classMapFromScheme(NonCubicBilinearAlgorithm alg) {
		int[][] ids = SerendipitousBudProduct.independentClassIds(alg);
		return compact(ids[2]);
	}

	private static int[] compact(int[] raw) {
		java.util.Map<Integer, Integer> remap = new java.util.LinkedHashMap<>();
		int[] out = new int[raw.length];
		for (int i = 0; i < raw.length; i++) {
			out[i] = remap.computeIfAbsent(raw[i], k -> remap.size());
		}
		return out;
	}

	/** Tied W init from a scheme: one column per class (first member's column). */
	private static double[][] tiedWFromScheme(NonCubicBilinearAlgorithm alg, int[] classOf) {
		double[][] w = alg.denseW();
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
