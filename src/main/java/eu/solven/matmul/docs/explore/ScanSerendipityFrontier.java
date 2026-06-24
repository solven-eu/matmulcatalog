package eu.solven.matmul.docs.explore;

import eu.solven.matmul.recombination.Recombination;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;

import eu.solven.matmul.NonCubicBilinearAlgorithm;
import eu.solven.matmul.Verifier;
import eu.solven.matmul.algebra.Field;
import eu.solven.matmul.catalog.FieldAwareLookup;
import eu.solven.matmul.catalog.SerendipitousBudProduct;
import eu.solven.matmul.catalog.SerendipitousBudProduct.BudType;
import eu.solven.matmul.catalog.SerendipityCeiling;
import eu.solven.matmul.search.als.StructuredWAls;
import lombok.extern.slf4j.Slf4j;

/**
 * GENERAL serendipity-frontier scan for a target ⟨N,M,P⟩ (user 2026-06-12:
 * "a procedure generally checking all these combinations — kronecker products,
 * bud-structure from rank+0 to rank+N"). For every Kronecker pairing
 * outer⊗inner (all divisor splits — both orientations arise naturally), every
 * bud axis, and every outer rank r₀+δ:
 *
 * <ol>
 *   <li>price σ per axis from the catalog ({@link SerendipityCeiling#sigmaTable});</li>
 *   <li>prune by the certified ceiling (productRankFloor ≥ catalog best);</li>
 *   <li>enumerate the cheapest class profiles reaching the needed savings;</li>
 *   <li>run the tied-slot ALS decision — WITH an unconstrained CONTROL at the
 *       same (shape, rank): a cold NO is labelled "evidence" only when the
 *       control solves at that rank, else "UNINFORMATIVE" (the optimizer has
 *       no power there; tight ranks typically need warm starts or SAT).</li>
 * </ol>
 *
 * <p>Honesty tiers: SOLVED results are constructive YES over ℝ (rationalize
 * before catalog entry); everything else is a bound/evidence, never a proof.</p>
 *
 * <p>Args: {@code N M P [maxDelta=2] [restarts=200] [iters=3000]
 * [maxOuterVolume=48]}.</p>
 */
@Slf4j
public class ScanSerendipityFrontier {

	public static void main(String[] args) {
		int idx = 0;
		int bigN = args.length > idx ? Integer.parseInt(args[idx++]) : 6;
		int bigM = args.length > idx ? Integer.parseInt(args[idx++]) : 8;
		int bigP = args.length > idx ? Integer.parseInt(args[idx++]) : 9;
		int maxDelta = args.length > idx ? Integer.parseInt(args[idx++]) : 2;
		int restarts = args.length > idx ? Integer.parseInt(args[idx++]) : 200;
		int iters = args.length > idx ? Integer.parseInt(args[idx++]) : 3_000;
		int maxOuterVolume = args.length > idx ? Integer.parseInt(args[idx++]) : 48;

		FieldAwareLookup lookup = new FieldAwareLookup(Field.Q);
		int best = lookup.findRank(bigN, bigM, bigP);
		log.info("frontier scan ⟨{},{},{}⟩ (catalog best {} over Q-chain): maxDelta={} "
				+ "restarts={} iters={} maxOuterVolume={}",
				bigN, bigM, bigP, best, maxDelta, restarts, iters, maxOuterVolume);

		long t0 = System.currentTimeMillis();
		int pairings = 0;
		int candidates = 0;
		int pruned = 0;
		int solvedCount = 0;
		for (int n1 : divisors(bigN)) {
			for (int m1 : divisors(bigM)) {
				for (int p1 : divisors(bigP)) {
					int vol = n1 * m1 * p1;
					int n2 = bigN / n1;
					int m2 = bigM / m1;
					int p2 = bigP / p1;
					if (vol == 1 || n2 * m2 * p2 == 1 || vol > maxOuterVolume) {
						continue;
					}
					int r0 = lookup.findRank(n1, m1, p1);
					int ri = lookup.findRank(n2, m2, p2);
					if (r0 >= Recombination.SotaResolver.UNKNOWN_RANK || ri >= Recombination.SotaResolver.UNKNOWN_RANK) {
						continue;
					}
					pairings++;
					log.info("pairing ⟨{},{},{}⟩(r0={}) ⊗ ⟨{},{},{}⟩(R={}): plain Kron {}",
							n1, m1, p1, r0, n2, m2, p2, ri, r0 * ri);
					for (int delta = 0; delta <= maxDelta; delta++) {
						int r = r0 + delta;
						Integer controlSolves = null;
						for (BudType axis : BudType.values()) {
							SerendipityCeiling.AxisCeiling ceil = SerendipityCeiling
									.forAxis(axis, n1, m1, p1, r, lookup, n2, m2, p2);
							if (ceil.productRankFloor() >= best) {
								pruned++;
								continue;
							}
							long needed = (long) r * ri - best + 1;
							long[] sigma = SerendipityCeiling.sigmaTable(axis, r, lookup,
									n2, m2, p2);
							int minClasses = SerendipityCeiling.minClasses(axis, n1, m1, p1, r);
							List<int[]> profiles = cheapestProfiles(sigma, r, minClasses,
									needed, 2);
							if (profiles.isEmpty()) {
								continue;
							}
							if (controlSolves == null) {
								controlSolves = control(n1, m1, p1, r, restarts / 2, iters);
								log.info("  control ⟨{},{},{}⟩ r={} unconstrained: {}/{} solved",
										n1, m1, p1, r, controlSolves, restarts / 2);
							}
							for (int[] sizes : profiles) {
								candidates++;
								long sig = Arrays.stream(sizes).mapToLong(k -> sigma[k]).sum();
								String label = String.format(
										"⟨%d,%d,%d⟩ r=%d %s-classes %s (σ=%d → ⟨%d,%d,%d⟩ ≤ %d)",
										n1, m1, p1, r, axis, Arrays.toString(sizes), sig,
										bigN, bigM, bigP, r * ri - sig);
								int[] classOf = toClassMap(sizes, r);
								boolean solved = decide(label, axis, n1, m1, p1, classOf,
										restarts, iters, lookup, n2, m2, p2, controlSolves);
								if (solved) {
									solvedCount++;
								}
							}
						}
					}
				}
			}
		}
		log.info("frontier scan done: {} pairings, {} candidates ({} ceiling-pruned), "
				+ "{} SOLVED, {}s", pairings, candidates, pruned, solvedCount,
				(System.currentTimeMillis() - t0) / 1000);
	}

	/** Unconstrained (all-singleton) solve count — the POWER calibration. */
	private static int control(int n, int m, int p, int r, int restarts, int iters) {
		int[] classOf = new int[r];
		for (int l = 0; l < r; l++) {
			classOf[l] = l;
		}
		int solves = 0;
		for (long seed = 0; seed < restarts; seed++) {
			if (StructuredWAls.solve(n, m, p, classOf, seed, iters, null, null, null).solved()) {
				solves++;
			}
		}
		return solves;
	}

	private static boolean decide(String label, BudType axis, int n, int m, int p,
			int[] classOf, int restarts, int iters, FieldAwareLookup lookup,
			int n2, int m2, int p2, int controlSolves) {
		double bestRes = Double.MAX_VALUE;
		for (long seed = 0; seed < restarts; seed++) {
			StructuredWAls.Result r = StructuredWAls.solveTied(axis, n, m, p, classOf,
					seed, iters);
			bestRes = Math.min(bestRes, r.residual());
			if (r.solved()) {
				NonCubicBilinearAlgorithm alg = StructuredWAls.expandTied(axis, n, m, p,
						classOf, r);
				boolean exact = Verifier.passesRandomMatmulSpotCheck(alg);
				long cost = SerendipitousBudProduct.serendipitousCost(alg, lookup, n2, m2, p2);
				log.info("  {} → SOLVED over ℝ (residual={} maxAbs={} spotCheck={} "
						+ "serendipitousCost={}) — rationalization pending",
						label, r.residual(), r.maxAbs(), exact, cost);
				return true;
			}
		}
		log.info("  {} → no solve in {} restarts (best residual {}) — {}",
				label, restarts, bestRes,
				controlSolves > 0
						? "EVIDENCE of infeasibility (control solved " + controlSolves + ")"
						: "UNINFORMATIVE (control solved 0 — no optimizer power at this rank)");
		return false;
	}

	/**
	 * Cheapest fused-class size multisets with {@code Σσ(k) ≥ needed}, total
	 * fused ≤ r, class count (fused + singletons) ≥ minClasses. Sorted by
	 * burden Σ(k−1) ascending then σ descending; at most {@code max} returned.
	 */
	static List<int[]> cheapestProfiles(long[] sigma, int r, int minClasses, long needed,
			int max) {
		List<int[]> all = new ArrayList<>();
		enumerate(sigma, r, minClasses, new ArrayList<>(), 2, all, needed);
		all.sort(java.util.Comparator
				.<int[]> comparingInt(s -> Arrays.stream(s).map(k -> k - 1).sum())
				.thenComparing(s -> -Arrays.stream(s).mapToLong(k -> sigma[k]).sum()));
		return all.subList(0, Math.min(max, all.size()));
	}

	private static void enumerate(long[] sigma, int r, int minClasses, List<Integer> cur,
			int minSize, List<int[]> out, long needed) {
		int used = cur.stream().mapToInt(Integer::intValue).sum();
		long sig = cur.stream().mapToLong(k -> sigma[k]).sum();
		int classes = cur.size() + (r - used);
		if (sig >= needed && classes >= minClasses && !cur.isEmpty()) {
			out.add(cur.stream().mapToInt(Integer::intValue).toArray());
		}
		if (cur.size() >= 8) {
			return;
		}
		for (int k = minSize; k < sigma.length && used + k <= r; k++) {
			if (sigma[k] <= 0) {
				continue;
			}
			// Class-count constraint is monotone-decreasing in added size: prune.
			if (cur.size() + 1 + (r - used - k) < minClasses) {
				continue;
			}
			cur.add(k);
			enumerate(sigma, r, minClasses, cur, k, out, needed);
			cur.remove(cur.size() - 1);
		}
	}

	/** Size multiset → class map (fused classes first, then singletons). */
	private static int[] toClassMap(int[] sizes, int rank) {
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

	private static List<Integer> divisors(int x) {
		List<Integer> out = new ArrayList<>();
		for (int d = 1; d <= x; d++) {
			if (x % d == 0) {
				out.add(d);
			}
		}
		return out;
	}
}
