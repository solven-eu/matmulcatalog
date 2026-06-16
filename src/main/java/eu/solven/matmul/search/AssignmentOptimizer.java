package eu.solven.matmul.search;

import java.util.ArrayList;
import java.util.Comparator;
import java.util.HashMap;
import java.util.HashSet;
import java.util.List;
import java.util.Set;

import eu.solven.matmul.NonCubicBilinearAlgorithm;
import eu.solven.matmul.catalog.Recombination.SotaResolver;

/**
 * Exact single-base allocation optimiser organised as <b>partition + assignment</b>
 * (2026-06-04). Where {@link AllocationOptimizer} branches over the flat grid of
 * {@code C(T-1,b-1)³} compositions, this decomposes each axis into:
 * <ol>
 *   <li>the <b>partition</b> — the multiset of block sizes (e.g. {@code {4,3,3,3,3}});
 *       only {@code p(T, b)} of them (37 for 16 into 5), enumerated balance-first;</li>
 *   <li>the <b>assignment</b> — which size sits at which block-position. Because the
 *       base is block-asymmetric the assignment matters, and the best one is found
 *       by an exact staged branch-and-bound over the (few) distinct arrangements,
 *       costed via {@link PerAxisFingerprint}.</li>
 * </ol>
 *
 * <p>Two pruning layers keep it exact yet cheap:</p>
 * <ul>
 *   <li><b>Partition-triple LB</b>: every product's axis sub-dim is ≥ the smallest
 *       part on that axis (a max over a non-empty support ≥ the min placed size), so
 *       {@code r · R(minA, minB, minC)} lower-bounds the whole triple — skip it
 *       outright if that already meets the incumbent.</li>
 *   <li><b>Staged arrangement LB</b>: after fixing axis A's arrangement, free axes
 *       sit at their partition minima; same after fixing A,B. Monotone {@code R} ⇒
 *       admissible.</li>
 * </ul>
 *
 * <p>Returns the SAME global optimum as {@link AllocationOptimizer} (verified by
 * test), reached through a much smaller, structurally-organised search.</p>
 */
public final class AssignmentOptimizer {

	private AssignmentOptimizer() {}

	/** {@code rank} = optimal recombination rank; {@code nodes} = arrangement
	 *  leaves evaluated; {@code triplesExamined}/{@code triplesTotal} = partition
	 *  triples that survived/were generated; {@code improvedOnBound} = best beats
	 *  the {@code upperBound} passed in; {@code exhaustive} = ran to proof (no node
	 *  budget hit). */
	public record Result(int[] allocA, int[] allocB, int[] allocC, long rank,
			long nodes, long triplesExamined, long triplesTotal,
			boolean improvedOnBound, boolean exhaustive) {}

	/** Unbounded exact optimisation. */
	public static Result optimize(NonCubicBilinearAlgorithm base, SotaResolver rawSota, int N, int M, int P) {
		return optimize(base, rawSota, N, M, P, SearchBudget.EXACT);
	}

	/**
	 * Exact optimisation under a {@link SearchBudget}. The {@code upperBound}
	 * tightens the partition-triple and arrangement pruning thresholds so the
	 * search collapses when a good bound is supplied (how
	 * {@code SchemeSweep}/{@code BlockSplitSearch} call it); {@code maxNodes} is a
	 * safety cap returning best-found ({@code exhaustive=false}). The stagnation
	 * knob is not used by the partition+assignment search (its enumeration is
	 * already partition-bounded).
	 */
	public static Result optimize(NonCubicBilinearAlgorithm base, SotaResolver rawSota,
			int N, int M, int P, SearchBudget budget) {
		long upperBound = budget.upperBound();
		long maxNodes = budget.maxNodes();
		HashMap<Long, Integer> memo = new HashMap<>();
		SotaResolver sota = (a, b, c) -> {
			long key = ((long) a * 1024 + b) * 1024 + c;
			Integer v = memo.get(key);
			if (v != null) return v;
			int rr = rawSota.getRank(a, b, c);
			memo.put(key, rr);
			return rr;
		};
		PerAxisFingerprint fp = PerAxisFingerprint.of(base);
		int r = base.r;

		// Partitions per axis, balance-first (sum-of-squares ascending = balanced first).
		List<int[]> pA = partitions(N, base.n);
		List<int[]> pB = partitions(M, base.m);
		List<int[]> pC = partitions(P, base.p);
		Comparator<int[]> byBalance = Comparator.comparingLong(AssignmentOptimizer::sumSq);
		pA.sort(byBalance); pB.sort(byBalance); pC.sort(byBalance);
		long triplesTotal = (long) pA.size() * pB.size() * pC.size();

		// Cache per-arrangement sub-dim vectors (arrangement → int[r]); the same
		// arrangement recurs across many partition triples.
		HashMap<String, List<int[]>> arrCache = new HashMap<>();   // partition-key → distinct arrangements
		HashMap<String, int[]> subACache = new HashMap<>();
		HashMap<String, int[]> subBCache = new HashMap<>();
		HashMap<String, int[]> subCCache = new HashMap<>();

		// Concrete incumbent = balanced partition, identity arrangement (the
		// partitions are already balance-first, so index 0 is the most balanced).
		int[] bestA = pA.get(0).clone(), bestB = pB.get(0).clone(), bestC = pC.get(0).clone();
		long bestCost = costOf(fp, sota, subACache, subBCache, subCCache, bestA, bestB, bestC, r);
		// Pruning threshold = tighter of (incumbent, supplied upper bound).
		long thr = Math.min(bestCost, upperBound);
		long nodes = 0;
		long triplesExamined = 0;
		boolean budgetHit = false;

		outer:
		for (int[] partA : pA) {
			int minA = partA[partA.length - 1]; // partitions are descending → last = min
			for (int[] partB : pB) {
				int minB = partB[partB.length - 1];
				for (int[] partC : pC) {
					int minC = partC[partC.length - 1];
					triplesExamined++;
					// Partition-triple lower bound: every sub-dim ≥ the axis min part.
					long tripleLB = (long) r * sota.getRank(minA, minB, minC);
					if (tripleLB >= thr) continue;

					List<int[]> arrsA = arrangements(arrCache, partA);
					List<int[]> arrsB = arrangements(arrCache, partB);
					List<int[]> arrsC = arrangements(arrCache, partC);

					for (int[] aA : arrsA) {
						int[] subA = subDims(subACache, fp.axisA, aA, r);
						long lbA = 0;
						for (int k = 0; k < r; k++) {
							lbA += sota.getRank(subA[k], minB, minC);
							if (lbA >= thr) break;
						}
						if (lbA >= thr) continue;

						for (int[] aB : arrsB) {
							int[] subB = subDims(subBCache, fp.axisB, aB, r);
							long lbAB = 0;
							for (int k = 0; k < r; k++) {
								lbAB += sota.getRank(subA[k], subB[k], minC);
								if (lbAB >= thr) break;
							}
							if (lbAB >= thr) continue;

							for (int[] aC : arrsC) {
								int[] subC = subDims(subCCache, fp.axisC, aC, r);
								nodes++;
								if (nodes >= maxNodes) { budgetHit = true; break outer; }
								long c = 0;
								for (int k = 0; k < r; k++) {
									c += sota.getRank(subA[k], subB[k], subC[k]);
									if (c >= thr) break;
								}
								if (c < bestCost) {
									bestCost = c; bestA = aA.clone(); bestB = aB.clone(); bestC = aC.clone();
									thr = Math.min(bestCost, upperBound);
								}
							}
						}
					}
				}
			}
		}
		return new Result(bestA.clone(), bestB.clone(), bestC.clone(), bestCost,
				nodes, triplesExamined, triplesTotal, bestCost < upperBound, !budgetHit);
	}

	private static long costOf(PerAxisFingerprint fp, SotaResolver sota,
			HashMap<String, int[]> subACache, HashMap<String, int[]> subBCache,
			HashMap<String, int[]> subCCache, int[] a, int[] b, int[] c, int r) {
		int[] subA = subDims(subACache, fp.axisA, a, r);
		int[] subB = subDims(subBCache, fp.axisB, b, r);
		int[] subC = subDims(subCCache, fp.axisC, c, r);
		long tot = 0;
		for (int k = 0; k < r; k++) tot += sota.getRank(subA[k], subB[k], subC[k]);
		return tot;
	}

	// ----- per-arrangement sub-dim caching -----

	private static int[] subDims(HashMap<String, int[]> cache,
			List<PerAxisFingerprint.AxisGroup> groups, int[] alloc, int r) {
		String key = java.util.Arrays.toString(alloc);
		int[] cached = cache.get(key);
		if (cached != null) return cached;
		int[] out = new int[r];
		for (PerAxisFingerprint.AxisGroup g : groups) {
			int d = Math.min(maxIndexed(alloc, g.s1()), maxIndexed(alloc, g.s2()));
			for (int k : g.members()) out[k] = d;
		}
		cache.put(key, out);
		return out;
	}

	private static int maxIndexed(int[] alloc, int[] indices) {
		int max = 0;
		for (int i : indices) if (alloc[i] > max) max = alloc[i];
		return max;
	}

	// ----- partitions & arrangements -----

	private static long sumSq(int[] a) {
		long s = 0;
		for (int v : a) s += (long) v * v;
		return s;
	}

	/** Partitions of {@code total} into exactly {@code parts} positive integers,
	 *  each returned in DESCENDING order (so {@code last = min part}). */
	static List<int[]> partitions(int total, int parts) {
		List<int[]> out = new ArrayList<>();
		buildPartition(total, parts, total, new int[parts], 0, out);
		return out;
	}

	private static void buildPartition(int remaining, int partsLeft, int maxPart,
			int[] cur, int idx, List<int[]> out) {
		if (partsLeft == 1) {
			if (remaining >= 1 && remaining <= maxPart) { cur[idx] = remaining; out.add(cur.clone()); }
			return;
		}
		// next part ≤ maxPart, and leave ≥1 for each remaining slot.
		int hi = Math.min(maxPart, remaining - (partsLeft - 1));
		int lo = (remaining + partsLeft - 1) / partsLeft; // ceil(remaining/partsLeft): keep descending feasible
		for (int v = hi; v >= lo; v--) {
			cur[idx] = v;
			buildPartition(remaining - v, partsLeft - 1, v, cur, idx + 1, out);
		}
	}

	private static List<int[]> arrangements(HashMap<String, List<int[]>> cache, int[] partition) {
		String key = java.util.Arrays.toString(partition);
		List<int[]> cached = cache.get(key);
		if (cached != null) return cached;
		List<int[]> out = distinctPermutations(partition);
		cache.put(key, out);
		return out;
	}

	/** All DISTINCT orderings of the multiset {@code parts}. */
	static List<int[]> distinctPermutations(int[] parts) {
		List<int[]> raw = new ArrayList<>();
		permute(parts.clone(), 0, raw);
		List<int[]> uniq = new ArrayList<>();
		Set<String> seen = new HashSet<>();
		for (int[] x : raw) if (seen.add(java.util.Arrays.toString(x))) uniq.add(x);
		return uniq;
	}

	private static void permute(int[] a, int k, List<int[]> out) {
		if (k == a.length) { out.add(a.clone()); return; }
		for (int i = k; i < a.length; i++) {
			int t = a[k]; a[k] = a[i]; a[i] = t;
			permute(a, k + 1, out);
			t = a[k]; a[k] = a[i]; a[i] = t;
		}
	}
}
