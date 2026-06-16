package eu.solven.matmul.search;

import java.util.ArrayList;
import java.util.List;

import eu.solven.matmul.NonCubicBilinearAlgorithm;
import eu.solven.matmul.catalog.Recombination.SotaResolver;
import eu.solven.matmul.search.AnalyticalMaskSearch.SchemeSupports;

/**
 * Finds the allocation (block-split of a target ⟨N,M,P⟩ along a base scheme
 * ⟨bn,bm,bp⟩) that MINIMISES the recombination rank — exactly, but without the
 * full {@code (N-1 choose bn-1)·…} exhaustive sweep.
 *
 * <p>Per the 2026-06-03 theory discussion: the objective
 * {@code f(alloc) = Σ_k R(shapeₖ(alloc))} is NOT unimodal (R is a
 * number-theoretic function with dips/bumps → many spurious local minima at
 * scale: ⟨24,24,24⟩ has 413 of them), so gradient/coordinate descent gives no
 * global guarantee. Instead we use <b>branch-and-bound</b>:</p>
 * <ol>
 *   <li><b>Incumbent seed</b> = the balanced split (round-robin into bn/bm/bp
 *       equal-ish parts). Empirically the global optimum sits at (or next to)
 *       the balanced split, so this is an excellent, usually-optimal upper
 *       bound that makes almost every subtree prune immediately.</li>
 *   <li><b>Admissible lower bound</b> at a partial node (allocA fixed, allocB/C
 *       free): {@code Σ_k R(aDimₖ, bForcedₖ, cForcedₖ)} where a free axis dim is
 *       its forced minimum — {@code ⌈M/bm⌉} when product k touches ALL base
 *       blocks on that axis (so its {@code max} is unavoidably that big),
 *       else 1. R is monotone, so this never overestimates → safe to prune
 *       any subtree whose LB ≥ incumbent.</li>
 * </ol>
 *
 * <p>Result: the SAME global optimum as exhaustive search (B&amp;B is exact),
 * while visiting far fewer allocations. {@link Result#nodes} vs
 * {@link Result#fullSpace} reports the reduction.</p>
 */
public final class AllocationOptimizer {

	private AllocationOptimizer() {}

	/** The optimal allocation, its recombination rank, and the induced
	 *  per-product shape multiset. {@code nodes} = allocations actually
	 *  evaluated/visited; {@code fullSpace} = size of the exhaustive grid;
	 *  {@code improvedOnBound} = the best recombination strictly beats the
	 *  {@code upperBound} passed in (e.g. the Kronecker rank). When false, the
	 *  base does not improve on the bound and the caller should keep the bound's
	 *  construction. {@code nodes == 0} signals an early drop (the root lower
	 *  bound already met the bound, so the allocation sweep was skipped).
	 *  {@code exhaustive} = the search ran to completion (proven optimum within
	 *  {@code upperBound}); {@code false} means a {@code maxNodes} budget was hit
	 *  and {@code rank} is only an upper bound (best-found, anytime). */
	public record Result(int[] allocA, int[] allocB, int[] allocC, long rank,
			int[][] multiset, long nodes, long fullSpace, boolean improvedOnBound,
			boolean exhaustive) {}

	/** Notified each time the branch-and-bound incumbent strictly improves —
	 *  the "run 1 / run 2 / …" trajectory. {@code nodes} = allocations visited
	 *  when the improvement landed; {@code cost} = the new incumbent rank. */
	public interface IncumbentTrace {
		void onIncumbent(long nodes, long cost, int[] a, int[] b, int[] c);
	}

	/** Unbounded optimisation (returns the exact rank-minimising allocation). */
	public static Result optimize(NonCubicBilinearAlgorithm base, SotaResolver sota, int N, int M, int P) {
		return optimize(base, sota, N, M, P, SearchBudget.EXACT, null);
	}

	/**
	 * Optimise with a known {@code upperBound} for the target — typically the
	 * Kronecker rank of a factorisation, so the search can prove a base useless
	 * and stop.
	 */
	public static Result optimize(NonCubicBilinearAlgorithm base, SotaResolver rawSota,
			int N, int M, int P, long upperBound) {
		return optimize(base, rawSota, N, M, P, SearchBudget.upTo(upperBound), null);
	}

	/**
	 * Anytime / instrumented optimisation under a {@link SearchBudget} (upper
	 * bound + absolute and stagnation node caps). Hitting a node cap returns the
	 * best allocation found so far (an upper bound — {@link Result#exhaustive} is
	 * {@code false}). {@code trace} (nullable) is notified on every incumbent
	 * improvement, exposing the "how fast does the incumbent drop" trajectory.
	 */
	public static Result optimize(NonCubicBilinearAlgorithm base, SotaResolver rawSota,
			int N, int M, int P, SearchBudget budget, IncumbentTrace trace) {
		long upperBound = budget.upperBound();
		long maxNodes = budget.maxNodes();
		long maxNodesWithoutImprovement = budget.maxNodesWithoutImprovement();
		// Memoise sub-shape ranks: across the millions of allocations explored in
		// one call there are only O(maxDim³) DISTINCT sub-shapes, so caching
		// getRank(a,b,c) turns repeated catalog lookups into a hash hit. This is
		// the "subproblems recur" sharing — safe because getRank is a pure
		// function of (a,b,c). Key packs a,b,c (each ≤ target dim) in base-1024.
		java.util.HashMap<Long, Integer> rankMemo = new java.util.HashMap<>();
		SotaResolver sota = (a, b, c) -> {
			long key = ((long) a * 1024 + b) * 1024 + c;
			Integer v = rankMemo.get(key);
			if (v != null) return v;
			int rr = rawSota.getRank(a, b, c);
			rankMemo.put(key, rr);
			return rr;
		};
		SchemeSupports sup = SchemeSupports.extract(base);
		int bn = base.n, bm = base.m, bp = base.p, r = base.r;

		// Per-product "touches every block on this axis" flags (drives the LB):
		// a free axis dim is forced ≥ ⌈T/b⌉ only when BOTH relevant supports are full.
		boolean[] aFull = new boolean[r], bFull = new boolean[r], cFull = new boolean[r];
		for (int k = 0; k < r; k++) {
			aFull[k] = sup.uRowSupport[k].length == bn && sup.wRowSupport[k].length == bn;
			bFull[k] = sup.uColSupport[k].length == bm && sup.vRowSupport[k].length == bm;
			cFull[k] = sup.vColSupport[k].length == bp && sup.wColSupport[k].length == bp;
		}
		int ceilA = (N + bn - 1) / bn, ceilB = (M + bm - 1) / bm, ceilC = (P + bp - 1) / bp;

		List<int[]> allocAs = compositions(N, bn);
		List<int[]> allocBs = compositions(M, bm);
		List<int[]> allocCs = compositions(P, bp);
		// Balance-first node ordering: visit near-balanced allocations first, so
		// the incumbent drops fast in the zone that almost always holds the
		// optimum. This makes pruning bite early (every later subtree is compared
		// against an already-tight bound) AND makes the anytime / stagnation cut
		// land on a good answer instead of a lex-prefix artefact. Pure ordering
		// change — the branch-and-bound stays exact. Key = Σ aᵢ² (minimised at
		// balanced for a fixed sum, monotone in imbalance).
		allocAs.sort(java.util.Comparator.comparingLong(AllocationOptimizer::imbalanceKey));
		allocBs.sort(java.util.Comparator.comparingLong(AllocationOptimizer::imbalanceKey));
		allocCs.sort(java.util.Comparator.comparingLong(AllocationOptimizer::imbalanceKey));
		long fullSpace = (long) allocAs.size() * allocBs.size() * allocCs.size();
		// Collapse each axis to its DISTINCT per-product effective-dim vectors. The
		// recombination cost depends ONLY on those vectors (cost = Σ R(aDimₖ,bDimₖ,cDimₖ)),
		// so two allocations with the same vector are interchangeable for ANY completion.
		// This exactly canonicalises the base's allocation symmetry + structural degeneracy
		// (the ⟨2,3,3⟩ 3!-orderings collapse) WITHOUT needing the base stabiliser, AND
		// precomputes the dim vectors so the hot loop never recomputes axisDims. Balance-first
		// order is preserved (first representative kept). Exact — same optimum.
		List<AxisClass> uA = uniqueAxis(allocAs, sup.uRowSupport, sup.wRowSupport);
		List<AxisClass> uB = uniqueAxis(allocBs, sup.uColSupport, sup.vRowSupport);
		List<AxisClass> uC = uniqueAxis(allocCs, sup.vColSupport, sup.wColSupport);
		long reducedSpace = (long) uA.size() * uB.size() * uC.size();

		// Best recombination known so far = balanced split (the incumbent seed).
		int[] balA = balanced(N, bn), balB = balanced(M, bm), balC = balanced(P, bp);
		long bestCost = cost(sup, sota, balA, balB, balC);
		int[] bestA = balA, bestB = balB, bestC = balC;

		// Root lower bound: the most optimistic recombination — every free axis
		// at its forced minimum (⌈T/b⌉ for a full-support product, else 1).
		long rootLB = 0;
		for (int k = 0; k < r; k++) {
			rootLB += sota.getRank(aFull[k] ? ceilA : 1, bFull[k] ? ceilB : 1, cFull[k] ? ceilC : 1);
		}
		// Kronecker-bound drop: if even the optimistic root LB cannot beat the
		// known bound, no allocation can — skip the sweep entirely (nodes = 0).
		if (rootLB >= upperBound) {
			int[][] ms0 = AnalyticalMaskSearch.shapesAt(sup, bestA, bestB, bestC);
			return new Result(bestA.clone(), bestB.clone(), bestC.clone(), bestCost, ms0,
					0, reducedSpace, bestCost < upperBound, true);
		}
		// Report the initial (balanced) incumbent as the first trace point.
		if (trace != null) trace.onIncumbent(0, bestCost, bestA, bestB, bestC);

		// Prune threshold = the tighter of (best recombination found, the bound).
		long thr = Math.min(bestCost, upperBound);
		long[] nodes = { 0 };
		long lastImproveNode = 0; // node index of the most recent incumbent drop
		boolean budgetHit = false;
		outer:
		for (AxisClass cA : uA) {
			nodes[0]++;
			if (nodes[0] >= maxNodes
					|| nodes[0] - lastImproveNode >= maxNodesWithoutImprovement) { budgetHit = true; break outer; }
			// aDim per product is fixed; b/c at their forced minima (precomputed vector).
			int[] aDim = cA.dim();
			long lbA = 0;
			for (int k = 0; k < r; k++) {
				lbA += sota.getRank(aDim[k], bFull[k] ? ceilB : 1, cFull[k] ? ceilC : 1);
				if (lbA >= thr) break;
			}
			if (lbA >= thr) continue; // prune ALL (B,C) completions under this A

			for (AxisClass cB : uB) {
				nodes[0]++;
				if (nodes[0] >= maxNodes
						|| nodes[0] - lastImproveNode >= maxNodesWithoutImprovement) { budgetHit = true; break outer; }
				int[] bDim = cB.dim();
				long lbAB = 0;
				for (int k = 0; k < r; k++) {
					lbAB += sota.getRank(aDim[k], bDim[k], cFull[k] ? ceilC : 1);
					if (lbAB >= thr) break;
				}
				if (lbAB >= thr) continue; // prune all C completions under this (A,B)

				for (AxisClass cC : uC) {
					nodes[0]++;
					if (nodes[0] >= maxNodes
							|| nodes[0] - lastImproveNode >= maxNodesWithoutImprovement) { budgetHit = true; break outer; }
					int[] cDim = cC.dim();
					long c = 0;
					for (int k = 0; k < r; k++) {
						c += sota.getRank(aDim[k], bDim[k], cDim[k]);
						if (c >= bestCost) break; // can't beat the incumbent — stop summing
					}
					if (c < bestCost) {
						bestCost = c; bestA = cA.alloc(); bestB = cB.alloc(); bestC = cC.alloc();
						thr = Math.min(bestCost, upperBound);
						lastImproveNode = nodes[0];
						if (trace != null) trace.onIncumbent(nodes[0], bestCost, bestA, bestB, bestC);
					}
				}
			}
		}
		int[][] ms = AnalyticalMaskSearch.shapesAt(sup, bestA, bestB, bestC);
		return new Result(bestA.clone(), bestB.clone(), bestC.clone(), bestCost, ms,
				nodes[0], reducedSpace, bestCost < upperBound, !budgetHit);
	}

	/** An axis allocation collapsed to its per-product effective-dim vector plus one
	 *  representative allocation (the first, balance-first, that induced it). Distinct
	 *  vectors are the only thing the recombination cost can distinguish. */
	private record AxisClass(int[] dim, int[] alloc) {}

	/** The DISTINCT per-product dim vectors over {@code allocs} (balance-first order
	 *  preserved; first representative kept). Collapses base symmetry + structural
	 *  degeneracy exactly — two allocations with the same vector are interchangeable
	 *  for any completion — and precomputes the vectors for the hot loop. */
	private static List<AxisClass> uniqueAxis(List<int[]> allocs, int[][] supA, int[][] supB) {
		List<AxisClass> out = new ArrayList<>();
		java.util.HashSet<String> seen = new java.util.HashSet<>();
		for (int[] alloc : allocs) {
			int[] dim = axisDims(alloc, supA, supB);
			if (seen.add(java.util.Arrays.toString(dim))) out.add(new AxisClass(dim, alloc.clone()));
		}
		return out;
	}

	/** Per-product effective dim on one axis = min over the two relevant
	 *  support views of (max block size in that support). Mirrors
	 *  {@link AnalyticalMaskSearch#shapesAt}. */
	private static int[] axisDims(int[] alloc, int[][] supA, int[][] supB) {
		int[] out = new int[supA.length];
		for (int k = 0; k < supA.length; k++) {
			out[k] = Math.min(maxIndexed(alloc, supA[k]), maxIndexed(alloc, supB[k]));
		}
		return out;
	}

	private static int maxIndexed(int[] alloc, int[] indices) {
		int max = 0;
		for (int i : indices) if (alloc[i] > max) max = alloc[i];
		return max;
	}

	private static long cost(SchemeSupports sup, SotaResolver sota, int[] aA, int[] aB, int[] aC) {
		int[][] shapes = AnalyticalMaskSearch.shapesAt(sup, aA, aB, aC);
		long tot = 0;
		for (int[] s : shapes) tot += sota.getRank(s[0], s[1], s[2]);
		return tot;
	}

	/** Round-robin split of {@code total} into {@code parts} positive, as-equal-
	 *  as-possible blocks (the first {@code total mod parts} get the extra +1). */
	static int[] balanced(int total, int parts) {
		int[] a = new int[parts];
		int base = total / parts, rem = total % parts;
		for (int i = 0; i < parts; i++) a[i] = base + (i < rem ? 1 : 0);
		return a;
	}

	/** Imbalance proxy: Σ aᵢ². For a fixed sum this is minimised exactly at the
	 *  balanced split and grows monotonically with imbalance — so sorting
	 *  compositions ascending by this key is "balance-first". */
	static long imbalanceKey(int[] a) {
		long s = 0;
		for (int v : a) s += (long) v * v;
		return s;
	}

	/** All compositions of {@code total} into exactly {@code parts} positive
	 *  integers (ordered, so (1,4) ≠ (4,1)). Count = C(total-1, parts-1). */
	static List<int[]> compositions(int total, int parts) {
		List<int[]> out = new ArrayList<>();
		compose(total, parts, new int[parts], 0, out);
		return out;
	}

	private static void compose(int remaining, int partsLeft, int[] cur, int idx, List<int[]> out) {
		if (partsLeft == 1) {
			if (remaining >= 1) { cur[idx] = remaining; out.add(cur.clone()); }
			return;
		}
		for (int v = 1; v <= remaining - (partsLeft - 1); v++) {
			cur[idx] = v;
			compose(remaining - v, partsLeft - 1, cur, idx + 1, out);
		}
	}
}
