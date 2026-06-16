package eu.solven.matmul.docs.explore;

import eu.solven.matmul.catalog.Recombination;

import java.io.File;
import java.util.ArrayList;
import java.util.List;

import eu.solven.matmul.NonCubicBilinearAlgorithm;
import eu.solven.matmul.catalog.FieldAwareLookup;
import eu.solven.matmul.catalog.Recombination.SotaResolver;
import eu.solven.matmul.catalog.SchemeIO;
import eu.solven.matmul.search.AllocationOptimizer;

/**
 * One-off investigation (2026-06-04, in response to the "deep-dive the
 * 5×5×5-base / 16×16×16-target explosion" request): instruments the
 * {@link AllocationOptimizer} branch-and-bound on a single (target, base) pair
 * and prints
 * <ul>
 *   <li>the exhaustive grid size {@code C(T-1,b-1)^3} (the brute-force estimate
 *       the user asked us to predict up front);</li>
 *   <li>the symmetry-reduced size: distinct <em>partitions</em> per axis (the
 *       lex/sorted-multiset reduction) — the ceiling on what block-permutation
 *       symmetry could buy IF the base scheme were block-symmetric;</li>
 *   <li>the initial (balanced) incumbent and the full incumbent-improvement
 *       trajectory ("run 1 / run 2 / …") under a node budget;</li>
 *   <li>nodes visited, and whether the search proved optimality or hit budget.</li>
 * </ul>
 *
 * <p>Run: {@code mvn -q -o exec:java
 * -Dexec.mainClass=eu.solven.matmul.docs.explore.AllocationDeepDive
 * -Dexec.args="16 16 16 src/main/resources/schemes/known/section5/alphaevolve-5x5x5_m93_a846.json 50000000"}</p>
 */
public final class AllocationDeepDive {

	private AllocationDeepDive() {}

	/** Catalog-backed SOTA over R, naive fallback for un-catalogued sub-shapes
	 *  (identical to TestAllocationOptimizer's resolver). */
	private static SotaResolver sota(FieldAwareLookup lk) {
		return (p, q, r) -> {
			if (p == 0 || q == 0 || r == 0) return 0;
			if (p == 1) return q * r;
			if (q == 1) return p * r;
			if (r == 1) return p * q;
			int v = lk.findRank(p, q, r);
			return v >= Recombination.SotaResolver.UNKNOWN_RANK ? p * q * r : v;
		};
	}

	/** C(n,k) as a long (n,k small here). */
	private static long binom(int n, int k) {
		if (k < 0 || k > n) return 0;
		k = Math.min(k, n - k);
		long r = 1;
		for (int i = 0; i < k; i++) r = r * (n - i) / (i + 1);
		return r;
	}

	/** Number of partitions of {@code total} into exactly {@code parts} positive
	 *  parts (= distinct sorted compositions). */
	private static long partitionsExact(int total, int parts) {
		// p(total, exactly parts) = p(total-parts, at most parts)
		return atMost(total - parts, parts);
	}

	/** Partitions of n into at most k parts (k,n small). */
	private static long atMost(int n, int k) {
		if (n < 0) return 0;
		if (n == 0) return 1;
		if (k == 0) return 0;
		long[][] dp = new long[k + 1][n + 1];
		for (int j = 0; j <= k; j++) dp[j][0] = 1;
		for (int j = 1; j <= k; j++)
			for (int s = 1; s <= n; s++)
				dp[j][s] = dp[j - 1][s] + (s >= j ? dp[j][s - j] : 0);
		return dp[k][n];
	}

	public static void main(String[] args) throws Exception {
		int N = args.length > 0 ? Integer.parseInt(args[0]) : 16;
		int M = args.length > 1 ? Integer.parseInt(args[1]) : 16;
		int P = args.length > 2 ? Integer.parseInt(args[2]) : 16;
		String basePath = args.length > 3 ? args[3]
				: "src/main/resources/schemes/known/section5/alphaevolve-5x5x5_m93_a846.json";
		long maxNodes = args.length > 4 ? Long.parseLong(args[4])
				: eu.solven.matmul.search.SearchHeuristics.DEFAULT_ALLOCATION_MAX_NODES;
		// Stagnation cap: stop after this many nodes with no incumbent improvement.
		long stagnation = args.length > 5 ? Long.parseLong(args[5]) : Long.MAX_VALUE;

		NonCubicBilinearAlgorithm base = SchemeIO.read(eu.solven.matmul.catalog.SchemeResolver.byHint(basePath));
		FieldAwareLookup lk = new FieldAwareLookup("R");
		SotaResolver sota = sota(lk);

		int bn = base.n, bm = base.m, bp = base.p;
		long compA = binom(N - 1, bn - 1), compB = binom(M - 1, bm - 1), compC = binom(P - 1, bp - 1);
		long fullSpace = compA * compB * compC;
		long partA = partitionsExact(N, bn), partB = partitionsExact(M, bm), partC = partitionsExact(P, bp);
		long partSpace = partA * partB * partC;

		System.out.println("=== AllocationOptimizer deep-dive ===");
		System.out.printf("target ⟨%d,%d,%d⟩  base ⟨%d,%d,%d⟩=%d  (%s)%n",
				N, M, P, bn, bm, bp, base.r, new File(basePath).getName());
		System.out.println();
		System.out.println("-- search-space estimate (the 'predict before brute-forcing' ask) --");
		System.out.printf("compositions/axis : C(%d,%d)=%d, C(%d,%d)=%d, C(%d,%d)=%d%n",
				N - 1, bn - 1, compA, M - 1, bm - 1, compB, P - 1, bp - 1, compC);
		System.out.printf("FULL grid         : %,d  (%.3e) allocations%n", fullSpace, (double) fullSpace);
		System.out.printf("partitions/axis   : %d, %d, %d   (sorted-multiset = lex reduction)%n",
				partA, partB, partC);
		System.out.printf("partition grid    : %,d  (%.3e)  → ×%.1f reduction *only if* base block-symmetric%n",
				partSpace, (double) partSpace, fullSpace / (double) Math.max(1, partSpace));
		System.out.printf("cubic-axis swap   : up to ×6 more (S3 on axes) *only if* base cyclic-symmetric%n");
		System.out.printf("  ⚠ this reduction is HYPOTHETICAL: it needs block-permutation to be a symmetry%n");
		System.out.printf("    of the base scheme itself. AlphaEvolve ⟨5,5,5⟩=93 has no such block%n");
		System.out.printf("    structure, so the FULL grid (%.3e) is the operative search space here —%n",
				(double) fullSpace);
		System.out.printf("    the optimizer branches over all of it, not the %,d partitions.%n", partSpace);
		System.out.println();

		System.out.printf("-- branch-and-bound run (maxNodes=%,d, stagnationCap=%s) --%n",
				maxNodes, stagnation == Long.MAX_VALUE ? "off" : String.format("%,d", stagnation));
		List<long[]> trail = new ArrayList<>();
		List<String> allocs = new ArrayList<>();
		AllocationOptimizer.IncumbentTrace trace = (nodes, cost, a, b, c) -> {
			trail.add(new long[] { nodes, cost });
			allocs.add(fmt(a) + " " + fmt(b) + " " + fmt(c));
		};

		long t0 = System.nanoTime();
		AllocationOptimizer.Result r = AllocationOptimizer.optimize(
				base, sota, N, M, P,
				new eu.solven.matmul.search.SearchBudget(Long.MAX_VALUE, maxNodes, stagnation), trace);
		long ms = (System.nanoTime() - t0) / 1_000_000;

		System.out.println("incumbent trajectory (run / nodes-when-found / incumbent-rank / allocation):");
		for (int i = 0; i < trail.size(); i++) {
			long nodes = trail.get(i)[0], cost = trail.get(i)[1];
			String tag = (i == 0) ? "init(balanced)" : "run " + i;
			long delta = (i == 0) ? 0 : trail.get(i - 1)[1] - cost;
			System.out.printf("  %-14s @ %,15d nodes  rank=%-6d  %s  %s%n",
					tag, nodes, cost, (i == 0 ? "" : "(−" + delta + ")"), allocs.get(i));
		}
		System.out.println();
		long lastImprove = trail.isEmpty() ? 0 : trail.get(trail.size() - 1)[0];
		System.out.printf("best rank      : %d%n", r.rank());
		System.out.printf("nodes visited  : %,d  (of %,d full grid = %.4f%%)%n",
				r.nodes(), fullSpace, 100.0 * r.nodes() / fullSpace);
		System.out.printf("last improve   : node %,d  → %,d nodes spent with NO improvement after it%n",
				lastImprove, r.nodes() - lastImprove);
		System.out.printf("exhaustive?    : %s%n", r.exhaustive()
				? "YES — proven optimum"
				: "NO — budget hit; rank is an upper bound (best-found, anytime)");
		System.out.printf("wall-clock     : %,d ms%n", ms);
		if (!r.exhaustive()) {
			double frac = r.nodes() / (double) fullSpace;
			System.out.printf("projection     : at this rate, exhaustive ≈ %,.0f ms (%.1f h) to scan the full grid%n",
					ms / Math.max(1e-9, frac), ms / Math.max(1e-9, frac) / 3_600_000.0);
		}

		// EXACT partition+assignment optimiser for comparison.
		System.out.println();
		System.out.println("-- exact partition+assignment B&B (AssignmentOptimizer) --");
		long t1 = System.nanoTime();
		// Prune against the anytime result (the representative wiring scenario:
		// SchemeSweep always supplies a cross-base bound). Safety node budget.
		eu.solven.matmul.search.AssignmentOptimizer.Result asg =
				eu.solven.matmul.search.AssignmentOptimizer.optimize(
						base, sota, N, M, P,
						new eu.solven.matmul.search.SearchBudget(r.rank() + 1, 300_000_000L, Long.MAX_VALUE));
		long ms2 = (System.nanoTime() - t1) / 1_000_000;
		System.out.printf("EXACT rank     : %d   %s%n", asg.rank(),
				asg.rank() <= r.rank() ? "(≤ anytime B&B result " + r.rank() + ")" : "(> anytime?!)");
		System.out.printf("arrangement leaves evaluated : %,d%n", asg.nodes());
		System.out.printf("partition triples examined   : %,d of %,d%n", asg.triplesExamined(), asg.triplesTotal());
		System.out.printf("best allocation: %s %s %s%n", fmt(asg.allocA()), fmt(asg.allocB()), fmt(asg.allocC()));
		System.out.printf("wall-clock     : %,d ms   (exhaustive=%s)%n", ms2, asg.exhaustive());
	}

	private static String fmt(int[] a) {
		StringBuilder sb = new StringBuilder("(");
		for (int i = 0; i < a.length; i++) { if (i > 0) sb.append(','); sb.append(a[i]); }
		return sb.append(')').toString();
	}
}
