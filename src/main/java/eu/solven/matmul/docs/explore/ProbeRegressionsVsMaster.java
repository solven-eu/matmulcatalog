package eu.solven.matmul.docs.explore;

import java.util.List;
import java.util.Optional;

import eu.solven.matmul.catalog.FieldAwareLookup;
import eu.solven.matmul.recombination.BlockSplitSearch;
import eu.solven.matmul.recombination.Recombination;

/**
 * Throwaway probe: for each shape where fullDerive regressed vs master, run the
 * TOP-LEVEL {@link BlockSplitSearch#findBestStrategy} ONCE, valuing leaves at the
 * CURRENT disk-catalog rank (sota = lookup.findRank). NO recursion, NO writes.
 *
 * <p>This isolates the single decision that matters: given the now-complete leaf
 * catalog, does the top-level search pick a strategy ≤ master? If YES → the
 * fullDerive regression was a single-pass DERIVE-ORDER staleness (the leaf wasn't
 * yet at its final rank when this shape was first derived); a closure/second pass
 * recovers it. If NO → a genuine top-level engine gap (scoring / enumeration /
 * missing base).</p>
 */
public class ProbeRegressionsVsMaster {

	// {n, m, p, fullDeriveRank, masterRank}
	private static final int[][] CASES = {
			// Group A — branch picked trivial AxisSplit<1,1,2>; master a real recomb/concat
			{ 13, 13, 32, 3408, 3291 },
			{ 7, 7, 30, 1016, 999 },
			{ 7, 7, 32, 1078, 1069 },
			{ 3, 15, 30, 1006, 1000 },
			{ 3, 14, 28, 876, 874 },
			{ 3, 13, 26, 756, 755 },
			// Group B — same base/multiset, split ORDER flipped
			{ 21, 21, 21, 5240, 5202 },
			{ 23, 23, 23, 6707, 6672 },
			// Group D — recombination degraded to Project
			{ 7, 14, 28, 1802, 1769 },
			// Group C — Pan-TA within recombination (master RecombinationTa / PeeledViaTa)
			{ 26, 29, 29, 11808, 11693 },
			{ 28, 31, 31, 14094, 14043 },
			// Group E — AE 5x5x5 alloc-order ties
			{ 14, 15, 29, 3592, 3591 },
			{ 14, 15, 31, 3840, 3839 },
	};

	public static void main(String[] args) {
		FieldAwareLookup lookup = new FieldAwareLookup("Q");
		List<BlockSplitSearch.NamedBase> pool = BlockSplitSearch.defaultPool();
		// FLAT sota: value every sub-shape at its CURRENT catalog rank (the live
		// catalog already holds master-quality leaves). balancedOnly=false for full
		// allocation freedom — the sweep default.
		Recombination.SotaResolver sota = (a, b, c) -> lookup.findRank(a, b, c);

		System.out.println("shape         fd  master   topPick  pickLabel                          verdict");
		for (int[] c : CASES) {
			int n = c[0], m = c[1], p = c[2], fd = c[3], master = c[4];
			long t0 = System.currentTimeMillis();
			Optional<BlockSplitSearch.NonCubicStrategy> best =
					BlockSplitSearch.findBestStrategy(n, m, p, pool, sota, false);
			long ms = System.currentTimeMillis() - t0;
			int got = best.map(s -> (int) s.rank()).orElse(-1);
			String label = best.map(BlockSplitSearch.NonCubicStrategy::label).orElse("(none)");
			String verdict;
			if (got < 0) {
				verdict = "NO-RESULT";
			} else if (got <= master) {
				verdict = "RECOVERED (<= master)  => stale derive-order";
			} else if (got < fd) {
				verdict = "improved-but-short";
			} else {
				verdict = "STILL-GAP (top-level engine gap)";
			}
			System.out.printf("%2d,%2d,%2d  %6d %6d  %7d  %-34s %s  (%dms)%n",
					n, m, p, fd, master, got, label, verdict, ms);
		}
	}
}
