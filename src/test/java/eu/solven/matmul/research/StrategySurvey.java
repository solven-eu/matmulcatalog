package eu.solven.matmul.research;

import eu.solven.matmul.search.RecursiveClosureSota;

import java.io.File;
import java.util.List;
import java.util.Optional;

import eu.solven.matmul.NonCubicBilinearAlgorithm;
import eu.solven.matmul.search.BlockSplitSearch;
import eu.solven.matmul.catalog.FieldAwareLookup;
import eu.solven.matmul.catalog.Recombination;
import eu.solven.matmul.catalog.SchemeIO;

/**
 * Survey: for each non-cubic gap shape, show which strategy
 * (recombination vs concat-right vs concat-below) the integrated
 * {@link BlockSplitSearch#findBestStrategy} picks, and the resulting
 * rank vs the catalog's current best.
 *
 * <p>Demonstrates that the integration correctly prefers concat for
 * narrow shapes (any axis = 2) and Strassen-recombination for thick
 * shapes.</p>
 */
public final class StrategySurvey {

	private StrategySurvey() {}

	public static void main(String[] args) throws Exception {
		FieldAwareLookup lookup = new FieldAwareLookup("Q");
		// Default pool: Strassen ⟨2,2,2⟩=7 + Laderman ⟨3,3,3⟩=23.
		List<BlockSplitSearch.NamedBase> pool = BlockSplitSearch.defaultPool();

		// One-shot vs recursive resolvers. Recursive closure (DIS09 §3.4)
		// memoises sub-shapes and is consulted by `findBestStrategy` for
		// every leaf, so deeper compositions become reachable.
		Recombination.SotaResolver flat = (a, b, c) ->
				lookup.find(a, b, c).map(alg -> alg.r).orElse(Integer.MAX_VALUE / 100);
		eu.solven.matmul.search.RecursiveClosureSota recursive =
				new eu.solven.matmul.search.RecursiveClosureSota(lookup, pool, true, true);
		Recombination.SotaResolver sota = recursive;

		// Probes — the 3 remaining FMM-better gaps + a few thick shapes for contrast +
		// larger shapes where direct catalog entries may be missing or weak.
		int[][] targets = {
				{2, 10, 15},
				{2, 10, 16},
				{2, 12, 16},
				{6, 8, 9},
				{12, 12, 12},
				{14, 14, 14},  // Pan pair-fusion known profitable (3 pairs + 1 solo, k=7)
				{18, 18, 18},  // beyond default scheme range; recursive should help
				{22, 22, 22},  // Pan pair-fusion profitable (k=11)
				{24, 24, 24},
				{32, 32, 32},
		};

		System.out.printf("%-12s  %-6s  %-32s  %-6s%n", "shape", "best", "strategy", "vs FMM");
		System.out.println("-".repeat(70));
		for (int[] s : targets) {
			int n = s[0], m = s[1], p = s[2];
			long t0 = System.currentTimeMillis();
			Optional<BlockSplitSearch.NonCubicStrategy> picked =
					BlockSplitSearch.findBestStrategy(n, m, p, pool, sota, true);
			long elapsed = System.currentTimeMillis() - t0;
			int local = lookup.find(n, m, p).map(a -> a.r).orElse(-1);
			if (picked.isEmpty()) {
				System.out.printf("⟨%d,%d,%d⟩      no strategy%n", n, m, p);
				continue;
			}
			String label = picked.get().label();
			long predictedRank = picked.get().rank();
			int flatRank = (int) Math.min(
					Long.MAX_VALUE / 4,
					BlockSplitSearch.findBestStrategy(n, m, p, pool, flat, true)
							.map(BlockSplitSearch.NonCubicStrategy::rank).orElse((long) (n * m * p)));
			System.out.printf("⟨%d,%d,%d⟩      %-6d  flat=%-6d  %-32s  catalog=%-6d  (%dms)%n",
					n, m, p, predictedRank, flatRank, label, local, elapsed);
		}
	}
}
