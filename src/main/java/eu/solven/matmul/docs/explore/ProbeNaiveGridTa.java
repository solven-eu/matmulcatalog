package eu.solven.matmul.docs.explore;

import java.util.List;
import java.util.Optional;

import eu.solven.matmul.catalog.FieldAwareLookup;
import eu.solven.matmul.recombination.BlockSplitSearch;
import eu.solven.matmul.search.RecombinationPoolConfig;
import eu.solven.matmul.recombination.Recombination;

/**
 * Throwaway probe: the fullDerive ran with {@code RecombinationPoolConfig.simple()} (cubicOnly=true),
 * which filters out EVERY non-cubic base (AxisSplits + the new Naive grids) at
 * {@link BlockSplitSearch}:233. So the in-recombination Pan-TA never fires and Group-A
 * shapes fall to a cubic Strassen. This compares the materialize-default pool against
 * the cubicOnly=false presets to find the lightest config that recovers the wins.
 */
public class ProbeNaiveGridTa {

	// {n, m, p, masterRank}
	private static final int[][] CASES = {
			{ 26, 29, 29, 11693 },
			{ 28, 31, 31, 14043 },
			{ 7, 14, 28, 1769 },
			{ 13, 13, 32, 3291 },
	};

	public static void main(String[] args) {
		FieldAwareLookup lookup = new FieldAwareLookup("Q");
		Recombination.SotaResolver sota = (a, b, c) -> lookup.findRank(a, b, c);

		record Pool(String name, List<BlockSplitSearch.NamedBase> pool) {}
		List<Pool> pools = List.of(
				new Pool("simple(cubicOnly)", BlockSplitSearch.buildPool(RecombinationPoolConfig.simple(), "Q")),
				new Pool("rectangular", BlockSplitSearch.buildPool(RecombinationPoolConfig.rectangular(), "Q")),
				new Pool("thorough", BlockSplitSearch.buildPool(RecombinationPoolConfig.thorough(), "Q")));

		for (Pool pl : pools) {
			long grids = pl.pool().stream()
					.filter(nb -> Recombination.isNaiveGrid(nb.base())).count();
			System.out.printf("%n=== pool %s (%d bases, %d naïve grids) ===%n",
					pl.name(), pl.pool().size(), grids);
			System.out.println("shape       master  topPick  pickLabel");
			for (int[] c : CASES) {
				int n = c[0], m = c[1], p = c[2], master = c[3];
				Optional<BlockSplitSearch.NonCubicStrategy> best =
						BlockSplitSearch.findBestStrategy(n, m, p, pl.pool(), sota, false);
				int got = best.map(s -> (int) s.rank()).orElse(-1);
				String label = best.map(BlockSplitSearch.NonCubicStrategy::label).orElse("(none)");
				String mark = got >= 0 && got <= master ? "<=master" : "GAP";
				System.out.printf("%2d,%2d,%2d  %6d  %7d  %-30s %s%n", n, m, p, master, got, label, mark);
			}
		}
	}
}
