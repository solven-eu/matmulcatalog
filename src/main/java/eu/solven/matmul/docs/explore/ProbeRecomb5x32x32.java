package eu.solven.matmul.docs.explore;

import java.util.List;
import java.util.Optional;

import eu.solven.matmul.NonCubicBilinearAlgorithm;
import eu.solven.matmul.algebra.Field;
import eu.solven.matmul.catalog.FieldAwareLookup;
import eu.solven.matmul.recombination.Recombination;
import eu.solven.matmul.recombination.BlockSplitSearch;
import eu.solven.matmul.search.CitedBound;
import eu.solven.matmul.search.RecombinationPoolConfig;
import lombok.extern.slf4j.Slf4j;

/**
 * Probe: can our recombination engine reach FMM's ⟨5,32,32⟩=3320?
 *
 * <p>FMM's recipe = recombination with base ⟨2,4,4⟩=26, allocA=[3,2] (n: 5=3+2),
 * allocB=[8,8,8,8], allocC=[8,8,8,8] (each 32=4·8) → 16×⟨3,8,8⟩=145 +
 * 10×⟨2,8,8⟩=100 = 3320. This probe (1) checks ⟨2,4,4⟩ is in the pool, (2) runs
 * the full {@link BlockSplitSearch#findBestStrategy} to see what WE pick, and
 * (3) forces the exact ⟨2,4,4⟩ allocation to confirm it scores 3320.</p>
 */
@Slf4j
public class ProbeRecomb5x32x32 {

	public static void main(String[] args) {
		FieldAwareLookup lookup = new FieldAwareLookup(Field.Q);
		CitedBound sota = new CitedBound(lookup);

		log.info("catalog ranks: ⟨2,4,4⟩={} ⟨3,8,8⟩={} ⟨2,8,8⟩={} | ours ⟨5,32,32⟩={} (FMM 3320)",
				lookup.findRank(2, 4, 4), lookup.findRank(3, 8, 8), lookup.findRank(2, 8, 8),
				lookup.findRank(5, 32, 32));

		// (1) Is ⟨2,4,4⟩ in the pool?
		List<BlockSplitSearch.NamedBase> pool = BlockSplitSearch.buildPool(RecombinationPoolConfig.includeDerived());
		long count244 = pool.stream().filter(nb -> {
			NonCubicBilinearAlgorithm a = nb.base();
			int[] d = { a.n, a.m, a.p };
			java.util.Arrays.sort(d);
			return d[0] == 2 && d[1] == 4 && d[2] == 4;
		}).count();
		log.info("pool size={} | ⟨2,4,4⟩-shaped bases in pool: {}", pool.size(), count244);

		// (2) What does the full search pick (unbalanced, uncapped)?
		Optional<BlockSplitSearch.NonCubicStrategy> best = BlockSplitSearch.findBestStrategy(
				5, 32, 32, pool, sota, false,
				RecombinationPoolConfig.UNBOUNDED_IMBALANCE, RecombinationPoolConfig.UNBOUNDED_COMBINATIONS, 0, Long.MAX_VALUE);
		if (best.isPresent()) {
			BlockSplitSearch.NonCubicStrategy s = best.get();
			log.info("FULL SEARCH best: rank={} label={}", s.rank(), s.label());
			if (s.recombination() != null) {
				var r = s.recombination();
				log.info("  recombination base={} allocA={} allocB={} allocC={}",
						r.baseLabel(), java.util.Arrays.toString(r.allocA()),
						java.util.Arrays.toString(r.allocB()), java.util.Arrays.toString(r.allocC()));
			}
		} else {
			log.info("FULL SEARCH: no strategy found");
		}

		// (3) Force FMM's exact ⟨2,4,4⟩ allocation.
		NonCubicBilinearAlgorithm base244 = lookup.findWithSource(2, 4, 4).orElseThrow().alg();
		log.info("forcing base ⟨{},{},{}⟩=r{} with allocA=[3,2] allocB=[8,8,8,8] allocC=[8,8,8,8]",
				base244.n, base244.m, base244.p, base244.r);
		// base is ⟨2,4,4⟩ in SOME orientation; map allocs to its actual axes.
		int[] allocA = axisAlloc(base244.n, new int[] { 3, 2 }, new int[] { 8, 8, 8, 8 });
		int[] allocB = axisAlloc(base244.m, new int[] { 3, 2 }, new int[] { 8, 8, 8, 8 });
		int[] allocC = axisAlloc(base244.p, new int[] { 3, 2 }, new int[] { 8, 8, 8, 8 });
		Recombination.Result res = Recombination.recombineWithAllocation(base244, sota, allocA, allocB, allocC);
		log.info("FORCED ⟨2,4,4⟩ recombination → rank={}  (FMM=3320, ours-best={})",
				res.totalRank, lookup.findRank(5, 32, 32));
		java.util.Map<String, Integer> leafCounts = new java.util.TreeMap<>();
		for (int[] sh : res.smallMatrixSizes) {
			leafCounts.merge(sh[0] + "x" + sh[1] + "x" + sh[2], 1, Integer::sum);
		}
		log.info("  leaf breakdown: {}", leafCounts);
	}

	/** Pick the right alloc for an axis of length 5 (→[3,2]) vs 32 (→[8,8,8,8]). */
	private static int[] axisAlloc(int dim, int[] forFive, int[] forThirtyTwo) {
		if (dim == 2) return forFive;       // base axis length 2 carries the n=5 split
		if (dim == 4) return forThirtyTwo;  // base axis length 4 carries a 32-split
		throw new IllegalStateException("unexpected base axis length " + dim);
	}
}
