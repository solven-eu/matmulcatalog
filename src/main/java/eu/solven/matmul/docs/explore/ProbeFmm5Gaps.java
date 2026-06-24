package eu.solven.matmul.docs.explore;

import java.util.List;
import java.util.Optional;

import eu.solven.matmul.algebra.Field;
import eu.solven.matmul.catalog.FieldAwareLookup;
import eu.solven.matmul.recombination.BlockSplitSearch;
import eu.solven.matmul.search.CitedBound;
import eu.solven.matmul.search.PoolConfig;
import lombok.extern.slf4j.Slf4j;

/** Probe: for each gapped ⟨5,m,p⟩, what allocation does the full search pick, and how far from FMM? */
@Slf4j
public class ProbeFmm5Gaps {
	public static void main(String[] args) {
		FieldAwareLookup lookup = new FieldAwareLookup(Field.Q);
		CitedBound sota = new CitedBound(lookup);
		List<BlockSplitSearch.NamedBase> pool = BlockSplitSearch.buildPool(PoolConfig.includeDerived());
		// {n,m,p,fmm}
		int[][] T = {{5,23,32,2394},{5,29,30,2771},{5,23,28,2093},{5,23,24,1783},{5,24,26,2028},{5,20,23,1511}};
		for (int[] t : T) {
			Optional<BlockSplitSearch.NonCubicStrategy> best = BlockSplitSearch.findBestStrategy(
				t[0], t[1], t[2], pool, sota, false,
				PoolConfig.UNBOUNDED_IMBALANCE, PoolConfig.UNBOUNDED_COMBINATIONS, 0, Long.MAX_VALUE);
			if (best.isPresent()) {
				var s = best.get();
				String alloc = "";
				if (s.recombination() != null) {
					var r = s.recombination();
					alloc = String.format(" base=%s A=%s B=%s C=%s", r.baseLabel(),
						java.util.Arrays.toString(r.allocA()), java.util.Arrays.toString(r.allocB()),
						java.util.Arrays.toString(r.allocC()));
				}
				log.info("⟨{},{},{}⟩ ours={} FMM={} gap={} | {}{}", t[0],t[1],t[2], s.rank(), t[3], s.rank()-t[3], s.label(), alloc);
			} else {
				log.info("⟨{},{},{}⟩ NO strategy", t[0],t[1],t[2]);
			}
		}
	}
}
