package eu.solven.matmul.docs.explore;

import java.util.List;
import java.util.Optional;

import eu.solven.matmul.algebra.Field;
import eu.solven.matmul.catalog.FieldAwareLookup;
import eu.solven.matmul.search.BlockSplitSearch;
import eu.solven.matmul.search.CitedBound;
import eu.solven.matmul.search.PoolConfig;
import lombok.extern.slf4j.Slf4j;

/**
 * Re-check the top "our-deriv vs FMM-deriv" behind cases with the
 * {@code includeDerived} pool (maxBaseDim=5, so ⟨2,4,4⟩-class bases are in
 * scope) — does {@link BlockSplitSearch#findBestStrategy} now predict a rank
 * below our committed catalog entry? A "yes" means the deficit is a stale
 * materialise (wrong pool), not a missing mechanism. Prediction-only (fast).
 */
@Slf4j
public class ProbeBucket1Recheck {

	// top ⟨5,·,·⟩ and a few others from the bucket-1 behind list (FMM rank in comment)
	static final int[][] SHAPES = {
			// (a) verification: prime-dim shapes FMM builds WITHOUT TA — do we reach them
			// with the unbalanced includeDerived recombination?  vs the lone TA holdout.
			{ 5, 23, 24 },   // FMM 1783 = ⟨2,4,4⟩ recombination, uneven 23={6,6,5,6}
			{ 5, 23, 32 },   // FMM 2394 = ⟨2,4,4⟩ recombination (recheck: we matched)
			{ 13, 20, 21 },  // FMM 3306 = Strassen recombination (we already match; Perminov 3165)
			{ 26, 29, 29 } };// FMM 11693 = peel + TA(cross-rotations) — the genuine TA case

	public static void main(String[] args) {
		FieldAwareLookup lookup = new FieldAwareLookup(Field.Q);
		CitedBound sota = new CitedBound(lookup);
		List<BlockSplitSearch.NamedBase> pool = BlockSplitSearch.buildPool(PoolConfig.includeDerived());
		log.info("pool size={} (includeDerived, maxBaseDim=5)", pool.size());

		int improved = 0;
		for (int[] s : SHAPES) {
			int n = s[0], m = s[1], p = s[2];
			long committed = lookup.findRank(n, m, p);
			long t0 = System.currentTimeMillis();
			Optional<BlockSplitSearch.NonCubicStrategy> best = BlockSplitSearch.findBestStrategy(
					n, m, p, pool, sota, false,
					PoolConfig.UNBOUNDED_IMBALANCE, PoolConfig.UNBOUNDED_COMBINATIONS, 0, Long.MAX_VALUE);
			long pred = best.map(BlockSplitSearch.NonCubicStrategy::rank).orElse(-1L);
			boolean win = pred > 0 && pred < committed;
			if (win) improved++;
			log.info("⟨{},{},{}⟩ committed={} predicted={} {} [{}] {}ms",
					n, m, p, committed, pred,
					win ? "← IMPROVES (" + (committed - pred) + ")" : "(no change)",
					best.map(b -> b.label()).orElse("none"),
					System.currentTimeMillis() - t0);
		}
		log.info("=== {}/{} shapes improve with the includeDerived pool ===", improved, SHAPES.length);
	}
}
