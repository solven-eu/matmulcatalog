package eu.solven.matmul.catalog;

import java.util.regex.Matcher;
import java.util.regex.Pattern;

/**
 * Multi-objective "best scheme" registration over the <strong>(rank, buds)</strong>
 * metric. Rank-minimality alone is too narrow: a slightly higher-rank but
 * <em>bud-richer</em> scheme can be the better building block, because its buds
 * give a cheaper serendipitous product downstream (paper §serendipitous —
 * Strassen ⟨2,2,2⟩=7 has 0 buds, the naive ⟨2,2,2⟩=8 has 4 U-buds; both belong
 * on the frontier).
 *
 * <p>So instead of collapsing each format to its single lowest rank, we keep the
 * <strong>Pareto frontier</strong> under (rank ↓, bud-richness ↑): a scheme is
 * registered as "best" iff no other scheme of the same format/field is at least
 * as good on <em>both</em> axes and strictly better on one. The rank-minimal
 * scheme is always on the frontier; a bud-richer higher-rank scheme joins it; a
 * higher-rank scheme that is <em>not</em> bud-richer is dominated and drops off.</p>
 *
 * <p>Bud-richness is scored by {@link #budScore(String)} — the total number of
 * rank-one terms that participate in a bud, summed over the three (independent)
 * factor partitions (Σ over buds of their size). It is a monotone scalar proxy
 * for serendipitous capacity; a finer multiset partial order is possible but the
 * scalar is enough to register the frontier.</p>
 */
public final class BudParetoSelection {

	private BudParetoSelection() {}

	/** Matches one bud term in a summary, e.g. {@code 4×U⟨1,1,2⟩} → (count=4, sizes=1,1,2). */
	private static final Pattern BUD_TERM =
			Pattern.compile("(\\d+)×[UVW]⟨(\\d+),(\\d+),(\\d+)⟩");

	/**
	 * Bud-richness score parsed from a {@link SerendipitousBudProduct.BudSummary}
	 * / {@link LineageBudInference.Profile} summary string. Each {@code c×T⟨…,k,…⟩}
	 * contributes {@code c·k} where {@code k ≥ 2} is the bud size (the other two
	 * triple entries are 1). Empty/null → 0.
	 */
	public static int budScore(String summary) {
		if (summary == null || summary.isEmpty()) return 0;
		Matcher m = BUD_TERM.matcher(summary);
		int score = 0;
		while (m.find()) {
			int count = Integer.parseInt(m.group(1));
			int k = Math.max(Integer.parseInt(m.group(2)),
					Math.max(Integer.parseInt(m.group(3)), Integer.parseInt(m.group(4))));
			score += count * k;
		}
		return score;
	}

	/**
	 * Does {@code (rankA, budA)} Pareto-dominate {@code (rankB, budB)}? — at least
	 * as good on both axes (lower rank, higher bud score) and strictly better on
	 * at least one.
	 */
	public static boolean dominates(int rankA, int budA, int rankB, int budB) {
		boolean atLeastAsGood = rankA <= rankB && budA >= budB;
		boolean strictlyBetter = rankA < rankB || budA > budB;
		return atLeastAsGood && strictlyBetter;
	}

	/**
	 * Pareto-frontier mask over (rank ↓, bud score ↑): {@code out[i]} is true iff
	 * candidate {@code i} is not dominated by any other candidate. Ties on both
	 * axes are co-optimal (both kept).
	 */
	public static boolean[] frontierMask(int[] ranks, int[] budScores) {
		int n = ranks.length;
		boolean[] onFrontier = new boolean[n];
		for (int i = 0; i < n; i++) {
			boolean dominated = false;
			for (int j = 0; j < n; j++) {
				if (i == j) continue;
				if (dominates(ranks[j], budScores[j], ranks[i], budScores[i])) {
					dominated = true;
					break;
				}
			}
			onFrontier[i] = !dominated;
		}
		return onFrontier;
	}
}
