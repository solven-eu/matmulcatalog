package eu.solven.matmul.search.flip;

import java.util.ArrayList;
import java.util.List;
import java.util.Random;

import lombok.extern.slf4j.Slf4j;

/**
 * Randomized flip-graph walk (Kauers–Moosbauer 2022; ternary discipline per
 * Perminov arXiv:2511.20317) with a pluggable {@link FlipObjective}. One step:
 * copy the working scheme, apply a random flip (random sign-sharing class,
 * random pair, random variant), reject it if the alphabet cap is exceeded,
 * reduce to fixpoint, then accept greedily-with-sideways
 * ({@code cost' ≤ cost}, plus an optional worse-accept probability). Plateaus
 * trigger a perturbation kick — a few unconditional flips, optionally preceded
 * by a plus-transition split when the rank budget allows.
 *
 * <p>The returned best is a <strong>bound</strong> (best found by a heuristic,
 * anytime walk) — never report it as an optimum. Exactness is structural: every
 * move rewrites two rank-one terms without changing their sum, in exact integer
 * arithmetic, so any walk vertex is a valid scheme; callers still re-verify at
 * the write boundary ({@code Verifier}) per catalog discipline.</p>
 *
 * <p>Phase 1 is fixed-format. The meta-moves (extend / project / merge /
 * product across formats, Kauers–Wood arXiv:2510.19787) are the follow-up —
 * they slot in as additional kick-time moves once the cross-format objective
 * bookkeeping (per-shape SOTA oracle) is wired.</p>
 */
@Slf4j
public final class FlipGraphWalk {

	private FlipGraphWalk() {}

	/**
	 * @param maxSteps        walk length (each step = one attempted flip)
	 * @param rngSeed         deterministic randomness — same config, same walk
	 * @param coefCap         max |coefficient| allowed; {@code 1} = ternary (ZT),
	 *                        {@code 0} = unbounded integer
	 * @param worseAcceptProb probability of accepting a cost-increasing move
	 *                        (0 = pure greedy-with-sideways)
	 * @param plateauKick     steps without a new best before a perturbation kick
	 * @param maxRankAbove    rank headroom over the seed for plus-transition
	 *                        splits during kicks (0 = never split)
	 * @param progressEvery   emit a {@code [progress]} line every N steps
	 * @param splitProb       probability a step proposes a plus-transition split
	 *                        instead of a flip (cost-gated like any move, within
	 *                        the {@code maxRankAbove} budget) — the lever that
	 *                        lets weighted objectives trade rank for structure
	 * @param fullReduce      {@code true}: reduce trials to fixpoint (zero-drops
	 *                        + merges) — right for rank-lexicographic walks.
	 *                        {@code false}: zero-drops only; merges would consume
	 *                        exactly the doubly-proportional pairs a structure
	 *                        walk builds, so they are left to the cost gate
	 */
	public record Config(long maxSteps, long rngSeed, int coefCap, double worseAcceptProb,
			long plateauKick, int maxRankAbove, long progressEvery,
			double splitProb, boolean fullReduce) {

		/** Back-compat: rank-lexicographic walk (no split moves, full reduce). */
		public Config(long maxSteps, long rngSeed, int coefCap, double worseAcceptProb,
				long plateauKick, int maxRankAbove, long progressEvery) {
			this(maxSteps, rngSeed, coefCap, worseAcceptProb, plateauKick, maxRankAbove,
					progressEvery, 0.0, true);
		}

		public static Config defaults(long rngSeed) {
			return new Config(100_000, rngSeed, 1, 0.0, 2_000, 0, 25_000);
		}
	}

	/** Best-found vertex of one walk. {@code bestCost < seedCost} ⇔ improved. */
	public record Result(FlipScheme best, long bestCost, long seedCost, long steps,
			long accepted, long kicks, long restarts) {

		public boolean improvedOverSeed() {
			return bestCost < seedCost;
		}
	}

	public static Result walk(FlipScheme seed, FlipObjective objective, Config cfg) {
		long start = System.currentTimeMillis();
		Random rng = new Random(cfg.rngSeed());
		FlipScheme origin = seed.copy();
		origin.reduce();
		int seedRank = origin.rank();
		long seedCost = objective.cost(origin);
		FlipScheme cur = origin.copy();
		long cost = seedCost;
		FlipScheme best = origin.copy();
		long bestCost = seedCost;
		long accepted = 0;
		long kicks = 0;
		long restarts = 0;
		long sinceBest = 0;
		long step = 0;
		for (; step < cfg.maxSteps(); step++) {
			FlipScheme trial = cur.copy();
			boolean moved;
			if (cfg.splitProb() > 0 && rng.nextDouble() < cfg.splitProb()
					&& trial.rank() < seedRank + cfg.maxRankAbove()) {
				// Plus-transition as a regular, cost-gated move: with a weighted
				// objective this is how the walk BUYS structure with rank.
				moved = trial.split(rng);
			} else {
				moved = randomFlip(trial, rng, cfg.coefCap());
			}
			if (!moved) {
				// No usable flip pivot: ternary walks dry up regularly (sign-shared
				// pairs are scarce over ZT). Split within the rank budget if we can;
				// otherwise RESTART from the seed — the literature's main tool (a
				// fresh random path through the component), never a dead stop.
				if (cfg.maxRankAbove() > 0 && cur.rank() < seedRank + cfg.maxRankAbove()
						&& cur.split(rng)) {
					cost = objective.cost(cur);
				} else {
					restarts++;
					cur = origin.copy();
					cost = seedCost;
				}
				continue;
			}
			if (cfg.fullReduce()) {
				trial.reduce();
			} else {
				trial.dropZeroProducts();
			}
			long c = objective.cost(trial);
			if (c <= cost || (cfg.worseAcceptProb() > 0 && rng.nextDouble() < cfg.worseAcceptProb())) {
				cur = trial;
				cost = c;
				accepted++;
				if (c < bestCost) {
					best = cur.copy();
					bestCost = c;
					sinceBest = 0;
					log.info("step {}: new best cost={} (rank={}, {})", step, c, cur.rank(),
							objective.describe());
					continue;
				}
			}
			sinceBest++;
			if (cfg.plateauKick() > 0 && sinceBest >= cfg.plateauKick()) {
				// Basin hop: perturb the best-so-far (a few unconditional flips,
				// optionally rank-budgeted splits) rather than drifting further.
				kicks++;
				sinceBest = 0;
				cur = best.copy();
				if (cfg.maxRankAbove() > 0 && cur.rank() < seedRank + cfg.maxRankAbove()) {
					cur.split(rng);
				}
				int forced = 2 + rng.nextInt(4);
				for (int k = 0; k < forced; k++) {
					randomFlip(cur, rng, cfg.coefCap());
				}
				if (cfg.fullReduce()) {
					cur.reduce();
				} else {
					cur.dropZeroProducts();
				}
				cost = objective.cost(cur);
			}
			if (cfg.progressEvery() > 0 && step > 0 && step % cfg.progressEvery() == 0) {
				log.info("[progress] {}/{} steps (cost={} best={} rank={} accepted={} kicks={}"
						+ " restarts={}) {}ms elapsed",
						step, cfg.maxSteps(), cost, bestCost, cur.rank(), accepted, kicks,
						restarts, System.currentTimeMillis() - start);
			}
		}
		return new Result(best, bestCost, seedCost, step, accepted, kicks, restarts);
	}

	/**
	 * Applies one random flip in place: uniform over sign-sharing classes (all
	 * slots pooled), then a random ordered pair and variant. Rejects (and leaves
	 * the scheme untouched) when the move would exceed {@code coefCap}; retries
	 * a few times before giving up. Returns whether a flip was applied.
	 */
	/** Public alias for drivers/probes outside this package. */
	public static boolean applyRandomFlip(FlipScheme s, Random rng, int coefCap) {
		return randomFlip(s, rng, coefCap);
	}

	static boolean randomFlip(FlipScheme s, Random rng, int coefCap) {
		record Pivot(FlipScheme.Slot slot, int[] members) {}
		List<Pivot> pivots = new ArrayList<>();
		for (FlipScheme.Slot slot : FlipScheme.Slot.values()) {
			for (int[] cls : s.signClasses(slot)) {
				pivots.add(new Pivot(slot, cls));
			}
		}
		if (pivots.isEmpty()) {
			return false;
		}
		for (int attempt = 0; attempt < 8; attempt++) {
			Pivot pick = pivots.get(rng.nextInt(pivots.size()));
			int[] cls = pick.members();
			int i = cls[rng.nextInt(cls.length)];
			int j = cls[rng.nextInt(cls.length)];
			if (i == j) {
				continue;
			}
			int sign = FlipScheme.signRatio(s.vec(pick.slot(), i), s.vec(pick.slot(), j));
			if (s.flipWithinCap(pick.slot(), i, j, sign, rng.nextBoolean(), coefCap)) {
				return true;
			}
		}
		return false;
	}
}
