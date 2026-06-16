package eu.solven.matmul.search.flip;

import java.util.LinkedHashMap;
import java.util.Map;
import java.util.Random;

import lombok.extern.slf4j.Slf4j;

/**
 * META flip-graph walk (Kauers–Wood, arXiv:2510.19787): a random walk that is
 * free to change the FORMAT, not just the scheme. On top of the fixed-format
 * moves (flips, splits, cost-gated reductions), two meta-moves fire with
 * configured probabilities:
 *
 * <ul>
 *   <li>{@link FlipScheme#extendAxis extend} — grow an axis by one, covering
 *       the new slice with naive products (rank +m·p / +n·p / +n·m);</li>
 *   <li>{@link FlipScheme#projectAxis project} — drop a random index on an
 *       axis and dead-code-eliminate (rank −μ at best).</li>
 * </ul>
 *
 * Meta-moves are exploratory (taken unconditionally, within the dim/rank
 * bounds) — structure discovered in one format leaks into its neighbours, the
 * Kauers–Wood mechanism. Within a format, moves stay cost-gated by the
 * objective. The walk tracks the BEST vertex per visited format; costs are
 * never compared across formats.
 *
 * <p>Every per-shape best is a <strong>bound</strong> (heuristic anytime
 * walk). Callers re-verify at the write boundary and must compare against the
 * catalog before claiming any improvement.</p>
 */
@Slf4j
public final class MetaFlipWalk {

	private MetaFlipWalk() {}

	/**
	 * @param maxSteps        walk length
	 * @param rngSeed         deterministic randomness
	 * @param coefCap         alphabet cap (1 = ternary/ZT, 0 = unbounded Z)
	 * @param worseAcceptProb probability of accepting a cost-increasing
	 *                        fixed-format move
	 * @param splitProb       probability a fixed-format step proposes a split
	 * @param extendProb      probability a step takes the extend meta-move
	 * @param projectProb     probability a step takes the project meta-move
	 * @param minDim          never project an axis below this
	 * @param maxDim          never extend an axis above this
	 * @param maxRank         never extend/split past this rank
	 * @param fullReduce      reduce trials to fixpoint vs zero-drops only
	 * @param progressEvery   emit a {@code [progress]} line every N steps
	 */
	public record Config(long maxSteps, long rngSeed, int coefCap, double worseAcceptProb,
			double splitProb, double extendProb, double projectProb,
			int minDim, int maxDim, int maxRank, boolean fullReduce, long progressEvery) {

		public static Config defaults(long rngSeed) {
			return new Config(200_000, rngSeed, 1, 0.01, 0.005, 0.002, 0.004,
					2, 6, 400, true, 50_000);
		}
	}

	/** Best vertex found at one format ({@code cost} under the walk's objective). */
	public record ShapeBest(FlipScheme scheme, long cost) {}

	public record Result(Map<String, ShapeBest> bestByShape, long steps, long accepted,
			long extendsTaken, long projectsTaken, long restarts) {}

	public static Result walk(FlipScheme seed, FlipObjective objective, Config cfg) {
		long start = System.currentTimeMillis();
		Random rng = new Random(cfg.rngSeed());
		FlipScheme origin = seed.copy();
		origin.reduce();
		FlipScheme cur = origin.copy();
		long cost = objective.cost(cur);
		Map<String, ShapeBest> best = new LinkedHashMap<>();
		record(best, cur, cost);
		long accepted = 0;
		long extendsTaken = 0;
		long projectsTaken = 0;
		long restarts = 0;
		long step = 0;
		for (; step < cfg.maxSteps(); step++) {
			double roll = rng.nextDouble();
			if (roll < cfg.extendProb()) {
				FlipScheme.Axis axis = pickExtendable(cur, rng, cfg);
				if (axis != null) {
					cur = cur.extendAxis(axis);
					cost = objective.cost(cur);
					record(best, cur, cost);
					extendsTaken++;
				}
				continue;
			}
			if (roll < cfg.extendProb() + cfg.projectProb()) {
				FlipScheme.Axis axis = pickProjectable(cur, rng, cfg);
				if (axis != null) {
					int dim = dimOf(cur, axis);
					cur = cur.projectAxis(axis, rng.nextInt(dim));
					cost = objective.cost(cur);
					record(best, cur, cost);
					projectsTaken++;
				}
				continue;
			}
			FlipScheme trial = cur.copy();
			boolean moved;
			if (cfg.splitProb() > 0 && rng.nextDouble() < cfg.splitProb()
					&& trial.rank() < cfg.maxRank()) {
				moved = trial.split(rng);
			} else {
				moved = FlipGraphWalk.randomFlip(trial, rng, cfg.coefCap());
			}
			if (!moved) {
				// Dry vertex (no sign-shared pair — e.g. a bud-free Strassen-class
				// scheme): escaping THROUGH a bigger format is the meta graph's
				// whole point, so prefer extend over a restart; restart only when
				// the dim/rank bounds forbid growing.
				FlipScheme.Axis axis = pickExtendable(cur, rng, cfg);
				if (axis != null) {
					cur = cur.extendAxis(axis);
					extendsTaken++;
				} else {
					restarts++;
					cur = origin.copy();
				}
				cost = objective.cost(cur);
				record(best, cur, cost);
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
				record(best, cur, cost);
			}
			if (cfg.progressEvery() > 0 && step > 0 && step % cfg.progressEvery() == 0) {
				log.info("[progress] {}/{} steps (shape=⟨{},{},{}⟩ rank={} cost={} shapes={}"
						+ " accepted={} extends={} projects={} restarts={}) {}ms elapsed",
						step, cfg.maxSteps(), cur.n, cur.m, cur.p, cur.rank(), cost, best.size(),
						accepted, extendsTaken, projectsTaken, restarts,
						System.currentTimeMillis() - start);
			}
		}
		return new Result(best, step, accepted, extendsTaken, projectsTaken, restarts);
	}

	private static void record(Map<String, ShapeBest> best, FlipScheme s, long cost) {
		String key = s.n + "x" + s.m + "x" + s.p;
		ShapeBest prev = best.get(key);
		if (prev == null || cost < prev.cost()) {
			best.put(key, new ShapeBest(s.copy(), cost));
			log.info("new best at ⟨{},{},{}⟩: cost={} rank={}", s.n, s.m, s.p, cost, s.rank());
		}
	}

	private static FlipScheme.Axis pickExtendable(FlipScheme s, Random rng, Config cfg) {
		FlipScheme.Axis[] axes = FlipScheme.Axis.values();
		int start = rng.nextInt(3);
		for (int d = 0; d < 3; d++) {
			FlipScheme.Axis ax = axes[(start + d) % 3];
			int added = switch (ax) {
				case N -> s.m * s.p;
				case M -> s.n * s.p;
				case P -> s.n * s.m;
			};
			if (dimOf(s, ax) < cfg.maxDim() && s.rank() + added <= cfg.maxRank()) {
				return ax;
			}
		}
		return null;
	}

	private static FlipScheme.Axis pickProjectable(FlipScheme s, Random rng, Config cfg) {
		FlipScheme.Axis[] axes = FlipScheme.Axis.values();
		int start = rng.nextInt(3);
		for (int d = 0; d < 3; d++) {
			FlipScheme.Axis ax = axes[(start + d) % 3];
			if (dimOf(s, ax) > cfg.minDim()) {
				return ax;
			}
		}
		return null;
	}

	private static int dimOf(FlipScheme s, FlipScheme.Axis ax) {
		return switch (ax) { case N -> s.n; case M -> s.m; case P -> s.p; };
	}
}
