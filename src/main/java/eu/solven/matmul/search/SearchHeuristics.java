package eu.solven.matmul.search;

/**
 * Single home for the <em>computation heuristics</em> — the bounded /
 * anytime knobs that trade search exhaustiveness for wall-clock. These are
 * NOT correctness parameters: every value here is a cost cap, so hitting one
 * yields an <strong>upper bound</strong> (best found so far), never a proven
 * optimum. Keep them here rather than scattered as loose magic numbers so a
 * single edit retunes the whole branch-and-bound pipeline.
 *
 * <p>Distinct from {@link eu.solven.matmul.catalog.CatalogLimits}, which holds
 * the on-disk JSON-format / dimension constants (max dim, materialise cap,
 * sparse threshold). This file is about <em>how hard we search</em>; that file
 * is about <em>how we store the result</em>.</p>
 *
 * <p>Related structured carriers (left where they are — they are types, not
 * loose constants): {@link SearchBudget} bundles the per-search
 * {@code (upperBound, maxNodes, maxNodesWithoutImprovement)} triple, and
 * {@link PoolConfig} carries the pool-shape knobs with their own
 * {@code UNBOUNDED_*} / {@code NO_PADDING} sentinels.</p>
 *
 * @see SearchBudget
 * @see PoolConfig
 */
public final class SearchHeuristics {

	private SearchHeuristics() {}

	/**
	 * Default anytime cap on allocations visited by the per-base
	 * <em>assignment</em> branch-and-bound ({@link AssignmentOptimizer}).
	 * Used to initialise {@link BlockSplitSearch#ASSIGNMENT_MAX_NODES}, which
	 * stays runtime-tunable (e.g. {@code SchemeSweep} can raise it). Hitting
	 * the cap returns the best assignment found (a bound, flagged
	 * non-exhaustive).
	 */
	public static final long DEFAULT_ASSIGNMENT_MAX_NODES = 20_000_000L;

	/**
	 * Default anytime cap on nodes visited by the <em>allocation</em>
	 * branch-and-bound ({@link AllocationOptimizer}); the default used by the
	 * {@code AllocationDeepDive} driver when none is passed on the command
	 * line. Looser than {@link #DEFAULT_ASSIGNMENT_MAX_NODES} because the
	 * allocation layer is the outer, cheaper-per-node loop.
	 */
	public static final long DEFAULT_ALLOCATION_MAX_NODES = 50_000_000L;
}
