package eu.solven.matmul.search;

/**
 * Stopping budget for a branch-and-bound allocation/assignment search,
 * bundling the three knobs that used to be loose method parameters:
 *
 * <ul>
 *   <li>{@code upperBound} — a known achievable rank (e.g. the Kronecker rank);
 *       the search prunes anything that cannot beat it and can drop a base
 *       outright when even its root lower bound meets it.</li>
 *   <li>{@code maxNodes} — absolute cap on allocations visited (anytime).</li>
 *   <li>{@code maxNodesWithoutImprovement} — <b>stagnation</b> cap: stop once
 *       this many nodes pass with no incumbent improvement. Keeps searching
 *       while progress is being made, bails when it spins.</li>
 * </ul>
 *
 * <p>Hitting either node cap returns the best allocation found so far (an upper
 * bound; the result is flagged non-exhaustive).</p>
 */
public record SearchBudget(long upperBound, long maxNodes, long maxNodesWithoutImprovement) {

	/** Unbounded, exact search (no upper bound, no node caps). */
	public static final SearchBudget EXACT =
			new SearchBudget(Long.MAX_VALUE, Long.MAX_VALUE, Long.MAX_VALUE);

	/** Exact search pruned by a known achievable {@code upperBound}. */
	public static SearchBudget upTo(long upperBound) {
		return new SearchBudget(upperBound, Long.MAX_VALUE, Long.MAX_VALUE);
	}

	public SearchBudget withUpperBound(long ub) {
		return new SearchBudget(ub, maxNodes, maxNodesWithoutImprovement);
	}

	public SearchBudget withMaxNodes(long mn) {
		return new SearchBudget(upperBound, mn, maxNodesWithoutImprovement);
	}

	public SearchBudget withStagnation(long stagnation) {
		return new SearchBudget(upperBound, maxNodes, stagnation);
	}
}
