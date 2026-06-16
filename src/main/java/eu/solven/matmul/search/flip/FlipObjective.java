package eu.solven.matmul.search.flip;

/**
 * Cost function a {@link FlipGraphWalk} MINIMIZES. Unlike the published flip
 * engines (Kauers–Moosbauer, Moosbauer–Poole, Perminov), which hard-code rank
 * descent, our walks target catalog structure metrics — bud-richness for the
 * upward serendipitous product, projection margin for the downward operator —
 * with rank as the leading lexicographic key so a rank reduction is never
 * traded away for structure.
 *
 * <p>Every walk result is a <strong>bound</strong> (best-found over a heuristic
 * walk), never an optimum — see the optimality discipline in CLAUDE.md.</p>
 */
public interface FlipObjective {

	long cost(FlipScheme s);

	/** Human label for logs/reports, e.g. {@code "rank↓ then budScore↑"}. */
	String describe();
}
