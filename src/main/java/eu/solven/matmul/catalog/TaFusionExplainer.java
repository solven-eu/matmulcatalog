package eu.solven.matmul.catalog;

import java.util.Optional;
import java.util.regex.Matcher;
import java.util.regex.Pattern;

import eu.solven.matmul.NonCubicBilinearAlgorithm;
import eu.solven.matmul.recombination.Recombination;
import eu.solven.matmul.recombination.Recombination.TaFusionBreakdown;

/**
 * Turns a {@link Lineage.RecombinationTaN} node into a {@link TaFusionBreakdown} so the
 * Pan-TA fusion inside a recombination can be HIGHLIGHTED — both as a one-line summary
 * for displays/logs and as a structured object stamped into {@code catalog.json}. This
 * is what makes a recombination's final rank explainable: it names which products TA
 * fused and how many multiplications that bought.
 *
 * <p>Pan TA is a saving WITHIN the recombination's multiplications, not a separate
 * strategy — so the highlight lives on the recombination node, recomputed from
 * {@code base + allocations} exactly as the build/replay does (no stored pair list to
 * drift). A non-TA recombination (block-combining base, or no profitable pair) yields
 * {@link Optional#empty()}.</p>
 */
public final class TaFusionExplainer {

	private TaFusionExplainer() {}

	/** {@code "naive-NxMxP"} — the base ref a {@link Lineage.RecombinationTaN} pins. */
	private static final Pattern NAIVE_REF = Pattern.compile("naive-(\\d+)x(\\d+)x(\\d+)");

	/**
	 * Breakdown for a TA-fused recombination node, valuing leaves at {@code sota}. Returns
	 * empty when the node's base is not a recoverable naïve grid (so the caller surfaces
	 * "no TA" rather than a fake zero-saving entry).
	 */
	public static Optional<TaFusionBreakdown> describe(
			Lineage.RecombinationTaN node, Recombination.SotaResolver sota) {
		NonCubicBilinearAlgorithm base = naiveBaseOf(node.base());
		if (base == null) {
			return Optional.empty();
		}
		TaFusionBreakdown bd = Recombination.describeTaFusion(
				base, sota, node.allocA(), node.allocB(), node.allocC());
		return Optional.ofNullable(bd);
	}

	/** The naïve grid a {@code "naive-NxMxP"} Atom names, or {@code null} if the base node
	 *  isn't such an Atom (we only highlight TA on the generic naïve-grid recombination). */
	private static NonCubicBilinearAlgorithm naiveBaseOf(Lineage.Node baseNode) {
		if (!(baseNode instanceof Lineage.Atom atom)) {
			return null;
		}
		Matcher mtch = NAIVE_REF.matcher(atom.ref());
		if (!mtch.matches()) {
			return null;
		}
		return NonCubicBilinearAlgorithm.naive(
				Integer.parseInt(mtch.group(1)),
				Integer.parseInt(mtch.group(2)),
				Integer.parseInt(mtch.group(3)));
	}
}
