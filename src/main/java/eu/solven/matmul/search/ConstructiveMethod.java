package eu.solven.matmul.search;

import java.util.Optional;

import eu.solven.matmul.NonCubicBilinearAlgorithm;
import eu.solven.matmul.catalog.Recombination.SotaResolver;

/**
 * Unified interface for every closed-form / formula-driven matmul
 * construction in the catalog (Strassen 1969, Hopcroft-Kerr 1971,
 * Waksman 1970, Rosowski 2019, Pan TA, Sedoglavic 2017 Prop 1, …).
 *
 * <h2>Two tiers</h2>
 *
 * <ol>
 *   <li><strong>Predictors</strong>: return only the rank upper bound the
 *       method gives at a target shape (when applicable). Cheap; queried
 *       by the search to pick the min candidate. Implements {@link #predict}.</li>
 *   <li><strong>Constructors</strong>: additionally materialise the actual
 *       bilinear algorithm. Optional; implements {@link #construct}. Some
 *       methods (e.g. {@code SedoglavicProp1} which composes catalog atoms)
 *       only have a predictor; the constructor falls back to recursive
 *       materialisation through {@link LineageReplayer}.</li>
 * </ol>
 *
 * <h2>Why this exists</h2>
 *
 * <p>Before this interface, every published construction lived in its own
 * package ({@code papers.strassen1969}, {@code papers.hopcroftkerr1971},
 * {@code papers.rosowski2019}, …) with bespoke entry-point signatures.
 * Wiring them all into the search required one ad-hoc hook per method.
 * After this interface, {@link MethodCatalog#all()} returns the full
 * registry; the search iterates uniformly.</p>
 *
 * <h2>Pattern</h2>
 *
 * <ul>
 *   <li>Each paper / family gets ONE class implementing {@code ConstructiveMethod}.</li>
 *   <li>The class delegates to the existing paper-specific implementation —
 *       no rewriting of underlying constructors.</li>
 *   <li>The class is registered in {@link MethodCatalog#all()}.</li>
 *   <li>{@link BlockSplitSearch#findBestStrategy} iterates the registry
 *       and adds each {@link Prediction} as a candidate strategy.</li>
 * </ul>
 *
 * @see MethodCatalog
 * @see SedoglavicProp1Method  initial adopter
 */
public interface ConstructiveMethod {

	/** Short identifier, e.g. {@code "SedoglavicProp1"}, {@code "HopcroftKerr1971"}. */
	String name();

	/** Bibliographic reference, e.g. {@code "Sedoglavic 2017 hal-01572046v2"}. */
	String paperRef();

	/**
	 * Predict the rank at {@code ⟨n,m,p⟩} (in axis-permutation orbit) using
	 * this method. Returns empty when the method doesn't apply (wrong
	 * shape family, missing sub-atoms, etc.).
	 *
	 * @param n,m,p target shape
	 * @param sota  catalog rank lookup
	 * @return the prediction with predicted rank + lineage formula, or
	 *         empty
	 */
	Optional<Prediction> predict(int n, int m, int p, SotaResolver sota);

	/**
	 * Optionally materialise the explicit bilinear algorithm. Default:
	 * {@code Optional.empty()} — most predictor-only methods don't ship
	 * their own materialiser; the recursive lineage replayer handles it
	 * by composing the sub-atoms named in the prediction's lineage.
	 *
	 * @param n,m,p target shape
	 * @param atoms catalog of on-disk schemes
	 * @return the algorithm, or empty
	 */
	default Optional<NonCubicBilinearAlgorithm> construct(int n, int m, int p,
			eu.solven.matmul.catalog.Recombination.AlgorithmLookup atoms) {
		return Optional.empty();
	}

	/**
	 * Output of {@link #predict}. Carries the rank, the constructive
	 * lineage (one-line canonical-key form, ready to embed as
	 * {@code lineage_compact} in a JSON), and a label for the search's
	 * candidate list.
	 *
	 * <p>{@code verified} = this method can actually <em>construct</em> a
	 * Verifier-passing scheme at the predicted rank (or it is backed by an
	 * on-disk scheme). When {@code false} the rank is a <em>theoretical /
	 * formula-only</em> bound (e.g. the Hopcroft-Kerr ⟨2,m,n⟩ closed form for
	 * parities/ranges we have not ported a construction for): it may be shown
	 * as a cited bound, but the search MUST NOT rely on it — for small matrices
	 * we require an actual realisable scheme, not an unproven extrapolation.
	 * (2026-06-04: HK1971 predicted ⟨2,10,15⟩=233, one below the buildable
	 * Perminov 234, which the materialiser could not realise.)</p>
	 */
	record Prediction(long predictedRank, String label, String lineageCompact,
			ConstructiveMethod method, boolean verified) {
		/** Back-compat: a prediction the method can construct (verified=true). */
		public Prediction(long predictedRank, String label, String lineageCompact, ConstructiveMethod method) {
			this(predictedRank, label, lineageCompact, method, true);
		}
	}
}
