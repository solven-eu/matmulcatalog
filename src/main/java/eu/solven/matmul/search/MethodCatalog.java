package eu.solven.matmul.search;

import java.util.List;
import java.util.Optional;

import eu.solven.matmul.papers.hopcroftkerr1971.HopcroftKerr1971Method;
import eu.solven.matmul.papers.waksman1970.Waksman1970Method;
import eu.solven.matmul.recombination.Recombination.SotaResolver;

/**
 * Registry of every {@link ConstructiveMethod} the catalog knows about.
 *
 * <h2>Current contents</h2>
 *
 * <ul>
 *   <li>{@link HopcroftKerr1971Method} — Hopcroft-Kerr 1971 {@code ⟨2,m,n⟩}
 *       closed form.</li>
 *   <li>{@link PanTrilinearAggregationMethod} — Pan 1980 / DIS09 §3 cubic
 *       odd/even closed form.</li>
 * </ul>
 *
 * <p>Note: Sedoglavic 2017 Proposition 1 is <em>not</em> a registered method.
 * Its bound — {@code ⟨u+v,u+v,u+v⟩ ≤ ⟨u,u,u⟩ + 3·⟨u,u,v⟩ + 3·⟨v,v,u⟩} — is the
 * recombination multiset of the Strassen {@code ⟨2,2,2⟩} base at a two-part
 * {@code [u,v]} block allocation, so it is already produced (and dominated when a
 * better split exists) by the generic recombination path in
 * {@link BlockSplitSearch#findBestStrategy}; its {@code u=v} doubling is the
 * {@link PairFusedRecombination} Pan-TA candidate. A standalone enumerator would
 * only re-derive what recombination already covers.</p>
 *
 * <h2>Adding a method</h2>
 *
 * <ol>
 *   <li>Implement {@link ConstructiveMethod} in the matching
 *       {@code papers.{author}{year}} sub-package (or {@code search/} if
 *       it's a general-purpose enumerator).</li>
 *   <li>Add a {@code new YourMethod()} entry to {@link #all()}.</li>
 *   <li>That's it — {@link BlockSplitSearch#findBestStrategy} picks it up
 *       on the next call.</li>
 * </ol>
 *
 * <h2>Future adopters (planned)</h2>
 *
 * <ul>
 *   <li>{@code HopcroftKerr1971Method} — wraps {@code HopcroftKerr2bcAsymmetric}
 *       for {@code ⟨2,b,c⟩} targets.</li>
 *   <li>{@code Waksman1970Method} — wraps {@code Waksman1970} for
 *       commutative cubic + {@code ⟨n,3,3⟩} + {@code ⟨2,2,n⟩}.</li>
 *   <li>{@code Rosowski2019Thm2Method} — wraps {@code RosowskiAlgorithm1}
 *       for non-bilinear commutative.</li>
 *   <li>{@code PanTrilinearAggregationMethod} — wraps
 *       {@code PanTrilinearAggregation.cubicBound} for cubic odd-n.</li>
 *   <li>{@code KnownTauIdentitiesMethod} — wraps the hand-extracted
 *       identity registry (Schönhage τ-style recipes from FMM-Lille).</li>
 * </ul>
 */
public final class MethodCatalog {

	private MethodCatalog() {}

	/**
	 * All registered <strong>non-commutative</strong> constructive methods.
	 * Iteration order ≈ method age (oldest paper first).
	 *
	 * <p>Commutative-only methods ({@link Waksman1970Method},
	 * Rosowski 2019 Thm 2/3, Makarov 1986, …) are <em>not</em>
	 * registered here — they don't compose under recursive matmul
	 * over non-commutative rings. A separate commutative-aware lookup
	 * will register them (see #65).</p>
	 */
	public static List<ConstructiveMethod> all() {
		return List.of(
				new HopcroftKerr1971Method(),
				new PanTrilinearAggregationMethod());
	}

	/**
	 * Predictions from every method that applies at {@code ⟨n,m,p⟩}, in
	 * registry order. Methods that don't apply (out-of-family, missing
	 * sub-atoms) are skipped silently.
	 */
	public static List<ConstructiveMethod.Prediction> predictAll(int n, int m, int p, SotaResolver sota) {
		List<ConstructiveMethod.Prediction> out = new java.util.ArrayList<>();
		for (ConstructiveMethod m1 : all()) {
			Optional<ConstructiveMethod.Prediction> pred = m1.predict(n, m, p, sota);
			pred.ifPresent(out::add);
		}
		return out;
	}
}
