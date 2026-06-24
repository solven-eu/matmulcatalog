package eu.solven.matmul.search;

import java.util.List;
import java.util.Optional;

import eu.solven.matmul.recombination.Recombination.SotaResolver;

/**
 * Parallel registry of {@link ConstructiveMethod} implementations that
 * produce <strong>commutative-only</strong> matmul algorithms (Waksman
 * 1970, Rosowski 2019 Thm 2/3, Makarov 1986, Islam 2009 generalised
 * Waksman).
 *
 * <h2>Why a separate registry</h2>
 *
 * <p>Commutative algorithms exploit scalar commutativity ({@code xy = yx})
 * and therefore <strong>do NOT lift to recursive matmul over a
 * non-commutative ring</strong> (matrix scalars). Using one as a
 * sub-product inside a Strassen-style recombination would silently
 * under-count the true rank — see the 2026-06-02 ⟨17,17,17⟩=2868
 * incident.</p>
 *
 * <p>The main {@link MethodCatalog} is therefore non-commutative-pure and
 * suitable as the SOTA resolver for any recombination. This separate
 * {@code CommutativeMethodCatalog} is consulted at the TOP LEVEL only —
 * when the caller accepts a commutative algorithm (e.g. computing a
 * single scalar matmul, not a recursive matmul over matrix scalars).</p>
 *
 * <h2>Pairs with</h2>
 *
 * <ul>
 *   <li>{@link eu.solven.matmul.catalog.FieldAwareLookup#findRankAllowCommutative}
 *       — catalog rank lookup that admits commutative entries.</li>
 *   <li>{@link eu.solven.matmul.catalog.FieldAwareLookup#findRankCommutativeOnly}
 *       — catalog rank lookup restricted to commutative entries.</li>
 * </ul>
 *
 * <h2>Current contents</h2>
 *
 * <ul>
 *   <li>{@link Waksman1970Method} — Waksman 1970 closed-form bound.</li>
 * </ul>
 *
 * <h2>Planned adopters</h2>
 *
 * <ul>
 *   <li>{@code Rosowski2019Thm2Method} — wraps
 *       {@code RosowskiAlgorithm1} / {@code RosowskiBound} for
 *       commutative {@code ⟨2,2,n⟩} and {@code ⟨3,3,n⟩} families.</li>
 *   <li>{@code Makarov1986Method} — wraps the Makarov {@code ⟨3,3,3⟩}=22
 *       commutative constructor (#149).</li>
 * </ul>
 */
public final class CommutativeMethodCatalog {

	private CommutativeMethodCatalog() {}

	/** All registered commutative constructive methods. */
	public static List<ConstructiveMethod> all() {
		return List.of(
				new Waksman1970Method());
	}

	/**
	 * Predictions from every commutative method applicable at
	 * {@code ⟨n,m,p⟩}, in registry order. Methods that don't apply
	 * (out-of-family, missing sub-atoms) are skipped silently.
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
