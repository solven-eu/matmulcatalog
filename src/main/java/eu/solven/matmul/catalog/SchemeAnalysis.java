package eu.solven.matmul.catalog;

import java.util.List;
import java.util.Map;

import eu.solven.matmul.NonCubicBilinearAlgorithm;
import eu.solven.matmul.Verifier;
import eu.solven.matmul.algebra.Field;

/**
 * An analysis that requires the fully-expanded factor matrices of a scheme —
 * the <em>costly</em> side of the catalog pipeline (#159). The catalog splits
 * into two op classes:
 *
 * <ul>
 *   <li><b>Metadata-only</b> (lineage + rank + shape): rank lookup,
 *       {@code atom}-ness, the FMM/Perminov comparison, dominance/shaving.
 *       No expansion — runs in the fast <em>detect</em> phase.</li>
 *   <li><b>Expanded-scheme</b> (needs the explicit {@code U/V/W}): symbolic
 *       verification, addition count, bud structure, field narrowing, … These
 *       implement {@code SchemeAnalysis}.</li>
 * </ul>
 *
 * <p>For {@code maxDim > 16} the matrices are not stored, so obtaining them
 * means replaying the lineage — expensive. The <em>validate-and-enrich</em>
 * (Phase 2) batch therefore expands each scheme <strong>once</strong> and runs
 * every registered {@code SchemeAnalysis} on it, stamping the results into the
 * scheme JSON. Implementing this interface <em>is</em> the declaration that an
 * operation needs expansion and belongs in Phase 2, never in the detect path
 * (which has no expanded scheme to hand it).</p>
 *
 * <p>{@link #analyse} returns the JSON key/value pairs to stamp, or
 * <strong>throws</strong> to signal the scheme is invalid (the batch deletes
 * it). The {@link Verify} analysis additionally returns {@code verified:false}
 * for a scheme that expands but does not compute matmul.</p>
 */
public interface SchemeAnalysis {

	/** Short identifier (logging / selective runs). */
	String name();

	/**
	 * Analyse the expanded scheme over algebra {@code field}; return the JSON
	 * fields to stamp into the scheme. Throw to mark the scheme INVALID (→
	 * pruned by the validate-and-enrich batch).
	 */
	Map<String, Object> analyse(NonCubicBilinearAlgorithm alg, Field field);

	/**
	 * The standard Phase-2 analysis set, in run order: {@link Verify} first so a
	 * scheme that fails verification is pruned before the (pointless) enrichment
	 * analyses run on it.
	 */
	static List<SchemeAnalysis> defaults() {
		return List.of(new Verify(), new Additions(), new Buds(), new ProjectionMargin(), new Ternary());
	}

	/**
	 * Verification → {@code verified}. Field-sensitive: F₂/F₃ use the exact
	 * modular checks; characteristic-0 uses the random matmul spot-check
	 * ({@link Verifier#passesRandomMatmulSpotCheck}), which is the operational
	 * bar the materialiser already applies and stays tractable for large
	 * (rank&nbsp;≫&nbsp;1000) composed schemes where full symbolic exactness is
	 * costly. {@code verified:false} ⇒ the batch deletes the scheme.
	 */
	final class Verify implements SchemeAnalysis {
		@Override
		public String name() {
			return "verify";
		}

		@Override
		public Map<String, Object> analyse(NonCubicBilinearAlgorithm alg, Field field) {
			boolean ok = switch (field) {
				case F2 -> Verifier.isExactNonCubicF2(alg);
				case F3 -> Verifier.isExactNonCubicF3(alg);
				default -> Verifier.passesRandomMatmulSpotCheck(alg);
			};
			return Map.of("verified", ok);
		}
	}

	/** Flat addition count → {@code additions}. Field-agnostic. */
	final class Additions implements SchemeAnalysis {
		@Override
		public String name() {
			return "additions";
		}

		@Override
		public Map<String, Object> analyse(NonCubicBilinearAlgorithm alg, Field field) {
			return Map.of("additions", Verifier.additionCount(alg));
		}
	}

	/**
	 * Bud (U/V/W shared-vector) structure → {@code has_buds} (+ {@code buds}
	 * summary when present). Needs the explicit u/v/w vectors to find
	 * proportional groups, so it is inherently an expanded-scheme op (this is
	 * why {@code maxDim>16} stubs carry no buds until validated).
	 */
	final class Buds implements SchemeAnalysis {
		@Override
		public String name() {
			return "buds";
		}

		@Override
		public Map<String, Object> analyse(NonCubicBilinearAlgorithm alg, Field field) {
			SerendipitousBudProduct.BudSummary b = SerendipitousBudProduct.summarise(alg);
			if (!b.hasBuds()) {
				return Map.of("has_buds", false);
			}
			return Map.of("has_buds", true, "buds", b.summary());
		}
	}

	/**
	 * Projection margin μ → {@code projection_margin}: the max products dead-code-
	 * eliminated when the single best index (over all three axes) is dropped, so
	 * {@code R_after = R − μ}. A high μ marks a strong <em>downward</em> (projection)
	 * parent even at higher rank — the (rank, μ) Pareto axis, dual to bud-richness
	 * for the (upward) serendipitous product. Needs the explicit U/V/W (per-product
	 * index supports), so it is computed in the SAME single-expansion pass as
	 * {@link Buds} — see {@link ProjectionSearch#projectionMargin}.
	 */
	final class ProjectionMargin implements SchemeAnalysis {
		@Override
		public String name() {
			return "projection_margin";
		}

		@Override
		public Map<String, Object> analyse(NonCubicBilinearAlgorithm alg, Field field) {
			return Map.of("projection_margin", ProjectionSearch.projectionMargin(alg));
		}
	}

	/**
	 * Ternary-integer sub-class flag → {@code zt}: true ⟺ every U/V/W coefficient
	 * is in {@code {-1,0,1}} ({@link SchemeIO#isTernary}). This is Perminov's
	 * "ternary integer" marker — a SUB-CLASS of {@code Z}, NOT a field (it has
	 * nothing to do with F₂/Z₂ characteristic 2, nor F₃'s ternary <i>modular</i>).
	 *
	 * <p>Gated on {@code field == Field.Z}: the enrichment driver's {@code fieldOf}
	 * returns {@code Z} exactly when {@code Z ∈ fields[]} (Z has top priority), so
	 * this stamps {@code zt} for — and only for — integer schemes, matching the
	 * {@code MaterialiseZT} pass. For a non-Z scheme the flag is meaningless, so we
	 * emit nothing (empty map) and leave any stale value to {@code MaterialiseZT}'s
	 * clear path. Computing it HERE — from the dense {@code alg} already expanded
	 * for the bud/projection metrics — means new schemes (notably the whole
	 * {@code derived/} tree, previously 0/7784 stamped) get {@code zt} in the same
	 * single-expansion enrichment pass, instead of only via a separate stamp run.</p>
	 */
	final class Ternary implements SchemeAnalysis {
		@Override
		public String name() {
			return "ternary";
		}

		@Override
		public Map<String, Object> analyse(NonCubicBilinearAlgorithm alg, Field field) {
			if (field != Field.Z) {
				return Map.of();
			}
			return Map.of("zt", SchemeIO.isTernary(alg));
		}
	}
}
