package eu.solven.matmul.catalog;

import eu.solven.matmul.algebra.Algebra;

/**
 * One known fast bilinear matmul algorithm: identifies the format and
 * algebra (field + commutativity), captures the rank attained, and
 * traces provenance (who, when, where).
 *
 * <p>Multiple entries for the same {@code (format, algebra)} are
 * expected and encouraged — they form the historical timeline of
 * improvements. {@link KnownAlgorithmCatalog#bestKnown} returns the
 * current minimum.</p>
 *
 * <p><strong>Every entry is a CONSTRUCTION</strong> — an algorithm that
 * attains {@code rank}, i.e. an upper bound. An impossibility result
 * (Winograd 1971, Hopcroft–Kerr 1971, Bläser, Wang …) is NOT an entry
 * here: it belongs in {@code docs/lower-bounds.json} / {@link LowerBoundRegistry}.
 * The catalog once carried {@code ⟨2,2,2⟩ commutative = 6} attributed to
 * Hopcroft–Kerr and Winograd — the two papers that <em>prove that rank
 * impossible</em>. A lower-bound paper was recorded as if it supplied an
 * algorithm, and nothing in the type could say otherwise; {@link Optimality}
 * plus the lower-bound cross-check ({@code TestBoundsVsLowerBounds}) is the
 * guard against a repeat.</p>
 *
 * <p>Mirrors the row schema in {@code SMALL_MATMUL_CATALOG.md} §3.</p>
 */
public final class KnownAlgorithm {

	/**
	 * Honesty tier of the {@code rank} value (CLAUDE.md optimality discipline).
	 * A construction always proves {@code R ≤ rank}; the tier says how much
	 * more than that we know.
	 */
	public enum Optimality {
		/**
		 * A published lower bound matches this rank, so it IS the optimum for
		 * this {@code (format, algebra)} — cross-checked against
		 * {@code docs/lower-bounds.json}. May be stated as "= optimum".
		 */
		PROVEN_OPTIMAL,
		/**
		 * Best construction known; no matching floor. An UPPER BOUND — never
		 * write "minimal" or "the minimum" for one of these.
		 */
		BOUND,
		/**
		 * Optimal over a restricted rule/atom set only. Unused by the curated
		 * catalog so far; kept so the repo-wide tier vocabulary stays complete.
		 */
		OPTIMAL_WITHIN_SCOPE
	}

	public final int n;
	public final int m;
	public final int p;
	public final Algebra algebra;
	/** Multiplication count (= tensor rank for the chosen algebra), attained by a construction. */
	public final int rank;
	/** Whether {@link #rank} is proven optimal or merely the best construction known. */
	public final Optimality optimality;
	/** First-publication year. */
	public final int year;
	/** Short algorithm/source name ("Strassen", "Laderman", "AlphaTensor"...). */
	public final String source;
	/** Optional URL / DOI / arXiv ID for the primary reference. */
	public final String link;
	/** Free-form note ("Strassen²", "RL-discovered", "complex-specific"...). */
	public final String notes;

	public KnownAlgorithm(int n, int m, int p, Algebra algebra, int rank, Optimality optimality,
			int year, String source, String link, String notes) {
		this.n = n;
		this.m = m;
		this.p = p;
		this.algebra = algebra;
		this.rank = rank;
		this.optimality = optimality;
		this.year = year;
		this.source = source;
		this.link = link;
		this.notes = notes;
	}

	public String formatLabel() {
		return String.format("⟨%d,%d,%d⟩", n, m, p);
	}

	@Override
	public String toString() {
		return String.format("%s/%s rank=%d [%s] (%s %d) %s",
				formatLabel(), algebra, rank,
				optimality == Optimality.PROVEN_OPTIMAL ? "optimum" : "upper bound",
				source, year, notes == null ? "" : notes);
	}
}
