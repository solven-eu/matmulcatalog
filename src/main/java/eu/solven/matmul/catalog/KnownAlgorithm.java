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
 * <p>Mirrors the row schema in {@code SMALL_MATMUL_CATALOG.md} §3.</p>
 */
public final class KnownAlgorithm {

	public final int n;
	public final int m;
	public final int p;
	public final Algebra algebra;
	/** Multiplication count (= tensor rank for the chosen algebra). */
	public final int rank;
	/** First-publication year. */
	public final int year;
	/** Short algorithm/source name ("Strassen", "Laderman", "AlphaTensor"...). */
	public final String source;
	/** Optional URL / DOI / arXiv ID for the primary reference. */
	public final String link;
	/** Free-form note ("Strassen²", "RL-discovered", "complex-specific"...). */
	public final String notes;

	public KnownAlgorithm(int n, int m, int p, Algebra algebra, int rank, int year,
			String source, String link, String notes) {
		this.n = n;
		this.m = m;
		this.p = p;
		this.algebra = algebra;
		this.rank = rank;
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
		return String.format("%s/%s rank=%d (%s %d) %s",
				formatLabel(), algebra, rank, source, year, notes == null ? "" : notes);
	}
}
