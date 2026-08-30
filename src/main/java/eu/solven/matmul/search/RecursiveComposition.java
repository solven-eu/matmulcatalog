package eu.solven.matmul.search;

import java.util.ArrayList;
import java.util.List;
import java.util.Optional;

import eu.solven.matmul.catalog.KnownAlgorithm;
import eu.solven.matmul.catalog.KnownAlgorithmCatalog;

/**
 * Compute upper bounds on matmul rank by composing smaller catalogued algorithms.
 *
 * <p><b>Theorem (folklore)</b>: if
 * {@code ⟨n, m, p⟩ = ⟨n_1·n_2, m_1·m_2, p_1·p_2⟩}, then
 * {@code R(⟨n, m, p⟩) ≤ R(⟨n_1, m_1, p_1⟩) · R(⟨n_2, m_2, p_2⟩)}, provided
 * both sub-algorithms are valid over the matrix-entry algebra (this is why
 * commutative-only algorithms don't recurse).</p>
 *
 * <p>This class evaluates such compositions over the
 * {@link KnownAlgorithmCatalog} — see {@code SMALL_MATMUL_CATALOG.md} §6.</p>
 */
public final class RecursiveComposition {

	private RecursiveComposition() {}

	/**
	 * Multiply best-known ranks down a chain of nested factorizations.
	 *
	 * <p>Each {@link Factor} represents one level of the recursive split. The
	 * returned rank is the product of best-known ranks at each level (or
	 * {@link Optional#empty()} if any level has no catalog entry).</p>
	 *
	 * <p>Example: for {@code ⟨32,32,32⟩} as {@code ⟨16,16,16⟩ × ⟨2,2,2⟩}, where
	 * {@code ⟨16,16,16⟩ = ⟨4,4,4⟩ × ⟨4,4,4⟩} (two levels of AlphaTensor),
	 * pass three {@code Factor} entries: {@code (2,2,2)}, {@code (4,4,4)},
	 * {@code (4,4,4)}.</p>
	 *
	 * @param algebra arithmetic setting — every catalog entry used must be valid here.
	 * @param factors ordered list of factor blocks (order doesn't affect the
	 *                product but is recorded in the breakdown).
	 */
	public static Optional<Result> evaluate(eu.solven.matmul.algebra.Algebra algebra, List<Factor> factors) {
		// Block recursion substitutes MATRICES for scalars, and matrices don't
		// commute — so a commutative-only algorithm cannot be a level of this
		// product (the theorem in the class javadoc requires each sub-algorithm
		// to be valid over the matrix-entry algebra). Nothing used to stop it:
		// with the old commutative ⟨2,2,2⟩=6 row this happily returned 6^k as a
		// bound for ⟨2^k,2^k,2^k⟩ — wrong twice over. Refuse the whole algebra
		// rather than hope the current curated rows happen to be ring-valid.
		if (algebra.commutative()) {
			throw new IllegalArgumentException(
					"RecursiveComposition is only valid for non-commutative algebras — "
					+ "commutative-only schemes (Waksman, Rosowski, Makarov) do not lift to "
					+ "block recursion. Got: " + algebra);
		}
		long product = 1;
		List<String> breakdown = new ArrayList<>();
		for (Factor f : factors) {
			Optional<KnownAlgorithm> best =
					KnownAlgorithmCatalog.bestKnown(f.n, f.m, f.p, algebra);
			if (best.isEmpty()) {
				return Optional.empty();
			}
			product *= best.get().rank;
			breakdown.add(String.format("R(⟨%d,%d,%d⟩/%s) = %d (%s %d)",
					f.n, f.m, f.p, algebra, best.get().rank,
					best.get().source, best.get().year));
		}
		return Optional.of(new Result(product, breakdown));
	}

	/** Convenience: rank of {@code n} levels of the same factor (Strassen, Laderman, ...). */
	public static Optional<Result> evaluatePower(eu.solven.matmul.algebra.Algebra algebra, Factor base, int levels) {
		List<Factor> factors = new ArrayList<>();
		for (int i = 0; i < levels; i++) factors.add(base);
		return evaluate(algebra, factors);
	}

	public static final class Factor {
		public final int n, m, p;

		public Factor(int n, int m, int p) {
			this.n = n;
			this.m = m;
			this.p = p;
		}

		public static Factor cube(int n) { return new Factor(n, n, n); }
	}

	public static final class Result {
		public final long rank;
		public final List<String> breakdown;

		public Result(long rank, List<String> breakdown) {
			this.rank = rank;
			this.breakdown = breakdown;
		}

		@Override
		public String toString() {
			return String.format("rank=%d via [%s]", rank, String.join(" × ", breakdown));
		}
	}
}
