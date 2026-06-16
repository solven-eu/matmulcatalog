package eu.solven.matmul.algebra;

/**
 * A fully-specified matmul rank question: <em>what is</em>
 * (or <em>what's the best upper bound on</em>) the rank of the
 * matmul tensor for an {@code n×m} times {@code m×p} product, in a
 * given {@link Algebra}?
 *
 * <p>Every API in the project that talks about "the rank of a matmul"
 * should accept a {@link MatmulProblem} rather than juggling
 * {@code (n, m, p, field, commutative)} as loose arguments — that's
 * how the field-confusion bug (claiming F₂ ranks over R/Q/Z) crept
 * in originally.</p>
 *
 * <p>Format is stored in declared order. Sorted/canonical access via
 * {@link #sortedShape()} is provided because matmul-tensor S₃ symmetry
 * means rank is the same for any permutation of {@code (n,m,p)}.</p>
 */
public record MatmulProblem(Algebra algebra, int n, int m, int p) {

	public MatmulProblem {
		if (algebra == null) throw new IllegalArgumentException("algebra is required");
		if (n < 1 || m < 1 || p < 1) {
			throw new IllegalArgumentException(
					"shape must be positive: (" + n + "," + m + "," + p + ")");
		}
	}

	/** Convenience: non-commutative problem over a given field. */
	public static MatmulProblem nc(Field field, int n, int m, int p) {
		return new MatmulProblem(Algebra.nonCommutative(field), n, m, p);
	}

	/** Convenience: commutative-only problem over a given field. */
	public static MatmulProblem cmt(Field field, int n, int m, int p) {
		return new MatmulProblem(Algebra.commutative(field), n, m, p);
	}

	/** Returns {@code (n, m, p)} as an int array (not aliased; safe to mutate). */
	public int[] shape() {
		return new int[] { n, m, p };
	}

	/**
	 * Returns {@code (n, m, p)} sorted ascending. Used for catalog
	 * lookups since matmul tensor rank is invariant under axis
	 * permutation (S₃ symmetry); a sorted key collapses all 6 orderings.
	 */
	public int[] sortedShape() {
		int[] s = shape();
		java.util.Arrays.sort(s);
		return s;
	}

	/** Compact format string {@code ⟨n,m,p⟩}. */
	public String formatTag() {
		return "⟨" + n + "," + m + "," + p + "⟩";
	}

	/** Returns max(n, m, p) — useful for catalog section bucketing. */
	public int maxDim() {
		return Math.max(n, Math.max(m, p));
	}

	/** True iff {@code n == m == p}. */
	public boolean isCubic() {
		return n == m && m == p;
	}

	@Override
	public String toString() {
		return formatTag() + " over " + algebra;
	}
}
