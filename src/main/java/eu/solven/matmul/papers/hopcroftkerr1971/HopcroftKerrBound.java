package eu.solven.matmul.papers.hopcroftkerr1971;

/**
 * Hopcroft-Kerr 1971 closed-form NON-COMMUTATIVE upper bound for
 * matmul where one dimension equals 2 — DIS09 Table 2 row 3.
 *
 * <p>For {@code ⟨a, 2, c⟩}: {@code R ≤ (3ac + max(a,c)) / 2}.</p>
 *
 * <p>By the {@code S₃} symmetry of the matmul tensor, the same bound
 * applies to {@code ⟨2, b, c⟩} (any permutation that puts a {@code 2}
 * in one slot): the rank is invariant under cyclic and transpose
 * symmetries.</p>
 *
 * <p>Reference: J.E. Hopcroft and L.R. Kerr, "On minimizing the
 * number of multiplications necessary for matrix multiplication",
 * SIAM J. Appl. Math. 20 (1971), pp. 30–36, DOI 10.1137/0120004.
 * DIS09 §1 (Table 2) cites this as the canonical NC bound for
 * matrices with a dim-2 axis.</p>
 *
 * <p>Field discipline: non-commutative — the construction lifts to
 * recursive matmul over any ring.</p>
 */
public final class HopcroftKerrBound {

	private HopcroftKerrBound() {}

	/**
	 * Non-commutative Hopcroft-Kerr bound for a shape with at least
	 * one dim equal to 2. Returns {@code -1} if no axis is 2.
	 *
	 * @param a first axis
	 * @param b second axis
	 * @param c third axis
	 * @return upper bound on {@code R(⟨a,b,c⟩)} if one of {a,b,c}=2, else -1
	 */
	public static long forShape(int a, int b, int c) {
		if (a < 1 || b < 1 || c < 1) throw new IllegalArgumentException("dims must be ≥ 1");
		// Rotate so the middle (b) carries the 2, then apply the formula.
		// Pick the permutation that minimises the result — all rotations
		// give the same rank by tensor symmetry, but the formula is
		// written for a specific axis layout.
		if (b == 2) return canonical(a, c);
		if (a == 2) return canonical(b, c);
		if (c == 2) return canonical(a, b);
		return -1;
	}

	/**
	 * Formula for the canonical {@code ⟨x, 2, y⟩} arrangement.
	 * The numerator {@code 3xy + max(x,y)} can be odd (e.g. x=2, y=3 → 21);
	 * DIS09 Table 2 row 3 specifies {@code ⌈(3xy + max(x,y))/2⌉}, so
	 * round up rather than truncate.
	 */
	private static long canonical(long x, long y) {
		long num = 3L * x * y + Math.max(x, y);
		return (num + 1L) / 2L;
	}
}
