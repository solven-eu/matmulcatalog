package eu.solven.matmul.papers.waksman1970;

/**
 * Waksman 1970 closed-form COMMUTATIVE upper bound for {@code ⟨a,b,c⟩}
 * matmul, per DIS09 Table 2:
 *
 * <ul>
 *   <li>{@code b} even: {@code b · (ac + a + c − 1) / 2}</li>
 *   <li>{@code b} odd:  {@code (b − 1) · (ac + a + c − 1) / 2 + ac}</li>
 * </ul>
 *
 * <p>This is a <strong>commutative-only</strong> bound: the algorithm
 * exploits commutativity of the underlying scalars and does NOT lift to
 * recursive matmul over non-commutative rings (i.e., does not contribute
 * to ω). It does provide tight upper bounds for matmul of fixed-size
 * matrices with commutative entries, which is what DIS09 Table 4
 * reports.</p>
 *
 * <p>Note: the formula is asymmetric in the three axes (the "b" axis
 * is special). Use {@link #bestCubic} to take the min over axis
 * orientations for cubic ⟨n,n,n⟩ targets.</p>
 */
public final class WaksmanBound {

	private WaksmanBound() {}

	/** Waksman bound for {@code ⟨a,b,c⟩} (b is the special middle axis). */
	public static long forShape(int a, int b, int c) {
		if (a < 1 || b < 1 || c < 1) {
			throw new IllegalArgumentException("dims must be ≥ 1");
		}
		long ac = (long) a * c;
		long term = ac + a + c - 1;
		if (b % 2 == 0) {
			return (long) b * term / 2L;
		}
		return (long) (b - 1) * term / 2L + ac;
	}

	/**
	 * Min over the three axis orientations of {@code ⟨n,n,n⟩} — Waksman
	 * with "b" set to whichever axis is even (if any), else odd-branch.
	 * For cubic this collapses to a single formula evaluation since
	 * all axes are equal.
	 */
	public static long bestCubic(int n) {
		return forShape(n, n, n);
	}
}
