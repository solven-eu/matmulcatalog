package eu.solven.matmul.papers.dis2009;

import eu.solven.matmul.NonCubicBilinearAlgorithm;

/**
 * Pan/Islam trilinear-aggregation (TA) cubic matmul: closed-form upper
 * bound {@link #cubicBound(int)} plus an explicit factor-matrix
 * constructor {@link #build(int)} that emits a {@link NonCubicBilinearAlgorithm}
 * with rank matching the formula. Verified against the Magma reference
 * implementation Islam shipped with DIS09 (recovered from the Wayback
 * Machine and archived at
 * {@code references/islam2009/magma/TA.mgm}) and against a sympy port
 * of that Magma at {@code references/islam2009/sympy/04c_*.py} and
 * {@code 05_*.py}.
 *
 * <p>The formulas are:</p>
 * <ul>
 *   <li>{@code n} even: {@code (n³ + 12n² + 11n) / 3}</li>
 *   <li>{@code n} odd:  {@code (n³ + 15n² + 14n − 6) / 3}</li>
 * </ul>
 *
 * <h2>Attribution</h2>
 *
 * <p>The TA <em>framework</em> is Pan's (Pan 1978, FOCS; Pan 1980,
 * <em>How to Multiply Matrices Faster</em>, LNCS 179). The specific
 * <strong>closed-form formulas above are Islam 2009 MSc thesis
 * Proposition 1</strong> — both the {@code n} even formula and (Islam's
 * introduction of) the {@code n} odd formula. The same material is
 * republished as DIS09 (Drevet–Islam–Schost 2009) §3, which is the
 * easier-to-cite shared paper. Whenever a third party writes "DIS09 §3"
 * about the TA bound, the underlying result is Islam's thesis.</p>
 *
 * <p>DIS09's introduction states the best <em>previous</em> result
 * (what Islam's Prop.~1 improves on), for {@code n} even:</p>
 *
 * <pre>{@code R(n,n,n) ≤ n³/3 + (1/3)·min{ 12n² + 17n , (45/4)n² + 32n + 27 }   (n even)}</pre>
 *
 * <p>The two branches are different families, each best in its own range:
 * the first is the earlier {@code 12n²} family; the second is
 * Hadas–Schwartz 1982 — exactly {@link #panHadasSchwartz1982Bound} =
 * {@code (4n³+45n²+128n+108)/12} — whose smaller {@code (45/4)n²}
 * quadratic wins for large {@code n}. Islam 2009 sharpens the first
 * branch's linear term ({@code 17n → 11n}), giving {@link #cubicBound};
 * Schwartz–Zwecher 2025 ({@link #schwartzZwecher2025Bound}) sharpens the
 * second. Crucially {@link #bestPanTaBound} returns the <strong>min over
 * all</strong> of these, so the {@code 12n²} (Islam) branch is used for
 * even {@code n ≤ ~26} and the {@code (45/4)n²} (Schwartz–Zwecher) branch
 * takes over for even {@code n ≥ 28} — the large-{@code n} regime is NOT
 * left on the table.</p>
 *
 * <p>Islam writes verbatim: <em>"we do not know of any previous mention
 * of the case of n odd"</em>. So the odd-case formula in our code is a
 * 2009 result, not a 1978 result.</p>
 *
 * <p>This is one of several cases where an attribution chain
 * (TA technique → Pan 1978 → Pan 1992 → Islam 2009 thesis → DIS09
 * publication → modern citations just saying "Pan TA") elides who
 * actually proved what. The thesis-to-paper transition with Drevet and
 * Schost as co-authors is normal academic flow; the cost is that the
 * thesis (single-author, less indexed) tends to fall out of citation,
 * while DIS09 keeps Islam as second-author of his own underlying work.
 * See <code>references/case_studies/makarov_1986_recovery.md</code> for
 * the same pattern playing out for Makarov 1986.</p>
 *
 * <h2>Practical range</h2>
 *
 * <p>These formulas exactly reproduce the "TA" column of DIS09 Table 3
 * (n=18→3306, n=23→6806, n=29→12468). They dominate naïve
 * Strassen-recursion only for {@code n ≥ ~18}; for smaller {@code n},
 * recursive Strassen with modern small-matmul bases is better.</p>
 *
 * <h2>Field discipline</h2>
 *
 * <p>Applies to <em>non-commutative</em> rings (R/Q/Z/C/F₂ all qualify).
 * The construction is genuinely non-commutative — does NOT exploit
 * commutativity — so it lifts to recursive matmul.</p>
 */
public final class PanTrilinearAggregation {

	private PanTrilinearAggregation() {}

	/**
	 * Closed-form Pan-style TA upper bound for {@code R(⟨n,n,n⟩)} over
	 * non-commutative rings. The formulas are <strong>Islam 2009
	 * Proposition 1</strong> (the odd case is Islam's introduction);
	 * see class-level Javadoc for the full attribution chain.
	 *
	 * @param n side length (must be ≥ 1)
	 * @return the upper bound from the appropriate parity branch
	 */
	public static long cubicBound(int n) {
		if (n < 1) throw new IllegalArgumentException("n must be ≥ 1, got " + n);
		long n3 = (long) n * n * n;
		long n2 = (long) n * n;
		if (n % 2 == 0) {
			// (n³ + 12n² + 11n) / 3
			return (n3 + 12L * n2 + 11L * n) / 3L;
		}
		// (n³ + 15n² + 14n − 6) / 3
		return (n3 + 15L * n2 + 14L * n - 6L) / 3L;
	}

	/**
	 * Which parity-branch formula was used (for attribution).
	 */
	/**
	 * Construct the explicit Pan-TA bilinear algorithm for cubic
	 * {@code ⟨n, n, n⟩} matrix multiplication.
	 *
	 * <p>Port of Islam's {@code TA.mgm} reference Magma implementation
	 * (Wayback snapshot in
	 * {@code references/islam2009/magma/TA.mgm}), verified against
	 * sympy in {@code references/islam2009/sympy/04c_*.py}
	 * (even case) and {@code 05_*.py} (odd case): both check the
	 * exact multiplication count matches {@link #cubicBound(int)}
	 * AND the construction symbolically equals direct {@code A·B}.</p>
	 *
	 * @param n input matrix side length (≥ 2)
	 * @return a {@link NonCubicBilinearAlgorithm} for {@code ⟨n,n,n⟩}
	 *         with rank {@link #cubicBound(int)} and rational
	 *         coefficients (denominators of {@code n/2 + 1} from the
	 *         zero-sum padding R matrix).
	 */
	public static NonCubicBilinearAlgorithm build(int n) {
		if (n < 2) throw new IllegalArgumentException("n must be ≥ 2, got " + n);
		return n % 2 == 0
				? PanTrilinearAggregationBuilder.buildEven(n)
				: PanTrilinearAggregationBuilder.buildOdd(n);
	}

	public static String branchLabel(int n) {
		// Even formula goes back to Pan 1992 (in Islam 2009's reading);
		// odd formula is Islam 2009 Prop 1 — credit accordingly.
		return n % 2 == 0
				? "Pan TA (n even, Islam 2009 Prop 1)"
				: "Pan TA (n odd, Islam 2009 Prop 1 — introduced)";
	}

	/**
	 * Pan 1980 (SIAM J. Comput. 9(2), <em>Strassen's algorithm is not
	 * optimal: trilinear technique of aggregating, uniting and cancelling
	 * for constructing fast algorithms for matrix operations</em>)
	 * closed-form upper bound on {@code R(⟨n,n,n⟩)} for {@code n} even.
	 *
	 * <p>Formula: {@code (n³ + 9n²/2 − 3n) / 2} (equivalently
	 * {@code (2n³ + 9n² − 6n) / 4}).</p>
	 *
	 * <p>Construction lives over Q-strict (divides by {@code n+1}); lifts
	 * to R/C/any characteristic-0 ring. NOT valid over F₂.</p>
	 *
	 * <p>Returns {@code -1} when {@code n} is odd (formula doesn't apply).</p>
	 *
	 * @param n side length (must be ≥ 1)
	 * @return the Pan 1980 bound, or {@code -1} if {@code n} is odd
	 */
	public static long panSiam1980Bound(int n) {
		if (n < 1) throw new IllegalArgumentException("n must be ≥ 1, got " + n);
		if (n % 2 != 0) return -1L;
		long n3 = (long) n * n * n;
		long n2 = (long) n * n;
		// (2n³ + 9n² − 6n) / 4
		return (2L * n3 + 9L * n2 - 6L * n) / 4L;
	}

	/**
	 * Pan 1982 / Hadas–Schwartz formalisation closed-form upper bound on
	 * {@code R(⟨n,n,n⟩)} for {@code n} even, {@code n ≠ 16}.
	 *
	 * <p>Formula: {@code n³/3 + (15/4)n² + (32/3)n + 9} (equivalently
	 * {@code (4n³ + 45n² + 128n + 108) / 12}).</p>
	 *
	 * <p>The {@code n = 16} exclusion comes from the construction's
	 * division by {@code (1 − 9/(n/2+1))}, which vanishes at
	 * {@code n/2+1 = 9}, i.e. {@code n = 16}. Construction lives over
	 * Q-strict; lifts to R/C/any char-0 ring. NOT valid over F₂.</p>
	 *
	 * <p>Returns {@code -1} when {@code n} is odd or {@code n == 16}.</p>
	 *
	 * @param n side length (must be ≥ 1)
	 * @return the Hadas–Schwartz 1982 bound, or {@code -1} if not applicable
	 */
	public static long panHadasSchwartz1982Bound(int n) {
		if (n < 1) throw new IllegalArgumentException("n must be ≥ 1, got " + n);
		if (n % 2 != 0) return -1L;
		if (n == 16) return -1L;
		long n3 = (long) n * n * n;
		long n2 = (long) n * n;
		// (4n³ + 45n² + 128n + 108) / 12
		return (4L * n3 + 45L * n2 + 128L * n + 108L) / 12L;
	}

	/**
	 * Schwartz–Zwecher 2025 closed-form upper bound on {@code R(⟨n,n,n⟩)}
	 * (arXiv:2508.01748, Theorem 3.4 — the {@code t^New25} bound).
	 *
	 * <p>Same applicability domain as {@link #panHadasSchwartz1982Bound(int)}:
	 * {@code n} even, {@code n ≠ 16}.</p>
	 *
	 * <p>Formula: {@code n³/3 + (15/4)n² + (61/6)n + 8} (equivalently
	 * {@code (4n³ + 45n² + 122n + 96) / 12}).</p>
	 *
	 * <p>STRICTLY DOMINATES {@link #panHadasSchwartz1982Bound(int)} for
	 * all valid {@code n} — the kin-row unification reduces the linear-term
	 * constant from {@code 32/3} to {@code 61/6} (i.e. {@code 64/6 → 61/6})
	 * and the absolute constant from {@code 9} to {@code 8}.</p>
	 *
	 * <p>Construction lives over Q-strict; lifts to R/C/any char-0 ring.
	 * NOT valid over F₂.</p>
	 *
	 * <p>Returns {@code -1} when {@code n} is odd or {@code n == 16}.</p>
	 *
	 * @param n side length (must be ≥ 1)
	 * @return the Schwartz–Zwecher 2025 bound, or {@code -1} if not applicable
	 */
	public static long schwartzZwecher2025Bound(int n) {
		if (n < 1) throw new IllegalArgumentException("n must be ≥ 1, got " + n);
		if (n % 2 != 0) return -1L;
		if (n == 16) return -1L;
		long n3 = (long) n * n * n;
		long n2 = (long) n * n;
		// (4n³ + 45n² + 122n + 96) / 12
		return (4L * n3 + 45L * n2 + 122L * n + 96L) / 12L;
	}

	/**
	 * Returns the MIN over all applicable Pan-TA closed-form formulas for
	 * a given {@code n}.
	 *
	 * <p>Dispatches by parity (Islam's odd-{@code n} formula in
	 * {@link #cubicBound(int)} is the only odd-{@code n} formula) and
	 * skips formulas whose exclusions apply (parity / {@code n = 16} /
	 * {@code n < 2}). For practical {@code n} (4..32) the MIN is usually
	 * Islam-even 2009 for small even {@code n}, Schwartz–Zwecher 2025 for
	 * large even {@code n}, and Islam-odd 2009 for any odd {@code n}.</p>
	 *
	 * <p>Returns {@code -1} if no Pan-TA formula applies (shouldn't happen
	 * for {@code n ≥ 2}).</p>
	 *
	 * @param n side length (must be ≥ 1)
	 * @return the minimum applicable Pan-TA bound
	 */
	public static long bestPanTaBound(int n) {
		if (n < 1) throw new IllegalArgumentException("n must be ≥ 1, got " + n);
		long best = -1L;
		// Islam 2009 (Prop 1) covers all n ≥ 1, both parities.
		long islam = cubicBound(n);
		if (islam > 0) best = islam;
		long pan80 = panSiam1980Bound(n);
		if (pan80 > 0 && (best < 0 || pan80 < best)) best = pan80;
		long hs82 = panHadasSchwartz1982Bound(n);
		if (hs82 > 0 && (best < 0 || hs82 < best)) best = hs82;
		long sz25 = schwartzZwecher2025Bound(n);
		if (sz25 > 0 && (best < 0 || sz25 < best)) best = sz25;
		return best;
	}
}
