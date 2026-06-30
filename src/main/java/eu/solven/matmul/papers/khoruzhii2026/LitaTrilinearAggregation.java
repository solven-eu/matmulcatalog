package eu.solven.matmul.papers.khoruzhii2026;

/**
 * Khoruzhii, Gelß &amp; Pokutta 2026 — <em>Local Improvements to Trilinear
 * Aggregation</em> (LITA) — closed-form NON-COMMUTATIVE rank upper bound for
 * the cubic matmul tensor {@code ⟨N,N,N⟩}.
 *
 * <p>LITA refines Pan-style trilinear aggregation with local re-pairings,
 * lowering the multiplicative constant. The per-parity closed forms are
 * (rank = number of active products):</p>
 *
 * <ul>
 *   <li><b>N even:</b> {@code Rₑ(N) = N³/3 + 15N²/4 + 29N/3 + 7}
 *       {@code = (4N³ + 45N² + 116N + 84) / 12}</li>
 *   <li><b>N odd:</b>  {@code Rₒ(N) = (4N³ + 57N² + 14N − 15)/12 − ⌊3(N−1)/8⌋}</li>
 * </ul>
 *
 * <p>Both branches are integer-valued and the construction is defined only for
 * {@code N ≥ 19} (see {@link #MIN_N}). The cubic constant beats dense / recursive
 * schemes from there up — for odd N it already wins at N=19, for even N around
 * N≥26 (below that a direct scheme wins and the LITA value, where defined, is a
 * correct-but-dominated upper bound). These are the values FMM-Lille's index
 * reports for the large cubic formats (e.g. {@code ⟨21,21,21⟩ = 5198},
 * {@code ⟨23,23,23⟩ = 6586}, {@code ⟨25,25,25⟩ = 8196}).</p>
 *
 * <p><b>Field.</b> The construction is RATIONAL — its factor matrices carry
 * denominators divisible by 12 (e.g. coefficients ±1/12 for ⟨21³⟩). So it is
 * Q-native and valid over {@code Q ⊂ R ⊂ C}, but does NOT reduce to F2 or F3
 * (gcd(12,2)=2, gcd(12,3)=3). It is non-commutative (lifts to recursive matmul,
 * so it bounds ω over those fields).</p>
 *
 * <p>Source: Kirill Khoruzhii, Patrick Gelß, Sebastian Pokutta,
 * <em>Local Improvements to Trilinear Aggregation</em>, 2026,
 * <a href="https://github.com/khoruzhii/lita">github.com/khoruzhii/lita</a>.</p>
 */
public final class LitaTrilinearAggregation {

	private LitaTrilinearAggregation() {}

	/** Smallest dimension the LITA construction is defined for ({@code N > 18}). */
	public static final int MIN_N = 19;

	/**
	 * LITA cubic rank for {@code ⟨N,N,N⟩}, exact integer arithmetic.
	 *
	 * <p>Only defined for {@code N ≥ 19} ({@link #MIN_N}): the construction has no
	 * realization below that (the Maple generators reject {@code N ≤ 18}), so the
	 * closed form must not be evaluated — and certainly not claimed as a bound —
	 * for smaller N. Conveniently, {@code N=19} is also where the trilinear
	 * constant first beats our catalog, so nothing is lost.</p>
	 *
	 * @param n cube side, {@code n ≥ 19}
	 * @return the LITA upper bound on the rank of {@code ⟨n,n,n⟩}
	 * @throws IllegalArgumentException if {@code n < 19}
	 */
	public static long cubicRank(int n) {
		if (n < MIN_N) {
			throw new IllegalArgumentException(
					"LITA is defined only for N ≥ " + MIN_N + " (N > 18), got " + n);
		}
		long N = n;
		long num;
		long sub;
		if ((n & 1) == 0) {
			// Even: (4N³ + 45N² + 116N + 84) / 12
			num = 4 * N * N * N + 45 * N * N + 116 * N + 84;
			sub = 0;
		} else {
			// Odd: (4N³ + 57N² + 14N − 15)/12 − ⌊3(N−1)/8⌋
			num = 4 * N * N * N + 57 * N * N + 14 * N - 15;
			sub = (3 * (N - 1)) / 8;
		}
		if (num % 12 != 0) {
			// The closed forms are integer-valued by construction; a non-exact
			// division means the formula was mis-transcribed — fail loud.
			throw new IllegalStateException(
					"LITA numerator not divisible by 12 for n=" + n + " (num=" + num + ")");
		}
		return num / 12 - sub;
	}
}
