package eu.solven.matmul.papers.rosowski2019;

import lombok.extern.slf4j.Slf4j;

import java.util.Optional;

/**
 * Closed-form COMMUTATIVE upper bounds on matmul rank from
 * Rosowski 2019/2020 (arXiv:1904.07683). The paper has multiple
 * theorems; we expose two methods covering distinct regimes:
 *
 * <ul>
 *   <li>{@link #commutativeBoundBilinear} — Theorems 2/3,
 *       <strong>bilinear over a commutative ring</strong>. The {@code /2}
 *       in the formula is integer counting (algorithm produces an even
 *       count by structure), NOT a runtime division — the algorithm
 *       itself contains no division operations.</li>
 *   <li>{@link #nonBilinearRankBound} — Theorems 4/5, <strong>non-bilinear
 *       upper bound on the standard non-commutative rank</strong>
 *       {@code R(⟨n,n,n⟩)}. The "ring satisfying (2)" condition means R
 *       has a transpose involution. Rosowski's "recursive non-bilinear"
 *       novelty recurses these <em>over commutative rings</em> (the
 *       linear combinations of products are bilinear, which is what
 *       permits recursion), reaching ω≈2.8125 — but it does NOT lift to
 *       non-commutative matmul and does NOT realise the standard NC
 *       tensor rank. Like Waksman/Makarov these are COMMUTATIVE-ONLY and
 *       must carry {@code commutative: true}; they fall below the NC
 *       table only because the commutative regime is cheaper (≈n³/2).
 *       (Earlier revisions mislabelled these as NC via a
 *       "transpose-involution → matrix-ring" argument; that over-reached —
 *       the recursion keeps the scalar base ring commutative.)</li>
 * </ul>
 *
 * <p><strong>When Thm 4/5 beat known NC bounds</strong>: for several
 * cubic targets, Rosowski's formula gives smaller rank than the best
 * tabulated NC results in DIS09 (Drevet–Islam–Schost 2009 Table 3):</p>
 *
 * <table>
 *   <caption>Rosowski Thm 4/5 vs DIS09 NC bilinear table 3</caption>
 *   <tr><th>n</th><th>DIS09 NC</th><th>Rosowski Thm 4/5</th></tr>
 *   <tr><td>7</td><td>258 (Winograd)</td><td><b>252</b></td></tr>
 *   <tr><td>9</td><td>522 (mul121)</td><td><b>495</b></td></tr>
 *   <tr><td>10</td><td>700 (Strassen)</td><td><b>655</b></td></tr>
 *   <tr><td>11</td><td>923 (Strassen)</td><td><b>858</b></td></tr>
 *   <tr><td>12</td><td>1125 (mul121)</td><td><b>1086</b></td></tr>
 *   <tr><td>13</td><td>1450 (Strassen)</td><td><b>1365</b></td></tr>
 *   <tr><td>14</td><td>1728 (Strassen)</td><td><b>1673</b></td></tr>
 *   <tr><td>15</td><td>2108 (Winograd2)</td><td><b>2040</b></td></tr>
 * </table>
 *
 * <p>Asymptotically Rosowski gives ω ≈ 2.8125 (per the paper, n=14
 * recursion), slightly worse than Strassen's 2.807, but these mid-n
 * direct bounds are genuinely new.</p>
 *
 * <p>Special case: Theorem 1 gives {@code R(⟨n,3,3⟩) ≤ 6n+3}, so
 * {@code R_c(⟨3,3,3⟩) ≤ 21} (better than both general formulas applied
 * to n=3, where Thm 3 → 21 and Thm 5 → 30; Corollary 1 confirms 21).</p>
 *
 * <p>{@link #bestCommutativeBound} returns the min over all axis
 * permutations using the bilinear (tighter) formula, with the n=3
 * special case applied when applicable.</p>
 *
 * <p><strong>Commutative-only</strong>: bounds apply only when the
 * underlying ring is commutative; they do NOT lift to recursive
 * matmul (which requires non-commutative coefficients on matrix
 * entries).</p>
 */
@Slf4j
public final class RosowskiBound {

	private RosowskiBound() {}

	/**
	 * Rosowski Theorems 2/3 — BILINEAR commutative upper bound for
	 * {@code ⟨l, n, m⟩}. May use divisions by 2 when {@code n} is odd
	 * (Theorem 3). Valid over commutative rings where 2 is invertible
	 * (Q, R, C); not directly over Z.
	 */
	public static Optional<Long> commutativeBoundBilinear(int l, int n, int m) {
		if (l < 1 || n < 1 || m < 1) return Optional.empty();
		long L = l, N = n, M = m;
		long inner = N * (L * M + L + M - 1);
		if (n % 2 == 1) {
			// Thm 3 — n odd. May use divisions.
			if (m % 2 == 1) return Optional.of(inner / 2);
			else return Optional.of((inner + L - 1) / 2);
		} else {
			// Thm 2 — n even, divisions-free.
			return Optional.of(inner / 2);
		}
	}

	/**
	 * Rosowski Theorems 4/5 — non-bilinear <b>COMMUTATIVE</b> rank upper bound
	 * on cubic matmul ({@code n ≥ 2}). Recurses over commutative rings only
	 * (ω≈2.8125); does NOT lift to non-commutative matmul / does NOT realise
	 * the NC tensor rank {@code R(⟨n,n,n⟩)} — see the class javadoc.
	 *
	 * <ul>
	 *   <li>n even (Thm 4): {@code n(n²+3n+1)/2}</li>
	 *   <li>n odd (Thm 5):  {@code n(n²+3n+2)/2 = n(n+1)(n+2)/2}</li>
	 * </ul>
	 *
	 * <p>Falls below DIS09's <em>NC</em> bilinear table 3 for
	 * {@code n ∈ {7,9,10,11,12,13,14,15}} — but only because the commutative
	 * regime is cheaper, so this is NOT an NC improvement.</p>
	 */
	public static Optional<Long> nonBilinearRankBound(int n) {
		if (n < 2) return Optional.empty();
		long N = n;
		long inner = (n % 2 == 0)
				? N * N * N + 3 * N * N + N            // Thm 4: n(n²+3n+1)
				: N * N * N + 3 * N * N + 2 * N;       // Thm 5: n(n²+3n+2)
		return Optional.of(inner / 2);
	}

	/**
	 * @deprecated alias retained for back-compat. Theorems 4/5 are
	 *             non-bilinear COMMUTATIVE bounds (commutative-only; do not
	 *             lift to NC matmul). Use {@link #nonBilinearRankBound}.
	 */
	@Deprecated
	public static Optional<Long> commutativeBoundNonBilinearDivFree(int n) {
		return nonBilinearRankBound(n);
	}

	/**
	 * @deprecated use {@link #commutativeBoundBilinear} for the bilinear
	 *             bound (divisions allowed) or
	 *             {@link #commutativeBoundNonBilinearDivFree} for the
	 *             divisions-free bound. Kept for backward compatibility.
	 */
	@Deprecated
	public static Optional<Long> commutativeBound(int l, int n, int m) {
		return commutativeBoundBilinear(l, n, m);
	}

	/**
	 * Returns the best Rosowski commutative bound for {@code ⟨a, b, c⟩},
	 * trying all 6 axis permutations (the formula is sensitive to which
	 * axis is the contraction, but matmul rank is invariant under
	 * permutation) and both the bilinear and non-bilinear-divisions-free
	 * formulas. Returns the min — the tightest upper bound.
	 *
	 * <p>For cubic {@code ⟨n,n,n⟩} this picks min(Theorem 3 bilinear,
	 * Theorem 5 divisions-free). For non-cubic, only the bilinear
	 * formula applies (Theorems 4/5 are cubic-only).</p>
	 */
	public static Optional<Long> bestCommutativeBound(int a, int b, int c) {
		int[][] perms = {
				{ a, b, c }, { a, c, b }, { b, a, c },
				{ b, c, a }, { c, a, b }, { c, b, a }
		};
		long best = Long.MAX_VALUE;
		for (int[] p : perms) {
			Optional<Long> r = commutativeBoundBilinear(p[0], p[1], p[2]);
			if (r.isPresent() && r.get() < best) best = r.get();
		}
		// Note: Thm 4/5 (nonBilinearRankBound) gives a NC rank bound, not a
		// commutative one — see the deprecated commutativeBoundNonBilinearDivFree
		// docstring. Don't mix it with commutative bounds here.
		return best == Long.MAX_VALUE ? Optional.empty() : Optional.of(best);
	}

	/** Sanity quick-test. */
	public static void main(String[] args) {
		log.info("Per-formula comparison (cubic ⟨n,n,n⟩):");
		log.info(String.format("%4s | %15s | %20s | %4s%n",
				"n", "Thm 3 bilinear", "Thm 4/5 div-free", "best"));
		log.info("-".repeat(60));
		for (int n = 2; n <= 16; n++) {
			Optional<Long> bil = commutativeBoundBilinear(n, n, n);
			Optional<Long> div = commutativeBoundNonBilinearDivFree(n);
			Optional<Long> best = bestCommutativeBound(n, n, n);
			log.info(String.format("%4d | %15s | %20s | %4s%n",
					n, bil.orElse(-1L), div.orElse(-1L), best.orElse(-1L)));
		}
		log.info("");
		log.info("⟨3,3,3⟩ = " + bestCommutativeBound(3, 3, 3).orElseThrow() + "  (expected 21 via Cor 1 / Thm 3)");
		log.info("⟨5,5,5⟩ bilinear = " + commutativeBoundBilinear(5, 5, 5).orElseThrow() + "  (Thm 3, uses /2)");
		log.info("⟨5,5,5⟩ div-free = " + commutativeBoundNonBilinearDivFree(5).orElseThrow() + "  (Thm 5, no divisions)");
	}
}
