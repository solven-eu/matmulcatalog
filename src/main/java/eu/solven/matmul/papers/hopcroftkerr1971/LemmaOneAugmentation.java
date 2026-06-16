package eu.solven.matmul.papers.hopcroftkerr1971;

import java.util.Arrays;

/**
 * Hopcroft-Kerr 1971 Lemma 1 augmentation: given a target row count
 * {@code n} and column count {@code p}, construct an {@code n × p}
 * integer matrix {@code M} such that:
 *
 * <ol>
 *   <li>The first {@code p} rows form the identity matrix.</li>
 *   <li>Any {@code p} cyclically consecutive rows form a non-singular
 *       {@code p × p} matrix.</li>
 * </ol>
 *
 * <p>This is the augmentation table used by the paper's asymmetric
 * {@code ⟨p, 2, n⟩} construction (page 9-12): the original input
 * matrix {@code A} (of size {@code p × 2}) gets extended via
 * {@code Ā = M · A} to size {@code n × 2}. Lemma 1's property
 * guarantees that for any cyclic window of {@code p} consecutive rows
 * of {@code Ā}, the corresponding sub-matrix of {@code M} is
 * invertible — which is what allows back-substitution to recover the
 * original {@code p} output rows from any {@code p} cyclic-consecutive
 * output rows of {@code Ā · X}.</p>
 *
 * <p>Construction strategy (per the paper's incremental proof):
 * iterate the {@code n - p} augmented rows; for each, pick the
 * smallest positive integer coefficient sequence that makes every
 * relevant determinant non-zero. In practice the rows often end up
 * with small coefficients (single-digit) for the {@code p ≤ 16,
 * n ≤ 32} range we care about — verified by the unit tests.</p>
 *
 * <p>Reference: Hopcroft-Kerr 1971, page 5-6 (Lemma 1 statement and
 * inductive proof).</p>
 */
public final class LemmaOneAugmentation {

	private LemmaOneAugmentation() {}

	/**
	 * Build an {@code n × p} augmentation matrix satisfying Lemma 1.
	 *
	 * @param p the column count; must be {@code ≥ 1}
	 * @param n the row count; must satisfy {@code n ≥ p}
	 * @return an {@code n × p} integer matrix with the first {@code p}
	 *         rows being the identity and every cyclic window of
	 *         {@code p} consecutive rows non-singular
	 * @throws IllegalArgumentException for bad dimensions
	 * @throws IllegalStateException    if no small-integer solution is
	 *                                  found within the search budget
	 *                                  (should not happen for {@code n ≤ 32}
	 *                                  empirically)
	 */
	public static int[][] build(int p, int n) {
		if (p < 1 || n < p) {
			throw new IllegalArgumentException("require p ≥ 1 and n ≥ p, got p=" + p + ", n=" + n);
		}
		int[][] M = new int[n][p];
		for (int i = 0; i < p; i++) M[i][i] = 1;
		for (int newRow = p; newRow < n; newRow++) {
			pickRow(M, newRow, p, n);
		}
		return M;
	}

	/**
	 * Long-precision variant of {@link #build(int, int)}. Necessary for
	 * {@code p ≥ 11} where Vandermonde powers {@code (i+1)^(p-1)} exceed
	 * {@link Integer#MAX_VALUE} (e.g. {@code 11^10 ≈ 2.6 × 10^10}).
	 * Returns the same structural matrix; the first {@code p} rows are
	 * identity and augmented rows hold {@code M[p+i][k] = (i+1)^k}.
	 *
	 * <p>Caveat: still {@code long}-bounded, so for very large
	 * {@code (p, n)} (combined exponents past {@code 9.2 × 10^18}) this
	 * also throws. Realistic targets up to {@code p ≤ 20, n ≤ 40} fit.</p>
	 */
	public static long[][] buildLong(int p, int n) {
		if (p < 1 || n < p) {
			throw new IllegalArgumentException("require p ≥ 1 and n ≥ p, got p=" + p + ", n=" + n);
		}
		long[][] M = new long[n][p];
		for (int i = 0; i < p; i++) M[i][i] = 1;
		for (int newRow = p; newRow < n; newRow++) {
			pickRowLong(M, newRow, p, n);
		}
		return M;
	}

	/**
	 * Fill {@code M[newRow][0..p-1]} with the smallest integer
	 * coefficient assignment (lexicographic over the search order)
	 * that makes every cyclic-window determinant non-zero given the
	 * already-fixed rows {@code M[0..newRow-1]} and the future
	 * (not-yet-fixed) zero-rows from {@code newRow+1 onward}.
	 */
	private static void pickRow(int[][] M, int newRow, int p, int n) {
		// Vandermonde: M[p+i][k] = (i+1)^k for the i-th augmented row.
		// Sub-matrices of p cyclically-consecutive Vandermonde rows are
		// Vandermonde matrices with distinct evaluation points, hence
		// non-singular by classical Vandermonde-determinant identity.
		//
		// Fixing each augmented row's COLUMN VALUES this way (rather than
		// "all ones except last") avoids the p≥3 trap where successive
		// augmented rows share the first p-1 columns and force later
		// determinants to be identically zero in the free variable.
		int i = newRow - p;  // augmented-row index, 0-based starting from i=0
		long pow = 1;
		long base = i + 1L;
		for (int k = 0; k < p; k++) {
			// Overflow check: (i+1)^(p-1) ≤ n^p for our (p ≤ 16, n ≤ 32) range,
			// max ≈ 33^16 ≈ 10^24 — overflows long. Bound the coefficient by clamping
			// powers to a safe Vandermonde-equivalent that still gives invertibility.
			// (For our actual targets — n ≤ 16, p ≤ 12 — pow stays within long.)
			if (pow > Long.MAX_VALUE / Math.max(1L, base) || pow > Integer.MAX_VALUE) {
				throw new IllegalStateException("Vandermonde coefficient overflow at "
						+ "newRow=" + newRow + ", k=" + k + "; need wider-precision construction");
			}
			M[newRow][k] = (int) pow;
			pow *= base;
		}
		// Sanity check: assert the cyclic windows involving this row are non-singular.
		// (Vandermonde guarantees this for non-wrapping windows; the wrap windows are
		// the special case checked when newRow == n-1.)
		if (!allWindowsInvertible(M, newRow, p, n)) {
			throw new IllegalStateException("Vandermonde failed at newRow=" + newRow
					+ " — cyclic-wrap window singular. This shouldn't happen for n ≤ 2p.");
		}
	}

	/**
	 * Check every cyclic window of {@code p} consecutive rows in the
	 * partial matrix {@code M[0..newRow]} for non-singularity. The
	 * undefined rows {@code newRow+1 .. n-1} are treated as zero and
	 * only windows that don't cross them (i.e. windows whose largest
	 * index is {@code ≤ newRow}) are checked. Cyclic windows that
	 * would wrap into the undefined region are skipped at this stage —
	 * they'll be checked when those rows are filled.
	 */
	private static boolean allWindowsInvertible(int[][] M, int newRow, int p, int n) {
		for (int start = 0; start <= newRow - p + 1 && start + p - 1 <= newRow; start++) {
			int end = start + p - 1;
			if (end != newRow) continue;  // only check windows whose last row is the just-set one
			if (det(M, start, p) == 0L) return false;
		}
		// Once newRow == n-1 (the last row), also check the cyclic windows
		// that wrap past the end into rows 0, 1, ..., p-2.
		if (newRow == n - 1) {
			for (int wrap = 1; wrap < p; wrap++) {
				int[][] window = new int[p][p];
				for (int r = 0; r < p; r++) {
					int sourceRow = (n - p + wrap + r) % n;
					System.arraycopy(M[sourceRow], 0, window[r], 0, p);
				}
				if (detMatrix(window) == 0L) return false;
			}
		}
		return true;
	}

	/** Long-precision row picker — mirrors {@link #pickRow} but uses long arithmetic. */
	private static void pickRowLong(long[][] M, int newRow, int p, int n) {
		int i = newRow - p;
		long pow = 1;
		long base = i + 1L;
		for (int k = 0; k < p; k++) {
			if (pow > Long.MAX_VALUE / Math.max(1L, base)) {
				throw new IllegalStateException("Vandermonde coefficient long-overflow at "
						+ "newRow=" + newRow + ", k=" + k + "; need BigInteger");
			}
			M[newRow][k] = pow;
			pow *= base;
		}
		if (!allWindowsInvertibleLong(M, newRow, p, n)) {
			throw new IllegalStateException("Vandermonde failed at newRow=" + newRow
					+ " (long path) — cyclic-wrap window singular");
		}
	}

	private static boolean allWindowsInvertibleLong(long[][] M, int newRow, int p, int n) {
		for (int start = 0; start <= newRow - p + 1 && start + p - 1 <= newRow; start++) {
			int end = start + p - 1;
			if (end != newRow) continue;
			long[][] T = new long[p][p];
			for (int r = 0; r < p; r++) System.arraycopy(M[start + r], 0, T[r], 0, p);
			if (bareissDet(T, p) == 0L) return false;
		}
		if (newRow == n - 1) {
			for (int wrap = 1; wrap < p; wrap++) {
				long[][] window = new long[p][p];
				for (int r = 0; r < p; r++) {
					int sourceRow = (n - p + wrap + r) % n;
					System.arraycopy(M[sourceRow], 0, window[r], 0, p);
				}
				if (bareissDet(window, p) == 0L) return false;
			}
		}
		return true;
	}

	/**
	 * Compute the determinant (as a {@code long}) of the {@code p × p}
	 * sub-matrix {@code M[start..start+p-1]}. Uses Gaussian
	 * elimination with integer arithmetic via Bareiss algorithm —
	 * avoids floating-point altogether to detect singularity exactly.
	 */
	static long det(int[][] M, int start, int p) {
		long[][] T = new long[p][p];
		for (int r = 0; r < p; r++)
			for (int c = 0; c < p; c++)
				T[r][c] = M[start + r][c];
		return bareissDet(T, p);
	}

	static long detMatrix(int[][] window) {
		int p = window.length;
		long[][] T = new long[p][p];
		for (int r = 0; r < p; r++)
			for (int c = 0; c < p; c++)
				T[r][c] = window[r][c];
		return bareissDet(T, p);
	}

	/**
	 * Small-coefficient Lemma-1 matrix with EXACT verification — the numerically
	 * robust replacement for the Vandermonde rows of {@link #buildLong} (task #7).
	 *
	 * <p>Why: Vandermonde augmented rows ({@code (i+1)^k}, up to {@code ~10⁹})
	 * (a) overflow the {@code long} Bareiss determinant in the wrap-window check —
	 * the spurious "cyclic-wrap window singular" failures — and (b) make the
	 * per-column window inverses so ill-conditioned that the double back-sub in
	 * {@code HopcroftKerr2bcAsymmetric.buildOdd} loses exactness from
	 * {@code n ≳ p+2}. Small entries fix both: window determinants stay tiny,
	 * the exact-fraction inverses have small numerators/denominators, and the
	 * final double W is exact to verifier tolerance.</p>
	 *
	 * <p>Construction: augmented rows drawn from {@code {-2..2}} by a
	 * deterministic splitmix64 stream (seeded by {@code (p, n, attempt)});
	 * every completed window is checked nonsingular with BigInteger Bareiss
	 * (no overflow), wrap windows included; on failure, retry with the next
	 * attempt seed. Nonsingularity is codim-1, so a few attempts suffice.</p>
	 */
	public static long[][] buildSmallLong(int p, int n) {
		if (p < 1 || n < p) {
			throw new IllegalArgumentException("require p ≥ 1 and n ≥ p, got p=" + p + ", n=" + n);
		}
		for (int attempt = 0; attempt < 400; attempt++) {
			long[][] M = new long[n][p];
			for (int i = 0; i < p; i++) M[i][i] = 1;
			long seed = 0x9E3779B97F4A7C15L * (((long) p << 32) ^ ((long) n << 8) ^ attempt) + 1;
			// NONZERO-TERNARY-FIRST (task #8): {-1,1} rows keep window
			// determinants — hence the exact back-substitution denominators —
			// small, AND satisfy the wrap-window requirement that every entry be
			// nonzero (for n−p ∈ {1,2} the wrap dets reduce to single entries /
			// 2×2 minors of the augmented rows, so a zero anywhere is fatal —
			// zero-allowing draws almost never succeed at large p). Widen to
			// {-2,-1,1,2} if the signs alone cannot break the singularities.
			boolean wide = attempt >= 300;
			boolean ok = true;
			for (int newRow = p; newRow < n && ok; newRow++) {
				boolean rowOk = false;
				// The LAST row's draw gets full wrap-window feedback (otherwise the
				// p−1 wrap determinants are only tested once per outer attempt and
				// the success probability collapses at large p).
				boolean last = (newRow == n - 1);
				for (int rowTry = 0; rowTry < 500 && !rowOk; rowTry++) {
					if (n - p == 2 && newRow == p + 1) {
						// CONSTRUCTIVE second row: the wrap windows omit consecutive
						// identity pairs (j, j+1); their dets are the 2×2 minors
						// a_j·b_{j+1} − a_{j+1}·b_j. The sign chain
						// b_{j+1} = −b_j·a_j·a_{j+1} makes every such minor ±2 ≠ 0.
						long[] a = M[p];
						M[newRow][0] = 1;
						for (int k = 1; k < p; k++) {
							M[newRow][k] = -M[newRow][k - 1] * a[k - 1] * a[k];
						}
					} else {
						for (int k = 0; k < p; k++) {
							seed = splitmix64(seed);
							if (wide) {
								int v = (int) Math.floorMod(seed, 4); // {0,1,2,3}
								M[newRow][k] = (v < 2) ? (v == 0 ? -1 : 1) : (v == 2 ? -2 : 2);
							} else {
								M[newRow][k] = (Math.floorMod(seed, 2) == 0) ? -1 : 1;
							}
						}
					}
					rowOk = windowsEndingAtNonsingular(M, newRow, p)
							&& (!last || wrapWindowsNonsingular(M, p, n));
				}
				ok = rowOk;
			}
			if (ok && wrapWindowsNonsingular(M, p, n)) {
				return M;
			}
		}
		throw new IllegalStateException("buildSmallLong: no small-coefficient Lemma-1 matrix found "
				+ "for p=" + p + ", n=" + n + " after 400 attempts (unexpected — codim-1 failure)");
	}

	private static long splitmix64(long z) {
		z += 0x9E3779B97F4A7C15L;
		z = (z ^ (z >>> 30)) * 0xBF58476D1CE4E5B9L;
		z = (z ^ (z >>> 27)) * 0x94D049BB133111EBL;
		return z ^ (z >>> 31);
	}

	/**
	 * <strong>All-unimodular</strong> Lemma-1 matrix (task #11): every cyclic
	 * {@code p}-window determinant is {@code ±1}, so the per-column window
	 * inverses — hence the back-substituted W — are INTEGER and the constructed
	 * schemes upgrade from Q to Z (gaining F₂/F₃ membership wholesale).
	 *
	 * <p>Window dets of {@code [I; B]} reduce to three minor families of the
	 * augmented block {@code B} (m = n−p rows): leading k×k (cols 0..k−1),
	 * trailing w×w (rows m−w.., cols p−w..), and sliding m×m over contiguous
	 * column windows. Construction — <em>Euclidean recursion</em>: a period-m
	 * comb body ({@code B[u][c] = 1 ⟺ c ≡ u (mod m)} on the first
	 * {@code p − (p mod m)} columns; sliding windows are then permutation
	 * matrices and leading minors identity blocks), and for {@code r = p mod
	 * m > 0} a tail of r columns equal to {@code B(m, r)}<sup>T</sup>:
	 * eliminating the comb's permutation part shows the crossing/trailing
	 * window minors of B are exactly the leading/trailing/sliding minors of
	 * the tail block transposed — the SAME problem at size {@code (m, r)}.
	 * The recursion descends like gcd and terminates at {@code r = 0} (pure
	 * comb). Pure 0/1 entries, no search; verified exhaustively over the band
	 * range in {@code references/hopcroftkerr1971/sympy/unimodular/euclidean.py}.</p>
	 *
	 * <p>Returns {@code null} if the (belt-and-braces) exact verification of
	 * the recursion's output fails — the caller falls back to
	 * {@link #buildSmallLong}, keeping that shape over Q (loud, never silent).</p>
	 */
	public static long[][] buildUnimodular(int p, int n) {
		if (p < 1 || n < p) {
			throw new IllegalArgumentException("require p ≥ 1 and n ≥ p, got p=" + p + ", n=" + n);
		}
		int m = n - p;
		long[][] M = new long[n][p];
		for (int i = 0; i < p; i++) M[i][i] = 1;
		if (m == 0) return M;
		long[][] B = euclideanBlock(p, m);
		for (int u = 0; u < m; u++) {
			System.arraycopy(B[u], 0, M[p + u], 0, p);
		}
		return allWindowsUnimodular(M, p, n) ? M : null;
	}

	/**
	 * The m×p augmented block of the Euclidean construction: comb body +
	 * transposed recursive tail. Terminates at {@code p % m == 0}.
	 */
	private static long[][] euclideanBlock(int p, int m) {
		int r = p % m;
		int body = p - r;
		long[][] B = new long[m][p];
		for (int u = 0; u < m; u++) {
			for (int c = u; c < body; c += m) {
				B[u][c] = 1;
			}
		}
		if (r > 0) {
			long[][] C = euclideanBlock(m, r); // r × m
			for (int u = 0; u < m; u++) {
				for (int j = 0; j < r; j++) {
					B[u][body + j] = C[j][u];
				}
			}
		}
		return B;
	}

	private static boolean isUnit(java.math.BigInteger d) {
		return d.abs().equals(java.math.BigInteger.ONE);
	}

	/** Exact check: every cyclic p-window of M has determinant ±1. */
	static boolean allWindowsUnimodular(long[][] M, int p, int n) {
		for (int j = 0; j < n; j++) {
			java.math.BigInteger[][] T = new java.math.BigInteger[p][p];
			for (int t = 0; t < p; t++)
				for (int c = 0; c < p; c++)
					T[t][c] = java.math.BigInteger.valueOf(M[(j + t) % n][c]);
			if (!isUnit(bigDet(T, p))) return false;
		}
		return true;
	}

	private static boolean windowsEndingAtNonsingular(long[][] M, int newRow, int p) {
		int start = newRow - p + 1;
		if (start < 0) return true; // window not complete yet
		java.math.BigInteger[][] T = new java.math.BigInteger[p][p];
		for (int r = 0; r < p; r++)
			for (int c = 0; c < p; c++)
				T[r][c] = java.math.BigInteger.valueOf(M[start + r][c]);
		return bigDet(T, p).signum() != 0;
	}

	private static boolean wrapWindowsNonsingular(long[][] M, int p, int n) {
		for (int wrap = 1; wrap < p; wrap++) {
			java.math.BigInteger[][] T = new java.math.BigInteger[p][p];
			for (int r = 0; r < p; r++) {
				int src = (n - p + wrap + r) % n;
				for (int c = 0; c < p; c++) T[r][c] = java.math.BigInteger.valueOf(M[src][c]);
			}
			if (bigDet(T, p).signum() == 0) return false;
		}
		return true;
	}

	/** Fraction-free Bareiss determinant in BigInteger — exact at any size. Mutates T. */
	static java.math.BigInteger bigDet(java.math.BigInteger[][] T, int p) {
		java.math.BigInteger prev = java.math.BigInteger.ONE;
		int sign = 1;
		for (int k = 0; k < p; k++) {
			if (T[k][k].signum() == 0) {
				int swap = -1;
				for (int r = k + 1; r < p; r++) if (T[r][k].signum() != 0) { swap = r; break; }
				if (swap < 0) return java.math.BigInteger.ZERO;
				java.math.BigInteger[] tmp = T[k]; T[k] = T[swap]; T[swap] = tmp;
				sign = -sign;
			}
			for (int r = k + 1; r < p; r++) {
				for (int c = k + 1; c < p; c++) {
					T[r][c] = T[k][k].multiply(T[r][c]).subtract(T[r][k].multiply(T[k][c])).divide(prev);
				}
				T[r][k] = java.math.BigInteger.ZERO;
			}
			prev = T[k][k];
		}
		return sign > 0 ? T[p - 1][p - 1] : T[p - 1][p - 1].negate();
	}

	/**
	 * Bareiss-style fraction-free determinant on a {@code p × p} long
	 * matrix. Mutates {@code T}.
	 */
	private static long bareissDet(long[][] T, int p) {
		long prev = 1;
		long sign = 1;
		for (int k = 0; k < p; k++) {
			// Partial-pivot for non-zero diagonal.
			if (T[k][k] == 0) {
				int swap = -1;
				for (int r = k + 1; r < p; r++) if (T[r][k] != 0) { swap = r; break; }
				if (swap < 0) return 0;
				long[] tmp = T[k]; T[k] = T[swap]; T[swap] = tmp;
				sign = -sign;
			}
			for (int r = k + 1; r < p; r++) {
				for (int c = k + 1; c < p; c++) {
					T[r][c] = (T[k][k] * T[r][c] - T[r][k] * T[k][c]) / prev;
				}
				T[r][k] = 0;
			}
			prev = T[k][k];
		}
		return sign * T[p - 1][p - 1];
	}
}
