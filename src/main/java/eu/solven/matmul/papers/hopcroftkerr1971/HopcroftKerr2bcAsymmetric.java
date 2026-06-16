package eu.solven.matmul.papers.hopcroftkerr1971;

import eu.solven.matmul.NonCubicBilinearAlgorithm;

/**
 * Asymmetric Hopcroft-Kerr ⟨p, 2, n⟩ via Lemma 1 augmentation. Two
 * implementations:
 *
 * <ol>
 *   <li>{@link #buildNaive(int, int)} — sub-optimal but
 *       <em>constructively complete</em>: augment {@code A} to
 *       {@code Ā = M·A} via {@link LemmaOneAugmentation}, run the
 *       full {@link HopcroftKerr2bc#buildSquare} on {@code ⟨n, 2, n⟩}
 *       computing all {@code n²} outputs of {@code ĀX}, then discard
 *       the augmented output rows (they're linear combinations of the
 *       genuine ones, but we don't need them). Rank:
 *       {@code (3n² + n)/2} — larger than the HK-optimal
 *       {@code (3pn + max(p,n))/2}. Use as a baseline / framework
 *       validator.</li>
 *   <li>{@link #buildBanded(int, int)} — TODO. Restrict the square
 *       HK to a column band of width {@code p}, then back-substitute
 *       via the cyclic-{@code p}-window Vandermonde inverse. Achieves
 *       the HK-optimal rank. The Vandermonde inverse introduces
 *       the HK-optimal rank; with the unimodular-first Lemma-1 matrix
 *       (task #11) the inverses are integer and the scheme is Z-tagged
 *       (Q only on loud fallback).</li>
 * </ol>
 *
 * <p>Both rely on {@link LemmaOneAugmentation} for the {@code n × p}
 * augmentation matrix {@code M} with first {@code p} rows = identity
 * and every cyclic {@code p}-window non-singular.</p>
 */
public final class HopcroftKerr2bcAsymmetric {

	private HopcroftKerr2bcAsymmetric() {}

	/**
	 * Naive asymmetric {@code ⟨p, 2, n⟩}: augment + square + discard.
	 * Rank: {@code R(⟨n, 2, n⟩) = (3n² + n)/2} (not HK-optimal).
	 *
	 * <p>For testing the framework: the U-matrix transformation
	 * {@code U_asym = M^T · U_sq} lets us go from "inputs in {@code Ā}
	 * coords" to "inputs in {@code A} coords" without changing the
	 * bilinear products themselves. Since the first {@code p} rows of
	 * {@code M} are the identity, rows {@code 1..p} of {@code ĀX}
	 * equal {@code AX} directly — so the discard step has no
	 * back-substitution cost.</p>
	 */
	public static NonCubicBilinearAlgorithm buildNaive(int p, int n) {
		if (p < 2 || n < p) {
			throw new IllegalArgumentException("require p ≥ 2 and n ≥ p, got p=" + p + ", n=" + n);
		}
		int[][] M = LemmaOneAugmentation.build(p, n);
		NonCubicBilinearAlgorithm sq = HopcroftKerr2bc.buildSquare(n);
		int r = sq.r;
		double[][] srcU = sq.denseU();
		double[][] srcV = sq.denseV();
		double[][] srcW = sq.denseW();

		// U_asym[j·2 + c][k] = Σ_i M[i, j] · U_sq[i·2 + c][k]
		double[][] Uasym = new double[p * 2][r];
		for (int j = 0; j < p; j++) {
			for (int c = 0; c < 2; c++) {
				for (int k = 0; k < r; k++) {
					double sum = 0;
					for (int i = 0; i < n; i++) {
						sum += M[i][j] * srcU[i * 2 + c][k];
					}
					Uasym[j * 2 + c][k] = sum;
				}
			}
		}

		// V_asym = V_sq unchanged.
		double[][] Vasym = new double[2 * n][r];
		for (int idx = 0; idx < 2 * n; idx++) {
			System.arraycopy(srcV[idx], 0, Vasym[idx], 0, r);
		}

		// W_asym[i·n + j][k] = W_sq[i·n + j][k] for i ∈ [0, p) (discard rows p..n-1)
		double[][] Wasym = new double[p * n][r];
		for (int i = 0; i < p; i++) {
			for (int j = 0; j < n; j++) {
				System.arraycopy(srcW[i * n + j], 0, Wasym[i * n + j], 0, r);
			}
		}

		return new NonCubicBilinearAlgorithm(p, 2, n, Uasym, Vasym, Wasym);
	}

	/**
	 * {@link #buildNaive} followed by <em>dead-code elimination</em>:
	 * drop every bilinear product {@code k} whose contribution to all
	 * kept output cells is zero (i.e. column {@code k} of {@code W_asym}
	 * is all zero). These products only existed to feed the discarded
	 * rows {@code p..n-1} of {@code ĀX}, so dropping them is exact.
	 *
	 * <p>Rank: between {@code R(⟨p, 2, n⟩)} (HK-optimal, if we hit
	 * every saveable product) and {@code R(⟨n, 2, n⟩)} (the naive
	 * baseline). Empirically falls in between; not optimal because the
	 * square HK construction shares products across "real" and "augmented"
	 * outputs, so some products remain even after DCE.</p>
	 *
	 * <p>This is a guaranteed improvement over {@link #buildNaive}
	 * without needing the full HK Lemma 1 band-restriction logic, which
	 * remains a TODO.</p>
	 */
	public static NonCubicBilinearAlgorithm buildNaiveDCE(int p, int n) {
		NonCubicBilinearAlgorithm naive = buildNaive(p, n);
		int r = naive.r;
		double[][] srcU = naive.denseU();
		double[][] srcV = naive.denseV();
		double[][] srcW = naive.denseW();
		boolean[] keep = new boolean[r];
		int kept = 0;
		for (int k = 0; k < r; k++) {
			boolean used = false;
			for (int row = 0; row < srcW.length; row++) {
				if (srcW[row][k] != 0.0) { used = true; break; }
			}
			if (used) {
				keep[k] = true;
				kept++;
			}
		}
		if (kept == r) return naive;

		double[][] U = new double[srcU.length][kept];
		double[][] V = new double[srcV.length][kept];
		double[][] W = new double[srcW.length][kept];
		int newIdx = 0;
		for (int k = 0; k < r; k++) {
			if (!keep[k]) continue;
			for (int row = 0; row < U.length; row++) U[row][newIdx] = srcU[row][k];
			for (int row = 0; row < V.length; row++) V[row][newIdx] = srcV[row][k];
			for (int row = 0; row < W.length; row++) W[row][newIdx] = srcW[row][k];
			newIdx++;
		}
		return new NonCubicBilinearAlgorithm(naive.n, naive.m, naive.p, U, V, W);
	}

	/**
	 * <strong>HK-optimal Case 1</strong> asymmetric ⟨p, 2, n⟩ via cyclic
	 * band restriction + Vandermonde back-sub. Implements HK 1971 §4
	 * Case 1 (paper p9, p11) — for {@code p = 2k+1 odd}, {@code p ≤ n ≤ 2p-1}.
	 *
	 * <p>Pipeline:</p>
	 * <ol>
	 *   <li>{@code M = LemmaOneAugmentation.build(p, n)} — first p rows
	 *       identity, augmented rows below.</li>
	 *   <li>{@code internal = HopcroftKerr2bc.buildOddBanded(n, k)} —
	 *       computes cells {@code (i', j) ∈ S} (cyclic band).</li>
	 *   <li>U-side translation: {@code U_asym = M^T · U_internal}.</li>
	 *   <li>W-side back-sub: per column j, with M_j = M[band(j), :],
	 *       {@code AX[:, j] = M_j⁻¹ · ĀX[band(j), j]}, materialised by
	 *       multiplying W_internal rows from the band against
	 *       M_j⁻¹.</li>
	 * </ol>
	 *
	 * <p>With the unimodular-first Lemma-1 matrix (task #11) the window
	 * inverses are integer, so the output scheme is over Z; only when the
	 * unimodular seam DFS fails (loud fallback to the ternary draw) do the
	 * window inverses introduce rational coefficients into {@code W_asym}
	 * (scheme over Q).</p>
	 *
	 * <p>Achieves the HK formula rank {@code (3pn + n)/2}.</p>
	 *
	 * @param p odd, ≥ 3
	 * @param n {@code p ≤ n ≤ 2p-1} so the band size {@code p} fits cyclically without overlap
	 */
	public static NonCubicBilinearAlgorithm buildOdd(int p, int n) {
		if ((p & 1) == 0) throw new IllegalArgumentException("buildOdd: p must be odd, got " + p);
		if (p < 3) throw new IllegalArgumentException("buildOdd: p must be ≥ 3, got " + p);
		if (n < p || n > 2 * p - 1) {
			throw new IllegalArgumentException("buildOdd: need p ≤ n ≤ 2p-1, got p=" + p + ", n=" + n);
		}
		int k = (p - 1) / 2;
		// Lemma-1 matrix: UNIMODULAR-first (task #11 — every window det ±1 ⇒
		// integer back-sub ⇒ Z scheme), falling back to the small-coefficient
		// ternary draw (task #7; scheme stays Q) when no seam is found.
		long[][] M = lemmaOneMatrix(p, n);
		NonCubicBilinearAlgorithm internal = HopcroftKerr2bc.buildOddBanded(n, k);
		int r = internal.r;
		double[][] srcU = internal.denseU();
		double[][] srcV = internal.denseV();
		double[][] srcW = internal.denseW();

		// U_asym[i*2+c][k] = Σ_{i'} M[i', i] · U_internal[i'*2+c][k]   (1-based: M[i'][i-1], etc.)
		double[][] U = new double[p * 2][r];
		for (int i = 0; i < p; i++) {
			for (int c = 0; c < 2; c++) {
				for (int kk = 0; kk < r; kk++) {
					double s = 0;
					for (int ip = 0; ip < n; ip++) {
						s += (double) M[ip][i] * srcU[ip * 2 + c][kk];
					}
					U[i * 2 + c][kk] = s;
				}
			}
		}

		// V_asym = V_internal.
		double[][] V = new double[2 * n][r];
		for (int row = 0; row < 2 * n; row++) {
			System.arraycopy(srcV[row], 0, V[row], 0, r);
		}

		// W back-sub: for each column j, compute M_j = M[band(j), :] (p×p)
		// and its inverse, then W_asym[(i-1)*n + (j-1)][k] =
		//   Σ_{i' ∈ band(j)} M_j^{-1}[i, i'] · W_internal[(i'-1)*n + (j-1)][k].
		// EXACT-FRACTION inversion (task #7): the double Gauss-Jordan lost
		// exactness for n ≳ p+2 even with small-coefficient M; with exact
		// rationals the only rounding is the final BigFraction→double, which is
		// exact to 1 ulp on the small fractions the small-coefficient M yields.
		double[][] W = new double[p * n][r];
		for (int j0 = 0; j0 < n; j0++) {
			int[] band = cyclicBand(n, j0, k);  // p row indices (0-based) for this column
			org.apache.commons.math3.fraction.BigFraction[][] Mj =
					new org.apache.commons.math3.fraction.BigFraction[p][p];
			for (int row = 0; row < p; row++) {
				for (int col = 0; col < p; col++) {
					Mj[row][col] = new org.apache.commons.math3.fraction.BigFraction(M[band[row]][col]);
				}
			}
			org.apache.commons.math3.fraction.BigFraction[][] Mjinv = exactInverse(Mj);
			for (int i = 0; i < p; i++) {
				for (int kk = 0; kk < r; kk++) {
					org.apache.commons.math3.fraction.BigFraction s =
							org.apache.commons.math3.fraction.BigFraction.ZERO;
					for (int bi = 0; bi < p; bi++) {
						int ip = band[bi];
						double w = srcW[ip * n + j0][kk];
						if (w == 0.0) continue;
						// Internal W entries are small integers (the HK emitters use ±1).
						long wl = Math.round(w);
						if (Math.abs(w - wl) > 1e-9) {
							throw new IllegalStateException("non-integer internal W entry " + w);
						}
						s = s.add(Mjinv[i][bi].multiply(new org.apache.commons.math3.fraction.BigFraction(wl)));
					}
					W[i * n + j0][kk] = s.doubleValue();
				}
			}
		}

		return new NonCubicBilinearAlgorithm(p, 2, n, U, V, W);
	}

	/**
	 * <strong>HK-optimal Case 2</strong> asymmetric ⟨p, 2, n⟩ for EVEN {@code p}
	 * (task #7). Same pipeline as {@link #buildOdd} but the internal computation
	 * is {@link HopcroftKerr2bc#buildEvenBanded}: the odd band of half-width
	 * {@code k = p/2 − 1} plus one extra distance-(k+1) cell per column, served
	 * by a circulant matching (full pairs, 3 products per 2 columns). Per-column
	 * windows are direction-aware: {@code [j−k, j+k+1]} when the extra cell is
	 * "up", {@code [j−k−1, j+k]} when "down" — both contiguous, so the Lemma-1
	 * window inverse applies unchanged.
	 *
	 * <p>Attains the HK formula exactly whenever the circulant {@code i ↔ i+k+1
	 * (mod n)} decomposes into even cycles (or a single odd cycle when n is odd —
	 * the ceiling absorbs one naive leftover); costs ~2 extra products per
	 * additional odd cycle otherwise (reported by the caller via rank).</p>
	 */
	public static NonCubicBilinearAlgorithm buildEven(int p, int n) {
		if ((p & 1) == 1) throw new IllegalArgumentException("buildEven: p must be even, got " + p);
		if (p < 4) throw new IllegalArgumentException("buildEven: p must be ≥ 4, got " + p);
		if (n < p || n > 2 * p - 1) {
			throw new IllegalArgumentException("buildEven: need p ≤ n ≤ 2p-1, got p=" + p + ", n=" + n);
		}
		int k = p / 2 - 1;
		long[][] M = lemmaOneMatrix(p, n);
		HopcroftKerr2bc.BandedEven banded = HopcroftKerr2bc.buildEvenBanded(n, k);
		NonCubicBilinearAlgorithm internal = banded.alg();
		boolean[] extraIsUp = banded.extraIsUp();
		int r = internal.r;
		double[][] srcU = internal.denseU();
		double[][] srcV = internal.denseV();
		double[][] srcW = internal.denseW();

		double[][] U = new double[p * 2][r];
		for (int i = 0; i < p; i++) {
			for (int c = 0; c < 2; c++) {
				for (int kk = 0; kk < r; kk++) {
					double s = 0;
					for (int ip = 0; ip < n; ip++) {
						s += (double) M[ip][i] * srcU[ip * 2 + c][kk];
					}
					U[i * 2 + c][kk] = s;
				}
			}
		}

		double[][] V = new double[2 * n][r];
		for (int row = 0; row < 2 * n; row++) {
			System.arraycopy(srcV[row], 0, V[row], 0, r);
		}

		// Precompute per-column windows + exact inverses (needed both for the
		// Z-pair corrections below and the final back-sub).
		var ZERO = org.apache.commons.math3.fraction.BigFraction.ZERO;
		int[][] bandPer = new int[n][];
		org.apache.commons.math3.fraction.BigFraction[][][] invPer =
				new org.apache.commons.math3.fraction.BigFraction[n][][];
		for (int j0 = 0; j0 < n; j0++) {
			int startOffset = extraIsUp[j0] ? -k : -(k + 1);
			int[] band = new int[p];
			for (int t = 0; t < p; t++) {
				band[t] = Math.floorMod(j0 + startOffset + t, n);
			}
			bandPer[j0] = band;
			org.apache.commons.math3.fraction.BigFraction[][] Mj =
					new org.apache.commons.math3.fraction.BigFraction[p][p];
			for (int row = 0; row < p; row++) {
				for (int col = 0; col < p; col++) {
					Mj[row][col] = new org.apache.commons.math3.fraction.BigFraction(M[band[row]][col]);
				}
			}
			invPer[j0] = exactInverse(Mj);
		}

		// Exact-fraction view of internal W rows, with a mutable overlay for the
		// Z-pair corrections (corrected rows leave ℤ).
		java.util.Map<Integer, org.apache.commons.math3.fraction.BigFraction[]> overlay =
				new java.util.HashMap<>();
		java.util.function.IntFunction<org.apache.commons.math3.fraction.BigFraction[]> rowOf = out -> {
			var cached = overlay.get(out);
			if (cached != null) return cached;
			var row = new org.apache.commons.math3.fraction.BigFraction[r];
			for (int kk = 0; kk < r; kk++) {
				double w = srcW[out][kk];
				long wl = Math.round(w);
				if (Math.abs(w - wl) > 1e-9) {
					throw new IllegalStateException("non-integer internal W entry " + w);
				}
				row[kk] = wl == 0 ? ZERO : new org.apache.commons.math3.fraction.BigFraction(wl);
			}
			return row;
		};

		// Z-pair corrections (repaired HK Step 3). The internal builder wrote
		//   W[i1,i2] = Z(β,α_x)   and   W[i3,i4] = −Z(α,β_x)
		// (on top of nothing — those cells had no other contribution). Complete:
		//   y_{i1,i2} = Z(β,α_x) − y_{i1,i3} + y_{i4,i2} + y_{i4,i3}
		//   y_{i3,i4} = −Z(α,β_x) + y_{i2,i1} − y_{i2,i4} + y_{i3,i1}
		// In-band cells are direct internal rows; out-of-band cells of COMPLETE
		// columns c are rational combinations: M[r,:]·Mwin(c)⁻¹·(window cells of c).
		for (HopcroftKerr2bc.ZPair z : banded.zPairs()) {
			int i1 = z.i1() - 1, i2 = z.i2() - 1, i3 = z.i3() - 1, i4 = z.i4() - 1;
			// y_{i1,i2} completion:
			var rowI1I2 = rowOf.apply(i1 * n + i2).clone();
			addScaled(rowI1I2, reconstructRow(i1, i3, bandPer, invPer, M, rowOf, n, p, r), -1);
			addScaled(rowI1I2, rowOf.apply(i4 * n + i2), +1);
			// y_{i4,i3} is at distance k+1 (i3 = i4−(k+1) by construction) — OUTSIDE
			// the band; recover it from complete column i3, not from a raw W row.
			addScaled(rowI1I2, reconstructRow(i4, i3, bandPer, invPer, M, rowOf, n, p, r), +1);
			overlay.put(i1 * n + i2, rowI1I2);
			// y_{i3,i4} completion:
			var rowI3I4 = rowOf.apply(i3 * n + i4).clone();
			addScaled(rowI3I4, reconstructRow(i2, i1, bandPer, invPer, M, rowOf, n, p, r), +1);
			addScaled(rowI3I4, rowOf.apply(i2 * n + i4), -1);
			addScaled(rowI3I4, reconstructRow(i3, i1, bandPer, invPer, M, rowOf, n, p, r), +1);
			overlay.put(i3 * n + i4, rowI3I4);
		}

		double[][] W = new double[p * n][r];
		for (int j0 = 0; j0 < n; j0++) {
			int[] band = bandPer[j0];
			var Mjinv = invPer[j0];
			// Fetch (possibly corrected) window rows once.
			org.apache.commons.math3.fraction.BigFraction[][] cellRows =
					new org.apache.commons.math3.fraction.BigFraction[p][];
			for (int bi = 0; bi < p; bi++) {
				cellRows[bi] = rowOf.apply(band[bi] * n + j0);
			}
			for (int i = 0; i < p; i++) {
				for (int kk = 0; kk < r; kk++) {
					var s = ZERO;
					for (int bi = 0; bi < p; bi++) {
						var w = cellRows[bi][kk];
						if (w.getNumerator().signum() == 0) continue;
						s = s.add(Mjinv[i][bi].multiply(w));
					}
					W[i * n + j0][kk] = s.doubleValue();
				}
			}
		}

		return new NonCubicBilinearAlgorithm(p, 2, n, U, V, W);
	}

	/** {@code dst += sign · src} over fraction rows. */
	private static void addScaled(org.apache.commons.math3.fraction.BigFraction[] dst,
			org.apache.commons.math3.fraction.BigFraction[] src, int sign) {
		for (int kk = 0; kk < dst.length; kk++) {
			if (src[kk].getNumerator().signum() == 0) continue;
			dst[kk] = (sign > 0) ? dst[kk].add(src[kk]) : dst[kk].subtract(src[kk]);
		}
	}

	/**
	 * Exact W row of the out-of-band cell {@code y_{row, col}} of a COMPLETE
	 * column {@code col}: {@code M[row,:] · Mwin(col)⁻¹ · (window cells of col)}.
	 * (0-based args.)
	 */
	private static org.apache.commons.math3.fraction.BigFraction[] reconstructRow(int row, int col,
			int[][] bandPer, org.apache.commons.math3.fraction.BigFraction[][][] invPer, long[][] M,
			java.util.function.IntFunction<org.apache.commons.math3.fraction.BigFraction[]> rowOf,
			int n, int p, int r) {
		var ZERO = org.apache.commons.math3.fraction.BigFraction.ZERO;
		int[] band = bandPer[col];
		var Mjinv = invPer[col];
		// coeff over window positions: v[bi] = Σ_c M[row][c] · Mjinv[c][bi]
		var v = new org.apache.commons.math3.fraction.BigFraction[p];
		for (int bi = 0; bi < p; bi++) {
			var s = ZERO;
			for (int c = 0; c < p; c++) {
				if (M[row][c] == 0) continue;
				s = s.add(Mjinv[c][bi].multiply(new org.apache.commons.math3.fraction.BigFraction(M[row][c])));
			}
			v[bi] = s;
		}
		var out = new org.apache.commons.math3.fraction.BigFraction[r];
		java.util.Arrays.fill(out, ZERO);
		for (int bi = 0; bi < p; bi++) {
			if (v[bi].getNumerator().signum() == 0) continue;
			var cell = rowOf.apply(band[bi] * n + col);
			for (int kk = 0; kk < r; kk++) {
				if (cell[kk].getNumerator().signum() == 0) continue;
				out[kk] = out[kk].add(v[bi].multiply(cell[kk]));
			}
		}
		return out;
	}

	/**
	 * <strong>Smart dispatcher</strong> for asymmetric {@code ⟨p, 2, n⟩}.
	 * Returns the best constructively-valid scheme this class can produce:
	 * <ul>
	 *   <li>odd {@code p ≥ 3}, {@code p ≤ n ≤ 2p−1}: {@link #buildOdd} —
	 *       at the HK formula {@code (3pn+n)/2} for every shape.</li>
	 *   <li>even {@code p ≥ 4}, {@code p ≤ n ≤ 2p−1}: {@link #buildEven} —
	 *       at the formula except the g ≥ 6 circulant combos (then
	 *       {@code +⌈(g−1)/2⌉}-ish; see CONSTRUCTIVE_METHOD.md).</li>
	 *   <li>{@code p ≥ 3}, {@code n > 2p−1}: {@link #buildChained} —
	 *       DP-optimal concatenation of band segments; at the formula
	 *       whenever a slack-free, non-degraded partition exists.</li>
	 *   <li>{@code p = 2}: {@link #buildNaiveDCE} — sub-optimal but
	 *       always returns a verified scheme.</li>
	 * </ul>
	 */
	public static NonCubicBilinearAlgorithm build(int p, int n) {
		if (p < 2 || n < p) {
			throw new IllegalArgumentException("require p ≥ 2 and n ≥ p, got p=" + p + ", n=" + n);
		}
		// HK-optimal path: cyclic-band + exact back-sub, ALL odd p (task #7).
		// The historical p ∈ {3, 5} cap guarded two since-fixed defects: the
		// linear-vs-cyclic bridge-position wrap (same-method emitters walked off
		// the index range) and the Vandermonde overflow/conditioning in the
		// Lemma-1 matrix + W back-sub. With arc-interior bridge SELECTION the
		// impossible (2,2,bridge-3) case never arises, and with the
		// small-coefficient exact pipeline the result verifies at every
		// p ≤ n ≤ 2p−1 (empirically swept p ≤ 13, all exact at formula).
		if ((p & 1) == 1 && p >= 3 && n <= 2 * p - 1) {
			return buildOdd(p, n);
		}
		// Case 2 (even p): banded + circulant matching + direction-aware back-sub.
		// Attains formula when the circulant cycles are even (or single-odd, n odd);
		// otherwise lands within a couple of products of it — still far better than
		// the DCE fallback, and exact-verified either way.
		if ((p & 1) == 0 && p >= 4 && n <= 2 * p - 1) {
			return buildEven(p, n);
		}
		// n > 2p−1: chained augmentation — partition the columns into segments
		// of size in [p, 2p−1] and concatenate independent band constructions.
		if (p >= 3 && n > 2 * p - 1) {
			return buildChained(p, n);
		}
		// Fallback (p = 2): augment + full-square + discard + DCE.
		return buildNaiveDCE(p, n);
	}

	/**
	 * <strong>Chained augmentation</strong> for {@code n > 2p−1} (HK 1971's
	 * own large-n regime): one Lemma-1 band covers at most {@code 2p−1}
	 * columns, so partition the {@code n} columns into segments of size
	 * {@code s ∈ [p, 2p−1]}, build each segment with the band construction,
	 * and concatenate along the column axis (ranks add — no sharing).
	 *
	 * <p>Since {@code max(p, s) = s} throughout the segment range, per-segment
	 * formula values telescope to the global formula {@code ⌈n(3p+1)/2⌉}
	 * provided the ceiling slack is killed: for odd p every segment cost
	 * {@code s(3p+1)/2} is an integer (no slack, any partition works); for
	 * even p at most one odd-size segment may be used (matching n's parity),
	 * and g ≥ 6 triangle-family segment sizes carry their +1..+3 penalty.
	 * Rather than hand-coding those rules, the partition is chosen by a DP
	 * over the <em>actual achieved ranks</em> of the segment builds — so
	 * parity slack and degraded segments are avoided exactly when avoidable,
	 * and the result is optimal over all band-segment partitions (a bound,
	 * not a proven optimum, for the shape itself).</p>
	 */
	public static NonCubicBilinearAlgorithm buildChained(int p, int n) {
		if (p < 3) throw new IllegalArgumentException("buildChained: p must be ≥ 3, got " + p);
		if (n <= 2 * p - 1) throw new IllegalArgumentException(
				"buildChained: need n > 2p−1 (use build for the band range), got p=" + p + ", n=" + n);
		int hi = 2 * p - 1;
		// Lazily-built segment schemes, indexed by size. Deterministic builders,
		// so one build per size suffices.
		NonCubicBilinearAlgorithm[] bySize = new NonCubicBilinearAlgorithm[hi + 1];
		java.util.function.IntFunction<NonCubicBilinearAlgorithm> segment = s -> {
			if (bySize[s] == null) bySize[s] = build(p, s);
			return bySize[s];
		};
		// DP over prefix length: best[c] = minimal total rank covering c columns.
		int[] best = new int[n + 1];
		int[] choice = new int[n + 1];
		java.util.Arrays.fill(best, Integer.MAX_VALUE);
		best[0] = 0;
		for (int c = p; c <= n; c++) {
			for (int s = p; s <= Math.min(hi, c); s++) {
				if (c - s != 0 && best[c - s] == Integer.MAX_VALUE) continue;
				if (c - s != 0 && c - s < p) continue; // remainder must host ≥ 1 full segment
				int cand = (c == s ? 0 : best[c - s]) + segment.apply(s).r;
				if (cand < best[c]) {
					best[c] = cand;
					choice[c] = s;
				}
			}
		}
		if (best[n] == Integer.MAX_VALUE) {
			throw new IllegalStateException("no [p,2p−1] partition of n=" + n + " for p=" + p);
		}
		// Recover the partition (order is immaterial for the rank).
		java.util.List<NonCubicBilinearAlgorithm> parts = new java.util.ArrayList<>();
		for (int c = n; c > 0; c -= choice[c]) {
			parts.add(segment.apply(choice[c]));
		}
		NonCubicBilinearAlgorithm out = parts.get(0);
		for (int i = 1; i < parts.size(); i++) {
			out = concatColumns(out, parts.get(i));
		}
		if (out.r != best[n]) {
			throw new IllegalStateException("chained rank " + out.r + " != DP optimum " + best[n]);
		}
		return out;
	}

	/**
	 * Direct sum of two {@code ⟨p,2,·⟩} schemes along the column axis:
	 * {@code ⟨p,2,na⟩ ⊕ ⟨p,2,nb⟩ → ⟨p,2,na+nb⟩}. U rows (the shared A side)
	 * are copied per-segment into disjoint product columns; V (2×n row-major)
	 * and W (p×n row-major) rows are re-indexed into the widened column space.
	 */
	static NonCubicBilinearAlgorithm concatColumns(NonCubicBilinearAlgorithm a, NonCubicBilinearAlgorithm b) {
		if (a.n != b.n || a.m != 2 || b.m != 2) {
			throw new IllegalArgumentException("concatColumns: incompatible shapes ⟨" + a.n + ",2," + a.p
					+ "⟩ + ⟨" + b.n + ",2," + b.p + "⟩");
		}
		int p = a.n, na = a.p, nb = b.p, n = na + nb, r = a.r + b.r;
		double[][] aU = a.denseU(), aV = a.denseV(), aW = a.denseW();
		double[][] bU = b.denseU(), bV = b.denseV(), bW = b.denseW();
		double[][] U = new double[p * 2][r];
		for (int row = 0; row < p * 2; row++) {
			System.arraycopy(aU[row], 0, U[row], 0, a.r);
			System.arraycopy(bU[row], 0, U[row], a.r, b.r);
		}
		double[][] V = new double[2 * n][r];
		double[][] W = new double[p * n][r];
		for (int row = 0; row < 2; row++) {
			for (int j = 0; j < na; j++) {
				System.arraycopy(aV[row * na + j], 0, V[row * n + j], 0, a.r);
			}
			for (int j = 0; j < nb; j++) {
				System.arraycopy(bV[row * nb + j], 0, V[row * n + na + j], a.r, b.r);
			}
		}
		for (int i = 0; i < p; i++) {
			for (int j = 0; j < na; j++) {
				System.arraycopy(aW[i * na + j], 0, W[i * n + j], 0, a.r);
			}
			for (int j = 0; j < nb; j++) {
				System.arraycopy(bW[i * nb + j], 0, W[i * n + na + j], a.r, b.r);
			}
		}
		return new NonCubicBilinearAlgorithm(p, 2, n, U, V, W);
	}

	/**
	 * Lemma-1 matrix selection (task #11): unimodular-first — when every cyclic
	 * window has det ±1, the back-substitution is integer and the output scheme
	 * is over Z (⇒ F₂/F₃ too) instead of Q. Falls back to the small-coefficient
	 * ternary draw when the seam DFS finds nothing (scheme stays Q; the field
	 * stamping reads the actual coefficients, so the claim is always honest).
	 */
	private static long[][] lemmaOneMatrix(int p, int n) {
		long[][] M = LemmaOneAugmentation.buildUnimodular(p, n);
		if (M == null) {
			// Loud, not silent (fail-loud discipline): this shape's scheme will
			// carry rational coefficients (Q) instead of integer (Z).
			System.err.println("[lemma1] no unimodular seam for p=" + p + ", n=" + n
					+ " within DFS budget — falling back to ternary draw (scheme stays Q)");
			M = LemmaOneAugmentation.buildSmallLong(p, n);
		}
		return M;
	}

	/**
	 * Cyclic band of p = 2k+1 consecutive row indices (0-based) centered
	 * on column {@code j0}. Returned in cyclic order starting at
	 * {@code j0 - k mod n}.
	 */
	private static int[] cyclicBand(int n, int j0, int k) {
		int p = 2 * k + 1;
		int[] out = new int[p];
		for (int t = 0; t < p; t++) {
			out[t] = Math.floorMod(j0 - k + t, n);
		}
		return out;
	}

	/** Exact Gauss-Jordan inverse over BigFraction for small p (≤ 16). */
	private static org.apache.commons.math3.fraction.BigFraction[][] exactInverse(
			org.apache.commons.math3.fraction.BigFraction[][] A) {
		int n = A.length;
		var ZERO = org.apache.commons.math3.fraction.BigFraction.ZERO;
		var ONE = org.apache.commons.math3.fraction.BigFraction.ONE;
		org.apache.commons.math3.fraction.BigFraction[][] aug =
				new org.apache.commons.math3.fraction.BigFraction[n][2 * n];
		for (int i = 0; i < n; i++) {
			for (int j = 0; j < n; j++) aug[i][j] = A[i][j];
			for (int j = n; j < 2 * n; j++) aug[i][j] = (j - n == i) ? ONE : ZERO;
		}
		for (int col = 0; col < n; col++) {
			int piv = -1;
			for (int r = col; r < n; r++) {
				if (aug[r][col].getNumerator().signum() != 0) { piv = r; break; }
			}
			if (piv < 0) throw new ArithmeticException("singular M-submatrix at col " + col);
			if (piv != col) { var tmp = aug[col]; aug[col] = aug[piv]; aug[piv] = tmp; }
			var pivVal = aug[col][col];
			for (int j = 0; j < 2 * n; j++) aug[col][j] = aug[col][j].divide(pivVal);
			for (int r = 0; r < n; r++) {
				if (r == col) continue;
				var factor = aug[r][col];
				if (factor.getNumerator().signum() == 0) continue;
				for (int j = 0; j < 2 * n; j++) {
					aug[r][j] = aug[r][j].subtract(factor.multiply(aug[col][j]));
				}
			}
		}
		org.apache.commons.math3.fraction.BigFraction[][] inv =
				new org.apache.commons.math3.fraction.BigFraction[n][n];
		for (int i = 0; i < n; i++) System.arraycopy(aug[i], n, inv[i], 0, n);
		return inv;
	}

	/** Gauss-Jordan inverse for small p (≤ 16). */
	@SuppressWarnings("unused")
	private static double[][] inverse(double[][] A) {
		int n = A.length;
		double[][] aug = new double[n][2 * n];
		for (int i = 0; i < n; i++) {
			for (int j = 0; j < n; j++) aug[i][j] = A[i][j];
			aug[i][n + i] = 1;
		}
		for (int col = 0; col < n; col++) {
			int piv = col;
			double best = Math.abs(aug[col][col]);
			for (int r = col + 1; r < n; r++) {
				double v = Math.abs(aug[r][col]);
				if (v > best) { best = v; piv = r; }
			}
			if (best < 1e-12) throw new ArithmeticException("singular M-submatrix at col " + col);
			if (piv != col) { double[] tmp = aug[col]; aug[col] = aug[piv]; aug[piv] = tmp; }
			double pivVal = aug[col][col];
			for (int j = 0; j < 2 * n; j++) aug[col][j] /= pivVal;
			for (int r = 0; r < n; r++) {
				if (r == col) continue;
				double factor = aug[r][col];
				if (factor == 0) continue;
				for (int j = 0; j < 2 * n; j++) aug[r][j] -= factor * aug[col][j];
			}
		}
		double[][] inv = new double[n][n];
		for (int i = 0; i < n; i++) {
			System.arraycopy(aug[i], n, inv[i], 0, n);
		}
		return inv;
	}
}
