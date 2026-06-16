package eu.solven.matmul.papers.rosowski2019;

import eu.solven.matmul.NonBilinearAlgorithm;

/**
 * Rosowski 2019 Theorem 2 — explicit DIVISIONS-FREE non-bilinear algorithm
 * for the product of an {@code l × n} matrix by an {@code n × m} matrix over a
 * commutative ring, valid when the contraction dimension {@code n} is EVEN.
 * Cost: {@code n(lm + l + m − 1)/2} multiplications.
 *
 * <p>The rank {@code n(lm+l+m−1)/2} is NOT a Rosowski discovery — it was
 * already known (Waksman 1970 [19], Islam 2009 [10]) but those constructions
 * use divisions by 2 (so they need 2 invertible: Q/R/C, not Z). Rosowski's
 * Theorem 2 contribution is a construction with the <em>same</em> rank that
 * uses NO divisions, so it is valid over ANY commutative ring including
 * {@code Z}, {@code F₂}, {@code F₃}. That divisions-free property is the
 * reason to materialise it separately from Waksman.</p>
 *
 * <p>Specialised to {@code l = n = 2}, {@code m = p} this is the
 * {@code ⟨2,2,p⟩ = 3p+1} family the catalog tracks (e.g. ⟨2,2,3⟩=10,
 * ⟨2,2,4⟩=13, …, ⟨2,2,16⟩=49).</p>
 *
 * <p>Paper formulas (1-indexed; {@code k = 1..n/2}). For {@code i = 1..l}:</p>
 * <pre>
 *   c_{i,1} = Σ_k a_{i,2k−1}(b_{2k−1,1} + a_{i,2k})
 *           + Σ_k a_{i,2k}  (b_{2k,1}   − a_{i,2k−1})
 * </pre>
 * <p>For {@code i = 1..l} and {@code j = 2..m}:</p>
 * <pre>
 *   c_{i,j} = Σ_k (a_{i,2k−1} + b_{2k,j})(a_{i,2k} + b_{2k−1,1} + b_{2k−1,j})
 *           − Σ_k a_{i,2k−1}(b_{2k−1,1} + a_{i,2k})
 *           − Σ_k b_{2k,j}  (b_{2k−1,1} + b_{2k−1,j})
 * </pre>
 *
 * <p>The four product families and their counts:</p>
 * <ul>
 *   <li>{@code P1(i,k) = a_{i,2k−1}(b_{2k−1,1}+a_{i,2k})} — {@code l·n/2},
 *       reused in {@code c_{i,1}} (add) and every {@code c_{i,j≥2}} (subtract).</li>
 *   <li>{@code P2(i,k) = a_{i,2k}(b_{2k,1}−a_{i,2k−1})} — {@code l·n/2}.</li>
 *   <li>{@code S(j,k)  = b_{2k,j}(b_{2k−1,1}+b_{2k−1,j})} — {@code (m−1)·n/2},
 *       B-only so shared across all rows {@code i}.</li>
 *   <li>{@code Q(i,j,k)= (a_{i,2k−1}+b_{2k,j})(a_{i,2k}+b_{2k−1,1}+b_{2k−1,j})}
 *       — {@code l(m−1)·n/2}.</li>
 * </ul>
 * <p>Total {@code = n/2·[2l + (l+1)(m−1)] = n(lm+l+m−1)/2}.</p>
 *
 * <p><strong>Commutative-only</strong>: the cancellations (e.g. the
 * {@code a_{i,2k−1}a_{i,2k}} cross-terms in {@code c_{i,1}}, the
 * {@code b·b} cross-terms in {@code c_{i,j}}) require a commutative ring; this
 * does NOT lift to recursive matmul over a non-commutative base. See
 * {@code references/rosowski-algorithms.md} and
 * {@link RosowskiBound#commutativeBoundBilinear}.</p>
 */
public final class RosowskiTheorem2 {

	private RosowskiTheorem2() {}

	/**
	 * Builds the explicit divisions-free non-bilinear scheme for an
	 * {@code l × nContr} by {@code nContr × m} product.
	 *
	 * @param l       rows of A (and of C); {@code ≥ 1}
	 * @param nContr  contraction dimension (cols of A, rows of B); must be EVEN and {@code ≥ 2}
	 * @param m       cols of B (and of C); {@code ≥ 1}
	 * @return a non-bilinear algorithm with {@code r = nContr(l·m+l+m−1)/2} products
	 */
	public static NonBilinearAlgorithm build(int l, int nContr, int m) {
		if (l < 1 || m < 1) throw new IllegalArgumentException("l,m must be ≥ 1");
		if (nContr < 2 || nContr % 2 != 0)
			throw new IllegalArgumentException("Theorem 2 requires an even contraction dim ≥ 2, got " + nContr);

		final int half = nContr / 2;
		final int r = nContr * (l * m + l + m - 1) / 2;

		// Product column layout (deterministic):
		//   P region : 2 columns per (i,k)            -> 2·l·half
		//   S region : 1 column  per (j≥2,k)          -> (m−1)·half   (B-only, shared over i)
		//   Q region : 1 column  per (i,j≥2,k)        -> l·(m−1)·half
		final int sBase = 2 * l * half;
		final int qBase = sBase + (m - 1) * half;

		double[][] Ua = new double[l * nContr][r];
		double[][] Ub = new double[nContr * m][r];
		double[][] Va = new double[l * nContr][r];
		double[][] Vb = new double[nContr * m][r];
		double[][] W  = new double[l * m][r];

		for (int i = 1; i <= l; i++) {
			for (int k = 1; k <= half; k++) {
				int p1 = ((i - 1) * half + (k - 1)) * 2;
				int p2 = p1 + 1;

				// P1 = a_{i,2k−1} · (b_{2k−1,1} + a_{i,2k})
				Ua[aRow(i, 2 * k - 1, nContr)][p1] = 1;
				Vb[bRow(2 * k - 1, 1, m)][p1] += 1;
				Va[aRow(i, 2 * k, nContr)][p1] += 1;

				// P2 = a_{i,2k} · (b_{2k,1} − a_{i,2k−1})
				Ua[aRow(i, 2 * k, nContr)][p2] = 1;
				Vb[bRow(2 * k, 1, m)][p2] += 1;
				Va[aRow(i, 2 * k - 1, nContr)][p2] += -1;

				// c_{i,1} += P1 + P2
				W[cRow(i, 1, m)][p1] += 1;
				W[cRow(i, 1, m)][p2] += 1;
			}
		}

		// S(j,k) = b_{2k,j} · (b_{2k−1,1} + b_{2k−1,j})  — B-only, shared across rows.
		for (int j = 2; j <= m; j++) {
			for (int k = 1; k <= half; k++) {
				int s = sBase + ((j - 2) * half + (k - 1));
				Ub[bRow(2 * k, j, m)][s] = 1;
				Vb[bRow(2 * k - 1, 1, m)][s] += 1;
				Vb[bRow(2 * k - 1, j, m)][s] += 1;
			}
		}

		// Q(i,j,k) and the c_{i,j≥2} output combination.
		for (int i = 1; i <= l; i++) {
			for (int j = 2; j <= m; j++) {
				for (int k = 1; k <= half; k++) {
					int q = qBase + (((i - 1) * (m - 1) + (j - 2)) * half + (k - 1));
					int p1 = ((i - 1) * half + (k - 1)) * 2;
					int s = sBase + ((j - 2) * half + (k - 1));

					// Q = (a_{i,2k−1} + b_{2k,j}) · (a_{i,2k} + b_{2k−1,1} + b_{2k−1,j})
					Ua[aRow(i, 2 * k - 1, nContr)][q] += 1;
					Ub[bRow(2 * k, j, m)][q] += 1;
					Va[aRow(i, 2 * k, nContr)][q] += 1;
					Vb[bRow(2 * k - 1, 1, m)][q] += 1;
					Vb[bRow(2 * k - 1, j, m)][q] += 1;

					// c_{i,j} += Q − P1(i,k) − S(j,k)
					W[cRow(i, j, m)][q] += 1;
					W[cRow(i, j, m)][p1] += -1;
					W[cRow(i, j, m)][s] += -1;
				}
			}
		}

		return new NonBilinearAlgorithm(l, nContr, m, Ua, Ub, Va, Vb, W);
	}

	/** Convenience: the {@code ⟨2,2,p⟩ = 3p+1} family (l = nContr = 2). */
	public static NonBilinearAlgorithm build22p(int p) {
		return build(2, 2, p);
	}

	// (1-indexed paper convention) → 0-indexed flat array index helpers.
	/** Row index in {@code Ua/Va} for {@code A[i,t]} (A is l×nContr). */
	private static int aRow(int i, int t, int nContr) { return (i - 1) * nContr + (t - 1); }
	/** Row index in {@code Ub/Vb} for {@code B[s,t]} (B is nContr×m). */
	private static int bRow(int s, int t, int m) { return (s - 1) * m + (t - 1); }
	/** Row index in {@code W} for {@code C[i,j]} (C is l×m). */
	private static int cRow(int i, int j, int m) { return (i - 1) * m + (j - 1); }
}
