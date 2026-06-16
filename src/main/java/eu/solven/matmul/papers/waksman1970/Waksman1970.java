package eu.solven.matmul.papers.waksman1970;

import eu.solven.matmul.NonBilinearAlgorithm;

/**
 * Waksman 1970 generic {@code ⟨n,n,n⟩} COMMUTATIVE matmul, as exposed in
 * Mezzarobba 2007 (MSc thesis, MPRI, Fig. 3 + general formula).
 *
 * <p>Identity exploited (requires {@code 1/2 ∈ K}, hence characteristic
 * ≠ 2): {@code (x+y)(u+v) - (x-y)(u-v) = 2(xv + yu)}. Each pair of
 * columns {@code 2t-1, 2t} of {@code A} is "folded" with a pair of rows
 * {@code 2t-1, 2t} of {@code B} into a constant number of products
 * sharable across rows.</p>
 *
 * <p>For each {@code t ∈ [1..⌊n/2⌋]}, the algorithm forms:</p>
 * <pre>
 *   A(i,j,t) = (a[i,2t-1] + b[2t,j])   · (a[i,2t]   + b[2t-1,j])     ∀ i,j
 *   B(1,j,t) = A(1,j,t) + (a[1,2t-1] - b[2t,j])(a[1,2t] - b[2t-1,j]) ∀ j
 *   B(i,i,t) = A(i,i,t) + (a[i,2t-1] - b[2t,i])(a[i,2t] - b[2t-1,i]) i≥2
 *   B(i,j,t) = B(i,i,t) - B(1,i,t) + B(1,j,t)                        i≥2, j≠i
 *
 *   c[i,j] = Σ_t [ A(i,j,t) - (1/2)·B(i,j,t) ]
 *            + (n odd ? a[i,n]·b[n,j] : 0)
 * </pre>
 *
 * <p>Multiplication count: {@code (n²+2n−1)·⌊n/2⌋ + (n mod 2)·n²} —
 * confirmed for {@code n = 2..32} via the exact-rational cross-check
 * {@link VerifyWaksman1970} and the unit test {@code TestWaksman1970}.</p>
 *
 * <p>Non-bilinear: each rank-1 product {@code (linA + linB)·(linA' + linB')}
 * carries coefficients on BOTH operands, so this only computes matmul
 * over a commutative ring (cross-terms {@code a·a} and {@code b·b}
 * cancel after the {@code -1/2} reduction). Hence it lives under
 * {@link NonBilinearAlgorithm}, not the bilinear
 * {@link eu.solven.matmul.NonCubicBilinearAlgorithm}.</p>
 *
 * <p>Cited as the commutative baseline by DIS09 Table 4 (via
 * Mezzarobba 2007). Strictly worse than later commutative schemes
 * (Makarov 1986 ⟨3,3,3⟩=22, Hopcroft-Kerr 1971 ⟨3,3,3⟩=21, DIS09
 * compositions ≥ 11) for many {@code n}, but the unique source of a
 * uniform construction across all {@code n} ≥ 2.</p>
 */
public final class Waksman1970 {

	private Waksman1970() {}

	/**
	 * Build the explicit Waksman algorithm for {@code ⟨n, n, n⟩}.
	 *
	 * @param n side length; must be ≥ 2
	 * @return a non-bilinear algorithm with rank
	 *         {@code (n²+2n−1)·⌊n/2⌋ + (n mod 2)·n²}
	 */
	public static NonBilinearAlgorithm build(final int n) {
		if (n < 2) throw new IllegalArgumentException("n must be ≥ 2, got " + n);
		final int tmax = n / 2;
		final int r = (n * n + 2 * n - 1) * tmax + (n % 2) * n * n;
		final java.util.function.IntBinaryOperator idx = (row, col) -> row * n + col;

		double[][] Ua = new double[n * n][r];
		double[][] Ub = new double[n * n][r];
		double[][] Va = new double[n * n][r];
		double[][] Vb = new double[n * n][r];
		double[][] W  = new double[n * n][r];

		// Layout (0-indexed):
		//   block-A: a(t,i,j) at offset (t)·n² + i·n + j      , size n²·tmax
		//   block-B1: b1(t,j) at offset off1 + t·n + j        , size n·tmax
		//   block-Bii: bii(t,i') at offset off2 + t·(n-1) + i' ,
		//              where i' = i-1, i in [1..n-1]          , size (n-1)·tmax
		//   block-odd: o(i,j) at offset off3 + i·n + j (if n odd), size n²
		final int blockA = n * n * tmax;
		final int blockB1 = n * tmax;
		final int blockBii = (n - 1) * tmax;
		final int off1 = blockA;
		final int off2 = off1 + blockB1;
		final int off3 = off2 + blockBii;

		for (int t = 0; t < tmax; t++) {
			// 1-indexed paper columns:  col 2t-1 → 0-idx 2t,  col 2t → 0-idx 2t+1
			int c1 = 2 * t;       // = (2t-1) - 1
			int c2 = 2 * t + 1;   // = (2t)   - 1

			// A(t, i, j) = (a[i, c1] + b[c2, j]) · (a[i, c2] + b[c1, j])
			for (int i = 0; i < n; i++) {
				for (int j = 0; j < n; j++) {
					int k = t * n * n + i * n + j;
					Ua[idx.applyAsInt(i, c1)][k] = 1; Ub[idx.applyAsInt(c2, j)][k] = 1;
					Va[idx.applyAsInt(i, c2)][k] = 1; Vb[idx.applyAsInt(c1, j)][k] = 1;
				}
			}
			// B1(t, j) = (a[0, c1] - b[c2, j]) · (a[0, c2] - b[c1, j])
			for (int j = 0; j < n; j++) {
				int k = off1 + t * n + j;
				Ua[idx.applyAsInt(0, c1)][k] =  1; Ub[idx.applyAsInt(c2, j)][k] = -1;
				Va[idx.applyAsInt(0, c2)][k] =  1; Vb[idx.applyAsInt(c1, j)][k] = -1;
			}
			// Bii(t, i) = (a[i, c1] - b[c2, i]) · (a[i, c2] - b[c1, i])    for i in [1..n-1]
			for (int i = 1; i < n; i++) {
				int k = off2 + t * (n - 1) + (i - 1);
				Ua[idx.applyAsInt(i, c1)][k] =  1; Ub[idx.applyAsInt(c2, i)][k] = -1;
				Va[idx.applyAsInt(i, c2)][k] =  1; Vb[idx.applyAsInt(c1, i)][k] = -1;
			}
		}

		// Odd-n correction column: o(i,j) = a[i, n-1] · b[n-1, j]
		if (n % 2 == 1) {
			int last = n - 1;
			for (int i = 0; i < n; i++) {
				for (int j = 0; j < n; j++) {
					int k = off3 + i * n + j;
					Ua[idx.applyAsInt(i, last)][k] = 1;
					Vb[idx.applyAsInt(last, j)][k] = 1;
				}
			}
		}

		// W assembly: c[i,j] = Σ_t [ A_sum(t,i,j) - (1/2)·Bfull(t,i,j) ]
		//                    + (n odd) · o(i,j)
		// where (in the paper's accounting) A_sum(t,i,j) is the SUM
		// A(t,i,j); but in our rank-1-atom encoding, kA[t,i,j] holds the
		// product A_kernel(t,i,j) = (a[i,c1]+b[c2,j])·(a[i,c2]+b[c1,j])
		// which is exactly A(t,i,j) in the paper. So A_sum = A_kernel.
		// On the other hand B(t,1,j), B(t,i,i) in the paper are SUMS
		//   B(t,1,j) = A(t,0,j) + ker_B1(t,j)
		//   B(t,i,i) = A(t,i,i) + ker_Bii(t,i)
		// where the *atoms* we encoded (and that carry their own product
		// slot) are ker_B1 and ker_Bii. So when assembling W we must
		// expand Bfull in terms of the A products + ker atoms.
		//
		//   Bfull(t,0,j)        = A(t,0,j) + ker_B1(t,j)
		//   Bfull(t,i,i)  i≥1   = A(t,i,i) + ker_Bii(t,i)
		//   Bfull(t,i,j)  i≥1,j≠i
		//                       = (A(t,i,i)+ker_Bii(t,i))
		//                       - (A(t,0,i)+ker_B1(t,i))
		//                       + (A(t,0,j)+ker_B1(t,j))
		for (int t = 0; t < tmax; t++) {
			for (int i = 0; i < n; i++) {
				for (int j = 0; j < n; j++) {
					int kAij = t * n * n + i * n + j;
					int outRow = idx.applyAsInt(i, j);

					// + A(t,i,j)
					W[outRow][kAij] += 1.0;

					if (i == 0) {
						// - 1/2 · (A(t,0,j) + ker_B1(t,j))
						int kA0j = t * n * n + 0 * n + j; // == kAij here since i==0
						int kB1j = off1 + t * n + j;
						W[outRow][kA0j] += -0.5;
						W[outRow][kB1j] += -0.5;
					} else if (j == i) {
						// - 1/2 · (A(t,i,i) + ker_Bii(t,i))
						int kAii = t * n * n + i * n + i; // == kAij here since j==i
						int kBii = off2 + t * (n - 1) + (i - 1);
						W[outRow][kAii] += -0.5;
						W[outRow][kBii] += -0.5;
					} else {
						// - 1/2 · ( Bfull(t,i,i) - Bfull(t,0,i) + Bfull(t,0,j) )
						// = -1/2·A(t,i,i) -1/2·ker_Bii(t,i)
						//   +1/2·A(t,0,i) +1/2·ker_B1(t,i)
						//   -1/2·A(t,0,j) -1/2·ker_B1(t,j)
						int kAii = t * n * n + i * n + i;
						int kA0i = t * n * n + 0 * n + i;
						int kA0j = t * n * n + 0 * n + j;
						int kBii = off2 + t * (n - 1) + (i - 1);
						int kB1i = off1 + t * n + i;
						int kB1j = off1 + t * n + j;
						W[outRow][kAii] += -0.5;
						W[outRow][kBii] += -0.5;
						W[outRow][kA0i] += +0.5;
						W[outRow][kB1i] += +0.5;
						W[outRow][kA0j] += -0.5;
						W[outRow][kB1j] += -0.5;
					}
				}
			}
		}
		if (n % 2 == 1) {
			for (int i = 0; i < n; i++) {
				for (int j = 0; j < n; j++) {
					int k = off3 + i * n + j;
					W[idx.applyAsInt(i, j)][k] += 1.0;
				}
			}
		}

		return new NonBilinearAlgorithm(n, n, n, Ua, Ub, Va, Vb, W);
	}

	/** Closed-form rank: {@code (n²+2n−1)·⌊n/2⌋ + (n mod 2)·n²}. */
	public static int rank(int n) {
		return (n * n + 2 * n - 1) * (n / 2) + (n % 2) * n * n;
	}

}
