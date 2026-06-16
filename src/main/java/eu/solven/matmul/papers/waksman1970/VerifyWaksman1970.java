package eu.solven.matmul.papers.waksman1970;

import java.util.Random;

import org.apache.commons.math3.fraction.BigFraction;

/**
 * Exact-rational cross-check of the Waksman 1970 generic {@code ⟨n,n,n⟩}
 * commutative matmul (Mezzarobba 2007 thesis, Fig.3 + general formula).
 *
 * <p>This is the Java port of the former {@code tools/verify_waksman.py}.
 * It re-implements the algorithm <em>independently</em> from
 * {@link Waksman1970} (the production constructor), so it stays a genuine
 * cross-check rather than a tautology:
 * {@code TestWaksman1970} exercises {@code Waksman1970.build}; this driver
 * re-derives the same identity from the paper's formula and confirms both
 * give the correct product over ℚ.</p>
 *
 * <p>For each {@code n = 2..32} it multiplies a random integer matrix pair
 * with <strong>exact {@link BigFraction} arithmetic</strong> (no float
 * round-off), compares against the naive product, and checks the
 * multiplication count against the closed form
 * {@code (n²+2n−1)·⌊n/2⌋ + (n mod 2)·n²}. Exact arithmetic on integer inputs
 * makes a mismatch a definite defect witness (a coefficient/index bug cannot
 * hide behind floating-point noise), though a single random point is still a
 * witness, not a full symbolic proof.</p>
 *
 * <p>Run:</p>
 * <pre>
 *   mvn -q -o exec:java \
 *       -Dexec.mainClass=eu.solven.matmul.papers.waksman1970.VerifyWaksman1970
 * </pre>
 */
public final class VerifyWaksman1970 {

	private VerifyWaksman1970() {}

	/** Result of one verification run: exact-correctness flag + multiplication count. */
	private record Result(boolean ok, int mulCount) {}

	/**
	 * Compute {@code C = A·B} with Waksman's commutative algorithm over ℚ,
	 * counting the rank-1 scalar multiplications. Mirrors
	 * {@code tools/verify_waksman.py::waksman_matmul}.
	 */
	private static Result waksmanMatmul(BigFraction[][] a, BigFraction[][] b, int n) {
		int mulCount = 0;
		BigFraction half = new BigFraction(1, 2);
		BigFraction[][] c = zeros(n);
		int tmax = n / 2;
		for (int t = 1; t <= tmax; t++) {
			// A(i,j,t) = (a[i,2t-1] + b[2t,j]) · (a[i,2t] + b[2t-1,j])
			BigFraction[][] aT = new BigFraction[n][n];
			for (int i = 0; i < n; i++) {
				for (int j = 0; j < n; j++) {
					aT[i][j] = a[i][2 * t - 2].add(b[2 * t - 1][j])
							.multiply(a[i][2 * t - 1].add(b[2 * t - 2][j]));
					mulCount++;
				}
			}
			// B(1,j,t) = A(1,j,t) + (a[1,2t-1] - b[2t,j])(a[1,2t] - b[2t-1,j])
			BigFraction[] b1 = new BigFraction[n];
			for (int j = 0; j < n; j++) {
				BigFraction extra = a[0][2 * t - 2].subtract(b[2 * t - 1][j])
						.multiply(a[0][2 * t - 1].subtract(b[2 * t - 2][j]));
				mulCount++;
				b1[j] = aT[0][j].add(extra);
			}
			// B(i,i,t) for i = 2..n
			BigFraction[] bii = new BigFraction[n];
			bii[0] = b1[0];
			for (int i = 1; i < n; i++) {
				BigFraction extra = a[i][2 * t - 2].subtract(b[2 * t - 1][i])
						.multiply(a[i][2 * t - 1].subtract(b[2 * t - 2][i]));
				mulCount++;
				bii[i] = aT[i][i].add(extra);
			}
			// B(i,j,t): row 1 = b1; i≥2,j=i → bii[i]; i≥2,j≠i → bii[i] - b1[i] + b1[j]
			BigFraction[][] bFull = new BigFraction[n][n];
			for (int j = 0; j < n; j++) {
				bFull[0][j] = b1[j];
			}
			for (int i = 1; i < n; i++) {
				for (int j = 0; j < n; j++) {
					bFull[i][j] = (i == j) ? bii[i] : bii[i].subtract(b1[i]).add(b1[j]);
				}
			}
			// c[i,j] += A(i,j,t) - (1/2)·B(i,j,t)
			for (int i = 0; i < n; i++) {
				for (int j = 0; j < n; j++) {
					c[i][j] = c[i][j].add(aT[i][j].subtract(half.multiply(bFull[i][j])));
				}
			}
		}
		// Odd-dimension correction column: a[i,n]·b[n,j]
		if (n % 2 == 1) {
			for (int i = 0; i < n; i++) {
				for (int j = 0; j < n; j++) {
					c[i][j] = c[i][j].add(a[i][n - 1].multiply(b[n - 1][j]));
					mulCount++;
				}
			}
		}
		boolean ok = matchesNaive(a, b, c, n);
		return new Result(ok, mulCount);
	}

	/** Closed-form multiplication count, mirroring {@code closed_form_rank}. */
	static int closedFormRank(int n) {
		return (n * n + 2 * n - 1) * (n / 2) + (n % 2) * n * n;
	}

	private static boolean matchesNaive(BigFraction[][] a, BigFraction[][] b, BigFraction[][] c, int n) {
		for (int i = 0; i < n; i++) {
			for (int j = 0; j < n; j++) {
				BigFraction acc = BigFraction.ZERO;
				for (int k = 0; k < n; k++) {
					acc = acc.add(a[i][k].multiply(b[k][j]));
				}
				if (!acc.equals(c[i][j])) {
					return false;
				}
			}
		}
		return true;
	}

	private static BigFraction[][] zeros(int n) {
		BigFraction[][] m = new BigFraction[n][n];
		for (int i = 0; i < n; i++) {
			for (int j = 0; j < n; j++) {
				m[i][j] = BigFraction.ZERO;
			}
		}
		return m;
	}

	private static BigFraction[][] randomMatrix(int n, Random rng) {
		BigFraction[][] m = new BigFraction[n][n];
		for (int i = 0; i < n; i++) {
			for (int j = 0; j < n; j++) {
				m[i][j] = new BigFraction(rng.nextInt(11) - 5);
			}
		}
		return m;
	}

	public static void main(String[] args) {
		boolean allPass = true;
		for (int n = 2; n <= 32; n++) {
			// Fixed per-n seed → deterministic, reproducible (parity with the
			// python's random.Random(seed=0) intent: stable across runs).
			Random rng = new Random(n);
			BigFraction[][] a = randomMatrix(n, rng);
			BigFraction[][] b = randomMatrix(n, rng);
			Result r = waksmanMatmul(a, b, n);
			int expected = closedFormRank(n);
			boolean pass = r.ok() && r.mulCount() == expected;
			allPass &= pass;
			System.out.printf("n=%2d  mul_count=%5d  closed_form=%5d  verify=%s%n",
					n, r.mulCount(), expected, pass ? "PASS" : "FAIL");
		}
		if (!allPass) {
			System.err.println("FAIL: at least one n did not verify");
			System.exit(1);
		}
		System.out.println("All n=2..32 verified exact over Q.");
	}
}
