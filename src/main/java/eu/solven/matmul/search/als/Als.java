package eu.solven.matmul.search.als;

import eu.solven.matmul.LinAlg;

import eu.solven.matmul.Verifier;

import eu.solven.matmul.BilinearAlgorithm;

import java.util.SplittableRandom;

/**
 * Alternating Least Squares for rank-r decomposition of the n×n matmul tensor
 * T ∈ R^{n²×n²×n²}.
 *
 * Iterates: hold two of (U, V, W) fixed and solve the third as a linear least-
 * squares problem. Each sub-problem is independent across the n² rows of the
 * factor being updated, so it reduces to a single r×r Gram-matrix solve shared
 * across rows.
 *
 * Update formula for U, holding V, W fixed:
 *
 *   M[k][bc] = V[b][k] · W[c][k]       (shape r × n⁴)
 *   gram = M Mᵀ                        (r × r)
 *   rhs[a][k] = Σ_{b,c} T[a][b][c] · M[k][bc]
 *   U = rhs · gram⁻¹
 *
 * Symmetric formulas update V and W. A small ridge {@code RIDGE * I} is added to
 * the Gram matrix to keep it invertible across iterations where one factor
 * happens to drop rank.
 *
 * The objective is non-convex; ALS converges only to a local minimum. Use
 * {@link #fitWithRestarts} with many random inits to navigate the landscape.
 */
public class Als {

	/** Fallback ridge only used when the un-regularized Gram is numerically singular. */
	private static final double FALLBACK_RIDGE = 1e-10;

	public static class Result {
		public final BilinearAlgorithm algorithm;
		public final double residual;
		public final int iterations;
		public final boolean converged;

		public Result(BilinearAlgorithm algorithm, double residual, int iterations, boolean converged) {
			this.algorithm = algorithm;
			this.residual = residual;
			this.iterations = iterations;
			this.converged = converged;
		}
	}

	private final int n;
	private final int n2;
	private final int r;
	private final double[][][] T;

	/** ALS targeting the n×n matmul tensor at rank r. */
	public Als(int n, int r) {
		this(n, r, Verifier.matmulTensor(n));
	}

	/** ALS targeting an arbitrary {@code n²×n²×n²} tensor — useful for tests with synthetic targets. */
	public Als(int n, int r, double[][][] target) {
		this.n = n;
		this.n2 = n * n;
		this.r = r;
		this.T = target;
	}

	public Result fit(int maxIters, double tol, long seed) {
		SplittableRandom rng = new SplittableRandom(seed);
		return fitFrom(randomFactor(rng), randomFactor(rng), randomFactor(rng), maxIters, tol);
	}

	/**
	 * Runs ALS starting from explicit factors. The arrays are copied, so the caller
	 * can reuse them. Useful for testing implementation correctness (initialize
	 * near a known solution → expect rapid convergence) and for warm-starting from
	 * structured inits.
	 */
	public Result fitFrom(double[][] U0, double[][] V0, double[][] W0, int maxIters, double tol) {
		double[][] U = copy(U0);
		double[][] V = copy(V0);
		double[][] W = copy(W0);
		double prev = Double.POSITIVE_INFINITY;
		for (int it = 0; it < maxIters; it++) {
			U = solveFactor(V, W, 0);
			V = solveFactor(U, W, 1);
			W = solveFactor(U, V, 2);
			double res = residual(U, V, W);
			if (res < tol) {
				return new Result(new BilinearAlgorithm(n, U, V, W), res, it + 1, true);
			}
			if (Math.abs(prev - res) < 1e-15 && it > 5) {
				return new Result(new BilinearAlgorithm(n, U, V, W), res, it + 1, false);
			}
			prev = res;
		}
		return new Result(new BilinearAlgorithm(n, U, V, W), prev, maxIters, false);
	}

	public Result fitWithRestarts(int nRestarts, int maxIters, double tol, long baseSeed) {
		Result best = null;
		for (int s = 0; s < nRestarts; s++) {
			Result candidate = fit(maxIters, tol, baseSeed + s);
			if (best == null || candidate.residual < best.residual) {
				best = candidate;
				if (best.converged) return best;
			}
		}
		return best;
	}

	private double[][] randomFactor(SplittableRandom rng) {
		double[][] F = new double[n2][r];
		for (int a = 0; a < n2; a++) {
			for (int k = 0; k < r; k++) {
				F[a][k] = rng.nextGaussian();
			}
		}
		return F;
	}

	private static double[][] copy(double[][] M) {
		double[][] out = new double[M.length][];
		for (int i = 0; i < M.length; i++) out[i] = M[i].clone();
		return out;
	}

	/**
	 * Solves for one factor (axis = 0 → U, 1 → V, 2 → W) given the other two.
	 * The two input factors are passed in their natural order — the {@code axis}
	 * parameter selects which slot we are solving for and tells the inner loops
	 * how to index into T.
	 */
	private double[][] solveFactor(double[][] F1, double[][] F2, int axis) {
		// gram[i][j] = Σ_{a,b} F1[a][i]·F2[b][i] · F1[a][j]·F2[b][j]
		//            = (F1ᵀ F1)[i][j] · (F2ᵀ F2)[i][j]   (Hadamard of two Gram mats)
		double[][] g1 = gram(F1);
		double[][] g2 = gram(F2);
		double[][] gram = new double[r][r];
		for (int i = 0; i < r; i++) {
			for (int j = 0; j < r; j++) {
				gram[i][j] = g1[i][j] * g2[i][j];
			}
		}
		double[][] gramInv;
		try {
			gramInv = LinAlg.invert(gram);
		} catch (ArithmeticException singular) {
			// Rank-deficient step — fall back to ridge regression. Biased but recoverable.
			for (int i = 0; i < r; i++) gram[i][i] += FALLBACK_RIDGE;
			gramInv = LinAlg.invert(gram);
		}

		// rhs[x][k] = Σ_{a,b} T(axis-slice)[x, a, b] · F1[a][k] · F2[b][k]
		double[][] rhs = new double[n2][r];
		for (int x = 0; x < n2; x++) {
			for (int a = 0; a < n2; a++) {
				for (int b = 0; b < n2; b++) {
					double tval = sliceT(axis, x, a, b);
					if (tval == 0.0) continue;
					double[] F1a = F1[a];
					double[] F2b = F2[b];
					double[] rhsX = rhs[x];
					for (int k = 0; k < r; k++) {
						rhsX[k] += tval * F1a[k] * F2b[k];
					}
				}
			}
		}
		return LinAlg.mul(rhs, gramInv);
	}

	private double sliceT(int axis, int x, int a, int b) {
		// axis=0: solving for U. T-slice has U-index x, V-index a, W-index b.
		// axis=1: solving for V. T-slice has V-index x — so T[a][x][b].
		// axis=2: solving for W. T-slice has W-index x — so T[a][b][x].
		switch (axis) {
			case 0: return T[x][a][b];
			case 1: return T[a][x][b];
			case 2: return T[a][b][x];
			default: throw new IllegalArgumentException("axis: " + axis);
		}
	}

	private double[][] gram(double[][] F) {
		double[][] g = new double[r][r];
		for (int a = 0; a < n2; a++) {
			double[] Fa = F[a];
			for (int i = 0; i < r; i++) {
				double fi = Fa[i];
				if (fi == 0.0) continue;
				for (int j = 0; j < r; j++) {
					g[i][j] += fi * Fa[j];
				}
			}
		}
		return g;
	}

	private double residual(double[][] U, double[][] V, double[][] W) {
		double sumSq = 0.0;
		for (int a = 0; a < n2; a++) {
			for (int b = 0; b < n2; b++) {
				for (int c = 0; c < n2; c++) {
					double approx = 0.0;
					for (int k = 0; k < r; k++) {
						approx += U[a][k] * V[b][k] * W[c][k];
					}
					double d = T[a][b][c] - approx;
					sumSq += d * d;
				}
			}
		}
		return Math.sqrt(sumSq);
	}
}
