package eu.solven.matmul.search.als;

import eu.solven.matmul.Verifier;

import eu.solven.matmul.BilinearAlgorithm;

import java.util.SplittableRandom;

import org.apache.commons.math3.linear.Array2DRowRealMatrix;
import org.apache.commons.math3.linear.ArrayRealVector;
import org.apache.commons.math3.linear.RealMatrix;
import org.apache.commons.math3.linear.RealVector;
import org.apache.commons.math3.linear.SingularValueDecomposition;

/**
 * Z/3-equivariant ALS supporting both fixed triples and orbit generators.
 *
 * Parametrization (f + 3g = r):
 * <ul>
 *   <li>{@code f} fixed triples — each a single vector u_fix[s] ∈ R^{n²}.
 *       Contribution: u_fix[s][a] · u_fix[s][b] · u_fix[s][c] (one rank-1 term per fixed triple).</li>
 *   <li>{@code g} orbit generators — each a triple (u_s, v_s, w_s) ∈ R^{n²}³.
 *       Unfolds under Z/3 to three rank-1 terms (see {@link #toAlgorithm}).</li>
 * </ul>
 *
 * From {@code RANK_3X3_SEARCH.md} §4.1, the allowed structures at the
 * interesting ranks for n=3 are:
 * <ul>
 *   <li>r=21: 0+7, 3+6, 6+5, 9+4, 12+3, 15+2, 18+1, 21+0</li>
 *   <li>r=22: 1+7, 4+6, 7+5, 10+4, 13+3, 16+2, 19+1 (always ≥1 fixed: 22 ≡ 1 mod 3)</li>
 *   <li>r=23: 2+7, 5+6, … (always ≥2 fixed)</li>
 * </ul>
 *
 * ALS scheme:
 * <ul>
 *   <li>Orbit generator update: solve n²-dim linear LS — generator appears
 *       linearly in three rank-1 terms (one in each slot under cyclic rotation).</li>
 *   <li>Fixed triple update: the vector appears cubically (u⊗u⊗u). Use
 *       Jacobi-averaged symmetric ALS: solve three single-slot LS problems
 *       (holding the other two slots at u_old) and average the three solutions.
 *       Converges to a local minimum of ‖R − u⊗u⊗u‖².</li>
 * </ul>
 */
public class Z3Als {

	private static final double MIN_DENOM = 1e-30;

	public static class Result {
		public final BilinearAlgorithm algorithm;
		public final double residual;
		public final int iterations;
		public final boolean converged;

		public Result(BilinearAlgorithm a, double res, int it, boolean ok) {
			this.algorithm = a;
			this.residual = res;
			this.iterations = it;
			this.converged = ok;
		}
	}

	private final int n;
	private final int n2;
	private final int f;
	private final int g;
	private final int r;
	private final double[][][] T;

	/**
	 * Targets {@link Verifier#trilinTensor(int)} — the cyclic-symmetric form
	 * trace(A·B·C). The raw Z/3 cyclic action {@code (u,v,w) → (v,w,u)} makes
	 * orbit contributions cyclic-symmetric in the (a, b, c) tensor indices,
	 * and only {@code trilinTensor} is itself cyclic-symmetric in those indices
	 * — {@code matmulTensor} is not. Targeting {@code matmulTensor} directly
	 * with this action would search a tensor that has no Z/3-equivariant
	 * decomposition under the raw action, no matter the rank.
	 *
	 * A rank-r decomposition of {@code trilinTensor(n)} converts to a rank-r
	 * decomposition of {@code matmulTensor(n)} by transposing the W factor;
	 * see {@link Verifier#transposeW(BilinearAlgorithm)}.
	 */
	public Z3Als(int n, int f, int g) {
		this(n, f, g, Verifier.trilinTensor(n));
	}

	public Z3Als(int n, int f, int g, double[][][] target) {
		this.n = n;
		this.n2 = n * n;
		this.f = f;
		this.g = g;
		this.r = f + 3 * g;
		this.T = target;
	}

	/** Orbit-only convenience constructor (f = 0). */
	public static Z3Als orbitOnly(int n, int g) {
		return new Z3Als(n, 0, g);
	}

	public Result fit(int maxIters, double tol, long seed) {
		SplittableRandom rng = new SplittableRandom(seed);
		return fitFrom(
				randomFactor(rng, f),
				randomFactor(rng, g),
				randomFactor(rng, g),
				randomFactor(rng, g),
				maxIters,
				tol);
	}

	public Result fitFrom(double[][] uFix0, double[][] uGen0, double[][] vGen0, double[][] wGen0,
			int maxIters, double tol) {
		double[][] uFix = copy(uFix0);
		double[][] u = copy(uGen0);
		double[][] v = copy(vGen0);
		double[][] w = copy(wGen0);

		double prev = Double.POSITIVE_INFINITY;
		for (int it = 0; it < maxIters; it++) {
			updateOrbitBlock(uFix, u, v, w, 0);
			updateOrbitBlock(uFix, u, v, w, 1);
			updateOrbitBlock(uFix, u, v, w, 2);
			updateFixedBlock(uFix, u, v, w);

			double res = residual(uFix, u, v, w);
			if (res < tol) return new Result(toAlgorithm(uFix, u, v, w), res, it + 1, true);
			if (Math.abs(prev - res) < 1e-15 && it > 5) {
				return new Result(toAlgorithm(uFix, u, v, w), res, it + 1, false);
			}
			prev = res;
		}
		return new Result(toAlgorithm(uFix, u, v, w), prev, maxIters, false);
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

	// ----- orbit generator update (linear LS, dimension n²) -----------------

	private void updateOrbitBlock(double[][] uFix, double[][] u, double[][] v, double[][] w, int which) {
		double[][][] residualWithAllTerms = computeResidual(uFix, u, v, w);
		for (int s = 0; s < g; s++) {
			double[][][] Rs = addOrbit(residualWithAllTerms, u[s], v[s], w[s]);
			double[] p, q;
			switch (which) {
				case 0: p = v[s]; q = w[s]; break;
				case 1: p = w[s]; q = u[s]; break;
				case 2: p = u[s]; q = v[s]; break;
				default: throw new IllegalStateException();
			}
			double[] newGen = solveOrbitGenerator(Rs, p, q);
			switch (which) {
				case 0: u[s] = newGen; break;
				case 1: v[s] = newGen; break;
				case 2: w[s] = newGen; break;
			}
			residualWithAllTerms = computeResidual(uFix, u, v, w);
		}
	}

	/**
	 * Solves the n²-dim LS problem for one orbit generator given the residual Rs
	 * (which still includes this orbit's contribution) and the two other generators
	 * p, q (per the cyclic-rotation convention in {@link #updateOrbitBlock}).
	 *
	 * grad_i[a,b,c] = δ(a=i)·p[b]·q[c] + δ(b=i)·q[a]·p[c] + δ(c=i)·p[a]·q[b]
	 */
	private double[] solveOrbitGenerator(double[][][] Rs, double[] p, double[] q) {
		double[][] A = new double[n2][n2];
		double[] rhs = new double[n2];
		double[][][][] grad = new double[n2][n2][n2][n2]; // grad[i][a][b][c]
		for (int i = 0; i < n2; i++) {
			for (int b = 0; b < n2; b++) {
				for (int c = 0; c < n2; c++) grad[i][i][b][c] += p[b] * q[c];
			}
			for (int a = 0; a < n2; a++) {
				for (int c = 0; c < n2; c++) grad[i][a][i][c] += q[a] * p[c];
			}
			for (int a = 0; a < n2; a++) {
				for (int b = 0; b < n2; b++) grad[i][a][b][i] += p[a] * q[b];
			}
		}
		for (int i = 0; i < n2; i++) {
			double rhsI = 0;
			for (int a = 0; a < n2; a++) {
				for (int b = 0; b < n2; b++) {
					for (int c = 0; c < n2; c++) rhsI += Rs[a][b][c] * grad[i][a][b][c];
				}
			}
			rhs[i] = rhsI;
			for (int j = i; j < n2; j++) {
				double aij = 0;
				for (int a = 0; a < n2; a++) {
					for (int b = 0; b < n2; b++) {
						for (int c = 0; c < n2; c++) aij += grad[i][a][b][c] * grad[j][a][b][c];
					}
				}
				A[i][j] = aij;
				if (i != j) A[j][i] = aij;
			}
		}
		return solveSpd(A, rhs);
	}

	// ----- fixed triple update (symmetric ALS, Jacobi-averaged) -------------

	/**
	 * Damped symmetric Newton step for each fixed triple. Linearizing the
	 * symmetric rank-1 contribution {@code u_new⊗u_new⊗u_new} around
	 * {@code u_old} gives
	 * <pre>
	 *     u_new⊗u_new⊗u_new ≈ u_old⊗u_old⊗u_old + 3 · sym(δu, u_old, u_old)
	 * </pre>
	 * where {@code δu = u_new − u_old}. Solving the linear LS for {@code δu}
	 * against the asymmetric residual gives
	 * {@code δu = (uJacobi − u_old) / 3}, where
	 * {@code uJacobi = (uU + uV + uW) / 3} is the unconstrained per-slot LS
	 * average. So the correct symmetric ALS step is
	 * <pre>
	 *     u_new = u_old + (uJacobi − u_old) / 3
	 *           = (2/3) · u_old + (1/3) · uJacobi.
	 * </pre>
	 *
	 * Previously this code applied {@code u_new = uJacobi} (3× over-shoot),
	 * which had a stable wrong attractor at higher rank because the symmetric
	 * tensor problem is non-convex.
	 */
	private void updateFixedBlock(double[][] uFix, double[][] u, double[][] v, double[][] w) {
		double[][][] residualWithAllTerms = computeResidual(uFix, u, v, w);
		for (int s = 0; s < f; s++) {
			double[] uOld = uFix[s];
			double[][][] Rs = addFixed(residualWithAllTerms, uOld);

			double norm2 = 0;
			for (int i = 0; i < n2; i++) norm2 += uOld[i] * uOld[i];
			double denom = norm2 * norm2;
			if (denom < MIN_DENOM) denom = MIN_DENOM;

			double[] uU = new double[n2];
			double[] uV = new double[n2];
			double[] uW = new double[n2];
			for (int i = 0; i < n2; i++) {
				double sU = 0, sV = 0, sW = 0;
				for (int a = 0; a < n2; a++) {
					for (int b = 0; b < n2; b++) {
						sU += Rs[i][a][b] * uOld[a] * uOld[b];
						sV += Rs[a][i][b] * uOld[a] * uOld[b];
						sW += Rs[a][b][i] * uOld[a] * uOld[b];
					}
				}
				uU[i] = sU / denom;
				uV[i] = sV / denom;
				uW[i] = sW / denom;
			}
			double[] uNew = new double[n2];
			for (int i = 0; i < n2; i++) {
				double uJacobi = (uU[i] + uV[i] + uW[i]) / 3.0;
				uNew[i] = uOld[i] + (uJacobi - uOld[i]) / 3.0;
			}
			uFix[s] = uNew;
			residualWithAllTerms = computeResidual(uFix, u, v, w);
		}
	}

	// ----- residual / contribution bookkeeping ------------------------------

	private double[][][] computeResidual(double[][] uFix, double[][] u, double[][] v, double[][] w) {
		double[][][] R = new double[n2][n2][n2];
		for (int a = 0; a < n2; a++) {
			for (int b = 0; b < n2; b++) {
				for (int c = 0; c < n2; c++) {
					double approx = 0;
					for (int s = 0; s < f; s++) {
						approx += uFix[s][a] * uFix[s][b] * uFix[s][c];
					}
					for (int s = 0; s < g; s++) {
						approx += u[s][a] * v[s][b] * w[s][c];
						approx += v[s][a] * w[s][b] * u[s][c];
						approx += w[s][a] * u[s][b] * v[s][c];
					}
					R[a][b][c] = T[a][b][c] - approx;
				}
			}
		}
		return R;
	}

	private double[][][] addOrbit(double[][][] R, double[] u, double[] v, double[] w) {
		double[][][] out = new double[n2][n2][n2];
		for (int a = 0; a < n2; a++) {
			for (int b = 0; b < n2; b++) {
				for (int c = 0; c < n2; c++) {
					out[a][b][c] = R[a][b][c]
							+ u[a] * v[b] * w[c]
							+ v[a] * w[b] * u[c]
							+ w[a] * u[b] * v[c];
				}
			}
		}
		return out;
	}

	private double[][][] addFixed(double[][][] R, double[] u) {
		double[][][] out = new double[n2][n2][n2];
		for (int a = 0; a < n2; a++) {
			for (int b = 0; b < n2; b++) {
				for (int c = 0; c < n2; c++) {
					out[a][b][c] = R[a][b][c] + u[a] * u[b] * u[c];
				}
			}
		}
		return out;
	}

	private double residual(double[][] uFix, double[][] u, double[][] v, double[][] w) {
		double sumSq = 0;
		for (int a = 0; a < n2; a++) {
			for (int b = 0; b < n2; b++) {
				for (int c = 0; c < n2; c++) {
					double approx = 0;
					for (int s = 0; s < f; s++) {
						approx += uFix[s][a] * uFix[s][b] * uFix[s][c];
					}
					for (int s = 0; s < g; s++) {
						approx += u[s][a] * v[s][b] * w[s][c];
						approx += v[s][a] * w[s][b] * u[s][c];
						approx += w[s][a] * u[s][b] * v[s][c];
					}
					double d = T[a][b][c] - approx;
					sumSq += d * d;
				}
			}
		}
		return Math.sqrt(sumSq);
	}

	/**
	 * Unfolds (uFix, u, v, w) into a length-r BilinearAlgorithm.
	 * Column layout: first f columns are the fixed triples, then 3g orbit columns
	 * grouped per orbit (cols f+3s, f+3s+1, f+3s+2 = three cyclic shifts of generator s).
	 */
	public BilinearAlgorithm toAlgorithm(double[][] uFix, double[][] u, double[][] v, double[][] w) {
		double[][] U = new double[n2][r];
		double[][] V = new double[n2][r];
		double[][] W = new double[n2][r];
		for (int s = 0; s < f; s++) {
			for (int a = 0; a < n2; a++) {
				U[a][s] = uFix[s][a];
				V[a][s] = uFix[s][a];
				W[a][s] = uFix[s][a];
			}
		}
		for (int s = 0; s < g; s++) {
			int base = f + 3 * s;
			for (int a = 0; a < n2; a++) {
				U[a][base + 0] = u[s][a]; V[a][base + 0] = v[s][a]; W[a][base + 0] = w[s][a];
				U[a][base + 1] = v[s][a]; V[a][base + 1] = w[s][a]; W[a][base + 1] = u[s][a];
				U[a][base + 2] = w[s][a]; V[a][base + 2] = u[s][a]; W[a][base + 2] = v[s][a];
			}
		}
		return new BilinearAlgorithm(n, U, V, W);
	}

	// ----- utilities --------------------------------------------------------

	/**
	 * Solves {@code A x = rhs} via the SVD-based Moore–Penrose pseudoinverse
	 * (Apache Commons Math). Numerically robust at any condition number with
	 * no LS bias — replaces the previous Gauss-Jordan inversion which diverged
	 * at higher ranks due to ill-conditioning amplification.
	 *
	 * The slight allocation overhead per call (~1 µs for an n²×n² = 9×9 system
	 * at n=3) is dwarfed by the rest of the ALS iteration.
	 */
	private double[] solveSpd(double[][] A, double[] rhs) {
		RealMatrix matrix = new Array2DRowRealMatrix(A, false);
		SingularValueDecomposition svd = new SingularValueDecomposition(matrix);
		RealVector solution = svd.getSolver().solve(new ArrayRealVector(rhs, false));
		return solution.toArray();
	}

	private double[][] randomFactor(SplittableRandom rng, int rows) {
		double[][] M = new double[rows][n2];
		for (int s = 0; s < rows; s++) {
			for (int i = 0; i < n2; i++) M[s][i] = rng.nextGaussian();
		}
		return M;
	}

	private static double[][] copy(double[][] M) {
		double[][] out = new double[M.length][];
		for (int i = 0; i < M.length; i++) out[i] = M[i].clone();
		return out;
	}
}
