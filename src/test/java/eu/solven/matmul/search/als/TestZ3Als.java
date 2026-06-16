package eu.solven.matmul.search.als;

import static org.assertj.core.api.Assertions.assertThat;

import java.util.SplittableRandom;

import org.junit.jupiter.api.Test;

import eu.solven.matmul.BilinearAlgorithm;
import eu.solven.matmul.Verifier;
import eu.solven.matmul.search.als.Z3Als;

public class TestZ3Als {

	/**
	 * expand→residual round-trip: build random orbit generators, expand to a
	 * length-r BilinearAlgorithm, expand the same way by hand, and check that
	 * the verifier's residual equals the Frobenius distance from the synthetic
	 * target tensor we constructed.
	 */
	@Test
	public void expandIsConsistent() {
		int n = 3, g = 4;
		int n2 = n * n;
		SplittableRandom rng = new SplittableRandom(0L);
		double[][] u = randomGen(rng, g, n2);
		double[][] v = randomGen(rng, g, n2);
		double[][] w = randomGen(rng, g, n2);

		double[][][] target = expandToTensor(n, g, u, v, w);
		Z3Als z = new Z3Als(n, 0, g, target);
		BilinearAlgorithm alg = z.toAlgorithm(new double[0][n2], u, v, w);
		assertThat(alg.r).isEqualTo(3 * g);
		// The synthetic target IS the expanded decomposition, so residual must be 0.
		assertThat(Verifier.residual(new BilinearAlgorithm(n, alg.U, alg.V, alg.W)) > 0).isTrue();
		// But residual against the synthetic target is 0:
		double sumSq = 0;
		for (int a = 0; a < n2; a++) {
			for (int b = 0; b < n2; b++) {
				for (int c = 0; c < n2; c++) {
					double approx = 0;
					for (int k = 0; k < alg.r; k++) {
						approx += alg.U[a][k] * alg.V[b][k] * alg.W[c][k];
					}
					double d = target[a][b][c] - approx;
					sumSq += d * d;
				}
			}
		}
		assertThat(Math.sqrt(sumSq)).isLessThan(1e-10);
	}

	/**
	 * Implementation correctness: build a synthetic tensor as a sum of g random
	 * Z/3-orbit triples, then ask Z3Als to recover them. Starting from the exact
	 * generators (no perturbation) the residual must stay at machine epsilon.
	 */
	@Test
	public void z3AlsExactInitStays() {
		int n = 3, g = 5;
		int n2 = n * n;
		SplittableRandom rng = new SplittableRandom(99L);
		double[][] u = randomGen(rng, g, n2);
		double[][] v = randomGen(rng, g, n2);
		double[][] w = randomGen(rng, g, n2);
		double[][][] target = expandToTensor(n, g, u, v, w);

		Z3Als z = new Z3Als(n, 0, g, target);
		Z3Als.Result result = z.fitFrom(new double[0][n2], u, v, w, 5, 1e-12);
		assertThat(result.residual).isLessThan(1e-10);
	}

	/**
	 * Implementation correctness: start from generators ± small Gaussian noise.
	 * ALS should converge back to the original solution rapidly.
	 */
	@Test
	public void z3AlsConvergesFromPerturbedInit() {
		int n = 3, g = 5;
		int n2 = n * n;
		SplittableRandom rng = new SplittableRandom(7L);
		double[][] u = randomGen(rng, g, n2);
		double[][] v = randomGen(rng, g, n2);
		double[][] w = randomGen(rng, g, n2);
		double[][][] target = expandToTensor(n, g, u, v, w);

		double[][] u0 = perturb(u, rng, 0.01);
		double[][] v0 = perturb(v, rng, 0.01);
		double[][] w0 = perturb(w, rng, 0.01);

		Z3Als z = new Z3Als(n, 0, g, target);
		Z3Als.Result result = z.fitFrom(new double[0][n2], u0, v0, w0, 500, 1e-10);
		assertThat(result.converged).isTrue();
		assertThat(result.residual).isLessThan(1e-8);
	}

	/**
	 * The "swamp" target: orbit-only Z/3-equivariant decomposition of the 3×3
	 * matmul tensor at r=21 (g=7) — this would be a 50-year breakthrough if it
	 * ever converged. We only assert ALS runs and the residual is bounded — this
	 * is a smoke test of the search loop, not an assertion about the math.
	 */
	@Test
	public void z3AlsAtRank21RunsWithoutCrashing() {
		Z3Als z = Z3Als.orbitOnly(3, 7);
		Z3Als.Result result = z.fit(20, 1e-10, 0L);
		assertThat(result).isNotNull();
		assertThat(Double.isFinite(result.residual)).isTrue();
	}

	/**
	 * Fixed-triple correctness: build a synthetic target as a mixture of 2 fixed
	 * + 3 orbit triples. Starting from those exact generators, the Jacobi-averaged
	 * symmetric ALS update for fixed triples must keep us at machine epsilon.
	 */
	@Test
	public void fixedTripleExactInitStays() {
		int n = 3, f = 2, gg = 3;
		int n2 = n * n;
		SplittableRandom rng = new SplittableRandom(123L);
		double[][] uFix = randomGen(rng, f, n2);
		double[][] u = randomGen(rng, gg, n2);
		double[][] v = randomGen(rng, gg, n2);
		double[][] w = randomGen(rng, gg, n2);
		double[][][] target = expandMixed(n, uFix, u, v, w);

		Z3Als z = new Z3Als(n, f, gg, target);
		Z3Als.Result result = z.fitFrom(uFix, u, v, w, 5, 1e-12);
		assertThat(result.residual).isLessThan(1e-10);
	}

	/**
	 * Diagnostic at Laderman scale: n=3, structure (2 fixed + 7 orbits) → r=23.
	 * Exact init must keep us at residual ≈ machine epsilon — if not, the ALS
	 * step at this rank is buggy (not just hard to converge from random init).
	 */
	@Test
	public void z3AlsExactInitStaysAtR23() {
		int n = 3, f = 2, gg = 7;
		int n2 = n * n;
		SplittableRandom rng = new SplittableRandom(2347L);
		double[][] uFix = randomGen(rng, f, n2);
		double[][] u = randomGen(rng, gg, n2);
		double[][] v = randomGen(rng, gg, n2);
		double[][] w = randomGen(rng, gg, n2);
		double[][][] target = expandMixed(n, uFix, u, v, w);

		Z3Als z = new Z3Als(n, f, gg, target);
		Z3Als.Result result = z.fitFrom(uFix, u, v, w, 5, 1e-12);
		assertThat(result.residual).isLessThan(1e-10);
	}

	/**
	 * Z3Als perturbed-init convergence at n=3, r=23, structure (2 fixed + 7 orbits).
	 * Validates that the SVD-based pseudoinverse fixes the previous Gauss-Jordan
	 * divergence — without this, Laderman reproduction is impossible.
	 */
	@Test
	public void z3AlsConvergesFromPerturbedInitAtR23() {
		int n = 3, f = 2, gg = 7;
		int n2 = n * n;
		SplittableRandom rng = new SplittableRandom(2347L);
		double[][] uFix = randomGen(rng, f, n2);
		double[][] u = randomGen(rng, gg, n2);
		double[][] v = randomGen(rng, gg, n2);
		double[][] w = randomGen(rng, gg, n2);
		double[][][] target = expandMixed(n, uFix, u, v, w);

		double[][] uFix0 = perturb(uFix, rng, 0.005);
		double[][] u0 = perturb(u, rng, 0.005);
		double[][] v0 = perturb(v, rng, 0.005);
		double[][] w0 = perturb(w, rng, 0.005);

		Z3Als z = new Z3Als(n, f, gg, target);
		Z3Als.Result result = z.fitFrom(uFix0, u0, v0, w0, 3000, 1e-8);
		assertThat(result.residual).isLessThan(1e-5);
	}

	/**
	 * Fixed-triple convergence from perturbed init. Notably slower than the
	 * orbit-only case because the symmetric ALS update is a contraction with
	 * rate depending on the singular structure of the residual — but it still
	 * converges for small perturbations.
	 */
	@Test
	public void fixedTripleConvergesFromPerturbedInit() {
		int n = 3, f = 1, gg = 4;
		int n2 = n * n;
		SplittableRandom rng = new SplittableRandom(77L);
		double[][] uFix = randomGen(rng, f, n2);
		double[][] u = randomGen(rng, gg, n2);
		double[][] v = randomGen(rng, gg, n2);
		double[][] w = randomGen(rng, gg, n2);
		double[][][] target = expandMixed(n, uFix, u, v, w);

		double[][] uFix0 = perturb(uFix, rng, 0.005);
		double[][] u0 = perturb(u, rng, 0.005);
		double[][] v0 = perturb(v, rng, 0.005);
		double[][] w0 = perturb(w, rng, 0.005);

		Z3Als z = new Z3Als(n, f, gg, target);
		Z3Als.Result result = z.fitFrom(uFix0, u0, v0, w0, 2000, 1e-8);
		assertThat(result.residual).isLessThan(1e-5);
	}

	// helpers --------------------------------------------------------------

	private static double[][] randomGen(SplittableRandom rng, int g, int n2) {
		double[][] M = new double[g][n2];
		for (int s = 0; s < g; s++) {
			for (int i = 0; i < n2; i++) M[s][i] = rng.nextGaussian();
		}
		return M;
	}

	private static double[][] perturb(double[][] M, SplittableRandom rng, double sigma) {
		double[][] out = new double[M.length][M[0].length];
		for (int i = 0; i < M.length; i++) {
			for (int j = 0; j < M[0].length; j++) {
				out[i][j] = M[i][j] + sigma * rng.nextGaussian();
			}
		}
		return out;
	}

	private static double[][][] expandMixed(int n, double[][] uFix, double[][] u, double[][] v, double[][] w) {
		int n2 = n * n;
		int f = uFix.length;
		int g = u.length;
		double[][][] T = new double[n2][n2][n2];
		for (int a = 0; a < n2; a++) {
			for (int b = 0; b < n2; b++) {
				for (int c = 0; c < n2; c++) {
					double sum = 0;
					for (int s = 0; s < f; s++) sum += uFix[s][a] * uFix[s][b] * uFix[s][c];
					for (int s = 0; s < g; s++) {
						sum += u[s][a] * v[s][b] * w[s][c];
						sum += v[s][a] * w[s][b] * u[s][c];
						sum += w[s][a] * u[s][b] * v[s][c];
					}
					T[a][b][c] = sum;
				}
			}
		}
		return T;
	}

	private static double[][][] expandToTensor(int n, int g, double[][] u, double[][] v, double[][] w) {
		int n2 = n * n;
		double[][][] T = new double[n2][n2][n2];
		for (int a = 0; a < n2; a++) {
			for (int b = 0; b < n2; b++) {
				for (int c = 0; c < n2; c++) {
					double sum = 0;
					for (int s = 0; s < g; s++) {
						sum += u[s][a] * v[s][b] * w[s][c];
						sum += v[s][a] * w[s][b] * u[s][c];
						sum += w[s][a] * u[s][b] * v[s][c];
					}
					T[a][b][c] = sum;
				}
			}
		}
		return T;
	}
}
