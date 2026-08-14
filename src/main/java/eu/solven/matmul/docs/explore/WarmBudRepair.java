package eu.solven.matmul.docs.explore;

import java.io.File;

import org.apache.commons.math3.linear.Array2DRowRealMatrix;
import org.apache.commons.math3.linear.RealMatrix;
import org.apache.commons.math3.linear.SingularValueDecomposition;

import eu.solven.matmul.FactorMatrix;
import eu.solven.matmul.NonCubicBilinearAlgorithm;
import eu.solven.matmul.algebra.Field;
import eu.solven.matmul.catalog.FieldAwareLookup;
import eu.solven.matmul.catalog.SchemeIO;
import eu.solven.matmul.catalog.SerendipitousBudProduct;
import eu.solven.matmul.verifiers.Verifier;
import lombok.extern.slf4j.Slf4j;

/**
 * Warm-repair tied-ALS solver: given a bud-rich base scheme, try to fuse one
 * more serendipity bud by tying two products' factor on a chosen axis and
 * re-solving the CP (matmul-tensor) decomposition by alternating least squares,
 * warm-started from the base. A converged, exactly-verified tied scheme has one
 * more bud → the serendipitous product against the inner drops by σ(k=2).
 *
 * <p>The doc's power caveat: COLD tied-ALS has no power at tight rank; this is
 * the WARM variant (seed = the incumbent, one tie perturbation) that retains
 * local power. Args: {@code baseFile  iN iM iP  targetSOTA}.</p>
 */
@Slf4j
public final class WarmBudRepair {
	private WarmBudRepair() {}

	public static void main(String[] args) throws Exception {
		NonCubicBilinearAlgorithm base = SchemeIO.read(new File(args[0]));
		int n2 = Integer.parseInt(args[1]), m2 = Integer.parseInt(args[2]), p2 = Integer.parseInt(args[3]);
		long sota = Long.parseLong(args[4]);
		FieldAwareLookup q = new FieldAwareLookup(Field.Q);
		int n = base.n, m = base.m, p = base.p, r = base.r;
		int NM = n * m, MP = m * p, NP = n * p;
		double[][] U = dense(base.u(), NM, r), V = dense(base.v(), MP, r), W = dense(base.w(), NP, r);
		double[][][] T = matmulTensor(n, m, p);
		long baseCost = SerendipitousBudProduct.serendipitousCost(base, q, n2, m2, p2);
		log.info("base ⟨{},{},{}⟩ r={}  serendip→{}  (SOTA {})  residual0={}",
				n, m, p, r, baseCost, sota, residual(U, V, W, T));

		long best = baseCost;
		NonCubicBilinearAlgorithm bestAlg = null;
		int tried = 0, converged = 0;
		// Try tying every pair of products on each of the three axes (U/V/W).
		for (int axis = 0; axis < 3; axis++) {
			for (int i = 0; i < r; i++) {
				for (int j = i + 1; j < r; j++) {
					tried++;
					double[][] u = copy(U), v = copy(V), w = copy(W);
					boolean ok = tiedAls(u, v, w, T, axis, i, j, 400);
					if (!ok) {
						continue;
					}
					converged++;
					NonCubicBilinearAlgorithm cand = rationaliseAndVerify(u, v, w, n, m, p, r);
					if (cand == null) {
						continue;
					}
					long cost = SerendipitousBudProduct.serendipitousCost(cand, q, n2, m2, p2);
					if (cost < best) {
						best = cost;
						bestAlg = cand;
						log.info("  tie axis={} ({},{}) → verified base, serendip predicts {} {}",
								"UVW".charAt(axis), i, j, cost, cost < sota ? "*** BEATS SOTA ***" : "");
					}
				}
			}
		}
		log.info("tried {} ties, {} converged (float), {} snapped to EXACT rational. "
				+ "best serendip cost = {} (SOTA {}) — {}",
				tried, converged, rationalised, best, sota,
				best < sota ? "*** IMPROVEMENT ***" : "no improvement");
		if (bestAlg != null && best < sota) {
			log.info("BUILT base rank={} exact={} — product would be ⟨{},{},{}⟩={}",
					bestAlg.r, Verifier.isExactNonCubic(bestAlg), n * n2, m * m2, p * p2, best);
		}
	}

	/** Tied ALS: enforce column i == column j on {@code axis} (0=U,1=V,2=W); returns true if
	 *  the reconstruction residual drops below tolerance. Warm-started from the passed factors. */
	private static boolean tiedAls(double[][] U, double[][] V, double[][] W, double[][][] T,
			int axis, int ci, int cj, int iters) {
		int r = U[0].length;
		// initialise the tie: average the two columns on the tied axis
		double[][] tiedF = axis == 0 ? U : axis == 1 ? V : W;
		for (int row = 0; row < tiedF.length; row++) {
			double avg = 0.5 * (tiedF[row][ci] + tiedF[row][cj]);
			tiedF[row][ci] = avg;
			tiedF[row][cj] = avg;
		}
		double prev = Double.MAX_VALUE;
		for (int it = 0; it < iters; it++) {
			// update U (tie on U handled by merged column)
			updateFactor(U, V, W, T, 0, axis == 0 ? ci : -1, axis == 0 ? cj : -1);
			updateFactor(U, V, W, T, 1, axis == 1 ? ci : -1, axis == 1 ? cj : -1);
			updateFactor(U, V, W, T, 2, axis == 2 ? ci : -1, axis == 2 ? cj : -1);
			double res = residual(U, V, W, T);
			if (res < 1e-11) {
				return true;
			}
			if (Math.abs(prev - res) < 1e-14 && it > 30) {
				return res < 1e-9;
			}
			prev = res;
		}
		return residual(U, V, W, T) < 1e-9;
	}

	/** ALS update of one factor (mode 0=U,1=V,2=W). If {@code ci>=0}, columns ci,cj of the
	 *  UPDATED factor are tied (share one variable via a merged design column). */
	private static void updateFactor(double[][] U, double[][] V, double[][] W, double[][][] T,
			int mode, int ci, int cj) {
		int r = U[0].length;
		double[][] A = mode == 0 ? V : U;      // the two OTHER factors form the Khatri-Rao design
		double[][] B = mode == 2 ? V : W;
		int dimTarget = mode == 0 ? U.length : mode == 1 ? V.length : W.length;
		int dimA = A.length, dimB = B.length;
		// design column l = kron(A[:,l], B[:,l]); merge ci,cj when tied.
		boolean tied = ci >= 0;
		int cols = tied ? r - 1 : r;
		int[] map = new int[r];              // product-index → design-column
		int nextCol = 0;
		for (int l = 0; l < r; l++) {
			if (tied && l == cj) {
				map[l] = map[ci];
			} else {
				map[l] = nextCol++;
			}
		}
		double[][] design = new double[dimA * dimB][cols];
		for (int l = 0; l < r; l++) {
			int c = map[l];
			for (int a = 0; a < dimA; a++) {
				double av = A[a][l];
				if (av == 0) {
					continue;
				}
				for (int b = 0; b < dimB; b++) {
					design[a * dimB + b][c] += av * B[b][l];   // += merges the tied pair
				}
			}
		}
		// RHS: mode-`mode` unfolding of T, dimTarget × (dimA*dimB)
		double[][] rhs = new double[dimA * dimB][dimTarget];
		for (int t = 0; t < dimTarget; t++) {
			for (int a = 0; a < dimA; a++) {
				for (int b = 0; b < dimB; b++) {
					rhs[a * dimB + b][t] = unfold(T, mode, t, a, b);
				}
			}
		}
		RealMatrix sol = new SingularValueDecomposition(new Array2DRowRealMatrix(design, false))
				.getSolver().solve(new Array2DRowRealMatrix(rhs, false));   // cols × dimTarget
		double[][] X = sol.getData();
		double[][] F = mode == 0 ? U : mode == 1 ? V : W;
		for (int l = 0; l < r; l++) {
			for (int t = 0; t < dimTarget; t++) {
				F[t][l] = X[map[l]][t];
			}
		}
	}

	private static double unfold(double[][][] T, int mode, int t, int a, int b) {
		// mode 0: T[t][a][b]; mode 1: A=U so design uses U(a)&W(b) → T[a][t][b]; mode 2: T[a][b][t]
		if (mode == 0) {
			return T[t][a][b];
		} else if (mode == 1) {
			return T[a][t][b];
		} else {
			return T[a][b][t];
		}
	}

	private static double residual(double[][] U, double[][] V, double[][] W, double[][][] T) {
		int NM = U.length, MP = V.length, NP = W.length, r = U[0].length;
		double s = 0;
		for (int a = 0; a < NM; a++) {
			for (int b = 0; b < MP; b++) {
				for (int c = 0; c < NP; c++) {
					double rec = 0;
					for (int l = 0; l < r; l++) {
						rec += U[a][l] * V[b][l] * W[c][l];
					}
					double d = rec - T[a][b][c];
					s += d * d;
				}
			}
		}
		return Math.sqrt(s);
	}

	static int rationalised = 0;

	/** Snap each factor entry to the nearest small-denominator rational, then accept only if
	 *  the SNAPPED scheme reconstructs the tensor exactly (residual≈0) and verifies. */
	private static NonCubicBilinearAlgorithm rationaliseAndVerify(double[][] U, double[][] V,
			double[][] W, int n, int m, int p, int r) {
		double[][] ru = round(U), rv = round(V), rw = round(W);
		double[][][] T = matmulTensor(n, m, p);
		if (residual(ru, rv, rw, T) > 1e-9) {
			return null;    // snapped point is off the (exact) matmul variety
		}
		rationalised++;
		try {
			NonCubicBilinearAlgorithm alg = new NonCubicBilinearAlgorithm(n, m, p, ru, rv, rw);
			return Verifier.isExactNonCubic(alg) ? alg : null;
		} catch (RuntimeException ex) {
			return null;
		}
	}

	private static final int[] DENOMS = { 1, 2, 3, 4, 6, 8, 12 };

	private static double[][] round(double[][] X) {
		double[][] out = new double[X.length][X[0].length];
		for (int i = 0; i < X.length; i++) {
			for (int j = 0; j < X[0].length; j++) {
				double v = X[i][j], bestv = Math.rint(v), bestErr = Math.abs(v - Math.rint(v));
				for (int d : DENOMS) {
					double num = Math.rint(v * d);
					double err = Math.abs(v - num / d);
					if (err < bestErr) {
						bestErr = err;
						bestv = num / d;
					}
				}
				out[i][j] = bestv;
			}
		}
		return out;
	}

	private static double[][][] matmulTensor(int n, int m, int p) {
		double[][][] T = new double[n * m][m * p][n * p];
		for (int i = 0; i < n; i++) {
			for (int j = 0; j < m; j++) {
				for (int k = 0; k < p; k++) {
					T[i * m + j][j * p + k][i * p + k] = 1;
				}
			}
		}
		return T;
	}

	private static double[][] dense(FactorMatrix f, int rows, int r) {
		double[][] out = new double[rows][r];
		for (int c = 0; c < r; c++) {
			final int col = c;
			f.forEachInColumn(c, (row, val) -> out[row][col] = val);
		}
		return out;
	}

	private static double[][] copy(double[][] a) {
		double[][] b = new double[a.length][];
		for (int i = 0; i < a.length; i++) {
			b[i] = a[i].clone();
		}
		return b;
	}
}
