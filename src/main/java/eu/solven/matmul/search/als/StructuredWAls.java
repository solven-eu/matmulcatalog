package eu.solven.matmul.search.als;

import java.util.Random;

import eu.solven.matmul.NonCubicBilinearAlgorithm;
import lombok.extern.slf4j.Slf4j;

/**
 * CONSTRUCTIVE decision procedure for "does a rank-r scheme with a prescribed
 * W-class structure exist?" (user 2026-06-11) — tied-W alternating least
 * squares. The class structure (which products share a W direction) is FIXED
 * by a class-assignment map; the W factor has one free column per CLASS, and
 * ALS alternates exact least-squares solves of U (r columns), V (r columns)
 * and W (d ≤ r tied columns) against the matmul tensor.
 *
 * <p>Outcome semantics (label honestly):</p>
 * <ul>
 *   <li>residual → 0 with bounded factors: <b>constructive YES over ℝ</b> —
 *       a scheme with this exact class structure exists (rationalize next).</li>
 *   <li>residual stalls or factors blow up while residual shrinks: typical
 *       border-rank trap — NOT a yes; report maxAbs alongside residual.</li>
 *   <li>persistent failure across restarts: <b>evidence</b> of infeasibility,
 *       never a proof — the certified NO needs an exact method (fixed-W SAT
 *       / Gröbner) on the same structure.</li>
 * </ul>
 */
@Slf4j
public final class StructuredWAls {

	private StructuredWAls() {}

	/** @param residual  final Frobenius residual ‖T − Σ u⊗v⊗w‖
	 *  @param maxAbs    largest |entry| across factors (border-rank tell-tale) */
	public record Result(double residual, double maxAbs, double[][] u, double[][] v,
			double[][] wByClass) {

		public boolean solved() {
			return residual < 1e-9 && maxAbs < 1e3;
		}
	}

	/**
	 * @param n,m,p     target matmul format
	 * @param classOf   length-r map product→class id (0..d−1); products sharing
	 *                  a class id share their W column (the tied structure)
	 * @param warmU/V/W optional warm start (U: nm×r, V: mp×r, W: np×d); null →
	 *                  random init
	 */
	public static Result solve(int n, int m, int p, int[] classOf, long rngSeed, int maxIters,
			double[][] warmU, double[][] warmV, double[][] warmW) {
		int r = classOf.length;
		int d = java.util.Arrays.stream(classOf).max().orElseThrow() + 1;
		int dimA = n * m;
		int dimB = m * p;
		int dimC = n * p;
		double[][][] t = tensor(n, m, p);
		Random rng = new Random(rngSeed);
		double[][] u = warmU != null ? deepCopy(warmU) : gaussian(dimA, r, rng);
		double[][] v = warmV != null ? deepCopy(warmV) : gaussian(dimB, r, rng);
		double[][] w = warmW != null ? deepCopy(warmW) : gaussian(dimC, d, rng);
		double ridge = 1e-8;
		double res = Double.MAX_VALUE;
		double prev = Double.MAX_VALUE;
		int slow = 0;
		int hops = 0;
		int maxHops = Math.max(2, maxIters / 500);
		double bestRes = Double.MAX_VALUE;
		double[][] bestU = null;
		double[][] bestV = null;
		double[][] bestW = null;
		for (int it = 0; it < maxIters; it++) {
			double[][] uOld = deepCopy(u);
			double[][] vOld = deepCopy(v);
			double[][] wOld = deepCopy(w);
			solveU(t, u, v, w, classOf, dimA, dimB, dimC, r, ridge);
			solveV(t, u, v, w, classOf, dimA, dimB, dimC, r, ridge);
			solveW(t, u, v, w, classOf, dimA, dimB, dimC, d, ridge);
			res = residual(t, u, v, w, classOf, dimA, dimB, dimC, r);
			if (res < 1e-12) {
				break;
			}
			// Bro line-search: ALS swamps crawl at <0.01%/iter; overshooting
			// along the sweep direction (s ≈ ∛it) escapes them — keep only if
			// the residual actually drops.
			double s = Math.min(6, Math.cbrt(it + 1.));
			if (s > 1.01) {
				double[][] uE = extrapolate(uOld, u, s);
				double[][] vE = extrapolate(vOld, v, s);
				double[][] wE = extrapolate(wOld, w, s);
				double resE = residual(t, uE, vE, wE, classOf, dimA, dimB, dimC, r);
				if (resE < res) {
					u = uE;
					v = vE;
					w = wE;
					res = resE;
				}
			}
			// Hops are destructive moves: remember the best visited point so a
			// restart never reports (or returns) worse than it achieved.
			if (res < bestRes) {
				bestRes = res;
				bestU = deepCopy(u);
				bestV = deepCopy(v);
				bestW = deepCopy(w);
			}
			// Stalled above 1e-6: basin-hop (perturb ∝ residual, continue) a few
			// times — plain ALS gets trapped in local minima it can escape from
			// nearby (the known-feasible 4-pair control solved 0/300 without
			// hops). Hops exhausted → exit so the budget flows to other seeds.
			// Below 1e-6 the run is in a convergent basin: let it finish.
			slow = res > prev * 0.9999 ? slow + 1 : 0;
			if (slow >= 25 && res > 1e-6) {
				if (hops >= maxHops) {
					break;
				}
				hops++;
				perturb(u, rng, res);
				perturb(v, rng, res);
				perturb(w, rng, res);
				slow = 0;
				prev = Double.MAX_VALUE;
				continue;
			}
			prev = res;
			// Decay the ridge with the residual: a FIXED ridge biases every LS
			// solve and floors the residual near sqrt(ridge)·‖factors‖ (~1e-7),
			// which the sanity warm-start exposed — a true solution would then
			// be misclassified as a failure by the 1e-9 solved() gate.
			ridge = Math.max(1e-13, Math.min(1e-8, res * 1e-6));
		}
		if (bestU != null && bestRes < res) {
			u = bestU;
			v = bestV;
			w = bestW;
			res = bestRes;
		}
		return new Result(res, maxAbs(u, v, w), u, v, w);
	}

	/**
	 * Tie any slot, not just W: tying U (resp. V) of ⟨n,m,p⟩ is tying W of the
	 * cyclically-rotated tensor ⟨m,p,n⟩ (resp. ⟨p,n,m⟩) — the matmul tensor is
	 * Z/3-symmetric and rotation is a row permutation of the factors, which
	 * preserves column proportionality. Pair with {@link #expandTied}.
	 */
	public static Result solveTied(eu.solven.matmul.catalog.SerendipitousBudProduct.BudType slot,
			int n, int m, int p, int[] classOf, long rngSeed, int maxIters) {
		return switch (slot) {
		case W -> solve(n, m, p, classOf, rngSeed, maxIters, null, null, null);
		case U -> solve(m, p, n, classOf, rngSeed, maxIters, null, null, null);
		case V -> solve(p, n, m, classOf, rngSeed, maxIters, null, null, null);
		};
	}

	/** Expanded ⟨n,m,p⟩ algorithm from a {@link #solveTied} result. */
	public static NonCubicBilinearAlgorithm expandTied(
			eu.solven.matmul.catalog.SerendipitousBudProduct.BudType slot,
			int n, int m, int p, int[] classOf, Result rr) {
		return switch (slot) {
		case W -> expand(n, m, p, classOf, rr);
		case U -> expand(m, p, n, classOf, rr).cyclicShift().cyclicShift();
		case V -> expand(p, n, m, classOf, rr).cyclicShift();
		};
	}

	/** Expanded ⟨n,m,p⟩ algorithm from a solved result (W columns un-tied). */
	public static NonCubicBilinearAlgorithm expand(int n, int m, int p, int[] classOf, Result rr) {
		int r = classOf.length;
		double[][] w = new double[n * p][r];
		for (int l = 0; l < r; l++) {
			for (int c = 0; c < n * p; c++) {
				w[c][l] = rr.wByClass()[c][classOf[l]];
			}
		}
		return new NonCubicBilinearAlgorithm(n, m, p, rr.u(), rr.v(), w);
	}

	// ── ALS sub-steps (normal equations + Gaussian solve; tiny systems) ─────

	private static void solveU(double[][][] t, double[][] u, double[][] v, double[][] w,
			int[] classOf, int dimA, int dimB, int dimC, int r, double ridge) {
		// Design K[bc][l] = v[b][l]·w[c][class(l)]; solve per row a.
		double[][] gram = new double[r][r];
		for (int l1 = 0; l1 < r; l1++) {
			for (int l2 = 0; l2 < r; l2++) {
				double sv = 0;
				for (int b = 0; b < dimB; b++) {
					sv += v[b][l1] * v[b][l2];
				}
				double sw = 0;
				for (int c = 0; c < dimC; c++) {
					sw += w[c][classOf[l1]] * w[c][classOf[l2]];
				}
				gram[l1][l2] = sv * sw + (l1 == l2 ? ridge : 0);
			}
		}
		for (int a = 0; a < dimA; a++) {
			double[] rhs = new double[r];
			for (int l = 0; l < r; l++) {
				double s = 0;
				for (int b = 0; b < dimB; b++) {
					for (int c = 0; c < dimC; c++) {
						s += t[a][b][c] * v[b][l] * w[c][classOf[l]];
					}
				}
				rhs[l] = s;
			}
			double[] x = gaussSolve(gram, rhs);
			for (int l = 0; l < r; l++) {
				u[a][l] = x[l];
			}
		}
	}

	private static void solveV(double[][][] t, double[][] u, double[][] v, double[][] w,
			int[] classOf, int dimA, int dimB, int dimC, int r, double ridge) {
		double[][] gram = new double[r][r];
		for (int l1 = 0; l1 < r; l1++) {
			for (int l2 = 0; l2 < r; l2++) {
				double su = 0;
				for (int a = 0; a < dimA; a++) {
					su += u[a][l1] * u[a][l2];
				}
				double sw = 0;
				for (int c = 0; c < dimC; c++) {
					sw += w[c][classOf[l1]] * w[c][classOf[l2]];
				}
				gram[l1][l2] = su * sw + (l1 == l2 ? ridge : 0);
			}
		}
		for (int b = 0; b < dimB; b++) {
			double[] rhs = new double[r];
			for (int l = 0; l < r; l++) {
				double s = 0;
				for (int a = 0; a < dimA; a++) {
					for (int c = 0; c < dimC; c++) {
						s += t[a][b][c] * u[a][l] * w[c][classOf[l]];
					}
				}
				rhs[l] = s;
			}
			double[] x = gaussSolve(gram, rhs);
			for (int l = 0; l < r; l++) {
				v[b][l] = x[l];
			}
		}
	}

	private static void solveW(double[][][] t, double[][] u, double[][] v, double[][] w,
			int[] classOf, int dimA, int dimB, int dimC, int d, double ridge) {
		// Design G[ab][cls] = Σ_{l∈cls} u[a][l]·v[b][l]; solve per output row c.
		int r = classOf.length;
		double[][] gram = new double[d][d];
		double[][][] g = new double[dimA][dimB][d];
		for (int a = 0; a < dimA; a++) {
			for (int b = 0; b < dimB; b++) {
				for (int l = 0; l < r; l++) {
					g[a][b][classOf[l]] += u[a][l] * v[b][l];
				}
			}
		}
		for (int c1 = 0; c1 < d; c1++) {
			for (int c2 = 0; c2 < d; c2++) {
				double s = 0;
				for (int a = 0; a < dimA; a++) {
					for (int b = 0; b < dimB; b++) {
						s += g[a][b][c1] * g[a][b][c2];
					}
				}
				gram[c1][c2] = s + (c1 == c2 ? ridge : 0);
			}
		}
		for (int c = 0; c < dimC; c++) {
			double[] rhs = new double[d];
			for (int cls = 0; cls < d; cls++) {
				double s = 0;
				for (int a = 0; a < dimA; a++) {
					for (int b = 0; b < dimB; b++) {
						s += t[a][b][c] * g[a][b][cls];
					}
				}
				rhs[cls] = s;
			}
			double[] x = gaussSolve(gram, rhs);
			for (int cls = 0; cls < d; cls++) {
				w[c][cls] = x[cls];
			}
		}
	}

	private static double residual(double[][][] t, double[][] u, double[][] v, double[][] w,
			int[] classOf, int dimA, int dimB, int dimC, int r) {
		double s = 0;
		for (int a = 0; a < dimA; a++) {
			for (int b = 0; b < dimB; b++) {
				for (int c = 0; c < dimC; c++) {
					double e = -t[a][b][c];
					for (int l = 0; l < r; l++) {
						e += u[a][l] * v[b][l] * w[c][classOf[l]];
					}
					s += e * e;
				}
			}
		}
		return Math.sqrt(s);
	}

	static double[][][] tensor(int n, int m, int p) {
		double[][][] t = new double[n * m][m * p][n * p];
		for (int i = 0; i < n; i++) {
			for (int j = 0; j < m; j++) {
				for (int k = 0; k < p; k++) {
					t[i * m + j][j * p + k][i * p + k] = 1;
				}
			}
		}
		return t;
	}

	private static double[] gaussSolve(double[][] aIn, double[] bIn) {
		int nn = bIn.length;
		double[][] a = deepCopy(aIn);
		double[] b = bIn.clone();
		for (int col = 0; col < nn; col++) {
			int piv = col;
			for (int row = col + 1; row < nn; row++) {
				if (Math.abs(a[row][col]) > Math.abs(a[piv][col])) {
					piv = row;
				}
			}
			double[] tmp = a[col];
			a[col] = a[piv];
			a[piv] = tmp;
			double tb = b[col];
			b[col] = b[piv];
			b[piv] = tb;
			double diag = a[col][col];
			if (Math.abs(diag) < 1e-14) {
				continue;  // ridge keeps this rare; skip degenerate pivot
			}
			for (int row = col + 1; row < nn; row++) {
				double f = a[row][col] / diag;
				for (int k = col; k < nn; k++) {
					a[row][k] -= f * a[col][k];
				}
				b[row] -= f * b[col];
			}
		}
		double[] x = new double[nn];
		for (int row = nn - 1; row >= 0; row--) {
			double s = b[row];
			for (int k = row + 1; k < nn; k++) {
				s -= a[row][k] * x[k];
			}
			x[row] = Math.abs(a[row][row]) < 1e-14 ? 0 : s / a[row][row];
		}
		return x;
	}

	private static double[][] gaussian(int rows, int cols, Random rng) {
		double[][] out = new double[rows][cols];
		for (int i = 0; i < rows; i++) {
			for (int j = 0; j < cols; j++) {
				out[i][j] = rng.nextGaussian();
			}
		}
		return out;
	}

	private static void perturb(double[][] x, Random rng, double res) {
		// Scale STRICTLY with the residual: a floor would blast nearly-converged
		// runs (a 4e-6 stall) out of their basin instead of nudging them.
		for (double[] row : x) {
			for (int j = 0; j < row.length; j++) {
				row[j] += res * rng.nextGaussian();
			}
		}
	}

	private static double[][] extrapolate(double[][] old, double[][] cur, double s) {
		double[][] out = new double[cur.length][cur[0].length];
		for (int i = 0; i < cur.length; i++) {
			for (int j = 0; j < cur[0].length; j++) {
				out[i][j] = old[i][j] + s * (cur[i][j] - old[i][j]);
			}
		}
		return out;
	}

	private static double[][] deepCopy(double[][] in) {
		double[][] out = new double[in.length][];
		for (int i = 0; i < in.length; i++) {
			out[i] = in[i].clone();
		}
		return out;
	}

	private static double maxAbs(double[][]... mats) {
		double max = 0;
		for (double[][] mat : mats) {
			for (double[] row : mat) {
				for (double x : row) {
					max = Math.max(max, Math.abs(x));
				}
			}
		}
		return max;
	}
}
