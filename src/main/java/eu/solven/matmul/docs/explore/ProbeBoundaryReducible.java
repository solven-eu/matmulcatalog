package eu.solven.matmul.docs.explore;

import java.util.Random;

import org.apache.commons.math3.linear.Array2DRowRealMatrix;
import org.apache.commons.math3.linear.RealMatrix;
import org.apache.commons.math3.linear.SingularValueDecomposition;

import eu.solven.matmul.NonCubicBilinearAlgorithm;
import eu.solven.matmul.papers.schwartzzwecher2025.TaNew25Construction;

/**
 * Boundary-reducibility test for the SZ aggregation {@code 2·N²} tail. The tail is
 * ENTIRELY the two-equal-index ("boundary") aggregation products (see
 * {@code ProbeAggTail}). Let {@code R_bd = Σ(boundary products)} — the tensor those
 * products must produce with everything else (generic + correction + diagonal) held
 * fixed, since {@code Σ(all products) = T_matmul} exactly. If {@code rank(R_bd) <
 * #boundary}, replacing the boundary block with a minimal decomposition yields a
 * VALID scheme of strictly lower rank — a direct cut to the frozen {@code 15/4} N²
 * coefficient. We CP-ALS {@code R_bd} at target ranks below {@code #boundary}; a
 * residual→0 at rank {@code < #boundary} proves reducibility (constructively).
 *
 * <p>Args: {@code n0 [minRankFraction]} (default scans down to 0.6·#boundary).</p>
 */
public final class ProbeBoundaryReducible {
	private ProbeBoundaryReducible() {}

	public static void main(String[] args) {
		int n0 = Integer.parseInt(args[0]);
		double floorFrac = args.length > 1 ? Double.parseDouble(args[1]) : 0.6;
		TaNew25Construction.Tagged tg = TaNew25Construction.buildTagged(n0);
		NonCubicBilinearAlgorithm alg = tg.alg();
		int[] cls = tg.productClass();
		int r = alg.r, dim = n0 * n0;
		double[][] U = alg.denseU(), V = alg.denseV(), W = alg.denseW();

		int[] hist = new int[4];
		for (int c : cls) hist[c]++;
		System.out.printf("%n=== SZ ⟨%d³⟩ boundary-reducibility : r=%d ===%n", n0, r);
		System.out.printf("classes: generic=%d  boundary=%d  correction=%d  diagonal=%d%n",
				hist[0], hist[1], hist[2], hist[3]);
		int B = hist[1];

		// R_bd = Σ over boundary products of u⊗v⊗w  (dim×dim×dim)
		double[][][] T = new double[dim][dim][dim];
		double nrm2 = 0;
		for (int l = 0; l < r; l++) {
			if (cls[l] != 1) continue;
			for (int a = 0; a < dim; a++) {
				if (U[a][l] == 0) continue;
				for (int b = 0; b < dim; b++) {
					if (V[b][l] == 0) continue;
					double uv = U[a][l] * V[b][l];
					for (int c = 0; c < dim; c++) T[a][b][c] += uv * W[c][l];
				}
			}
		}
		for (double[][] Ta : T) for (double[] Tab : Ta) for (double t : Tab) nrm2 += t * t;
		System.out.printf("#boundary products B=%d ; ||R_bd||=%.4f ; scanning CP rank ↓%n", B, Math.sqrt(nrm2));

		// Rank continuation: start from the KNOWN exact rank-B decomposition (the boundary
		// product columns), then repeatedly drop the smallest-norm rank-1 and ALS-refine.
		// Far more reliable than cold ALS in this high-rank (B > dim) regime.
		double[][] X = new double[dim][B], Y = new double[dim][B], Z = new double[dim][B];
		int col = 0;
		for (int l = 0; l < r; l++) {
			if (cls[l] != 1) continue;
			for (int a = 0; a < dim; a++) X[a][col] = U[a][l];
			for (int b = 0; b < dim; b++) Y[b][col] = V[b][l];
			for (int c = 0; c < dim; c++) Z[c][col] = W[c][l];
			col++;
		}
		System.out.printf("init rel-residual at R=B (exact decomposition): %.2e%n",
				Math.sqrt(residual2(T, X, Y, Z, dim, B) / nrm2));

		int floor = Math.max(1, (int) Math.floor(floorFrac * B));
		int firstWall = -1;
		double[][] cX = X, cY = Y, cZ = Z;
		for (int R = B - 1; R >= floor; R--) {
			// drop the smallest-norm rank-1 component from the current (R+1)-solution
			double[][][] dropped = dropSmallestTriple(cX, cY, cZ, dim, R + 1);
			double[][][] best = alsRefine(T, dropped[0], dropped[1], dropped[2], dim, R, nrm2, 4000);
			double bestRel = Math.sqrt(residual2(T, best[0], best[1], best[2], dim, R) / nrm2);
			// a couple of extra random restarts in case continuation got stuck
			for (int rs = 0; rs < 2 && bestRel > 1e-8; rs++) {
				double[][][] rr = alsRefine(T, randn(dim, R, new Random(7 * R + rs)),
						randn(dim, R, new Random(31 * R + rs)), randn(dim, R, new Random(53 * R + rs)),
						dim, R, nrm2, 4000);
				double rel = Math.sqrt(residual2(T, rr[0], rr[1], rr[2], dim, R) / nrm2);
				if (rel < bestRel) { bestRel = rel; best = rr; }
			}
			double maxc = Math.max(maxAbs(best[0]), Math.max(maxAbs(best[1]), maxAbs(best[2])));
			String flag = bestRel < 1e-8 ? "  EXACT" : bestRel < 1e-4 ? "  ~approx (border?)" : "  ** WALL **";
			System.out.printf("  rank %3d : rel-residual %.2e   max|coef|=%.1e%s%n", R, bestRel, maxc, flag);
			cX = best[0]; cY = best[1]; cZ = best[2];
			if (bestRel >= 1e-8) { firstWall = R; break; }
		}
		System.out.println("--------");
		if (firstWall < 0) {
			System.out.printf("EXACT at every rank down to %d — R_bd deeply reducible (B=%d).%n", floor, B);
		} else if (firstWall == B - 1) {
			System.out.printf("WALL at rank %d = B−1: R_bd is NOT reducible; the %d boundary products are "
					+ "an essentially minimal decomposition of their sum. The 2·N² tail is tight for THIS "
					+ "aggregation skeleton.%n", B - 1, B);
		} else {
			System.out.printf("EXACT down to rank %d, WALL at %d ⟹ R_bd reduces B=%d → %d  = −%d products, "
					+ "a strictly better scheme (N² coefficient drops).%n",
					firstWall + 1, firstWall, B, firstWall + 1, B - (firstWall + 1));
		}
	}

	/** Drop the rank-1 column with smallest ||x||·||y||·||z|| → (R-1)-column init triple. */
	private static double[][][] dropSmallestTriple(double[][] X, double[][] Y, double[][] Z, int dim, int R) {
		int drop = 0;
		double min = Double.MAX_VALUE;
		for (int rr = 0; rr < R; rr++) {
			double nx = 0, ny = 0, nz = 0;
			for (int i = 0; i < dim; i++) { nx += X[i][rr] * X[i][rr]; ny += Y[i][rr] * Y[i][rr]; nz += Z[i][rr] * Z[i][rr]; }
			double m = Math.sqrt(nx * ny * nz);
			if (m < min) { min = m; drop = rr; }
		}
		double[][] nX = new double[dim][R - 1], nY = new double[dim][R - 1], nZ = new double[dim][R - 1];
		int c = 0;
		for (int rr = 0; rr < R; rr++) {
			if (rr == drop) continue;
			for (int i = 0; i < dim; i++) { nX[i][c] = X[i][rr]; nY[i][c] = Y[i][rr]; nZ[i][c] = Z[i][rr]; }
			c++;
		}
		return new double[][][] { nX, nY, nZ };
	}

	/** ALS from a given init; returns {X,Y,Z}. */
	private static double[][][] alsRefine(double[][][] T, double[][] X0, double[][] Y0, double[][] Z0,
			int dim, int R, double nrm2, int iters) {
		double[][] X = X0, Y = Y0, Z = Z0;
		double prev = Double.MAX_VALUE;
		for (int it = 0; it < iters; it++) {
			X = update(T, X, Y, Z, dim, R, 0);
			Y = update(T, X, Y, Z, dim, R, 1);
			Z = update(T, X, Y, Z, dim, R, 2);
			if ((it & 63) == 0 || it == iters - 1) {
				double rel = Math.sqrt(residual2(T, X, Y, Z, dim, R) / nrm2);
				if (rel < 1e-11) break;
				if (Math.abs(prev - rel) < 1e-13 && it > 200) break;
				prev = rel;
			}
		}
		return new double[][][] { X, Y, Z };
	}

	/** One ALS factor update (mode 0=X,1=Y,2=Z) via MTTKRP · pinv(Hadamard of Gramians). */
	private static double[][] update(double[][][] T, double[][] X, double[][] Y, double[][] Z,
			int dim, int R, int mode) {
		double[][] A = mode == 0 ? Y : X;                 // the two OTHER factors
		double[][] Bm = mode == 2 ? Y : Z;
		// MTTKRP: M[i][r] = Σ_{p,q} T_(mode)[i][p][q] A[p][r] B[q][r]
		double[][] M = new double[dim][R];
		for (int i = 0; i < dim; i++) {
			for (int p = 0; p < dim; p++) {
				for (int q = 0; q < dim; q++) {
					double t = get(T, mode, i, p, q);
					if (t == 0) continue;
					for (int rr = 0; rr < R; rr++) M[i][rr] += t * A[p][rr] * Bm[q][rr];
				}
			}
		}
		double[][] H = hadamard(gram(A, R), gram(Bm, R), R);   // (AᵀA)∘(BᵀB)
		RealMatrix Hpinv = pinv(new Array2DRowRealMatrix(H, false));
		return new Array2DRowRealMatrix(M, false).multiply(Hpinv).getData();
	}

	private static double get(double[][][] T, int mode, int i, int p, int q) {
		// mode 0: T[i][p][q]; mode 1 (Y): factors X(p),Z(q) → T[p][i][q]; mode 2 (Z): X(p),Y(q) → T[p][q][i]
		return mode == 0 ? T[i][p][q] : mode == 1 ? T[p][i][q] : T[p][q][i];
	}

	private static double[][] gram(double[][] F, int R) {
		int n = F.length;
		double[][] g = new double[R][R];
		for (int a = 0; a < R; a++) for (int b = a; b < R; b++) {
			double s = 0;
			for (int i = 0; i < n; i++) s += F[i][a] * F[i][b];
			g[a][b] = s; g[b][a] = s;
		}
		return g;
	}

	private static double[][] hadamard(double[][] A, double[][] B, int R) {
		double[][] h = new double[R][R];
		for (int i = 0; i < R; i++) for (int j = 0; j < R; j++) h[i][j] = A[i][j] * B[i][j];
		return h;
	}

	private static RealMatrix pinv(RealMatrix m) {
		SingularValueDecomposition svd = new SingularValueDecomposition(m);
		double[] s = svd.getSingularValues();
		double tol = Math.max(m.getRowDimension(), m.getColumnDimension()) * s[0] * 1e-12;
		RealMatrix sInv = new Array2DRowRealMatrix(m.getColumnDimension(), m.getRowDimension());
		for (int i = 0; i < s.length; i++) if (s[i] > tol) sInv.setEntry(i, i, 1.0 / s[i]);
		return svd.getV().multiply(sInv).multiply(svd.getUT());
	}

	private static double residual2(double[][][] T, double[][] X, double[][] Y, double[][] Z, int dim, int R) {
		double s = 0;
		for (int a = 0; a < dim; a++) for (int b = 0; b < dim; b++) for (int c = 0; c < dim; c++) {
			double rec = 0;
			for (int rr = 0; rr < R; rr++) rec += X[a][rr] * Y[b][rr] * Z[c][rr];
			double d = rec - T[a][b][c];
			s += d * d;
		}
		return s;
	}

	private static double maxAbs(double[][] m) {
		double x = 0;
		for (double[] row : m) for (double v : row) x = Math.max(x, Math.abs(v));
		return x;
	}

	private static double[][] randn(int rows, int R, Random rng) {
		double[][] m = new double[rows][R];
		for (int i = 0; i < rows; i++) for (int j = 0; j < R; j++) m[i][j] = rng.nextGaussian();
		return m;
	}
}
