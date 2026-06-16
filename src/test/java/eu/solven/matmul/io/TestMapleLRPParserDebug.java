package eu.solven.matmul.io;

import java.io.File;

import org.junit.jupiter.api.Disabled;
import org.junit.jupiter.api.Test;
import org.junit.jupiter.api.condition.EnabledIf;

import eu.solven.matmul.NonCubicBilinearAlgorithm;
import eu.solven.matmul.Verifier;

/**
 * Debug-only: probes 4 flatten conventions (row-major vs col-major on each
 * matrix factor) to find which one yields a verifying ⟨17,17,17⟩ scheme.
 */
public class TestMapleLRPParserDebug {

	private static final File LRP = new File("references/fmm-lille/17x17x17/17x17x17_LRP.mpl");

	static boolean lrpAvailable() { return LRP.isFile(); }

	@Test
	@Disabled("Diagnostic-only — re-enable to re-probe flatten conventions for a new LRP file.")
	@EnabledIf("lrpAvailable")
	public void probe_index_conventions() throws Exception {
		MapleLRPParser.LRPMatrices mats = MapleLRPParser.parseMatrices(LRP);
		int rank = mats.rank();
		int n = 17, m = 17, p = 17;
		int dimU = n * m, dimV = m * p, dimW = n * p;

		// 8 convention combinations: each of {U, V, W} can be row-major (i*M+j)
		// or col-major (i + j*N). Test all.
		for (int conv = 0; conv < 8; conv++) {
			boolean uCol = (conv & 1) != 0;
			boolean vCol = (conv & 2) != 0;
			boolean wCol = (conv & 4) != 0;

			double[][] U = new double[dimU][rank];
			double[][] V = new double[dimV][rank];
			double[][] W = new double[dimW][rank];

			for (int k = 0; k < rank; k++) {
				double[] row = mats.L()[k];
				for (int idx = 0; idx < dimU; idx++) {
					int i, j;
					if (uCol) { j = idx / n; i = idx % n; }
					else      { i = idx / m; j = idx % m; }
					U[i * m + j][k] = row[idx];
				}
				double[] rrow = mats.R()[k];
				for (int idx = 0; idx < dimV; idx++) {
					int j, l;
					if (vCol) { l = idx / m; j = idx % m; }
					else      { j = idx / p; l = idx % p; }
					V[j * p + l][k] = rrow[idx];
				}
			}
			for (int idx = 0; idx < dimW; idx++) {
				double[] prow = mats.P()[idx];
				int i, l;
				if (wCol) { l = idx / n; i = idx % n; }
				else      { i = idx / p; l = idx % p; }
				for (int k = 0; k < rank; k++) W[i * p + l][k] = prow[k];
			}

			try {
				NonCubicBilinearAlgorithm alg = new NonCubicBilinearAlgorithm(n, m, p, U, V, W);
				boolean ok = Verifier.passesRandomMatmulSpotCheck(alg, 2, 1e-6);
				double resid = computeApproxResidual(alg);
				System.out.printf("conv uCol=%s vCol=%s wCol=%s : numeric_ok=%s residual~%.3f%n",
						uCol, vCol, wCol, ok, resid);
			} catch (Exception ex) {
				System.out.printf("conv uCol=%s vCol=%s wCol=%s : EXCEPTION %s%n",
						uCol, vCol, wCol, ex.getMessage());
			}
		}
	}

	/** Sum of |C_algo - C_naive| for ONE random sample (cheap probe). */
	private static double computeApproxResidual(NonCubicBilinearAlgorithm alg) {
		double[][] srcU = alg.denseU();
		double[][] srcV = alg.denseV();
		double[][] srcW = alg.denseW();
		java.util.Random rng = new java.util.Random(42);
		double[][] A = new double[alg.n][alg.m];
		double[][] B = new double[alg.m][alg.p];
		for (int i = 0; i < alg.n; i++) for (int j = 0; j < alg.m; j++) A[i][j] = rng.nextGaussian();
		for (int j = 0; j < alg.m; j++) for (int l = 0; l < alg.p; l++) B[j][l] = rng.nextGaussian();
		double[] alpha = new double[alg.r];
		double[] beta = new double[alg.r];
		for (int k = 0; k < alg.r; k++) {
			double a = 0, b = 0;
			for (int i = 0; i < alg.n; i++) for (int j = 0; j < alg.m; j++)
				a += srcU[i * alg.m + j][k] * A[i][j];
			for (int j = 0; j < alg.m; j++) for (int l = 0; l < alg.p; l++)
				b += srcV[j * alg.p + l][k] * B[j][l];
			alpha[k] = a;
			beta[k] = b;
		}
		double[][] cAlgo = new double[alg.n][alg.p];
		for (int k = 0; k < alg.r; k++) {
			double m = alpha[k] * beta[k];
			for (int i = 0; i < alg.n; i++) for (int l = 0; l < alg.p; l++)
				cAlgo[i][l] += srcW[i * alg.p + l][k] * m;
		}
		double tot = 0;
		for (int i = 0; i < alg.n; i++) {
			for (int l = 0; l < alg.p; l++) {
				double s = 0;
				for (int j = 0; j < alg.m; j++) s += A[i][j] * B[j][l];
				tot += Math.abs(cAlgo[i][l] - s);
			}
		}
		return tot;
	}
}
