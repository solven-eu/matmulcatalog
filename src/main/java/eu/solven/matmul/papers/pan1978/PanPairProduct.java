package eu.solven.matmul.papers.pan1978;

/**
 * Pan's pair-product trilinear-aggregation algorithm for **two**
 * cyclically-related matrix multiplications computed jointly:
 * {@code C = A·B (⟨a,b,c⟩)} and {@code C' = U·V (⟨b,c,a⟩)}.
 *
 * <p>Rank: {@code abc + ab + bc + ca} — strictly less than
 * {@code 2·R_naive(⟨a,b,c⟩) = 2abc} for {@code a,b,c ≥ 3}.</p>
 *
 * <p>Per Pan 2014 (arXiv:1412.1145) §5–6, derived from Pan 1972's
 * trilinear-aggregation generating table:</p>
 *
 * <pre>
 *   Aggregates (abc):   Q(i,j,h) = (A[i,j] + U[j,h]) · (B[j,h] + V[h,i])
 *   v-corrections (ac): P1(h,i) = V[h,i] · (Σ_j A[i,j] + Σ_j U[j,h])
 *   u·b shared (bc):    P2(j,h) = U[j,h] · B[j,h]
 *   a-corrections (ab): P3(i,j) = A[i,j] · (Σ_h B[j,h] + Σ_h V[h,i])
 *
 *   C[i,h]  = Σ_j Q(i,j,h) − P1(h,i) − Σ_j P2(j,h)
 *   C'[j,i] = Σ_h Q(i,j,h) − P3(i,j) − Σ_h P2(j,h)
 * </pre>
 *
 * <p>Symbolic verification (8 outputs for n=2, all PASS) lives in
 * the project's sympy script {@code tools/verify_pan_pair_222.py}.</p>
 *
 * <p><b>Encoding note</b>: a pair scheme is NOT a standard
 * {@link eu.solven.matmul.NonCubicBilinearAlgorithm} — it has two
 * distinct output blocks of incompatible shapes
 * ({@code C} is {@code a×c}, {@code C'} is {@code b×a}). This class
 * exposes the raw product list and the W-style output assembly
 * matrices instead, so callers wanting to fuse a pair of Strassen
 * sub-products can use the construction directly without forcing
 * the result into the single-output bilinear format.</p>
 *
 * <p>Integration into {@link Recombination} is left as a follow-up:
 * the recombination framework would need to detect cyclically-paired
 * sub-products and route them through this constructor instead of
 * looking up two independent {@code ⟨a,b,c⟩} schemes.</p>
 */
public final class PanPairProduct {

	private PanPairProduct() {}

	/** Closed-form rank: {@code abc + ab + bc + ca}. */
	public static int rank(int a, int b, int c) {
		if (a < 1 || b < 1 || c < 1) throw new IllegalArgumentException("dims must be ≥ 1");
		return a * b * c + a * b + b * c + c * a;
	}

	/**
	 * A pair-product scheme: two matmul outputs ({@code C = A·B} of
	 * shape {@code a×c} and {@code C' = U·V} of shape {@code b×a})
	 * produced jointly by a single set of {@code rank} bilinear products.
	 *
	 * <p>Each product {@code k ∈ [0, rank)} is a pair
	 * {@code (α_k, β_k)} where:</p>
	 * <ul>
	 *   <li>{@code α_k} is a linear form in entries of {@code A} and
	 *       {@code U} (combined "first-side" inputs)</li>
	 *   <li>{@code β_k} is a linear form in entries of {@code B} and
	 *       {@code V} (combined "second-side" inputs)</li>
	 * </ul>
	 *
	 * <p>The output assembly is given by {@code W_C[a*c][rank]} for
	 * {@code C} and {@code W_Cp[b*a][rank]} for {@code C'}. Each
	 * output entry is a {@code Σ_k W·γ_k} where {@code γ_k = α_k·β_k}.</p>
	 */
	public record PairScheme(
			int a, int b, int c,
			double[][] alphaA,  // [a*b][rank] coef of A[i,j] in α_k
			double[][] alphaU,  // [b*c][rank] coef of U[j,h] in α_k
			double[][] betaB,   // [b*c][rank] coef of B[j,h] in β_k
			double[][] betaV,   // [c*a][rank] coef of V[h,i] in β_k
			double[][] W_C,     // [a*c][rank] assembly for C[i,h]
			double[][] W_Cp     // [b*a][rank] assembly for C'[j,i]
	) {
		public int rank() { return alphaA[0].length; }
	}

	/**
	 * Build the pair scheme for the joint computation of
	 * {@code C = A·B (⟨a,b,c⟩)} and {@code C' = U·V (⟨b,c,a⟩)}.
	 */
	public static PairScheme build(final int a, final int b, final int c) {
		if (a < 1 || b < 1 || c < 1) {
			throw new IllegalArgumentException("dims must be ≥ 1, got (" + a + "," + b + "," + c + ")");
		}
		final int rTotal = a * b * c + a * b + b * c + c * a;

		double[][] alphaA = new double[a * b][rTotal];
		double[][] alphaU = new double[b * c][rTotal];
		double[][] betaB  = new double[b * c][rTotal];
		double[][] betaV  = new double[c * a][rTotal];
		double[][] W_C    = new double[a * c][rTotal];
		double[][] W_Cp   = new double[b * a][rTotal];

		// Product slot layout:
		//   [0,        a*b*c)        — Q(i,j,h) aggregates  (index = i·b·c + j·c + h)
		//   [a*b*c,    a*b*c + c*a)  — P1(h, i)             (index = base1 + h·a + i)
		//   [base2,    base2 + b*c)  — P2(j, h)             (index = base2 + j·c + h)
		//   [base3,    base3 + a*b)  — P3(i, j)             (index = base3 + i·b + j)
		final int OFF_Q  = 0;
		final int OFF_P1 = a * b * c;
		final int OFF_P2 = OFF_P1 + c * a;
		final int OFF_P3 = OFF_P2 + b * c;

		// ────────────────────── Q(i, j, h) ──────────────────────
		// α_k = A[i,j] + U[j,h]   β_k = B[j,h] + V[h,i]
		for (int i = 0; i < a; i++)
			for (int j = 0; j < b; j++)
				for (int h = 0; h < c; h++) {
					int k = OFF_Q + i * b * c + j * c + h;
					alphaA[i * b + j][k]  = 1;
					alphaU[j * c + h][k]  = 1;
					betaB[j * c + h][k]   = 1;
					betaV[h * a + i][k]   = 1;
				}

		// ────────────────────── P1(h, i): V[h,i] · (Σ_j A[i,j] + Σ_j U[j,h]) ──────────────────────
		for (int h = 0; h < c; h++)
			for (int i = 0; i < a; i++) {
				int k = OFF_P1 + h * a + i;
				// α_k = Σ_j A[i,j] + Σ_j U[j,h]
				for (int j = 0; j < b; j++) {
					alphaA[i * b + j][k] = 1;
					alphaU[j * c + h][k] = 1;
				}
				// β_k = V[h,i]
				betaV[h * a + i][k] = 1;
			}

		// ────────────────────── P2(j, h): U[j,h] · B[j,h] ──────────────────────
		for (int j = 0; j < b; j++)
			for (int h = 0; h < c; h++) {
				int k = OFF_P2 + j * c + h;
				alphaU[j * c + h][k] = 1;
				betaB[j * c + h][k]  = 1;
			}

		// ────────────────────── P3(i, j): A[i,j] · (Σ_h B[j,h] + Σ_h V[h,i]) ──────────────────────
		for (int i = 0; i < a; i++)
			for (int j = 0; j < b; j++) {
				int k = OFF_P3 + i * b + j;
				alphaA[i * b + j][k] = 1;
				for (int h = 0; h < c; h++) {
					betaB[j * c + h][k] = 1;
					betaV[h * a + i][k] = 1;
				}
			}

		// ────────────────────── W_C[i, h] = Σ_j Q(i,j,h) − P1(h,i) − Σ_j P2(j,h) ──────────────────────
		for (int i = 0; i < a; i++)
			for (int h = 0; h < c; h++) {
				int outRow = i * c + h;
				for (int j = 0; j < b; j++) {
					int kQ = OFF_Q + i * b * c + j * c + h;
					W_C[outRow][kQ] = +1;
				}
				int kP1 = OFF_P1 + h * a + i;
				W_C[outRow][kP1] = -1;
				for (int j = 0; j < b; j++) {
					int kP2 = OFF_P2 + j * c + h;
					W_C[outRow][kP2] = -1;
				}
			}

		// ────────────────────── W_C'[j, i] = Σ_h Q(i,j,h) − P3(i,j) − Σ_h P2(j,h) ──────────────────────
		for (int j = 0; j < b; j++)
			for (int i = 0; i < a; i++) {
				int outRow = j * a + i;
				for (int h = 0; h < c; h++) {
					int kQ = OFF_Q + i * b * c + j * c + h;
					W_Cp[outRow][kQ] = +1;
				}
				int kP3 = OFF_P3 + i * b + j;
				W_Cp[outRow][kP3] = -1;
				for (int h = 0; h < c; h++) {
					int kP2 = OFF_P2 + j * c + h;
					W_Cp[outRow][kP2] = -1;
				}
			}

		return new PairScheme(a, b, c, alphaA, alphaU, betaB, betaV, W_C, W_Cp);
	}

	/**
	 * Random spot-check: evaluate the pair scheme on a uniformly-sampled
	 * {@code (A, B, U, V)} input quadruple and compare each output
	 * against naive matmul. Returns {@code true} iff both outputs match
	 * within {@code eps}.
	 */
	public static boolean spotCheck(PairScheme s, int samples, double eps, long seed) {
		java.util.Random rng = new java.util.Random(seed);
		int a = s.a(), b = s.b(), c = s.c();
		double[][] alphaA = s.alphaA(), alphaU = s.alphaU();
		double[][] betaB = s.betaB(), betaV = s.betaV();
		double[][] W_C = s.W_C(), W_Cp = s.W_Cp();
		int r = s.rank();
		for (int t = 0; t < samples; t++) {
			double[] A = new double[a * b], B = new double[b * c];
			double[] U = new double[b * c], V = new double[c * a];
			for (int i = 0; i < A.length; i++) A[i] = rng.nextGaussian();
			for (int i = 0; i < B.length; i++) B[i] = rng.nextGaussian();
			for (int i = 0; i < U.length; i++) U[i] = rng.nextGaussian();
			for (int i = 0; i < V.length; i++) V[i] = rng.nextGaussian();

			// Compute γ_k = α_k · β_k
			double[] gamma = new double[r];
			for (int k = 0; k < r; k++) {
				double alpha = 0, beta = 0;
				for (int idx = 0; idx < a * b; idx++) alpha += alphaA[idx][k] * A[idx];
				for (int idx = 0; idx < b * c; idx++) alpha += alphaU[idx][k] * U[idx];
				for (int idx = 0; idx < b * c; idx++) beta  += betaB[idx][k] * B[idx];
				for (int idx = 0; idx < c * a; idx++) beta  += betaV[idx][k] * V[idx];
				gamma[k] = alpha * beta;
			}

			// Naive AB
			double[] cExpect = new double[a * c];
			for (int i = 0; i < a; i++)
				for (int h = 0; h < c; h++) {
					double sum = 0;
					for (int j = 0; j < b; j++) sum += A[i * b + j] * B[j * c + h];
					cExpect[i * c + h] = sum;
				}
			// Naive UV
			double[] cpExpect = new double[b * a];
			for (int j = 0; j < b; j++)
				for (int i = 0; i < a; i++) {
					double sum = 0;
					for (int h = 0; h < c; h++) sum += U[j * c + h] * V[h * a + i];
					cpExpect[j * a + i] = sum;
				}

			// Algorithm C
			for (int idx = 0; idx < a * c; idx++) {
				double sum = 0;
				for (int k = 0; k < r; k++) sum += W_C[idx][k] * gamma[k];
				if (Math.abs(sum - cExpect[idx]) > eps) return false;
			}
			// Algorithm C'
			for (int idx = 0; idx < b * a; idx++) {
				double sum = 0;
				for (int k = 0; k < r; k++) sum += W_Cp[idx][k] * gamma[k];
				if (Math.abs(sum - cpExpect[idx]) > eps) return false;
			}
		}
		return true;
	}
}
