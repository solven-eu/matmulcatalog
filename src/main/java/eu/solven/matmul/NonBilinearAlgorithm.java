package eu.solven.matmul;

/**
 * A <strong>non-bilinear</strong> ("quadratic") algorithm for matrix
 * multiplication {@code ⟨n, m, p⟩} over a COMMUTATIVE ring, per
 * DIS09 §1.2: each rank-1 product factor may carry coefficients on
 * BOTH input matrices.
 *
 * <p>Concretely, for {@code k = 0..r-1}:</p>
 * <pre>
 *   α_k = Σ_{i,j} Ua[i·m + j][k] · A[i][j]  +  Σ_{j,l} Ub[j·p + l][k] · B[j][l]
 *   β_k = Σ_{i,j} Va[i·m + j][k] · A[i][j]  +  Σ_{j,l} Vb[j·p + l][k] · B[j][l]
 *   γ_k = α_k · β_k
 *   C[i][l] = Σ_k W[i·p + l][k] · γ_k
 * </pre>
 *
 * <p>The product {@code γ_k} requires <em>commutativity</em>: when we
 * expand {@code (Ua·A + Ub·B)·(Va·A + Vb·B)} we get cross-terms
 * involving {@code A[i][j]·A[i'][j']}, {@code A·B}, {@code B·A}, and
 * {@code B·B} — over a non-commutative ring {@code A·B ≠ B·A}, so the
 * cancellation that makes the matmul output correct only works when
 * scalar entries commute.</p>
 *
 * <p>This format is the one needed to encode e.g.
 * <strong>Rosowski 2019</strong> Algorithm 1 ({@code ⟨n,3,3⟩ = 6n+3}
 * commutative) and Corollary 1 ({@code ⟨3,3,3⟩ = 21} commutative),
 * which are non-bilinear and so don't fit the standard bilinear
 * {@link NonCubicBilinearAlgorithm} format. See
 * {@code references/rosowski-algorithms.md} for the explicit products.</p>
 *
 * <p><strong>Setting bilinear = true:</strong> when {@code Ub} and
 * {@code Va} are both all-zero, this reduces to a standard bilinear
 * algorithm (the {@code α} factor is {@code Ua·A}, the {@code β}
 * factor is {@code Vb·B}). A bilinear algorithm can be losslessly
 * embedded in this format by setting those two matrices to zero.</p>
 */
public class NonBilinearAlgorithm {

	public final int n;
	public final int m;
	public final int p;
	public final int r;

	/** A-coefficients of the first factor (α_k). Shape {@code [n·m][r]}. */
	public final double[][] Ua;
	/** B-coefficients of the first factor (α_k). Shape {@code [m·p][r]}. */
	public final double[][] Ub;
	/** A-coefficients of the second factor (β_k). Shape {@code [n·m][r]}. */
	public final double[][] Va;
	/** B-coefficients of the second factor (β_k). Shape {@code [m·p][r]}. */
	public final double[][] Vb;
	/** Output combination matrix. Shape {@code [n·p][r]}. */
	public final double[][] W;

	public NonBilinearAlgorithm(int n, int m, int p,
			double[][] Ua, double[][] Ub,
			double[][] Va, double[][] Vb,
			double[][] W) {
		int dimA = n * m, dimB = m * p, dimC = n * p;
		check2D(Ua, dimA, "Ua");
		check2D(Ub, dimB, "Ub");
		check2D(Va, dimA, "Va");
		check2D(Vb, dimB, "Vb");
		check2D(W, dimC, "W");
		int rank = Ua[0].length;
		if (Ub[0].length != rank || Va[0].length != rank || Vb[0].length != rank || W[0].length != rank) {
			throw new IllegalArgumentException("inconsistent rank across factor matrices");
		}
		this.n = n;
		this.m = m;
		this.p = p;
		this.r = rank;
		this.Ua = Ua;
		this.Ub = Ub;
		this.Va = Va;
		this.Vb = Vb;
		this.W = W;
	}

	private static void check2D(double[][] M, int expectedRows, String name) {
		if (M.length != expectedRows) {
			throw new IllegalArgumentException(name + " rows must be " + expectedRows + ", got " + M.length);
		}
		int rank = M[0].length;
		for (int i = 1; i < M.length; i++) {
			if (M[i].length != rank) {
				throw new IllegalArgumentException(name + " row " + i + " has rank " + M[i].length + " ≠ " + rank);
			}
		}
	}

	/** Lift a bilinear algorithm into the non-bilinear representation (Ub = Va = 0). */
	public static NonBilinearAlgorithm fromBilinear(NonCubicBilinearAlgorithm alg) {
		double[][] zeroUb = new double[alg.m * alg.p][alg.r];
		double[][] zeroVa = new double[alg.n * alg.m][alg.r];
		return new NonBilinearAlgorithm(alg.n, alg.m, alg.p, alg.denseU(), zeroUb, zeroVa, alg.denseV(), alg.denseW());
	}

	/**
	 * Returns {@code true} when this non-bilinear algorithm is in fact bilinear
	 * — Ub and Va are both all-zero — and can be losslessly downcast to
	 * {@link NonCubicBilinearAlgorithm}.
	 */
	public boolean isPurelyBilinear() {
		for (double[] row : Ub) for (double v : row) if (v != 0) return false;
		for (double[] row : Va) for (double v : row) if (v != 0) return false;
		return true;
	}
}
