package eu.solven.matmul;

/**
 * Complex-valued counterpart to {@link NonCubicBilinearAlgorithm}: factors
 * live in {@code C}, stored as parallel real/imaginary {@code double[][]}
 * arrays. Needed for schemes like AlphaEvolve's {@code ⟨4,4,4⟩ = 48} where
 * the coefficient set is {@code 0.5·Z[i]} (half-integer Gaussian rationals).
 *
 * <p>For {@code k = 0..r-1}, the trilinear identity holds with complex
 * arithmetic:
 * {@code Σ_k U_k · V_k · W_k = T_{n,m,p}} (the real-valued matmul tensor —
 * the imaginary part of the sum must vanish).</p>
 */
public class ComplexNonCubicBilinearAlgorithm {

	public final int n;
	public final int m;
	public final int p;
	public final int r;
	public final double[][] uRe;  // [n·m][r]
	public final double[][] uIm;
	public final double[][] vRe;  // [m·p][r]
	public final double[][] vIm;
	public final double[][] wRe;  // [n·p][r]
	public final double[][] wIm;

	public ComplexNonCubicBilinearAlgorithm(int n, int m, int p,
			double[][] uRe, double[][] uIm,
			double[][] vRe, double[][] vIm,
			double[][] wRe, double[][] wIm) {
		int dimU = n * m, dimV = m * p, dimW = n * p;
		checkShape(uRe, dimU, "uRe");
		checkShape(uIm, dimU, "uIm");
		checkShape(vRe, dimV, "vRe");
		checkShape(vIm, dimV, "vIm");
		checkShape(wRe, dimW, "wRe");
		checkShape(wIm, dimW, "wIm");
		int rank = uRe[0].length;
		if (uIm[0].length != rank || vRe[0].length != rank || vIm[0].length != rank
				|| wRe[0].length != rank || wIm[0].length != rank) {
			throw new IllegalArgumentException("re/im factor matrices must share rank");
		}
		this.n = n;
		this.m = m;
		this.p = p;
		this.r = rank;
		this.uRe = uRe;
		this.uIm = uIm;
		this.vRe = vRe;
		this.vIm = vIm;
		this.wRe = wRe;
		this.wIm = wIm;
	}

	private static void checkShape(double[][] arr, int expectedRows, String label) {
		if (arr.length != expectedRows) {
			throw new IllegalArgumentException(
					label + " rows = " + arr.length + ", expected " + expectedRows);
		}
	}

	public boolean isCubic() {
		return n == m && m == p;
	}

	public int dimU() { return n * m; }
	public int dimV() { return m * p; }
	public int dimW() { return n * p; }
}
