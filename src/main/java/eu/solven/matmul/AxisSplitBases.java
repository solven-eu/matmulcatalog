package eu.solven.matmul;

/**
 * Trivial "axis-split" base algorithms — the building blocks DIS09 calls
 * mul211, mul121, mul112. Each is the rank-2 multiplication of a thin
 * vector by a vector/scalar; useful as a Kronecker base for compositions
 * that split a single axis.
 *
 * <ul>
 *   <li>{@link #mul211} = {@code ⟨2,1,1⟩=2} — splits the first axis</li>
 *   <li>{@link #mul121} = {@code ⟨1,2,1⟩=2} — splits the middle axis (dot product)</li>
 *   <li>{@link #mul112} = {@code ⟨1,1,2⟩=2} — splits the third axis</li>
 * </ul>
 *
 * <p>These are not "improvements" in any sense — naive matmul of an
 * {@code ⟨n,m,p⟩} where one axis is degenerate (=1 or =2) is the standard
 * dot product. DIS09 includes them in their pattern pool because using
 * one as a recursive base produces the standard axis-split composition
 * recipe, expressed uniformly with the rest of the framework.</p>
 */
public final class AxisSplitBases {

	private AxisSplitBases() {}

	/** {@code ⟨2,1,1⟩=2}. A is 2×1, B is 1×1, C is 2×1. C[i,0] = A[i,0]·B[0,0]. */
	public static NonCubicBilinearAlgorithm mul211() {
		double[][] U = new double[2][2];
		double[][] V = new double[1][2];
		double[][] W = new double[2][2];
		// mult k=0: c[0,0] = a[0,0]·b[0,0]
		U[0][0] = 1; V[0][0] = 1; W[0][0] = 1;
		// mult k=1: c[1,0] = a[1,0]·b[0,0]
		U[1][1] = 1; V[0][1] = 1; W[1][1] = 1;
		return new NonCubicBilinearAlgorithm(2, 1, 1, U, V, W);
	}

	/** {@code ⟨1,2,1⟩=2}. A is 1×2, B is 2×1, C is 1×1. Standard dot product. */
	public static NonCubicBilinearAlgorithm mul121() {
		double[][] U = new double[2][2];
		double[][] V = new double[2][2];
		double[][] W = new double[1][2];
		// mult k=0: a[0,0]·b[0,0]
		U[0][0] = 1; V[0][0] = 1; W[0][0] = 1;
		// mult k=1: a[0,1]·b[1,0]
		U[1][1] = 1; V[1][1] = 1; W[0][1] = 1;
		return new NonCubicBilinearAlgorithm(1, 2, 1, U, V, W);
	}

	/** {@code ⟨1,1,2⟩=2}. A is 1×1, B is 1×2, C is 1×2. C[0,j] = A[0,0]·B[0,j]. */
	public static NonCubicBilinearAlgorithm mul112() {
		double[][] U = new double[1][2];
		double[][] V = new double[2][2];
		double[][] W = new double[2][2];
		// mult k=0: c[0,0] = a[0,0]·b[0,0]
		U[0][0] = 1; V[0][0] = 1; W[0][0] = 1;
		// mult k=1: c[0,1] = a[0,0]·b[0,1]
		U[0][1] = 1; V[1][1] = 1; W[1][1] = 1;
		return new NonCubicBilinearAlgorithm(1, 1, 2, U, V, W);
	}
}
