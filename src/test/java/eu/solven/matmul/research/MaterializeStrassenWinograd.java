package eu.solven.matmul.research;

import java.io.File;

import eu.solven.matmul.NonCubicBilinearAlgorithm;
import eu.solven.matmul.Verifier;
import eu.solven.matmul.catalog.Lineage;
import eu.solven.matmul.catalog.SchemeIO;

/**
 * Materialise Winograd's 7-multiplication ⟨2,2,2⟩ variant (1971),
 * sometimes called "Strassen-Winograd". Same rank as Strassen 1969 (=7)
 * but a different bilinear form, designed so the computation has
 * <strong>15 scheduled additions</strong> when intermediates
 * {@code s1, s2, s3, s4, t1, t2, t3, t4, u1, u2, u3} are shared
 * (vs. Strassen 1969's 18). The matrix form U/V/W in the JSON
 * captures only the per-product coefficient vectors, so
 * {@link Verifier#additionCount} reports the naive non-CSE total
 * (24 here) — the 15 figure is a property of the schedule, not the
 * bilinear matrices.
 *
 * <p>Products (canonical Winograd schedule):</p>
 * <pre>
 *   s1 = a21 + a22;  s2 = s1 - a11;  s3 = a11 - a21;  s4 = a12 - s2;
 *   t1 = b12 - b11;  t2 = b22 - t1;  t3 = b22 - b12;  t4 = t2 - b21;
 *
 *   m1 = s2 * t2     = (a21+a22-a11) * (b11-b12+b22)
 *   m2 = a11 * b11
 *   m3 = a12 * b21
 *   m4 = s1 * t1     = (a21+a22) * (b12-b11)
 *   m5 = s3 * t3     = (a11-a21) * (b22-b12)
 *   m6 = s4 * b22    = (a11+a12-a21-a22) * b22
 *   m7 = a22 * t4    = a22 * (b11-b12-b21+b22)
 *
 *   c11 = m2 + m3
 *   c12 = m1 + m2 + m4 + m6
 *   c21 = m1 + m2 + m5 - m7
 *   c22 = m1 + m2 + m4 + m5
 * </pre>
 *
 * <p>Reference: S. Winograd (1971), "On multiplication of 2×2
 * matrices", Linear Algebra and its Applications 4 (1971), 381-388.
 * The schedule analysis is from Probert 1976 (Heun 1994 surveys).</p>
 */
public final class MaterializeStrassenWinograd {

	public static void main(String[] args) throws Exception {
		// U[i*2 + j][k] — a_{i,j} coefficient in product k. Row order: a11, a12, a21, a22.
		double[][] U = {
				// k=  1     2     3     4     5     6     7
				/* a11 */ {-1,  1,  0,  0,  1,  1,  0},
				/* a12 */ { 0,  0,  1,  0,  0,  1,  0},
				/* a21 */ { 1,  0,  0,  1, -1, -1,  0},
				/* a22 */ { 1,  0,  0,  1,  0, -1,  1},
		};
		// V[j*2 + l][k] — b_{j,l} coefficient. Row order: b11, b12, b21, b22.
		double[][] V = {
				/* b11 */ { 1,  1,  0, -1,  0,  0,  1},
				/* b12 */ {-1,  0,  0,  1, -1,  0, -1},
				/* b21 */ { 0,  0,  1,  0,  0,  0, -1},
				/* b22 */ { 1,  0,  0,  0,  1,  1,  1},
		};
		// W[i*2 + l][k] — c_{i,l} = sum over k of W[i*2+l][k] * m_k. Row order: c11, c12, c21, c22.
		double[][] W = {
				/* c11 */ { 0,  1,  1,  0,  0,  0,  0},
				/* c12 */ { 1,  1,  0,  1,  0,  1,  0},
				/* c21 */ { 1,  1,  0,  0,  1,  0, -1},
				/* c22 */ { 1,  1,  0,  1,  1,  0,  0},
		};
		NonCubicBilinearAlgorithm alg = new NonCubicBilinearAlgorithm(2, 2, 2, U, V, W);

		if (!Verifier.passesRandomMatmulSpotCheck(alg)) {
			System.err.println("Winograd 7-mult spot-check FAILED — derivation bug");
			return;
		}
		int naiveAdds = Verifier.additionCount(alg);
		System.out.printf("Winograd ⟨2,2,2⟩=7 verified.%n");
		System.out.printf("  matrix-form naive additions: %d%n", naiveAdds);
		System.out.printf("  scheduled additions (Winograd CSE): 15%n");

		Lineage.Node lineage = new Lineage.Atom("winograd-1971-7mult");
		File out = new File(String.format(
				"src/main/resources/schemes/known/section2/winograd-1971_2x2x2_m7_a%d.json", naiveAdds));
		out.getParentFile().mkdirs();
		SchemeIO.write(alg, out, lineage);
		System.out.println("wrote " + out);
	}
}
