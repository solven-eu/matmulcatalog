package eu.solven.matmul.catalog;

import java.util.ArrayList;
import java.util.List;

import eu.solven.matmul.NonCubicBilinearAlgorithm;

/**
 * Human-friendly textual rendering of a bilinear scheme (#188), mirroring the
 * Perminov FastMatrixMultiplication {@code multiplications} / {@code elements}
 * representation (e.g. {@code schemes/results/ZT/2x4x11_m70_ZT.json}):
 *
 * <ul>
 *   <li><b>multiplications</b> — one line per product
 *       {@code m1 = (a11 + a22)*(b11 + b22)}, built from the {@code U}/{@code V}
 *       factor columns.</li>
 *   <li><b>elements</b> — one line per output cell
 *       {@code c11 = m1 + m4 - m5 + m7}, built from the {@code W} factor.</li>
 * </ul>
 *
 * <p>Index convention follows {@link NonCubicBilinearAlgorithm}: {@code U} row
 * {@code i·m+j} is {@code a_{i,j}} (A is n×m), {@code V} row {@code j·p+l} is
 * {@code b_{j,l}} (B is m×p), {@code W} row {@code i·p+l} is {@code c_{i,l}}
 * (C is n×p), all 1-indexed in the output. Subscripts are concatenated for
 * dims ≤ 9 ({@code a23}) and comma-separated above ({@code a2,11}).</p>
 *
 * <p>Precomputed at manifest/JSON-generation time so the SPA only ever renders
 * the strings (user principle: the SPA should barely ever compute).</p>
 */
public final class HumanScheme {

	private HumanScheme() {}

	public record Readable(List<String> multiplications, List<String> elements) {}

	public static Readable of(NonCubicBilinearAlgorithm alg) {
		int n = alg.n, m = alg.m, p = alg.p, r = alg.r;
		List<String> mults = new ArrayList<>(r);
		for (int k = 0; k < r; k++) {
			String a = linearForm(alg.denseU(), k, n * m, idx -> cell("a", idx / m, idx % m, n, m));
			String b = linearForm(alg.denseV(), k, m * p, idx -> cell("b", idx / p, idx % p, m, p));
			mults.add("m" + (k + 1) + " = (" + a + ")*(" + b + ")");
		}
		List<String> elems = new ArrayList<>(n * p);
		for (int row = 0; row < n * p; row++) {
			// Element row index = i·p + l → c_{i,l}.
			String lhs = cell("c", row / p, row % p, n, p);
			String rhs = productCombo(alg.denseW(), row, r);
			elems.add(lhs + " = " + rhs);
		}
		return new Readable(mults, elems);
	}

	/** Linear form over A/B entries for product column {@code k}. */
	private static String linearForm(double[][] factor, int k, int dim, java.util.function.IntFunction<String> name) {
		StringBuilder sb = new StringBuilder();
		for (int idx = 0; idx < dim; idx++) {
			double c = factor[idx][k];
			if (c == 0) continue;
			appendTerm(sb, c, name.apply(idx));
		}
		return sb.length() == 0 ? "0" : sb.toString();
	}

	/** Combination of products m_k for one output cell (row of W). */
	private static String productCombo(double[][] w, int row, int r) {
		StringBuilder sb = new StringBuilder();
		for (int k = 0; k < r; k++) {
			double c = w[row][k];
			if (c == 0) continue;
			appendTerm(sb, c, "m" + (k + 1));
		}
		return sb.length() == 0 ? "0" : sb.toString();
	}

	/** Append "± [coeff*]token" with tidy signs (drops a leading "+ "). */
	private static void appendTerm(StringBuilder sb, double c, String token) {
		boolean neg = c < 0;
		double abs = Math.abs(c);
		String coeff = (abs == 1.0) ? "" : (trim(abs) + "*");
		if (sb.length() == 0) {
			sb.append(neg ? "-" : "").append(coeff).append(token);
		} else {
			sb.append(neg ? " - " : " + ").append(coeff).append(token);
		}
	}

	/** "a" + 1-indexed (row,col); concatenated for dims ≤ 9, comma'd above. */
	private static String cell(String base, int row0, int col0, int rows, int cols) {
		int i = row0 + 1, j = col0 + 1;
		if (rows <= 9 && cols <= 9) return base + i + "" + j;
		return base + i + "," + j;
	}

	private static String trim(double v) {
		if (v == Math.rint(v)) return Long.toString((long) v);
		return Double.toString(v);
	}
}
