package eu.solven.matmul.io;

import java.io.File;
import java.io.IOException;
import java.nio.charset.StandardCharsets;
import java.nio.file.Files;
import java.util.ArrayList;
import java.util.List;

import eu.solven.matmul.NonCubicBilinearAlgorithm;

/**
 * Parser for FMM-Lille's <em>tensor</em> artifacts
 * ({@code https://fmm.univ-lille.fr/{n}x{m}x{p}_tensor.mpl.bz2}, once bunzip2'ed) —
 * the Java counterpart of {@code tools/import_fmm_maple.py}, so an artifact can be
 * pulled without Python. ({@link MapleSchemeParser} reads the other FMM format, the
 * {@code *_raw.mpl} MUL/ADD listing.)
 *
 * <pre>
 *   Tensor:=TriadSet([
 *     Triad([Matrix(n, m, [[…]]), Matrix(m, p, [[…]]), Matrix(p, n, [[…]])]),
 *     …
 *   ]):
 *   map(expand, A.B - add(Trace(Transpose(T[i][1]).A) · Trace(Transpose(T[i][2]).B) · Transpose(T[i][3]), i=1..R));
 * </pre>
 *
 * <p>Per the file's own verification line: product {@code i} multiplies
 * {@code Σ U[a][b]·A[a][b]} by {@code Σ V[b][c]·B[b][c]}, and contributes to
 * {@code C[a][c]} with coefficient {@code W[c][a]} — the third matrix is
 * <b>p×n</b>, i.e. transposed with respect to this code base's {@code n×p}
 * row-major W. Entries are integers or {@code p/q} rationals. A recipe-only
 * placeholder (symbolic matrices and an empty {@code Tensor:=Tensor:}) has no
 * {@code Triad} and is reported as such.</p>
 */
public final class MapleTensorParser {

	private MapleTensorParser() {}

	/**
	 * An artifact this parser cannot turn into a scheme: a placeholder without tensor
	 * data (FMM publishes those for recipe-only ranks), or a PARAMETRIC family whose
	 * entries carry a free symbol (e.g. {@code 1/l} in {@code 7x7x17_tensor.mpl}).
	 */
	public static final class UnsupportedArtifactException extends IOException {
		private static final long serialVersionUID = 1L;

		public UnsupportedArtifactException(String message) {
			super(message);
		}
	}

	public static NonCubicBilinearAlgorithm parse(File f) throws IOException {
		return parse(Files.readString(f.toPath(), StandardCharsets.UTF_8));
	}

	public static NonCubicBilinearAlgorithm parse(String text) throws IOException {
		List<double[][][]> triads = new ArrayList<>();
		int at = 0;
		while (true) {
			int triad = text.indexOf("Triad(", at);
			if (triad < 0) break;
			double[][][] uvw = new double[3][][];
			int cursor = triad;
			for (int k = 0; k < 3; k++) {
				int matrix = text.indexOf("Matrix(", cursor);
				if (matrix < 0) throw new IOException("Triad #" + (triads.size() + 1) + ": matrix " + (k + 1) + "/3 missing");
				int[] end = new int[1];
				uvw[k] = parseMatrix(text, matrix + "Matrix(".length(), end);
				cursor = end[0];
			}
			triads.add(uvw);
			at = cursor;
		}
		if (triads.isEmpty()) {
			throw new UnsupportedArtifactException("no Triad in the artifact — a recipe-only placeholder, not a tensor");
		}
		double[][][] first = triads.get(0);
		int n = first[0].length, m = first[0][0].length, p = first[1][0].length;
		int r = triads.size();
		double[][] u = new double[n * m][r], v = new double[m * p][r], w = new double[n * p][r];
		for (int q = 0; q < r; q++) {
			double[][][] t = triads.get(q);
			requireShape(t[0], n, m, q, "U");
			requireShape(t[1], m, p, q, "V");
			requireShape(t[2], p, n, q, "W (p×n)");
			for (int a = 0; a < n; a++) for (int b = 0; b < m; b++) u[a * m + b][q] = t[0][a][b];
			for (int b = 0; b < m; b++) for (int c = 0; c < p; c++) v[b * p + c][q] = t[1][b][c];
			for (int c = 0; c < p; c++) for (int a = 0; a < n; a++) w[a * p + c][q] = t[2][c][a];
		}
		return new NonCubicBilinearAlgorithm(n, m, p, u, v, w);
	}

	private static void requireShape(double[][] mat, int rows, int cols, int triad, String which) throws IOException {
		if (mat.length != rows || mat[0].length != cols) {
			throw new IOException("Triad #" + (triad + 1) + ": " + which + " is " + mat.length + "×" + mat[0].length
					+ ", expected " + rows + "×" + cols);
		}
	}

	/** Parses {@code rows, cols, [[…],[…]])} starting right after {@code Matrix(}; {@code end[0]} = index after it. */
	private static double[][] parseMatrix(String s, int from, int[] end) throws IOException {
		int comma1 = s.indexOf(',', from);
		int comma2 = s.indexOf(',', comma1 + 1);
		int rows = Integer.parseInt(s.substring(from, comma1).trim());
		int cols = Integer.parseInt(s.substring(comma1 + 1, comma2).trim());
		double[][] out = new double[rows][cols];
		int i = s.indexOf('[', comma2);
		if (i < 0) throw new IOException("Matrix without data near offset " + from);
		i++; // inside the outer [
		for (int row = 0; row < rows; row++) {
			i = s.indexOf('[', i) + 1;
			for (int col = 0; col < cols; col++) {
				int stop = i;
				while (s.charAt(stop) != ',' && s.charAt(stop) != ']') stop++;
				out[row][col] = parseCoefficient(s.substring(i, stop).trim());
				i = stop + 1;
			}
		}
		int close = s.indexOf(')', i);
		end[0] = close + 1;
		return out;
	}

	private static double parseCoefficient(String token) throws IOException {
		try {
			int slash = token.indexOf('/');
			if (slash < 0) return Long.parseLong(token);
			return (double) Long.parseLong(token.substring(0, slash)) / Long.parseLong(token.substring(slash + 1));
		} catch (NumberFormatException e) {
			throw new UnsupportedArtifactException("entry '" + token + "' is not an integer or p/q — a parametric (symbolic) family?");
		}
	}
}
