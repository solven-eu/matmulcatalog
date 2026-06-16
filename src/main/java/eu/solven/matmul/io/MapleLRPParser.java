package eu.solven.matmul.io;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileReader;
import java.io.IOException;
import java.io.PushbackReader;
import java.io.Reader;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.List;
import java.util.Map;

import eu.solven.matmul.NonCubicBilinearAlgorithm;

/**
 * Parser for FMM-Lille's <em>LRP</em> Maple matmul-scheme files
 * ({@code *_LRP.mpl}, downloadable from
 * https://fmm.univ-lille.fr/{n}x{m}x{p}.html).
 *
 * <p>The file is a single Maple expression of the form:
 * <pre>
 *   LRP := [Matrix(r, n·m, [[row], [row], ...]),
 *           Matrix(r, m·p, [[row], [row], ...]),
 *           Matrix(n·p, r, [[row], [row], ...])]:
 * </pre>
 * encoding the bilinear scheme as three matrices {@code (L, R, P)}.
 * Convention (per the file's published shape; confirmed by parsing
 * {@code 17x17x17_LRP.mpl} — three matrices, dims {@code (r, n·m)},
 * {@code (r, m·p)}, {@code (n·p, r)} with {@code r=2931}, {@code n·m=289}):</p>
 * <ul>
 *   <li>{@code L} is {@code r × n·m} (rows = products, cols = A-entries)</li>
 *   <li>{@code R} is {@code r × m·p} (rows = products, cols = B-entries)</li>
 *   <li>{@code P} is {@code n·p × r} (rows = C-entries, cols = products) —
 *       already in our {@link NonCubicBilinearAlgorithm#W W} layout</li>
 * </ul>
 *
 * <p>So the mapping to our {@link NonCubicBilinearAlgorithm} is:
 * {@code U[?][k] = L[k][?]} (transpose), {@code V[?][k] = R[k][?]} (transpose),
 * {@code W[?][k] = P[?][k]} (identity).</p>
 *
 * <p>Coefficients seen in {@code 17x17x17_LRP.mpl} are
 * {@code {0, ±1, ±1/8}} — pure Q-rational. The parser handles
 * {@code N/D} fractions and integer literals; final storage is
 * {@code double[][]} (the JSON catalog re-emits fractions exactly via
 * {@link eu.solven.matmul.catalog.SchemeIO}).</p>
 *
 * <p>The file is one ~5 MB single-token line — the parser streams the
 * source through a {@link PushbackReader} rather than slurping it into
 * memory and regex-walking it (much smaller peak heap).</p>
 */
public final class MapleLRPParser {

	private MapleLRPParser() {}

	/** Raw {@code (L, R, P)} matrices as published in the LRP file, with their dims. */
	public record LRPMatrices(
			int rank,
			int dimA,
			int dimB,
			int dimC,
			double[][] L, // [rank][dimA]
			double[][] R, // [rank][dimB]
			double[][] P  // [dimC][rank]
	) {}

	/** Parse a Maple LRP file into the {@code (L, R, P)} matrices, no shape check. */
	public static LRPMatrices parseMatrices(File f) throws IOException {
		try (Reader r = new BufferedReader(new FileReader(f))) {
			return parseMatrices(r);
		}
	}

	public static LRPMatrices parseMatrices(Reader r) throws IOException {
		PushbackReader in = new PushbackReader(r, 16);
		// Skip header preamble (HTML-ish warning line, which also contains the
		// substring "LRP"). Look for the assignment token "LRP:=" specifically.
		seek(in, "LRP:=");
		skipWhitespace(in);
		expectChar(in, '[');

		double[][] L = readMatrix(in);
		skipWhitespace(in);
		expectChar(in, ',');
		double[][] R = readMatrix(in);
		skipWhitespace(in);
		expectChar(in, ',');
		double[][] P = readMatrix(in);
		skipWhitespace(in);
		expectChar(in, ']');

		int rank = L.length;
		if (R.length != rank) {
			throw new IOException("LRP: R has " + R.length + " rows but L has " + rank);
		}
		int dimA = L[0].length;
		int dimB = R[0].length;
		int dimC = P.length;
		if (P[0].length != rank) {
			throw new IOException("LRP: P has " + P[0].length + " cols but rank is " + rank);
		}
		return new LRPMatrices(rank, dimA, dimB, dimC, L, R, P);
	}

	/** Parse the LRP file and convert directly to a {@link NonCubicBilinearAlgorithm}. */
	public static NonCubicBilinearAlgorithm parse(File f, int n, int m, int p) throws IOException {
		LRPMatrices mats = parseMatrices(f);
		return toAlgorithm(mats, n, m, p);
	}

	/**
	 * Map the raw {@code (L, R, P)} into our {@code (U, V, W)}. Empirically
	 * verified index conventions for FMM-Lille's {@code _LRP.mpl} (confirmed
	 * against {@code 17x17x17_LRP.mpl} by an 8-way flatten-orientation probe):
	 * <ul>
	 *   <li>{@code L[k][j·n + i] = U[i·m + j][k]} — L is column-major over A
	 *       (i.e. {@code A[i][j]} is at column {@code j*n + i})</li>
	 *   <li>{@code R[k][l·m + j] = V[j·p + l][k]} — R is column-major over B</li>
	 *   <li>{@code P[i·p + l][k] = W[i·p + l][k]} — P is row-major over C (same
	 *       layout as our W)</li>
	 * </ul>
	 *
	 * <p>This is also what one would expect from Maple's default column-major
	 * convention for {@code Matrix(rows, cols, …)} when the source flattening
	 * is "Maple-style column-major"; FMM-Lille's choice of row-major-C for P
	 * is the asymmetry. The probe is replayed in {@code TestMapleLRPParserDebug}
	 * for any future LRP file with potentially different conventions.</p>
	 */
	public static NonCubicBilinearAlgorithm toAlgorithm(LRPMatrices mats, int n, int m, int p) throws IOException {
		int dimU = n * m, dimV = m * p, dimW = n * p;
		if (mats.dimA != dimU) {
			throw new IOException("LRP: L cols = " + mats.dimA + " but n·m = " + dimU);
		}
		if (mats.dimB != dimV) {
			throw new IOException("LRP: R cols = " + mats.dimB + " but m·p = " + dimV);
		}
		if (mats.dimC != dimW) {
			throw new IOException("LRP: P rows = " + mats.dimC + " but n·p = " + dimW);
		}
		int rank = mats.rank;
		double[][] U = new double[dimU][rank];
		double[][] V = new double[dimV][rank];
		double[][] W = new double[dimW][rank];
		// L: column-major over A → U[i·m+j][k] = L[k][j·n+i]
		for (int k = 0; k < rank; k++) {
			double[] row = mats.L[k];
			for (int idx = 0; idx < dimU; idx++) {
				int j = idx / n;
				int i = idx % n;
				U[i * m + j][k] = row[idx];
			}
		}
		// R: column-major over B → V[j·p+l][k] = R[k][l·m+j]
		for (int k = 0; k < rank; k++) {
			double[] row = mats.R[k];
			for (int idx = 0; idx < dimV; idx++) {
				int l = idx / m;
				int j = idx % m;
				V[j * p + l][k] = row[idx];
			}
		}
		// P: row-major over C (same as W) → W[i·p+l][k] = P[i·p+l][k]
		for (int c = 0; c < dimW; c++) {
			double[] row = mats.P[c];
			for (int k = 0; k < rank; k++) W[c][k] = row[k];
		}
		return new NonCubicBilinearAlgorithm(n, m, p, U, V, W);
	}

	// ──────────────────────────────────────────────────────────────────────
	// Streaming Maple-Matrix reader
	// ──────────────────────────────────────────────────────────────────────

	/** Read {@code Matrix(rows, cols, [[...], [...], ...])} into a dense {@code double[rows][cols]}. */
	private static double[][] readMatrix(PushbackReader in) throws IOException {
		skipWhitespace(in);
		expectKeyword(in, "Matrix");
		skipWhitespace(in);
		expectChar(in, '(');
		int rows = readInt(in);
		skipWhitespace(in);
		expectChar(in, ',');
		int cols = readInt(in);
		skipWhitespace(in);
		expectChar(in, ',');
		skipWhitespace(in);
		expectChar(in, '[');

		double[][] out = new double[rows][cols];
		for (int i = 0; i < rows; i++) {
			skipWhitespace(in);
			expectChar(in, '[');
			for (int j = 0; j < cols; j++) {
				skipWhitespace(in);
				out[i][j] = readRational(in);
				skipWhitespace(in);
				if (j < cols - 1) expectChar(in, ',');
			}
			expectChar(in, ']');
			skipWhitespace(in);
			if (i < rows - 1) expectChar(in, ',');
		}
		skipWhitespace(in);
		expectChar(in, ']');
		skipWhitespace(in);
		expectChar(in, ')');
		return out;
	}

	private static int readInt(PushbackReader in) throws IOException {
		skipWhitespace(in);
		StringBuilder sb = new StringBuilder();
		int c = in.read();
		if (c == '+' || c == '-') {
			sb.append((char) c);
			c = in.read();
		}
		while (c >= '0' && c <= '9') {
			sb.append((char) c);
			c = in.read();
		}
		if (c != -1) in.unread(c);
		if (sb.length() == 0 || (sb.length() == 1 && (sb.charAt(0) == '+' || sb.charAt(0) == '-'))) {
			throw new IOException("expected integer, got: " + sb);
		}
		return Integer.parseInt(sb.toString());
	}

	/**
	 * Read a Maple rational: signed integer, possibly followed by
	 * {@code /<integer>}. Returns the value as a {@code double}.
	 */
	private static double readRational(PushbackReader in) throws IOException {
		StringBuilder sb = new StringBuilder();
		int c = in.read();
		if (c == '+' || c == '-') {
			sb.append((char) c);
			c = in.read();
		}
		while (c >= '0' && c <= '9') {
			sb.append((char) c);
			c = in.read();
		}
		if (sb.length() == 0 || (sb.length() == 1 && (sb.charAt(0) == '+' || sb.charAt(0) == '-'))) {
			throw new IOException("expected rational numerator at char: " + (c == -1 ? "EOF" : (char) c));
		}
		double num = Double.parseDouble(sb.toString());
		if (c == '/') {
			StringBuilder den = new StringBuilder();
			c = in.read();
			while (c >= '0' && c <= '9') {
				den.append((char) c);
				c = in.read();
			}
			if (c != -1) in.unread(c);
			if (den.length() == 0) throw new IOException("expected rational denominator");
			return num / Double.parseDouble(den.toString());
		}
		if (c != -1) in.unread(c);
		return num;
	}

	private static void skipWhitespace(PushbackReader in) throws IOException {
		int c;
		while ((c = in.read()) != -1) {
			if (!Character.isWhitespace(c)) {
				in.unread(c);
				return;
			}
		}
	}

	private static void expectChar(PushbackReader in, char ch) throws IOException {
		int c = in.read();
		if (c != ch) {
			throw new IOException("expected '" + ch + "', got '"
					+ (c == -1 ? "EOF" : (char) c) + "'");
		}
	}

	/** Read the exact characters of {@code kw}, throwing if they don't match. */
	private static void expectKeyword(PushbackReader in, String kw) throws IOException {
		for (int i = 0; i < kw.length(); i++) {
			int c = in.read();
			if (c != kw.charAt(i)) {
				throw new IOException("expected keyword '" + kw + "', mismatch at index " + i
						+ " (got '" + (c == -1 ? "EOF" : (char) c) + "')");
			}
		}
	}

	/**
	 * Scan forward, character by character, until {@code keyword} has been
	 * fully matched at the current read position. Used to skip the leading
	 * HTML-ish warning line before {@code LRP:=}.
	 */
	private static void seek(PushbackReader in, String keyword) throws IOException {
		StringBuilder buf = new StringBuilder();
		int c;
		while ((c = in.read()) != -1) {
			buf.append((char) c);
			if (buf.length() > keyword.length()) buf.deleteCharAt(0);
			if (buf.toString().equals(keyword)) return;
		}
		throw new IOException("never found '" + keyword + "'");
	}

	/* ------------------------------------------------------------------ */
	/* Diagnostic helpers — exposed for investigation tooling.            */
	/* ------------------------------------------------------------------ */

	/**
	 * Count how many pairs {@code (i < j)} of product columns have identical
	 * {@code U[:,i] == U[:,j]} AND identical {@code V[:,i] == V[:,j]} —
	 * i.e. duplicated A·B combinations. Each such pair is a candidate for
	 * "shared multiplication" merging.
	 */
	public static int countDuplicateUVPairs(NonCubicBilinearAlgorithm alg) {
		double[][] srcU = alg.denseU();
		double[][] srcV = alg.denseV();
		int r = alg.r;
		int dimU = alg.dimU(), dimV = alg.dimV();
		List<String> sigs = new ArrayList<>(r);
		for (int k = 0; k < r; k++) {
			StringBuilder sb = new StringBuilder();
			for (int a = 0; a < dimU; a++) {
				double v = srcU[a][k];
				if (v != 0.0) sb.append(a).append(':').append(v).append(';');
			}
			sb.append('|');
			for (int b = 0; b < dimV; b++) {
				double v = srcV[b][k];
				if (v != 0.0) sb.append(b).append(':').append(v).append(';');
			}
			sigs.add(sb.toString());
		}
		Map<String, List<Integer>> grouped = new HashMap<>();
		for (int k = 0; k < r; k++) {
			grouped.computeIfAbsent(sigs.get(k), s -> new ArrayList<>()).add(k);
		}
		int pairs = 0;
		for (List<Integer> grp : grouped.values()) {
			int s = grp.size();
			if (s >= 2) pairs += s * (s - 1) / 2;
		}
		return pairs;
	}

	/**
	 * For each product k, return {@code (rowExtent, colExtent)} on the A side:
	 * how many distinct rows / cols of A the column {@code U[:,k]} touches.
	 * Used to cross-check the "polynomial decomposition" hypothesis
	 * (e.g. {@code 95×⟨4,4,4⟩ + 336×⟨3,3,3⟩ + ...}).
	 */
	public static int[][] rowExtentsA(NonCubicBilinearAlgorithm alg) {
		double[][] srcU = alg.denseU();
		int r = alg.r;
		int[][] out = new int[r][2];
		for (int k = 0; k < r; k++) {
			boolean[] rowSeen = new boolean[alg.n];
			boolean[] colSeen = new boolean[alg.m];
			for (int i = 0; i < alg.n; i++) {
				for (int j = 0; j < alg.m; j++) {
					if (srcU[i * alg.m + j][k] != 0.0) {
						rowSeen[i] = true;
						colSeen[j] = true;
					}
				}
			}
			int nr = 0, nc = 0;
			for (boolean b : rowSeen) if (b) nr++;
			for (boolean b : colSeen) if (b) nc++;
			out[k][0] = nr;
			out[k][1] = nc;
		}
		return out;
	}

	/** Same as {@link #rowExtentsA} but for the B (V) side. */
	public static int[][] rowExtentsB(NonCubicBilinearAlgorithm alg) {
		double[][] srcV = alg.denseV();
		int r = alg.r;
		int[][] out = new int[r][2];
		for (int k = 0; k < r; k++) {
			boolean[] rowSeen = new boolean[alg.m];
			boolean[] colSeen = new boolean[alg.p];
			for (int j = 0; j < alg.m; j++) {
				for (int l = 0; l < alg.p; l++) {
					if (srcV[j * alg.p + l][k] != 0.0) {
						rowSeen[j] = true;
						colSeen[l] = true;
					}
				}
			}
			int nr = 0, nc = 0;
			for (boolean b : rowSeen) if (b) nr++;
			for (boolean b : colSeen) if (b) nc++;
			out[k][0] = nr;
			out[k][1] = nc;
		}
		return out;
	}
}
