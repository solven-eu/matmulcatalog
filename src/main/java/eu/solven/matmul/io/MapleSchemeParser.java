package eu.solven.matmul.io;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileReader;
import java.io.IOException;
import java.io.Reader;
import java.util.regex.Matcher;
import java.util.regex.Pattern;

import eu.solven.matmul.NonCubicBilinearAlgorithm;

/**
 * Parser for FMM-Lille's <em>raw</em> Maple matmul-scheme files
 * (https://fmm.univ-lille.fr/{n}x{m}x{p}.html, {@code *_raw.mpl}).
 *
 * <p>The file declares three symbolic matrices {@code A} (n×m),
 * {@code B} (m×p), {@code C} (n×p), then a {@code MUL:=[...]:} block
 * holding {@code rank} equations of the form
 * <pre>
 *   m_K = (linear_combo_of_A_ij) * (linear_combo_of_B_jl),
 * </pre>
 * and finally an {@code ADD:=[...]:} block holding {@code n·p} equations
 * <pre>
 *   C_i_l = linear_combo_of_m_K,
 * </pre>
 * one per output cell. The trailing {@code 0=0} pads the block.
 *
 * <p>Coefficients seen in the {@code 17x17x17_raw.mpl} are
 * {@code {0, ±1, ±1/8}} — pure {@link eu.solven.matmul.algebra.Algebra Q-rational}.
 * The parser handles {@code N/D} fractions and integer literals; it produces
 * {@code double[][]} matrices because that's the in-memory representation of
 * {@link NonCubicBilinearAlgorithm} (the JSON catalog encodes fractions
 * exactly via string formatting, so no precision loss on the way out).</p>
 *
 * <p>Index convention: the file uses 1-based indices; the parser converts
 * to 0-based and uses row-major flatten {@code i·dim + j} matching
 * {@link NonCubicBilinearAlgorithm}.</p>
 */
public final class MapleSchemeParser {

	private MapleSchemeParser() {}

	/** Parse a raw Maple FMM-Lille scheme file. */
	public static NonCubicBilinearAlgorithm parseRawFmmLille(File f, int n, int m, int p) throws IOException {
		try (Reader r = new BufferedReader(new FileReader(f))) {
			return parseRawFmmLille(r, n, m, p);
		}
	}

	public static NonCubicBilinearAlgorithm parseRawFmmLille(Reader r, int n, int m, int p) throws IOException {
		StringBuilder full = new StringBuilder();
		try (BufferedReader br = (r instanceof BufferedReader b) ? b : new BufferedReader(r)) {
			String line;
			while ((line = br.readLine()) != null) {
				full.append(line).append('\n');
			}
		}
		String src = full.toString();

		// 1. Locate the MUL and ADD blocks.
		int mulStart = src.indexOf("MUL:=[");
		if (mulStart < 0) throw new IOException("missing MUL:=[ block");
		int mulEnd = src.indexOf("]:", mulStart);
		if (mulEnd < 0) throw new IOException("missing ]: terminator for MUL block");
		String mulBlock = src.substring(mulStart + "MUL:=[".length(), mulEnd);

		int addStart = src.indexOf("ADD:=[", mulEnd);
		if (addStart < 0) throw new IOException("missing ADD:=[ block");
		int addEnd = src.indexOf("]:", addStart);
		if (addEnd < 0) throw new IOException("missing ]: terminator for ADD block");
		String addBlock = src.substring(addStart + "ADD:=[".length(), addEnd);

		// 2. Split the MUL block into per-product entries.
		java.util.List<String> mulEntries = splitTopLevel(mulBlock);
		int rank = 0;
		for (String e : mulEntries) {
			String stripped = e.strip();
			if (stripped.isEmpty()) continue;
			if (stripped.startsWith("m_")) rank++;
		}
		if (rank == 0) throw new IOException("no m_K product entries found");

		int dimU = n * m, dimV = m * p, dimW = n * p;
		double[][] U = new double[dimU][rank];
		double[][] V = new double[dimV][rank];
		double[][] W = new double[dimW][rank];

		// 3. Parse each m_K = (A_combo) * (B_combo).
		int seen = 0;
		for (String e : mulEntries) {
			String stripped = e.strip();
			if (stripped.isEmpty() || !stripped.startsWith("m_")) continue;
			int eq = stripped.indexOf('=');
			if (eq < 0) throw new IOException("malformed MUL entry (no =): " + brief(stripped));
			int k = parseMIndex(stripped.substring(0, eq).strip()) - 1; // 1-based → 0-based
			if (k < 0 || k >= rank) throw new IOException("m_K index out of range: " + (k + 1));
			String rhs = stripped.substring(eq + 1).strip();
			String[] factors = splitProductFactors(rhs);
			if (factors.length != 2) throw new IOException("expected (U)*(V) for m_" + (k + 1) + ", got "
					+ factors.length + " factors: " + brief(rhs));
			fillFactorAB(factors[0], k, U, V, n, m, p, /*isA=*/true);
			fillFactorAB(factors[1], k, U, V, n, m, p, /*isA=*/false);
			seen++;
		}
		if (seen != rank) throw new IOException("expected " + rank + " products, parsed " + seen);

		// 4. Parse each C_i_l = sum of ± [coef·]m_K.
		java.util.List<String> addEntries = splitTopLevel(addBlock);
		int outputCount = 0;
		for (String e : addEntries) {
			String stripped = e.strip();
			if (stripped.isEmpty()) continue;
			if (stripped.startsWith("0")) continue; // padding line "0=0"
			if (!stripped.startsWith("C_")) continue;
			int eq = stripped.indexOf('=');
			if (eq < 0) throw new IOException("malformed ADD entry: " + brief(stripped));
			String head = stripped.substring(0, eq).strip();
			int[] il = parseCIndex(head);
			int i = il[0] - 1, l = il[1] - 1; // 1-based → 0-based
			if (i < 0 || i >= n || l < 0 || l >= p) {
				throw new IOException("C_i_l indices out of range: " + head);
			}
			String rhs = stripped.substring(eq + 1).strip();
			int row = i * p + l;
			for (Term t : parseLinearCombo(rhs)) {
				int kk = parseMIndex(t.var) - 1;
				if (kk < 0 || kk >= rank) {
					throw new IOException("W references unknown m_" + (kk + 1));
				}
				W[row][kk] += t.coef;
			}
			outputCount++;
		}
		if (outputCount != dimW) {
			throw new IOException("expected " + dimW + " C-output assignments, parsed " + outputCount);
		}

		return new NonCubicBilinearAlgorithm(n, m, p, U, V, W);
	}

	/* ------------------------------------------------------------------ */
	/* Splitting & parsing helpers                                         */
	/* ------------------------------------------------------------------ */

	/**
	 * Split a comma-separated Maple list at top level (depth 0 parentheses).
	 * Newlines are treated as whitespace.
	 */
	private static java.util.List<String> splitTopLevel(String src) {
		java.util.List<String> out = new java.util.ArrayList<>();
		int depth = 0;
		int last = 0;
		for (int i = 0; i < src.length(); i++) {
			char c = src.charAt(i);
			if (c == '(' || c == '[') depth++;
			else if (c == ')' || c == ']') depth--;
			else if (c == ',' && depth == 0) {
				out.add(src.substring(last, i));
				last = i + 1;
			}
		}
		out.add(src.substring(last));
		return out;
	}

	/**
	 * Split {@code (combo_A)*(combo_B)} into exactly two factors. The RHS of
	 * every MUL entry must be a product of two parenthesised linear combos;
	 * the boundary is the {@code )*(} between them.
	 */
	private static String[] splitProductFactors(String rhs) {
		// Find the two outermost balanced parens.
		int firstOpen = rhs.indexOf('(');
		if (firstOpen < 0) throw new IllegalArgumentException("no opening paren: " + brief(rhs));
		int depth = 0;
		int firstClose = -1;
		for (int i = firstOpen; i < rhs.length(); i++) {
			char c = rhs.charAt(i);
			if (c == '(') depth++;
			else if (c == ')') {
				depth--;
				if (depth == 0) { firstClose = i; break; }
			}
		}
		if (firstClose < 0) throw new IllegalArgumentException("unbalanced parens: " + brief(rhs));
		// Skip over the '*' between the two factors.
		int starIdx = firstClose + 1;
		while (starIdx < rhs.length() && Character.isWhitespace(rhs.charAt(starIdx))) starIdx++;
		if (starIdx >= rhs.length() || rhs.charAt(starIdx) != '*') {
			throw new IllegalArgumentException("expected '*' after first factor: " + brief(rhs));
		}
		int secondOpen = rhs.indexOf('(', starIdx);
		if (secondOpen < 0) throw new IllegalArgumentException("no second factor: " + brief(rhs));
		depth = 0;
		int secondClose = -1;
		for (int i = secondOpen; i < rhs.length(); i++) {
			char c = rhs.charAt(i);
			if (c == '(') depth++;
			else if (c == ')') {
				depth--;
				if (depth == 0) { secondClose = i; break; }
			}
		}
		if (secondClose < 0) throw new IllegalArgumentException("unbalanced second-factor parens");
		String left = rhs.substring(firstOpen + 1, firstClose);
		String right = rhs.substring(secondOpen + 1, secondClose);
		return new String[] { left, right };
	}

	private static int parseMIndex(String s) {
		if (!s.startsWith("m_")) throw new IllegalArgumentException("expected m_K, got: " + s);
		return Integer.parseInt(s.substring(2).trim());
	}

	private static int[] parseCIndex(String s) {
		Matcher m = C_RE.matcher(s);
		if (!m.matches()) throw new IllegalArgumentException("expected C_i_l, got: " + s);
		return new int[] { Integer.parseInt(m.group(1)), Integer.parseInt(m.group(2)) };
	}

	private static final Pattern C_RE = Pattern.compile("C_(\\d+)_(\\d+)");

	/**
	 * A single signed monomial: {@code ± [num/den ·] var}, where {@code var}
	 * is one of {@code A_i_j}, {@code B_j_l}, {@code m_K}.
	 */
	private record Term(double coef, String var) {}

	private static final Pattern TERM_RE = Pattern.compile(
			"\\s*([+\\-])?\\s*" +                       // optional sign
			"(?:(\\d+(?:/\\d+)?)\\s*\\*\\s*)?" +         // optional rational coef like "1/8*"
			"([ABm]_\\d+(?:_\\d+)?)" +                   // var
			"\\s*");

	/**
	 * Parse a linear combination like {@code -A_1_2 + 1/8*A_15_5 - B_3_4 + ...}
	 * (or {@code -1/8*m_42 + m_17}) into a list of {@link Term}s.
	 */
	private static java.util.List<Term> parseLinearCombo(String s) {
		java.util.List<Term> out = new java.util.ArrayList<>();
		Matcher m = TERM_RE.matcher(s);
		int pos = 0;
		while (pos < s.length()) {
			// skip whitespace and trailing commas/semicolons
			while (pos < s.length() && Character.isWhitespace(s.charAt(pos))) pos++;
			if (pos >= s.length()) break;
			m.region(pos, s.length());
			if (!m.lookingAt()) {
				throw new IllegalArgumentException("can't parse term at offset " + pos + ": "
						+ brief(s.substring(pos)));
			}
			String sign = m.group(1);
			String coefStr = m.group(2);
			String var = m.group(3);
			double coef = 1.0;
			if (coefStr != null) {
				int slash = coefStr.indexOf('/');
				if (slash >= 0) {
					double num = Double.parseDouble(coefStr.substring(0, slash));
					double den = Double.parseDouble(coefStr.substring(slash + 1));
					coef = num / den;
				} else {
					coef = Double.parseDouble(coefStr);
				}
			}
			if ("-".equals(sign)) coef = -coef;
			out.add(new Term(coef, var));
			pos = m.end();
		}
		return out;
	}

	/**
	 * Fill U or V from one of the two factor linear combinations of a MUL
	 * entry. {@code isA = true} means the factor references {@code A_i_j}
	 * (fills U); {@code false} means {@code B_j_l} (fills V).
	 */
	private static void fillFactorAB(String factor, int k, double[][] U, double[][] V,
			int n, int m, int p, boolean isA) {
		for (Term t : parseLinearCombo(factor)) {
			String[] parts = t.var.split("_");
			if (parts.length != 3) throw new IllegalArgumentException("malformed var: " + t.var);
			int row = Integer.parseInt(parts[1]) - 1;
			int col = Integer.parseInt(parts[2]) - 1;
			if (isA) {
				if (!"A".equals(parts[0])) {
					throw new IllegalArgumentException("expected A_*_* in U factor, got: " + t.var);
				}
				if (row < 0 || row >= n || col < 0 || col >= m) {
					throw new IllegalArgumentException("A index out of range: " + t.var);
				}
				U[row * m + col][k] += t.coef;
			} else {
				if (!"B".equals(parts[0])) {
					throw new IllegalArgumentException("expected B_*_* in V factor, got: " + t.var);
				}
				if (row < 0 || row >= m || col < 0 || col >= p) {
					throw new IllegalArgumentException("B index out of range: " + t.var);
				}
				V[row * p + col][k] += t.coef;
			}
		}
	}

	private static String brief(String s) {
		String t = s.replace('\n', ' ');
		return t.length() > 120 ? t.substring(0, 120) + "..." : t;
	}
}
