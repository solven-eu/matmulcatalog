package eu.solven.matmul.catalog;

import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.File;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.IOException;
import java.io.PrintWriter;
import java.io.Reader;
import java.util.List;
import java.util.Locale;
import java.util.Map;
import java.util.Optional;
import java.util.Map;

import eu.solven.matmul.BilinearAlgorithm;
import eu.solven.matmul.ComplexNonCubicBilinearAlgorithm;
import eu.solven.matmul.NonBilinearAlgorithm;
import eu.solven.matmul.NonCubicBilinearAlgorithm;
import tools.jackson.databind.JsonNode;
import tools.jackson.databind.json.JsonMapper;
import tools.jackson.databind.node.ArrayNode;
import tools.jackson.databind.node.JsonNodeFactory;
import tools.jackson.databind.node.ObjectNode;

/**
 * Read / write matmul algorithms in the
 * <a href="https://github.com/dronperminov/FastMatrixMultiplication#example">dronperminov
 * scheme format</a>.
 *
 * <p>The format is a JSON object with these fields:</p>
 * <ul>
 *   <li>{@code "n": [n_1, n_2, n_3]} — dimensions of {@code ⟨n_1, n_2, n_3⟩}.</li>
 *   <li>{@code "m": int} — rank (number of multiplications).</li>
 *   <li>{@code "fields": string[]} — the algebras the scheme is valid over
 *       (e.g. {@code ["F2","F3","Z","Q","R","C"]}); {@code "fields_not": string[]}
 *       lists explicitly-excluded algebras. These supersede the retired
 *       {@code "z2"} boolean — field membership is read from {@code fields[]},
 *       never from a per-field flag.</li>
 *   <li>{@code "u": int[m][n_1·n_2]} — coefficients of A-entries in each multiplication (row-major: position {@code i·n_2 + j}).</li>
 *   <li>{@code "v": int[m][n_2·n_3]} — coefficients of B-entries (row-major: {@code j·n_3 + l}).</li>
 *   <li>{@code "w": int[m][n_3·n_1]} — coefficients of C-entries in <b>column-major</b>
 *       flatten ({@code position = j·n_1 + i}), per dronperminov's spec. Our internal
 *       {@link BilinearAlgorithm} uses row-major C-flatten ({@code i·n_3 + j}), so
 *       {@link #read}/{@link #write} transpose W between the two conventions.</li>
 *   <li>{@code "multiplications"}, {@code "elements"} — human-readable strings (ignored on read).</li>
 * </ul>
 *
 * <p>Note: dronperminov's format stores {@code U/V/W} as {@code rank × dim}
 * (row-major per multiplication index), whereas our {@link BilinearAlgorithm}
 * uses {@code dim × rank} (one row per flatten position). {@link #read} and
 * {@link #write} transpose accordingly.</p>
 */
public final class SchemeIO {

	private static final JsonMapper MAPPER = JsonMapper.builder().build();

	private SchemeIO() {}

	/**
	 * Parse a JSON scheme file into a {@link JsonNode}. Use this once per
	 * file when you need to inspect the scheme kind ({@link #isNonBilinear},
	 * {@link #isComplex}, …) and then dispatch to the right reader — all
	 * the kind-checks and {@code read*} methods have JsonNode-accepting
	 * overloads, so loading once is sufficient instead of once per check.
	 */
	public static JsonNode parseJson(File f) throws IOException {
		try (Reader r = new BufferedReader(new FileReader(resolveSchemeFile(f)))) {
			return MAPPER.readTree(r);
		}
	}

	/** Parse a scheme from a JSON file (standard sparse/dense bilinear). */
	public static NonCubicBilinearAlgorithm read(File f) throws IOException {
		try (Reader r = new BufferedReader(new FileReader(resolveSchemeFile(f)))) {
			return read(r);
		}
	}

	/**
	 * Resolve a scheme-file path tolerantly. If {@code f} exists, return it.
	 * Otherwise look for a sibling whose stem extends {@code f}'s by a
	 * metric/field token — e.g. a stale {@code …_m7_a18.json} reference resolves
	 * to {@code …_m7_a18_b0.json} after the {@code _b{bud_score}} migration
	 * (2026-06-06). Returns {@code f} unchanged when nothing matches, so a genuine
	 * typo still surfaces as a {@code FileNotFoundException}. The {@code "_"}
	 * boundary stops {@code …_m7} from matching {@code …_m70}.
	 */
	static File resolveSchemeFile(File f) {
		if (f == null || f.exists()) return f;
		File dir = f.getParentFile();
		if (dir == null || !dir.isDirectory()) return f;
		String name = f.getName();
		String stem = name.endsWith(".json") ? name.substring(0, name.length() - 5) : name;
		File[] sibs = dir.listFiles((d, n) -> n.endsWith(".json") && n.startsWith(stem + "_"));
		if (sibs != null && sibs.length > 0) {
			java.util.Arrays.sort(sibs);
			return sibs[0];
		}
		return f;
	}

	/**
	 * Load a bilinear scheme file with one parse, auto-dispatching between
	 * standard and "reduced" sparse-list format. Use this in place of the
	 * idiomatic 2-parse pattern {@code isReduced(f) ? readReduced(f) : read(f)}
	 * — that pattern parses the JSON twice, this one parses it once.
	 */
	public static NonCubicBilinearAlgorithm readBilinear(File f) throws IOException {
		JsonNode root = parseJson(f);
		return isReduced(root) ? readReduced(root) : read(root);
	}

	public static NonCubicBilinearAlgorithm read(String json) throws IOException {
		return read(MAPPER.readTree(json));
	}

	public static NonCubicBilinearAlgorithm read(Reader r) throws IOException {
		return read(MAPPER.readTree(r));
	}

	/** Parse a pre-loaded JsonNode into a bilinear scheme. Auto-dispatches across the three
	 *  on-disk encodings: Perminov "reduced" sparse-list ({@code u}/{@code v}/{@code w} of
	 *  {@code {index,value}} objects, possibly with {@code *_fresh} shared intermediates),
	 *  {@code u_sparse} sparse, and dense. Previously {@code read} silently omitted the reduced
	 *  branch (only {@code readBilinear} dispatched it), so the 168 {@code *_reduced} files —
	 *  including the best NC integer bases like perminov ⟨2,4,4⟩=26 / ⟨4,4,4⟩=49 — threw
	 *  "row 0 has 1 cols, expected …" when loaded via the common {@code read} path. */
	public static NonCubicBilinearAlgorithm read(JsonNode root) throws IOException {
		if (isReduced(root)) return readReduced(root);
		return root.has("u_sparse") ? fromJsonSparse(root) : fromJson(root);
	}

	/**
	 * The unified {@code fields[]} tag list (task #174), or an empty list when
	 * the JSON predates the migration (a few complex / non-bilinear files still
	 * carry the legacy {@code field}/{@code complex} tags).
	 */
	public static java.util.List<String> fieldTags(JsonNode root) {
		JsonNode arr = root.get("fields");
		if (arr == null || !arr.isArray()) return java.util.List.of();
		java.util.List<String> out = new java.util.ArrayList<>(arr.size());
		for (JsonNode e : arr) {
			if (e.isTextual()) out.add(e.asString());
		}
		return out;
	}

	/**
	 * Cheap, EXACT necessary-condition gate: returns a (possibly empty) list of
	 * human-readable reasons for every {@code fields[]} tag that the scheme's
	 * U/V/W coefficient <em>denominators</em> make impossible — independent of
	 * any matmul-identity verification.
	 *
	 * <ul>
	 *   <li>{@code Z}  — every coefficient must be an integer (denominator 1);</li>
	 *   <li>{@code F2} — every reduced denominator must be odd (2 must be
	 *       invertible mod 2): a {@code 1/8} coefficient is NOT representable
	 *       in F₂;</li>
	 *   <li>{@code F3} — every reduced denominator must be coprime to 3.</li>
	 * </ul>
	 *
	 * <p>These are <b>necessary</b> conditions (a violated denominator is
	 * definitely unrepresentable), unlike {@link eu.solven.matmul.algebra.FieldCompliance}'s
	 * {@code {0,±1}}/{@code {0,±1,±2}} membership sets, which are too strict to be
	 * necessary for F₂/F₃ (the integer {@code 2} reduces cleanly mod 2 yet is
	 * outside {@code {0,±1}}). {@code Q}/{@code R}/{@code C} impose no denominator
	 * constraint and are never reported; complex schemes (no real denominator to
	 * analyse) and schemes carrying no coefficient arrays are skipped.</p>
	 *
	 * <p>This is the exact over-claim that has bitten this catalog repeatedly:
	 * a {@code 1/8}-coefficient derived scheme mechanically stamped
	 * {@code fields=[F2,F3,Z,Q,R,C]} from an "integer base ⇒ all fields" rule.
	 * A non-empty result is a field-discipline violation the catalog must never
	 * persist.</p>
	 */
	public static java.util.List<String> fieldsContradictedByCoefficients(JsonNode root) {
		java.util.List<String> tags = fieldTags(root);
		if (tags.isEmpty() || isComplex(root)) {
			return java.util.List.of();
		}
		java.util.List<JsonNode> tokens = collectCoefficientTokens(root);
		if (tokens == null || tokens.isEmpty()) {
			return java.util.List.of();  // stub / non-bilinear / reduced — no raw coefficients
		}
		boolean anyNonInteger = false, anyEvenDen = false, anyDen3 = false;
		for (JsonNode t : tokens) {
			long d = reducedDenominator(t);
			if (d != 1L) anyNonInteger = true;
			if (d % 2L == 0L) anyEvenDen = true;
			if (d % 3L == 0L) anyDen3 = true;
		}
		return contradictions(tags, anyNonInteger, anyEvenDen, anyDen3);
	}

	/**
	 * Same EXACT necessary-condition gate as {@link #fieldsContradictedByCoefficients(JsonNode)},
	 * but reading the coefficient denominators off an in-memory {@link NonCubicBilinearAlgorithm}'s
	 * dense U/V/W — for use at <em>write time</em>, before any JSON exists.
	 *
	 * <p>This is the production-side guard: the materialiser infers a derived
	 * scheme's {@code fields[]} from its lineage leaves (coefficient-BLIND — it
	 * assumes a composition is valid over a field iff every leaf is), which
	 * over-claims when a recombination / ½-polarization step introduces division
	 * (a {@code 1/8} coefficient from integer leaves). Intersecting the inferred
	 * set against this result drops exactly the contradicted tags, so the
	 * materialiser never emits the over-claim in the first place.</p>
	 *
	 * @param tags the candidate field tags (e.g. the lineage-inferred set)
	 */
	public static java.util.List<String> fieldsContradictedByCoefficients(
			NonCubicBilinearAlgorithm alg, java.util.List<String> tags) {
		if (alg == null || tags == null || tags.isEmpty()) {
			return java.util.List.of();
		}
		boolean anyNonInteger = false, anyEvenDen = false, anyDen3 = false;
		for (double[][] mat : new double[][][] { alg.denseU(), alg.denseV(), alg.denseW() }) {
			for (double[] row : mat) {
				for (double v : row) {
					long d = reducedDenominatorOfDouble(v);
					if (d != 1L) anyNonInteger = true;
					if (d % 2L == 0L) anyEvenDen = true;
					if (d % 3L == 0L) anyDen3 = true;
				}
			}
		}
		return contradictions(tags, anyNonInteger, anyEvenDen, anyDen3);
	}

	/**
	 * The candidate {@code tags} with every coefficient-contradicted field removed
	 * (drops {@code Z} if any coefficient is non-integer, {@code F2} if any
	 * denominator is even, {@code F3} if any is divisible by 3). Order-preserving;
	 * returns the same list instance when nothing is contradicted.
	 *
	 * <p>The write-time companion to {@link #fieldsContradictedByCoefficients(NonCubicBilinearAlgorithm, java.util.List)}:
	 * narrows a (coefficient-blind) lineage-inferred set down to what the
	 * materialized coefficients actually support, so the materialiser emits the
	 * correct {@code fields[]} rather than an over-claim a later gate would reject.</p>
	 */
	public static java.util.List<String> narrowFieldsToCoefficients(
			NonCubicBilinearAlgorithm alg, java.util.List<String> tags) {
		if (alg == null || tags == null || tags.isEmpty()) {
			return tags;
		}
		boolean anyNonInteger = false, anyEvenDen = false, anyDen3 = false;
		for (double[][] mat : new double[][][] { alg.denseU(), alg.denseV(), alg.denseW() }) {
			for (double[] row : mat) {
				for (double v : row) {
					long d = reducedDenominatorOfDouble(v);
					if (d != 1L) anyNonInteger = true;
					if (d % 2L == 0L) anyEvenDen = true;
					if (d % 3L == 0L) anyDen3 = true;
				}
			}
		}
		java.util.Set<String> drop = contradictedTagSet(tags, anyNonInteger, anyEvenDen, anyDen3);
		if (drop.isEmpty()) {
			return tags;
		}
		java.util.List<String> out = new java.util.ArrayList<>(tags.size());
		for (String t : tags) {
			if (!drop.contains(t)) out.add(t);
		}
		return out;
	}

	/** The contradicted field tags as a set (drop list), given the three denominator flags. */
	private static java.util.Set<String> contradictedTagSet(java.util.List<String> tags,
			boolean anyNonInteger, boolean anyEvenDen, boolean anyDen3) {
		java.util.Set<String> drop = new java.util.LinkedHashSet<>();
		if (anyNonInteger && tags.contains("Z")) drop.add("Z");
		if (anyEvenDen && tags.contains("F2")) drop.add("F2");
		if (anyDen3 && tags.contains("F3")) drop.add("F3");
		return drop;
	}

	/** Shared core: turn the three denominator flags into over-claim reasons for the claimed tags. */
	private static java.util.List<String> contradictions(java.util.List<String> tags,
			boolean anyNonInteger, boolean anyEvenDen, boolean anyDen3) {
		java.util.List<String> bad = new java.util.ArrayList<>();
		for (String t : contradictedTagSet(tags, anyNonInteger, anyEvenDen, anyDen3)) {
			switch (t) {
				case "Z" -> bad.add("declares Z but coefficients are not all integers (Z over-claim)");
				case "F2" -> bad.add("declares F2 but a coefficient denominator is even"
						+ " — not representable mod 2 (F2 over-claim)");
				case "F3" -> bad.add("declares F3 but a coefficient denominator is divisible by 3"
						+ " — not representable mod 3 (F3 over-claim)");
				default -> { }
			}
		}
		return bad;
	}

	/**
	 * Collect every raw U/V/W coefficient {@link JsonNode} (sparse {@code *_sparse}
	 * {@code c[]} arrays and dense {@code u/v/w} / {@code U/V/W} rows). Returns
	 * {@code null} when the scheme bears no coefficient arrays at all (stub).
	 */
	private static java.util.List<JsonNode> collectCoefficientTokens(JsonNode root) {
		java.util.List<JsonNode> out = new java.util.ArrayList<>();
		boolean found = false;
		for (String key : new String[] { "u_sparse", "v_sparse", "w_sparse" }) {
			JsonNode sp = root.get(key);
			if (sp != null && sp.isObject()) {
				found = true;
				for (JsonNode entry : sp) {
					JsonNode c = entry.get("c");
					if (c != null && c.isArray()) {
						c.forEach(out::add);
					}
				}
			}
		}
		for (String key : new String[] { "u", "v", "w", "U", "V", "W" }) {
			JsonNode arr = root.get(key);
			if (arr != null && arr.isArray()) {
				found = true;
				for (JsonNode row : arr) {
					if (row.isArray()) {
						row.forEach(out::add);
					}
				}
			}
		}
		return found ? out : null;
	}

	/**
	 * Reduced (lowest-terms) positive denominator of a single coefficient token:
	 * exact for integer / {@code "p/q"}-string tokens, rationalized for decimal
	 * tokens. Returns 1 (i.e. "treat as integer, do not flag") for anything that
	 * cannot be resolved — the gate stays a NECESSARY condition with no false
	 * positives.
	 */
	private static long reducedDenominator(JsonNode t) {
		if (t == null || t.isNull()) return 1L;
		if (t.isIntegralNumber()) return 1L;
		if (t.isTextual()) {
			String s = t.asString().trim();
			int slash = s.indexOf('/');
			if (slash < 0) {
				try {
					return reducedDenominatorOfDouble(Double.parseDouble(s));
				} catch (NumberFormatException e) {
					return 1L;
				}
			}
			try {
				long num = Long.parseLong(s.substring(0, slash).trim());
				long den = Long.parseLong(s.substring(slash + 1).trim());
				if (den == 0L) return 1L;
				long g = gcd(Math.abs(num), Math.abs(den));
				return g == 0L ? 1L : Math.abs(den) / g;
			} catch (NumberFormatException e) {
				return 1L;
			}
		}
		if (t.isNumber()) {
			return reducedDenominatorOfDouble(t.asDouble());
		}
		return 1L;
	}

	private static long reducedDenominatorOfDouble(double v) {
		if (Math.abs(v - Math.rint(v)) < INTEGER_TOL) return 1L;
		String frac = rationalize(v);  // "p/q" in lowest terms, or null when not resolvable
		if (frac == null) return 1L;
		int slash = frac.indexOf('/');
		try {
			return Long.parseLong(frac.substring(slash + 1).trim());
		} catch (NumberFormatException e) {
			return 1L;
		}
	}

	private static long gcd(long a, long b) {
		while (b != 0L) {
			long t = a % b;
			a = b;
			b = t;
		}
		return a;
	}

	/**
	 * Returns true if the scheme operates over GF(2) ONLY (characteristic-2
	 * universe — e.g. AlphaTensor ⟨4,4,4⟩=47). An integer scheme valid over
	 * many fields (Strassen: {@code fields=[F2,F3,Z,Q,R,C]}) is NOT "z2": its
	 * coefficients are integers, not mod-2 residues.
	 *
	 * <p>Decided solely from the unified {@code fields[]} (task #174): F2 is the
	 * sole field. The legacy {@code "z2"} boolean was retired (every scheme
	 * carried {@code z2:false}, even F2-valid integer schemes — it was both
	 * redundant with {@code fields[]} and routinely wrong), so a file with no
	 * {@code fields[]} is treated as not-F2-only.</p>
	 */
	public static boolean isZ2(JsonNode root) {
		java.util.List<String> tags = fieldTags(root);
		return tags.contains("F2")
				&& !tags.contains("Z") && !tags.contains("Q")
				&& !tags.contains("R") && !tags.contains("C") && !tags.contains("F3");
	}

	/**
	 * Returns true if the scheme uses genuine complex coefficients. Reads the
	 * unified {@code fields[]} (C is the sole characteristic-0 field, i.e. not
	 * also tagged R/Q/Z which would mean a real scheme merely liftable to C)
	 * when present; falls back to the legacy {@code "complex": true} flag.
	 */
	public static boolean isComplex(JsonNode root) {
		java.util.List<String> tags = fieldTags(root);
		if (!tags.isEmpty()) {
			return tags.contains("C")
					&& !tags.contains("R") && !tags.contains("Q") && !tags.contains("Z");
		}
		JsonNode v = root.get("complex");
		return v != null && v.isBoolean() && v.asBoolean();
	}

	/**
	 * Returns true if EVERY factor-matrix coefficient (U, V, W) is in the
	 * ternary set {@code {-1, 0, +1}}. This is the defining property of a
	 * "ZT" scheme in Perminov's terminology: a {@link eu.solven.matmul.algebra.Field#Z}
	 * (integer) scheme whose coefficients are further restricted to {-1,0,1}.
	 *
	 * <p><b>ZT is NOT a field</b> (and emphatically not F₂/Z₂ — that was a
	 * historical conflation). It is a sub-class of Z: every ZT scheme is a Z
	 * scheme, but not vice-versa. Callers should only treat the result as
	 * meaningful when the scheme is integer-valued (Z ∈ its fields); a
	 * rational/real scheme that happens to have all-ternary entries is not
	 * what "ZT" denotes.</p>
	 *
	 * <p>A small tolerance absorbs the {@code double} round-trip; coefficients
	 * are integers in practice so this only guards against e.g. {@code 1.0}
	 * vs {@code 0.9999999}.</p>
	 */
	public static boolean isTernary(NonCubicBilinearAlgorithm alg) {
		if (alg == null) return false;
		return isTernary(alg.denseU()) && isTernary(alg.denseV()) && isTernary(alg.denseW());
	}

	private static boolean isTernary(double[][] mat) {
		for (double[] row : mat) {
			for (double c : row) {
				double rounded = Math.rint(c);
				if (Math.abs(c - rounded) > 1e-9) return false;  // non-integer
				if (rounded < -1.0 || rounded > 1.0) return false;  // outside {-1,0,1}
			}
		}
		return true;
	}

	/**
	 * Reads the STORED {@code "zt"} boolean (the ternary-integer sub-class flag)
	 * from a scheme JSON, or {@code null} when absent.
	 *
	 * <p>The flag is computed once — from {@link #isTernary} — and stamped into
	 * each scheme JSON by the {@code MaterialiseZT} sanitization procedure.
	 * Consumers (e.g. {@code GenerateCatalogManifest}) READ it here rather than
	 * recomputing, so the definition lives in exactly one place
	 * (user 2026-06-04: "the catalog would rely on the given field").</p>
	 */
	public static Boolean readZT(JsonNode root) {
		JsonNode v = root.get("zt");
		return (v != null && v.isBoolean()) ? v.booleanValue() : null;
	}

	/** Read a complex-valued scheme (`"complex": true` in the JSON). */
	public static ComplexNonCubicBilinearAlgorithm readComplex(File f) throws IOException {
		try (Reader r = new BufferedReader(new FileReader(resolveSchemeFile(f)))) {
			return readComplex(r);
		}
	}

	public static ComplexNonCubicBilinearAlgorithm readComplex(String json) throws IOException {
		return readComplex(MAPPER.readTree(json));
	}

	public static ComplexNonCubicBilinearAlgorithm readComplex(Reader r) throws IOException {
		return readComplex(MAPPER.readTree(r));
	}

	public static ComplexNonCubicBilinearAlgorithm readComplex(JsonNode root) throws IOException {
		return fromJsonComplex(root);
	}

	private static ComplexNonCubicBilinearAlgorithm fromJsonComplex(JsonNode root) throws IOException {
		JsonNode nArr = expectArray(root, "n");
		if (nArr.size() != 3) {
			throw new IOException("'n' must have 3 entries, got " + nArr.size());
		}
		int n1 = nArr.get(0).asInt();
		int n2 = nArr.get(1).asInt();
		int n3 = nArr.get(2).asInt();
		int rank = expectField(root, "m").asInt();

		int dimU = n1 * n2, dimV = n2 * n3, dimW = n1 * n3;
		// U, V are stored row-major as [rank][dim] of [re, im] pairs.
		double[][][] u = readComplexTransposed(expectArray(root, "u"), rank, dimU);
		double[][][] v = readComplexTransposed(expectArray(root, "v"), rank, dimV);
		// W is col-major in dronperminov-extended; re-index to our row-major.
		double[][][] w = readComplexWColMajor(expectArray(root, "w"), rank, n1, n3);

		return new ComplexNonCubicBilinearAlgorithm(n1, n2, n3,
				u[0], u[1], v[0], v[1], w[0], w[1]);
	}

	/**
	 * Read rank × dim complex pairs into {@code {re[dim][rank], im[dim][rank]}}.
	 */
	private static double[][][] readComplexTransposed(JsonNode rows, int rank, int dim)
			throws IOException {
		if (rows.size() != rank) {
			throw new IOException("expected " + rank + " complex rows, got " + rows.size());
		}
		double[][] re = new double[dim][rank];
		double[][] im = new double[dim][rank];
		for (int k = 0; k < rank; k++) {
			JsonNode row = rows.get(k);
			if (row.size() != dim) {
				throw new IOException("complex row " + k + " has " + row.size()
						+ " cols, expected " + dim);
			}
			for (int j = 0; j < dim; j++) {
				JsonNode pair = row.get(j);
				if (pair.size() != 2) {
					throw new IOException("complex entry must be [re, im] pair, got size "
							+ pair.size());
				}
				re[j][k] = pair.get(0).asDouble();
				im[j][k] = pair.get(1).asDouble();
			}
		}
		return new double[][][] { re, im };
	}

	private static double[][][] readComplexWColMajor(JsonNode rows, int rank, int n1, int n3)
			throws IOException {
		int dim = n1 * n3;
		if (rows.size() != rank) {
			throw new IOException("expected " + rank + " w rows, got " + rows.size());
		}
		double[][] re = new double[dim][rank];
		double[][] im = new double[dim][rank];
		for (int k = 0; k < rank; k++) {
			JsonNode row = rows.get(k);
			if (row.size() != dim) {
				throw new IOException("w row " + k + " has " + row.size()
						+ " cols, expected " + dim);
			}
			for (int pDron = 0; pDron < dim; pDron++) {
				int i = pDron % n1;
				int j = pDron / n1;
				int pUs = i * n3 + j;
				JsonNode pair = row.get(pDron);
				re[pUs][k] = pair.get(0).asDouble();
				im[pUs][k] = pair.get(1).asDouble();
			}
		}
		return new double[][][] { re, im };
	}

	/**
	 * Returns true if the JSON declares {@code "scheme_type": "non_bilinear"}.
	 * Non-bilinear schemes (Waksman 1970, Rosowski 2019 Algorithm 1, …) have
	 * each rank-1 product mix entries from A and B inside both factors, so
	 * they don't fit {@link NonCubicBilinearAlgorithm}; use
	 * {@link #readNonBilinear} and verify via
	 * {@code Verifier.passesRandomMatmulSpotCheckNB} / {@code residualNonBilinear}.
	 */
	public static boolean isNonBilinear(JsonNode root) {
		JsonNode v = root.get("scheme_type");
		return v != null && v.isTextual() && "non_bilinear".equals(v.asString());
	}

	/**
	 * Returns true if the JSON is a lineage-only <strong>stub</strong> (no
	 * {@code u}/{@code v}/{@code w} arrays). Stubs carry just the
	 * {@code (n, m, p)} shape, the declared {@code (rank, additions)}
	 * invariants, attribution, and the {@code lineage} tree — consumers
	 * are expected to materialise the full tensors on demand by replaying
	 * the lineage through {@code RecursiveMaterialiser} (or equivalent in
	 * JS / notebook ports).
	 *
	 * <p>Used to compress the catalog: schemes with {@code maxDim > 16} are
	 * stripped to stubs because they're reproducible from their lineage,
	 * trimming ~9 GB. The {@code (rank, additions)} invariants are
	 * verified on every replay so drift is caught immediately.</p>
	 */
	public static boolean isStub(JsonNode root) {
		JsonNode v = root.get("scheme_type");
		return v != null && v.isTextual() && "stub".equals(v.asString());
	}

	/**
	 * Whether a scheme is flagged {@code "corrupted": true} — a lineage-only stub
	 * whose recipe can no longer be replayed into explicit matrices (legacy stubs
	 * that encode their concat operands in a permuted axis frame the replayer
	 * can't globally reconcile). Stamped by {@code EnrichSchemeMetrics} when replay
	 * fails and cleared automatically once a future replay succeeds. Such a scheme
	 * carries a <em>claimed</em> rank but is NOT currently reproducible, so it is
	 * treated as ABSENT by the search-gating lookup (it must not shadow a fresh,
	 * verifiable discovery) while still being shown — flagged — in the catalog.
	 */
	public static boolean isCorrupted(JsonNode root) {
		JsonNode v = root.get("corrupted");
		return v != null && v.isBoolean() && v.asBoolean();
	}

	/**
	 * Whether a scheme has been formally verified (its expanded matrices passed
	 * the symbolic / spot-check). Defaults to {@code true} when the field is
	 * absent: every legacy on-disk scheme was built and verified before being
	 * written, so only an explicit {@code "verified": false} (a detect-only stub
	 * recorded from a search prediction, not yet machine-checked) is unverified.
	 * The Phase-2 validate-and-prune batch flips it to {@code true} or deletes
	 * the scheme.
	 */
	public static boolean isVerified(JsonNode root) {
		JsonNode v = root.get("verified");
		return v == null || v.isNull() || v.asBoolean(true);
	}

	/**
	 * Read a non-bilinear scheme written by
	 * {@code MaterializeRosowskiAlgorithm1} / {@code MaterializeWaksman1970}.
	 * Format: sparse {@code [[row, value], …]} arrays for
	 * {@code u_a, u_b, v_a, v_b, w}, one per product (length = rank);
	 * {@code "n": [n, m, p]}, {@code "rank": r}.
	 */
	public static NonBilinearAlgorithm readNonBilinear(File f) throws IOException {
		try (Reader r = new BufferedReader(new FileReader(resolveSchemeFile(f)))) {
			return readNonBilinear(r);
		}
	}

	public static NonBilinearAlgorithm readNonBilinear(Reader r) throws IOException {
		return readNonBilinear(MAPPER.readTree(r));
	}

	public static NonBilinearAlgorithm readNonBilinear(JsonNode root) throws IOException {
		JsonNode nArr = expectArray(root, "n");
		if (nArr.size() != 3) {
			throw new IOException("'n' must have 3 entries, got " + nArr.size());
		}
		int n = nArr.get(0).asInt();
		int m = nArr.get(1).asInt();
		int p = nArr.get(2).asInt();
		// Tolerate either "rank" (our convention) or "m" (dronperminov-ish).
		int rank;
		if (root.has("rank")) rank = root.get("rank").asInt();
		else if (root.has("m")) rank = root.get("m").asInt();
		else throw new IOException("missing 'rank' (or 'm') field");

		// NonBilinearAlgorithm field dimensions:
		//   {Ua, Va} hold A-entry coefficients (one for each factor) → dimA = n·m
		//   {Ub, Vb} hold B-entry coefficients                       → dimB = m·p
		//   W holds output coefficients                              → dimC = n·p
		// For cubic shapes (n=m=p) all three are equal; for Rosowski ⟨n,3,3⟩
		// they differ (dimA = 3n, dimB = 9, dimC = 3n).
		int dimA = n * m, dimB = m * p, dimC = n * p;
		double[][] Ua = readSparseFactorNB(expectArray(root, "u_a"), rank, dimA);
		double[][] Ub = readSparseFactorNB(expectArray(root, "u_b"), rank, dimB);
		double[][] Va = readSparseFactorNB(expectArray(root, "v_a"), rank, dimA);
		double[][] Vb = readSparseFactorNB(expectArray(root, "v_b"), rank, dimB);
		double[][] W  = readSparseFactorNB(expectArray(root, "w"),   rank, dimC);
		return new NonBilinearAlgorithm(n, m, p, Ua, Ub, Va, Vb, W);
	}

	/**
	 * Parse {@code rank} arrays each containing {@code [row, value]} pairs into
	 * a dense {@code dim × rank} matrix (column-per-product layout matches
	 * {@link NonBilinearAlgorithm}'s field convention).
	 */
	private static double[][] readSparseFactorNB(JsonNode rows, int rank, int dim) throws IOException {
		if (rows.size() != rank) {
			throw new IOException("sparse non-bilinear factor: rows=" + rows.size()
					+ ", expected rank=" + rank);
		}
		double[][] dst = new double[dim][rank];
		for (int k = 0; k < rank; k++) {
			JsonNode col = rows.get(k);
			for (JsonNode pair : col) {
				if (pair.size() != 2) {
					throw new IOException("non-bilinear sparse entry must be [row, value], got size "
							+ pair.size());
				}
				int row = pair.get(0).asInt();
				if (row < 0 || row >= dim) {
					throw new IOException("non-bilinear row out of range: " + row + " / " + dim);
				}
				dst[row][k] += pair.get(1).asDouble();
			}
		}
		return dst;
	}

	/**
	 * Returns true if the JSON uses the "reduced" sparse-list format (entries
	 * in U/V/W are lists of {@code {"index": …, "value": …}} pairs instead of
	 * dense int arrays). Found in many dronperminov `*_reduced` files.
	 */
	public static boolean isReduced(JsonNode root) {
		return isReducedJson(root);
	}

	private static boolean isReducedJson(JsonNode root) {
		JsonNode u = root.get("u");
		if (u == null || !u.isArray() || u.isEmpty()) return false;
		JsonNode firstRow = u.get(0);
		if (!firstRow.isArray() || firstRow.isEmpty()) return false;
		// Reduced rows hold JSON objects; standard rows hold numbers.
		return firstRow.get(0).isObject();
	}

	/**
	 * Reads a "reduced" sparse-list scheme.
	 *
	 * <p>Layout (per file inspection):</p>
	 * <ul>
	 *   <li>{@code u}: rank rows; each row a list of {@code {"index": dim, "value": coef}}
	 *       (the non-zero positions of column k of U).</li>
	 *   <li>{@code v}: same shape (rank rows × sparse dim entries).</li>
	 *   <li>{@code w}: <b>dim rows × sparse rank entries</b> (note: transposed
	 *       relative to u/v). Each row corresponds to a C-position in col-major
	 *       flatten ({@code p_dron = j·n_1 + i}); entries are
	 *       {@code {"index": mult, "value": coef}}.</li>
	 * </ul>
	 *
	 * <p>Converted to our internal row-major C-flatten W on the way in.</p>
	 */
	public static NonCubicBilinearAlgorithm readReduced(File f) throws IOException {
		try (Reader r = new BufferedReader(new FileReader(resolveSchemeFile(f)))) {
			return readReduced(r);
		}
	}

	public static NonCubicBilinearAlgorithm readReduced(JsonNode root) throws IOException {
		return readReducedFromNode(root);
	}

	/**
	 * Parse just the {@code lineage} field from a scheme JSON, returning
	 * the {@link Lineage.Node} tree if present, else empty. Used by
	 * {@code RecursiveMaterialiser} to splice catalog leaves' own deep
	 * lineage into a composed scheme's lineage tree (avoiding the
	 * "shallow leaf" bug).
	 */
	public static Optional<Lineage.Node> readLineage(File f) throws IOException {
		return readLineage(parseJson(f));
	}

	public static Optional<Lineage.Node> readLineage(JsonNode root) {
		JsonNode lineage = root.get("lineage");
		if (lineage == null || lineage.isNull()) return Optional.empty();
		return Optional.of(parseLineageNode(lineage, new java.util.HashMap<>()));
	}

	private static Lineage.Node parseLineageNode(JsonNode node, Map<String, Lineage.Node> idMap) {
		String op = node.path("op").asText("");
		Lineage.Node built = switch (op) {
			case "Leaf", "Atom" -> new Lineage.Atom(node.path("ref").asText("?"));
			case "@ref" -> {
				Lineage.Node ref = idMap.get(node.path("id").asText());
				yield ref != null ? ref : new Lineage.Atom("@ref?:" + node.path("id").asText());
			}
			case "KronProduct" -> new Lineage.KronProduct(
					parseLineageNode(node.get("outer"), idMap),
					parseLineageNode(node.get("inner"), idMap));
			// "ConcatRight"/"ConcatBelow" are the legacy op names (pre-2026-06
			// rename); kept as read-compat aliases for on-disk catalog files.
			case "ConcatCols", "ConcatRight" -> new Lineage.ConcatCols(
					parseLineageNode(node.get("left"), idMap),
					parseLineageNode(node.get("right"), idMap));
			case "ConcatRows", "ConcatBelow" -> new Lineage.ConcatRows(
					parseLineageNode(node.get("top"), idMap),
					parseLineageNode(node.get("bottom"), idMap));
			case "SumInner" -> new Lineage.SumInner(
					parseLineageNode(node.get("left"), idMap),
					parseLineageNode(node.get("right"), idMap));
			case "Transpose" -> new Lineage.Transpose(
					parseLineageNode(node.get("child"), idMap),
					node.path("perm").asText(""));
			case "OrientAs" -> new Lineage.OrientAs(
					parseLineageNode(node.get("child"), idMap),
					node.path("n").asInt(0), node.path("m").asInt(0), node.path("p").asInt(0),
					node.hasNonNull("axisMap") ? Lineage.axisMapFromStr(node.get("axisMap").asString()) : null);
			case "DCE" -> new Lineage.Dce(parseLineageNode(node.get("child"), idMap));
			case "PeeledViaTa" -> new Lineage.PeeledViaTa(
					node.path("n").asInt(0), node.path("s").asInt(0),
					parseLineageNode(node.get("cube"), idMap),
					parseLineageNode(node.get("corner"), idMap));
			case "KronChain" -> {
				List<Lineage.Node> factors = new java.util.ArrayList<>();
				JsonNode f = node.get("factors");
				if (f != null && f.isArray()) {
					for (JsonNode child : f) factors.add(parseLineageNode(child, idMap));
				}
				yield new Lineage.KronChain(factors);
			}
			// "RecombinationN" is the legacy op name (pre-rename); read-compat alias.
			case "Recombination", "RecombinationN" -> {
				int[] aA = parseIntArray(node.get("allocA"));
				int[] aB = parseIntArray(node.get("allocB"));
				int[] aC = parseIntArray(node.get("allocC"));
				List<Lineage.Node> leaves = new java.util.ArrayList<>();
				JsonNode lvs = node.get("leaves");
				if (lvs != null && lvs.isArray()) {
					for (JsonNode child : lvs) leaves.add(parseLineageNode(child, idMap));
				}
				yield new Lineage.RecombinationN(
						parseLineageNode(node.get("base"), idMap), aA, aB, aC, leaves);
			}
			case "RecombinationTa", "RecombinationTaN" -> {
				int[] aA = parseIntArray(node.get("allocA"));
				int[] aB = parseIntArray(node.get("allocB"));
				int[] aC = parseIntArray(node.get("allocC"));
				List<Lineage.Node> leaves = new java.util.ArrayList<>();
				JsonNode lvs = node.get("leaves");
				if (lvs != null && lvs.isArray()) {
					for (JsonNode child : lvs) leaves.add(parseLineageNode(child, idMap));
				}
				yield new Lineage.RecombinationTaN(
						parseLineageNode(node.get("base"), idMap), aA, aB, aC, leaves);
			}
			case "RecombinationWithPair", "RecombinationWithPairN" -> {
				int[][] pairs = parseInt2DArray(node.get("pairs"));
				int[] solo = parseIntArray(node.get("solo"));
				List<Lineage.Node> leaves = new java.util.ArrayList<>();
				JsonNode lvs = node.get("leaves");
				if (lvs != null && lvs.isArray()) {
					for (JsonNode child : lvs) leaves.add(parseLineageNode(child, idMap));
				}
				yield new Lineage.RecombinationWithPairN(
						parseLineageNode(node.get("base"), idMap), pairs, solo, leaves);
			}
			case "AugmentSquareDiscard" -> new Lineage.AugmentSquareDiscard(
					node.path("p").asInt(0), node.path("n").asInt(0),
					parseLineageNode(node.get("square"), idMap));
			case "AxisFlip" -> {
				// Accept two shapes: {child: Node, mask: int} (current) and
				// the older {ref: String, mask: int} where the child was
				// inlined as a leaf reference.
				Lineage.Node child;
				if (node.has("child") && !node.get("child").isNull()) {
					child = parseLineageNode(node.get("child"), idMap);
				} else if (node.has("ref")) {
					child = new Lineage.Atom(node.path("ref").asText("?"));
				} else {
					child = new Lineage.Atom("axisflip-missing-child");
				}
				yield new Lineage.AxisFlip(child, node.path("mask").asInt(0));
			}
			case "AxisPermute" -> new Lineage.AxisPermute(
					parseLineageNode(node.get("child"), idMap),
					parseIntArray(node.get("permA")),
					parseIntArray(node.get("permB")),
					parseIntArray(node.get("permC")));
			case "DisjointSum" -> {
				List<Lineage.Node> children = new java.util.ArrayList<>();
				JsonNode ch = node.get("children");
				if (ch != null && ch.isArray()) {
					for (JsonNode c : ch) children.add(parseLineageNode(c, idMap));
				}
				List<List<Integer>> taLegs = new java.util.ArrayList<>();
				JsonNode tl = node.get("taLegs");
				if (tl != null && tl.isArray()) {
					for (JsonNode leg : tl) {
						List<Integer> g = new java.util.ArrayList<>();
						if (leg.isArray()) for (JsonNode i : leg) g.add(i.asInt());
						taLegs.add(g);
					}
				}
				yield new Lineage.DisjointSum(children, taLegs);
			}
			case "SerendipitousProduct" -> new Lineage.SerendipitousProduct(
					parseLineageNode(node.get("base"), idMap),
					node.path("n2").asInt(0), node.path("m2").asInt(0), node.path("p2").asInt(0));
			case "Project" -> new Lineage.Project(
					parseLineageNode(node.get("child"), idMap),
					parseIntArray(node.get("keepN")),
					parseIntArray(node.get("keepM")),
					parseIntArray(node.get("keepP")));
			default -> new Lineage.Atom("unknown-op:" + op);
		};
		// Register if this node has an id (for later @ref resolution).
		String id = node.path("id").asText("");
		if (!id.isEmpty()) idMap.put(id, built);
		return built;
	}

	private static int[] parseIntArray(JsonNode n) {
		if (n == null || !n.isArray()) return new int[0];
		int[] out = new int[n.size()];
		for (int i = 0; i < n.size(); i++) out[i] = n.get(i).asInt();
		return out;
	}

	private static int[][] parseInt2DArray(JsonNode n) {
		if (n == null || !n.isArray()) return new int[0][];
		int[][] out = new int[n.size()][];
		for (int i = 0; i < n.size(); i++) out[i] = parseIntArray(n.get(i));
		return out;
	}

	public static NonCubicBilinearAlgorithm readReduced(Reader r) throws IOException {
		return readReducedFromNode(MAPPER.readTree(r));
	}

	private static NonCubicBilinearAlgorithm readReducedFromNode(JsonNode root) throws IOException {
		JsonNode nArr = expectArray(root, "n");
		int n1 = nArr.get(0).asInt();
		int n2 = nArr.get(1).asInt();
		int n3 = nArr.get(2).asInt();
		int rank = expectField(root, "m").asInt();
		int dimU = n1 * n2, dimV = n2 * n3, dimW = n1 * n3;

		JsonNode uFresh = root.has("u_fresh") ? root.get("u_fresh") : null;
		JsonNode vFresh = root.has("v_fresh") ? root.get("v_fresh") : null;
		JsonNode wFresh = root.has("w_fresh") ? root.get("w_fresh") : null;

		// u_fresh/v_fresh combine A/B-entries (dim axis); w_fresh combines
		// MULTIPLICATIONS (rank axis) — different from u/v_fresh. Pass the
		// right "original space" size to expandFresh.
		double[][] freshU = expandFresh(uFresh, dimU);
		double[][] freshV = expandFresh(vFresh, dimV);
		double[][] freshW = expandFresh(wFresh, rank);

		double[][] U = readSparseFactorWithFresh(expectArray(root, "u"), rank, dimU, freshU, /*rowIsRank=*/ true);
		double[][] V = readSparseFactorWithFresh(expectArray(root, "v"), rank, dimV, freshV, /*rowIsRank=*/ true);
		double[][] W = readSparseWWithFresh(expectArray(root, "w"), rank, n1, n3, freshW);
		return new NonCubicBilinearAlgorithm(n1, n2, n3, U, V, W);
	}

	/**
	 * Expand a {@code u_fresh} / {@code v_fresh} / {@code w_fresh} block into a
	 * dense {@code dim × freshCount} matrix giving the linear combination of
	 * original positions that each fresh intermediate represents. Handles
	 * nested fresh references (fresh entry k may reference fresh[j] for j < k)
	 * by processing chronologically and expanding through the in-progress matrix.
	 */
	private static double[][] expandFresh(JsonNode freshRows, int dim) {
		int n = freshRows == null ? 0 : freshRows.size();
		double[][] expansion = new double[dim][n];
		for (int k = 0; k < n; k++) {
			JsonNode row = freshRows.get(k);
			for (JsonNode obj : row) {
				int idx = obj.get("index").asInt();
				double val = obj.get("value").asDouble();
				if (idx < dim) {
					expansion[idx][k] += val;
				} else {
					int freshIdx = idx - dim;
					// Nested ref: substitute the already-expanded fresh entry.
					for (int p = 0; p < dim; p++) {
						double f = expansion[p][freshIdx];
						if (f == 0.0) continue;
						expansion[p][k] += val * f;
					}
				}
			}
		}
		return expansion;
	}

	private static double[][] readSparseFactorWithFresh(JsonNode rows, int rank, int dim,
			double[][] fresh, boolean rowIsRank) throws IOException {
		double[][] dst = new double[dim][rank];
		int outerLen = rowIsRank ? rank : dim;
		if (rows.size() != outerLen) {
			throw new IOException("sparse rows = " + rows.size() + ", expected " + outerLen);
		}
		int freshCount = fresh[0].length;
		for (int outer = 0; outer < outerLen; outer++) {
			JsonNode row = rows.get(outer);
			for (JsonNode obj : row) {
				int idx = obj.get("index").asInt();
				double val = obj.get("value").asDouble();
				if (idx < dim) {
					if (rowIsRank) dst[idx][outer] += val;
					else dst[outer][idx] += val;
				} else {
					int freshIdx = idx - dim;
					if (freshIdx >= freshCount) {
						throw new IOException("fresh index out of range: " + idx
								+ " (dim=" + dim + ", fresh=" + freshCount + ")");
					}
					// Expand fresh[freshIdx] = Σ_p expansion[p][freshIdx] · a[p]
					for (int p = 0; p < dim; p++) {
						double f = fresh[p][freshIdx];
						if (f == 0.0) continue;
						double contrib = val * f;
						if (rowIsRank) dst[p][outer] += contrib;
						else dst[outer][p] += contrib;
					}
				}
			}
		}
		return dst;
	}

	private static double[][] readSparseWWithFresh(JsonNode rows, int rank, int n1, int n3,
			double[][] fresh) throws IOException {
		int dim = n1 * n3;
		double[][] dst = new double[dim][rank];
		if (rows.size() != dim) {
			throw new IOException("sparse w rows = " + rows.size() + ", expected " + dim);
		}
		int freshCount = fresh[0].length;
		for (int pDron = 0; pDron < dim; pDron++) {
			JsonNode row = rows.get(pDron);
			int i = pDron % n1;
			int j = pDron / n1;
			int pUs = i * n3 + j;
			for (JsonNode obj : row) {
				int k = obj.get("index").asInt();
				double val = obj.get("value").asDouble();
				if (k < rank) {
					dst[pUs][k] += val;
				} else {
					int freshIdx = k - rank;
					if (freshIdx >= freshCount) {
						throw new IOException("w fresh index out of range: " + k);
					}
					// w_fresh is keyed over the dim axis; same expansion idea.
					for (int q = 0; q < rank; q++) {
						double f = fresh[q][freshIdx];
						if (f == 0.0) continue;
						dst[pUs][q] += val * f;
					}
				}
			}
		}
		return dst;
	}

	private static NonCubicBilinearAlgorithm fromJson(JsonNode root) throws IOException {
		JsonNode nArr = expectArray(root, "n");
		if (nArr.size() != 3) {
			throw new IOException("'n' must have 3 entries, got " + nArr.size());
		}
		int n1 = nArr.get(0).asInt();
		int n2 = nArr.get(1).asInt();
		int n3 = nArr.get(2).asInt();
		int rank = expectField(root, "m").asInt();

		int dimU = n1 * n2, dimV = n2 * n3, dimW = n1 * n3;
		// dronperminov stores rank-rows × dim-cols; we transpose into dim-rows × rank-cols.
		double[][] U = readTransposed(expectArray(root, "u"), rank, dimU);
		double[][] V = readTransposed(expectArray(root, "v"), rank, dimV);
		// W uses col-major C-flatten in their format → re-index to our row-major.
		double[][] W = readWColMajor(expectArray(root, "w"), rank, n1, n3);

		return new NonCubicBilinearAlgorithm(n1, n2, n3, U, V, W);
	}

	private static double[][] readTransposed(JsonNode rows, int expectedRows, int expectedCols)
			throws IOException {
		if (rows.size() != expectedRows) {
			throw new IOException("expected " + expectedRows + " rows, got " + rows.size());
		}
		double[][] dst = new double[expectedCols][expectedRows];
		for (int k = 0; k < expectedRows; k++) {
			JsonNode row = rows.get(k);
			if (row.size() != expectedCols) {
				throw new IOException("row " + k + " has " + row.size()
						+ " cols, expected " + expectedCols);
			}
			for (int j = 0; j < expectedCols; j++) {
				dst[j][k] = parseCoef(row.get(j));
			}
		}
		return dst;
	}

	// ───────────────────────────────────────────────────────────────────────────
	// Sparse format (our extension for max-dim ≥ 9)
	//
	// Layout: alongside (or instead of) the dense `u`/`v`/`w` triple, the JSON
	// may carry `u_sparse`/`v_sparse`/`w_sparse`. Each is a rank-length array;
	// each element is an OBJECT mapping the position-index (as a string) to
	// the non-zero coefficient. Zero entries are omitted entirely:
	//
	//   "u_sparse": [
	//     {"0": 1, "40": 1},
	//     {"1": -1, "41": 1},
	//     ...
	//   ]
	//
	// W keeps dronperminov's col-major C-flatten convention (positions are
	// `j·n_1 + i`), so the sparse map for W uses those indices too.
	// ───────────────────────────────────────────────────────────────────────────

	private static NonCubicBilinearAlgorithm fromJsonSparse(JsonNode root) throws IOException {
		JsonNode nArr = expectArray(root, "n");
		if (nArr.size() != 3) {
			throw new IOException("'n' must have 3 entries, got " + nArr.size());
		}
		int n1 = nArr.get(0).asInt();
		int n2 = nArr.get(1).asInt();
		int n3 = nArr.get(2).asInt();
		int rank = expectField(root, "m").asInt();

		int dimU = n1 * n2, dimV = n2 * n3, dimW = n1 * n3;
		double[][] U = readSparseFactorAny(expectField(root, "u_sparse"), rank, dimU);
		double[][] V = readSparseFactorAny(expectField(root, "v_sparse"), rank, dimV);
		// W: convert dronperminov col-major positions to our row-major.
		double[][] W = readSparseFactorWAny(expectField(root, "w_sparse"), rank, n1, n3);

		return new NonCubicBilinearAlgorithm(n1, n2, n3, U, V, W);
	}

	/**
	 * Dispatch between the legacy array-of-per-key-maps format and the new
	 * row-oriented map format keyed by stringified product index with
	 * parallel {@code i} (indices) and {@code c} (coefficients) arrays.
	 * See class header for format details.
	 */
	private static double[][] readSparseFactorAny(JsonNode node, int rank, int dim) throws IOException {
		if (node.isArray()) {
			return readSparseObjectFactor(node, rank, dim);
		}
		if (node.isObject()) {
			return readSparseRowMapFactor(node, rank, dim);
		}
		throw new IOException("u_sparse/v_sparse must be array (legacy) or object (row-oriented), got "
				+ node.getNodeType());
	}

	private static double[][] readSparseFactorWAny(JsonNode node, int rank, int n1, int n3) throws IOException {
		if (node.isArray()) {
			return readSparseObjectFactorW(node, rank, n1, n3);
		}
		if (node.isObject()) {
			return readSparseRowMapFactorW(node, rank, n1, n3);
		}
		throw new IOException("w_sparse must be array (legacy) or object (row-oriented), got "
				+ node.getNodeType());
	}

	private static double[][] readSparseObjectFactor(JsonNode rows, int rank, int dim)
			throws IOException {
		if (rows.size() != rank) {
			throw new IOException("expected " + rank + " sparse rows, got " + rows.size());
		}
		double[][] dst = new double[dim][rank];
		for (int k = 0; k < rank; k++) {
			JsonNode obj = rows.get(k);
			for (var e : obj.properties()) {
				int pos = Integer.parseInt(e.getKey());
				double val = parseCoef(e.getValue());
				dst[pos][k] = val;
			}
		}
		return dst;
	}

	private static double[][] readSparseObjectFactorW(JsonNode rows, int rank, int n1, int n3)
			throws IOException {
		int dim = n1 * n3;
		if (rows.size() != rank) {
			throw new IOException("expected " + rank + " w_sparse rows, got " + rows.size());
		}
		double[][] dst = new double[dim][rank];
		for (int k = 0; k < rank; k++) {
			JsonNode obj = rows.get(k);
			for (var e : obj.properties()) {
				int pDron = Integer.parseInt(e.getKey());
				int i = pDron % n1;
				int j = pDron / n1;
				int pUs = i * n3 + j;
				dst[pUs][k] = parseCoef(e.getValue());
			}
		}
		return dst;
	}

	/**
	 * Reads the new row-oriented sparse format:
	 * <pre>
	 * "u_sparse": {
	 *   "0": {"i": [0, 4, 8], "c": [1, -1, 1]},
	 *   "1": {"i": [...],     "c": [...]},
	 *   ...
	 * }
	 * </pre>
	 * Keys are stringified product indices in [0, rank). Iteration order follows
	 * Jackson's insertion order — we look up each {@code "k"} explicitly so any
	 * key reordering by upstream tooling stays correct.
	 */
	private static double[][] readSparseRowMapFactor(JsonNode obj, int rank, int dim) throws IOException {
		double[][] dst = new double[dim][rank];
		for (int k = 0; k < rank; k++) {
			JsonNode row = obj.get(Integer.toString(k));
			if (row == null) {
				throw new IOException("row-oriented sparse: missing key \"" + k + "\"");
			}
			JsonNode iArr = expectArray(row, "i");
			JsonNode cArr = expectArray(row, "c");
			if (iArr.size() != cArr.size()) {
				throw new IOException("row-oriented sparse row " + k + ": i.size=" + iArr.size()
						+ " ≠ c.size=" + cArr.size());
			}
			for (int t = 0; t < iArr.size(); t++) {
				int pos = iArr.get(t).asInt();
				dst[pos][k] = parseCoef(cArr.get(t));
			}
		}
		return dst;
	}

	private static double[][] readSparseRowMapFactorW(JsonNode obj, int rank, int n1, int n3) throws IOException {
		int dim = n1 * n3;
		double[][] dst = new double[dim][rank];
		for (int k = 0; k < rank; k++) {
			JsonNode row = obj.get(Integer.toString(k));
			if (row == null) {
				throw new IOException("row-oriented w_sparse: missing key \"" + k + "\"");
			}
			JsonNode iArr = expectArray(row, "i");
			JsonNode cArr = expectArray(row, "c");
			if (iArr.size() != cArr.size()) {
				throw new IOException("row-oriented w_sparse row " + k + ": i.size=" + iArr.size()
						+ " ≠ c.size=" + cArr.size());
			}
			for (int t = 0; t < iArr.size(); t++) {
				int pDron = iArr.get(t).asInt();
				int i = pDron % n1;
				int j = pDron / n1;
				int pUs = i * n3 + j;
				dst[pUs][k] = parseCoef(cArr.get(t));
			}
		}
		return dst;
	}

	/**
	 * Parse a coefficient that may be a JSON number or a string fraction
	 * (dronperminov's `Q` schemes use strings like {@code "1/2"}, {@code "-1/4"}).
	 */
	private static double parseCoef(Object v) {
		if (v instanceof JsonNode node) {
			if (node.isNumber()) return node.asDouble();
			if (node.isTextual()) return parseCoefString(node.asText());
			throw new IllegalArgumentException("unexpected JsonNode coefficient: " + node);
		}
		if (v instanceof Number n) return n.doubleValue();
		if (v instanceof String s) return parseCoefString(s);
		throw new IllegalArgumentException("unexpected coefficient: " + v
				+ " (" + (v == null ? "null" : v.getClass().getSimpleName()) + ")");
	}

	private static double parseCoefString(String s) {
		String t = s.trim();
		int slash = t.indexOf('/');
		if (slash >= 0) {
			double num = Double.parseDouble(t.substring(0, slash));
			double den = Double.parseDouble(t.substring(slash + 1));
			return num / den;
		}
		return Double.parseDouble(t);
	}

	/**
	 * Read W with dronperminov's col-major C-flatten and convert to our row-major.
	 * Their position {@code p_dron = j·n_1 + i}; ours {@code p_us = i·n_3 + j}.
	 */
	private static double[][] readWColMajor(JsonNode rows, int rank, int n1, int n3)
			throws IOException {
		int dim = n1 * n3;
		if (rows.size() != rank) {
			throw new IOException("expected " + rank + " rows in w, got " + rows.size());
		}
		double[][] dst = new double[dim][rank];
		for (int k = 0; k < rank; k++) {
			JsonNode row = rows.get(k);
			if (row.size() != dim) {
				throw new IOException("w row " + k + " has " + row.size()
						+ " cols, expected " + dim);
			}
			for (int pDron = 0; pDron < dim; pDron++) {
				int i = pDron % n1;
				int j = pDron / n1;
				int pUs = i * n3 + j;
				dst[pUs][k] = parseCoef(row.get(pDron));
			}
		}
		return dst;
	}

	private static JsonNode expectArray(JsonNode obj, String key) throws IOException {
		JsonNode v = obj.get(key);
		if (v == null || !v.isArray()) {
			throw new IOException("expected array at key '" + key + "', got " + v);
		}
		return v;
	}

	private static JsonNode expectField(JsonNode obj, String key) throws IOException {
		JsonNode v = obj.get(key);
		if (v == null || v.isNull()) {
			throw new IOException("missing field '" + key + "'");
		}
		return v;
	}

	// ───────────────────────────────────────────────────────────────────────────
	// Writers
	// ───────────────────────────────────────────────────────────────────────────

	/**
	 * The {@code _b{bud_score}} filename token — bud-richness as a third
	 * optimisation axis alongside {@code _m{rank}} and {@code _a{additions}}
	 * (user 2026-06-06). {@code bud_score} is the total rank-one terms living in
	 * buds ({@link BudParetoSelection#budScore}); a materialised scheme is always
	 * checkable, so {@code _b0} legitimately means "checked, no buds" (vs a
	 * missing token = buds not computed, e.g. an un-validated stub). Insert it
	 * right after {@code _m…[_a…]} and before any field suffix / {@code .json}.
	 */
	public static String budFilenameToken(NonCubicBilinearAlgorithm alg) {
		SerendipitousBudProduct.BudSummary bs = SerendipitousBudProduct.summarise(alg);
		return "_b" + BudParetoSelection.budScore(bs.summary());
	}

	/**
	 * Switching threshold: dense {@code u}/{@code v}/{@code w} for
	 * {@code max(n,m,p) < this}, sparse {@code *_sparse} at or above (task #187,
	 * user 2026-06-03). Small schemes stay human-readable as dense matrices;
	 * large schemes use the compact sparse map.
	 *
	 * <p><strong>Canonical value lives in
	 * {@link CatalogLimits#SPARSE_DIM_THRESHOLD}</strong> (single home for all
	 * on-disk JSON-format constants); this is a delegating alias. Lowered from
	 * {@code 16} to {@code 9} in 2026-06 per the dense-vs-sparse size study —
	 * dense up to {@code 2³ = 8}, sparse from {@code 9}. Existing dim 9–15 dense
	 * files keep loading (the reader auto-detects {@code u_sparse}); only new
	 * writes switch to sparse until a migration re-serialises them.</p>
	 */
	public static final int SPARSE_DIM_THRESHOLD = CatalogLimits.SPARSE_DIM_THRESHOLD;

	public static void write(NonCubicBilinearAlgorithm alg, File f) throws IOException {
		write(alg, f, null);
	}

	/**
	 * Write a lineage-only <strong>stub</strong> (no factor matrices) — the
	 * on-disk form for derived schemes above {@link CatalogPolicy#MATERIALISE_MAX_DIM}.
	 * Carries only {@code {n=[n,m,p], m=rank, lineage_compact, lineage,
	 * scheme_type="stub"}}; the explicit U/V/W are reproduced on demand by
	 * {@link eu.solven.matmul.search.LineageReplayer}. Keeps the catalog small
	 * while still recording the construction + rank. {@code lineage} is required
	 * (a stub with no lineage would be irreproducible).
	 */
	// ── content hash: a unique, representation-invariant scheme reference ──
	/**
	 * Canonical content hash: SHA-256 over {@code (shape ; nonzero entries of
	 * U/V/W with fixed-precision coefficients)}. <strong>Representation
	 * invariant</strong> — the sparse and dense JSON of the same matrices hash
	 * identically — and cheap (cost ∝ nonzeros, not dense size). Distinguishes two
	 * schemes of the same {@code (shape, rank, additions)} that differ in content
	 * (e.g. a bud-rich vs bud-poor ⟨3,3,3⟩=23), which the canonical filename key
	 * cannot. The full hex is the authoritative unique id; {@link #shortHash} is
	 * the filename {@code _h} token (may collide → caller disambiguates by
	 * property). Coefficients are encoded at 1e-6 precision (same convention as the
	 * bud direction canonicaliser).
	 */
	public static String contentHash(NonCubicBilinearAlgorithm a) {
		// A scheme with any non-integer coefficient is hashed over the EXACT
		// canonical token ("1/4", …) — the same symbolic form written to disk —
		// so the hash is stable under the double⇄string round-trip (a truncated
		// decimal and the exact rational no longer hash apart). Pure-integer
		// schemes keep the legacy round(v·1e6) basis byte-for-byte, so their
		// hashes (and filenames) are unchanged.
		boolean exact = hasNonInteger(a.denseU()) || hasNonInteger(a.denseV()) || hasNonInteger(a.denseW());
		StringBuilder sb = new StringBuilder();
		sb.append(a.n).append('x').append(a.m).append('x').append(a.p).append('r').append(a.r).append(';');
		appendFactor(sb, 'U', a.denseU(), a.r, exact);
		appendFactor(sb, 'V', a.denseV(), a.r, exact);
		appendFactor(sb, 'W', a.denseW(), a.r, exact);
		return sha256Hex(sb.toString());
	}

	/**
	 * The canonical scheme filename — {@code {n}x{m}x{p}-r{rank}-{note}-{hash7}.json},
	 * the 2026-06 convention {@link eu.solven.matmul.docs.migrate.RenameSchemes}
	 * normalises to. <b>Every writer must use this</b> so the catalog never drifts
	 * back to legacy {@code {note}-{shape}_m{rank}_a{adds}} names. Filenames are pure
	 * labels: addition count, bud-score, field tag, etc. live in JSON content, NOT the
	 * name. {@code note} is the source/author for imports (e.g. {@code perminov_2025})
	 * or our-output tag for derivations ({@code derived}, {@code derived_recursive},
	 * {@code metaflip}, …). The 7-hex content-hash suffix dedups schemes sharing a
	 * {@code (shape, rank)} and matches a stub's stored {@code "hash"} field.
	 */
	public static String canonicalName(NonCubicBilinearAlgorithm alg, String note) {
		return alg.n + "x" + alg.m + "x" + alg.p + "-r" + alg.r + "-" + note + "-"
				+ contentHash(alg).substring(0, 7) + ".json";
	}

	/**
	 * Content hash for a COMPLEX scheme, whose factors are stored as {@code [re,im]}
	 * pairs and so cannot expand to a {@code double[][]}. Hashes the shape + the
	 * verbatim JSON list representation of each factor ("hashing the list is good
	 * enough" — it is deterministic and exact for the stored coefficients).
	 */
	public static String contentHashComplexJson(tools.jackson.databind.JsonNode root) {
		StringBuilder sb = new StringBuilder();
		tools.jackson.databind.JsonNode n = root.get("n");
		int rank = root.has("m") ? root.get("m").asInt() : root.path("rank").asInt(-1);
		sb.append(n.get(0).asInt()).append('x').append(n.get(1).asInt()).append('x').append(n.get(2).asInt())
				.append('r').append(rank).append(';');
		for (String key : new String[] { "u", "v", "w", "u_sparse", "v_sparse", "w_sparse" }) {
			tools.jackson.databind.JsonNode mat = root.get(key);
			if (mat != null) {
				sb.append(key).append('=').append(mat.toString()).append(';');
			}
		}
		return sha256Hex(sb.toString());
	}

	private static boolean hasNonInteger(double[][] m) {
		for (double[] row : m) {
			for (double v : row) {
				if (Math.abs(v - Math.rint(v)) > INTEGER_TOL) {
					return true;
				}
			}
		}
		return false;
	}

	/** Bare canonical coefficient token: {@code "n"} for integers, {@code "p/q"} for
	 *  small rationals, else the lossless decimal — the same symbolic form
	 *  {@link #formatCoef} writes to disk (without the surrounding quotes). */
	private static String canonicalToken(double v) {
		long rounded = Math.round(v);
		if (Math.abs(v - rounded) < INTEGER_TOL) {
			return Long.toString(rounded);
		}
		String frac = rationalize(v);
		return frac != null ? frac : Double.toString(v);
	}

	/**
	 * Short content-hash prefix for the filename's last token, i.e.
	 * {@code …_{4hex}.json} (no leading letter — avoids confusion with the field
	 * tokens). 4 hex (65536 values) is ample to disambiguate the handful of
	 * schemes that ever share one shape, and {@link FieldAwareLookup#findByHash}
	 * is shape-scoped so collisions only matter within a shape. The full
	 * {@link #contentHash} stays the authoritative id (stored in JSON / lineage).
	 */
	public static String shortHash(NonCubicBilinearAlgorithm a) {
		return contentHash(a).substring(0, 4);
	}

	/** The stamped {@code "hash"} field of a scheme JSON, or {@code null} if absent. */
	public static String readHash(tools.jackson.databind.JsonNode root) {
		tools.jackson.databind.JsonNode h = root.get("hash");
		return (h != null && h.isTextual()) ? h.asText() : null;
	}

	private static void appendFactor(StringBuilder sb, char tag, double[][] m, int r, boolean exact) {
		int rows = m.length;
		for (int c = 0; c < r; c++) {
			for (int row = 0; row < rows; row++) {
				double v = m[row][c];
				if (v == 0.0) {
					continue;
				}
				if (exact) {
					// Exact symbolic token (matches the on-disk "1/4" string).
					sb.append(tag).append(row).append(':').append(c).append(':').append(canonicalToken(v)).append(';');
				} else {
					long sv = Math.round(v * 1_000_000.0);
					if (sv == 0) {
						continue;
					}
					sb.append(tag).append(row).append(':').append(c).append(':').append(sv).append(';');
				}
			}
		}
	}

	private static String sha256Hex(String s) {
		try {
			java.security.MessageDigest md = java.security.MessageDigest.getInstance("SHA-256");
			byte[] d = md.digest(s.getBytes(java.nio.charset.StandardCharsets.UTF_8));
			StringBuilder hx = new StringBuilder(64);
			for (byte b : d) {
				hx.append(Character.forDigit((b >> 4) & 0xF, 16));
				hx.append(Character.forDigit(b & 0xF, 16));
			}
			return hx.toString();
		} catch (java.security.NoSuchAlgorithmException e) {
			throw new IllegalStateException("SHA-256 unavailable", e);
		}
	}

	public static void writeStub(NonCubicBilinearAlgorithm alg, File f, Lineage.Node lineage)
			throws IOException {
		writeStub(alg, f, lineage, null);
	}

	/**
	 * Write a lineage-only stub, optionally born-stamped with {@code fields[]}.
	 * A matrix-less stub cannot be content-stamped later ({@code
	 * BackfillMissingFields} skips it), and without {@code fields[]} the lookup
	 * treats it as ABSENT — so the materialiser passes the lineage-inferred field
	 * set here ({@link FieldAwareLookup#fieldNamesFromLineage}) to avoid writing a
	 * silently-invisible stub. {@code fields == null/empty} omits the property
	 * (back-compat path).
	 */
	public static void writeStub(NonCubicBilinearAlgorithm alg, File f, Lineage.Node lineage,
			java.util.List<String> fields) throws IOException {
		if (lineage == null) {
			throw new IllegalArgumentException("writeStub requires a lineage (else irreproducible)");
		}
		String fieldsLine = (fields == null || fields.isEmpty()) ? ""
				: "  \"fields\": [" + fields.stream().map(s -> "\"" + s + "\"")
						.collect(java.util.stream.Collectors.joining(", ")) + "],\n";
		// Stamp the content hash of the (built) scheme even though we don't store
		// its matrices — so a stub is still content-addressable and other lineages
		// can reference it precisely by hash.
		String base = "{\n"
				+ "  \"n\": [" + alg.n + ", " + alg.m + ", " + alg.p + "],\n"
				+ "  \"m\": " + alg.r + ",\n"
				+ "  \"hash\": \"" + contentHash(alg) + "\",\n"
				+ fieldsLine
				+ "  \"scheme_type\": \"stub\"\n"
				+ "}";
		String json = injectLineage(base, lineage);
		// Canonical formatting via the single shared formatter (never raw) so a
		// stub is written in the same style every other writer produces.
		try (PrintWriter pw = new PrintWriter(new BufferedWriter(new FileWriter(f)))) {
			pw.print(MatrixJsonFormatter.format(json));
		}
		FieldAwareLookup.onSchemeWritten(f);
	}

	/**
	 * Write a scheme with an optional {@link Lineage} expression
	 * injected as {@code lineage_str} (pretty form) + {@code lineage}
	 * (JSON DAG) header fields. The body (n, m, u, v, w) is
	 * unchanged. Pass {@code lineage = null} to omit.
	 */
	public static void write(NonCubicBilinearAlgorithm alg, File f,
			Lineage.Node lineage) throws IOException {
		int maxDim = Math.max(alg.n, Math.max(alg.m, alg.p));
		String json = (maxDim >= SPARSE_DIM_THRESHOLD)
				? toJsonSparse(alg)
				: toJson(alg, null, null);
		if (lineage != null) {
			json = injectLineage(json, lineage);
		}
		// Canonical formatting via the single shared formatter (never raw), so
		// the dense `toJson` path and the lineage-injected header agree with the
		// style every other writer (and ReformatSchemes) produces.
		try (PrintWriter pw = new PrintWriter(new BufferedWriter(new FileWriter(f)))) {
			pw.print(MatrixJsonFormatter.format(json));
		}
		// Surgical cache notification — keeps any FieldAwareLookup index that
		// was built before this write in sync. No-op for roots that haven't
		// been indexed yet.
		FieldAwareLookup.onSchemeWritten(f);
	}

	/**
	 * Insert {@code "lineage_str": ..., "lineage": ...} after the
	 * opening brace of an existing scheme JSON. Idempotent if called
	 * twice (the second call replaces only what it parses, but in
	 * practice {@link #write(NonCubicBilinearAlgorithm, File,
	 * Lineage.Node)} is the only caller and always passes raw JSON).
	 */
	private static String injectLineage(String json, Lineage.Node lineage) {
		String prettyStr = Lineage.prettyString(lineage);
		String compactStr = Lineage.prettyCompact(lineage);
		String lineageJson = Lineage.toJson(lineage);
		int openBrace = json.indexOf('{');
		if (openBrace < 0) return json;  // malformed
		String prefix = json.substring(0, openBrace + 1);
		String suffix = json.substring(openBrace + 1);
		String escapedPretty = prettyStr.replace("\\", "\\\\").replace("\"", "\\\"");
		String escapedCompact = compactStr.replace("\\", "\\\\").replace("\"", "\\\"");
		boolean dense = suffix.startsWith("\n");
		String sep = dense ? "\n  " : "";
		String comma = dense ? ",\n" : ",";
		String indent = dense ? "  " : "";
		return prefix
				+ sep + "\"lineage_str\": \"" + escapedPretty + "\"" + comma
				+ indent + "\"lineage_compact\": \"" + escapedCompact + "\"" + comma
				+ indent + "\"lineage\": " + lineageJson + comma
				+ suffix.substring(dense ? 1 : 0);
	}

	public static void write(BilinearAlgorithm alg, File f) throws IOException {
		write(NonCubicBilinearAlgorithm.fromCubic(alg), f);
	}

	/**
	 * Add top-level metadata fields to an <em>existing</em> scheme file, writing
	 * through the canonical {@link MatrixJsonFormatter} (never a hand-built
	 * string — see the contract on {@link MatrixJsonFormatter#write}). Existing
	 * keys are left untouched (no clobber). The whole file is re-emitted in the
	 * one canonical style, so a file whose formatting had drifted (e.g. a
	 * previously hand-appended field) is normalised in the same pass — passing an
	 * empty map performs a pure normalise.
	 *
	 * <p>Returns whether the on-disk bytes <em>would</em> change; writes (and
	 * notifies {@link FieldAwareLookup}) only when {@code apply} is true. This is
	 * THE entry point for stamping derived metrics ({@code verified, additions,
	 * has_buds, buds, projection_margin}) onto catalog files.</p>
	 */
	public static boolean addFields(File f, Map<String, Object> fields, boolean apply)
			throws IOException {
		return updateFields(f, fields, java.util.List.of(), apply);
	}

	/**
	 * Generalised {@link #addFields}: add {@code fields} (no clobber) AND remove
	 * {@code remove} keys, in one canonical re-emit. Used for self-healing flags
	 * (e.g. clear {@code corrupted} once a stub replays again). Returns whether the
	 * on-disk bytes would change; writes only when {@code apply} is true.
	 */
	public static boolean updateFields(File f, Map<String, Object> fields,
			java.util.Collection<String> remove, boolean apply) throws IOException {
		String existing = java.nio.file.Files.readString(f.toPath());
		ObjectNode root = (ObjectNode) MAPPER.readTree(existing);
		return updateFields(f, root, existing, fields, remove, apply);
	}

	/**
	 * Same as {@link #updateFields(File, Map, java.util.Collection, boolean)}, but
	 * reusing an <em>already-parsed</em> {@code root} and its original on-disk text
	 * {@code existing} — so a caller that has just read the file (e.g.
	 * {@code EnrichSchemeMetrics}) avoids a redundant re-read and re-parse. The
	 * byte-level {@code existing} is still required for the formatting-drift diff.
	 *
	 * <p><b>Mutates {@code root}</b> (applies {@code remove} / {@code fields}); the
	 * caller must not reuse it afterwards. {@code existing} must be the verbatim
	 * current on-disk content of {@code f} (no write may have intervened).</p>
	 */
	public static boolean updateFields(File f, ObjectNode root, String existing,
			Map<String, Object> fields, java.util.Collection<String> remove, boolean apply)
			throws IOException {
		for (String k : remove) {
			root.remove(k);
		}
		for (Map.Entry<String, Object> e : fields.entrySet()) {
			if (root.has(e.getKey())) {
				continue;  // never clobber an existing value
			}
			putValue(root, e.getKey(), e.getValue());
		}
		String formatted = MatrixJsonFormatter.format(root);
		boolean changed = !formatted.equals(existing);
		if (changed && apply) {
			java.nio.file.Files.writeString(f.toPath(), formatted);
			FieldAwareLookup.onSchemeWritten(f);
		}
		return changed;
	}

	/** {@link #addFields(File, Map, boolean)} reusing an already-parsed node + text
	 *  (no re-read). Mutates {@code root}; see the note on the {@code updateFields}
	 *  node overload. */
	public static boolean addFields(File f, ObjectNode root, String existing,
			Map<String, Object> fields, boolean apply) throws IOException {
		return updateFields(f, root, existing, fields, java.util.List.of(), apply);
	}

	/** Type-dispatching {@link ObjectNode} setter for {@link #addFields}. A
	 *  {@link java.util.Collection} becomes a JSON ARRAY (element-wise typed) —
	 *  it used to be silently stringified ({@code "[F2, F3]"}), which made a
	 *  stamped {@code fields} list invisible to the content-driven index. */
	private static void putValue(ObjectNode o, String k, Object v) {
		if (v instanceof Boolean b) {
			o.put(k, b);
		} else if (v instanceof Integer i) {
			o.put(k, i);
		} else if (v instanceof Long l) {
			o.put(k, l);
		} else if (v instanceof Double d) {
			o.put(k, d);
		} else if (v instanceof java.util.Collection<?> c) {
			ArrayNode arr = o.putArray(k);
			for (Object e : c) {
				if (e instanceof Boolean b) {
					arr.add(b);
				} else if (e instanceof Integer i) {
					arr.add(i);
				} else if (e instanceof Long l) {
					arr.add(l);
				} else if (e instanceof Double d) {
					arr.add(d);
				} else {
					arr.add(String.valueOf(e));
				}
			}
		} else {
			o.put(k, String.valueOf(v));
		}
	}

	/**
	 * Write a complex-valued scheme using the dense {@code [re, im]} pair
	 * encoding the AlphaEvolve corpus uses. W is converted to dronperminov's
	 * col-major C-flatten on the way out.
	 *
	 * <p>For large complex schemes the file size grows quickly — {@code ⟨16,16,16⟩}
	 * dense complex is ≈ 28 MB. A sparse-complex format (object of
	 * {@code "pos": [re, im]} per multiplication) is a TODO (see roadmap).</p>
	 */
	public static void write(ComplexNonCubicBilinearAlgorithm alg, File f, String fieldLabel)
			throws IOException {
		String json = toJsonComplex(alg, fieldLabel);
		try (PrintWriter pw = new PrintWriter(new BufferedWriter(new FileWriter(f)))) {
			pw.print(json);
		}
	}

	public static String toJsonComplex(ComplexNonCubicBilinearAlgorithm alg, String fieldLabel) {
		StringBuilder sb = new StringBuilder();
		sb.append('{');
		sb.append("\"n\":[").append(alg.n).append(',').append(alg.m).append(',').append(alg.p).append("],");
		sb.append("\"m\":").append(alg.r).append(',');
		sb.append("\"field\":\"").append(fieldLabel == null ? "C" : fieldLabel).append("\",");
		sb.append("\"complex\":true,");
		sb.append("\"u\":").append(writeComplexFactor(alg.uRe, alg.uIm, alg.r, alg.dimU(), false, alg.n, alg.p)).append(',');
		sb.append("\"v\":").append(writeComplexFactor(alg.vRe, alg.vIm, alg.r, alg.dimV(), false, alg.n, alg.p)).append(',');
		sb.append("\"w\":").append(writeComplexFactor(alg.wRe, alg.wIm, alg.r, alg.dimW(), true, alg.n, alg.p));
		sb.append("}");
		try {
			return MatrixJsonFormatter.format(sb.toString());
		} catch (IOException e) {
			throw new IllegalStateException("self-emitted complex JSON failed to parse", e);
		}
	}

	/**
	 * Build the dense complex factor as a {@code [rank][dim][2]} nested array
	 * of {@code [re, im]} pairs. W ({@code isW = true}) is re-indexed from our
	 * internal row-major C-flatten to dronperminov's col-major.
	 */
	private static String writeComplexFactor(double[][] re, double[][] im, int rank, int dim,
			boolean isW, int n1, int n3) {
		StringBuilder sb = new StringBuilder();
		sb.append('[');
		for (int k = 0; k < rank; k++) {
			if (k > 0) sb.append(',');
			sb.append('[');
			for (int p = 0; p < dim; p++) {
				int srcPos;
				if (isW) {
					int i = p / n3;
					int j = p % n3;
					srcPos = i * n3 + j;
					// We need the entry at dronperminov-position p (col-major
					// j*n1+i). Convert that back to our row-major (i*n3+j).
					int dronI = p % n1;
					int dronJ = p / n1;
					srcPos = dronI * n3 + dronJ;
				} else {
					srcPos = p;
				}
				if (p > 0) sb.append(',');
				sb.append('[').append(formatCoef(re[srcPos][k])).append(',')
						.append(formatCoef(im[srcPos][k])).append(']');
			}
			sb.append(']');
		}
		sb.append(']');
		return sb.toString();
	}

	/**
	 * Serialize to a JSON string. {@code multiplications} and {@code elements}
	 * are optional human-readable companions — pass {@code null} to omit them.
	 */
	public static String toJson(NonCubicBilinearAlgorithm alg,
			List<String> multiplications, List<String> elements) {
		int dimU = alg.dimU(), dimV = alg.dimV(), dimW = alg.dimW();
		StringBuilder sb = new StringBuilder();
		sb.append("{\n");
		sb.append("  \"n\": [").append(alg.n).append(", ").append(alg.m)
				.append(", ").append(alg.p).append("],\n");
		sb.append("  \"m\": ").append(alg.r).append(",\n");
		if (multiplications != null) {
			sb.append("  \"multiplications\": ").append(jsonStringArray(multiplications)).append(",\n");
		}
		if (elements != null) {
			sb.append("  \"elements\": ").append(jsonStringArray(elements)).append(",\n");
		}
		sb.append("  \"u\": ").append(writeTransposed(alg.denseU(), alg.r, dimU)).append(",\n");
		sb.append("  \"v\": ").append(writeTransposed(alg.denseV(), alg.r, dimV)).append(",\n");
		sb.append("  \"w\": ").append(writeWColMajor(alg.denseW(), alg.r, alg.n, alg.p)).append("\n");
		sb.append("}\n");
		return sb.toString();
	}

	/**
	 * Sparse-format JSON output (our extension for max-dim ≥ 9). Each
	 * {@code *_sparse} factor is a rank-length array of objects mapping
	 * non-zero position-index → coefficient.
	 *
	 * <p>Output is intentionally minified — pass through {@link MatrixJsonFormatter}
	 * (e.g. via {@code ReformatSchemes}) for matrix-friendly indentation.</p>
	 */
	public static String toJsonSparse(NonCubicBilinearAlgorithm alg) {
		StringBuilder sb = new StringBuilder();
		sb.append('{');
		sb.append("\"n\":[").append(alg.n).append(',').append(alg.m).append(',').append(alg.p).append("],");
		sb.append("\"m\":").append(alg.r).append(',');
		sb.append("\"u_sparse\":").append(writeSparseFactor(alg.denseU(), alg.r, alg.dimU(), false, alg.n, alg.p)).append(',');
		sb.append("\"v_sparse\":").append(writeSparseFactor(alg.denseV(), alg.r, alg.dimV(), false, alg.n, alg.p)).append(',');
		sb.append("\"w_sparse\":").append(writeSparseFactor(alg.denseW(), alg.r, alg.dimW(), true, alg.n, alg.p));
		sb.append("}");
		try {
			return MatrixJsonFormatter.format(sb.toString());
		} catch (IOException e) {
			throw new IllegalStateException("self-emitted JSON failed to parse", e);
		}
	}

	/**
	 * Build the sparse {@code u_sparse}/{@code v_sparse}/{@code w_sparse}
	 * array. For W ({@code isW = true}) the position index is converted from
	 * our internal row-major C-flatten ({@code i·n_3 + j}) to dronperminov's
	 * col-major ({@code j·n_1 + i}) so the stored indices match the dense
	 * format's convention.
	 */
	private static String writeSparseFactor(double[][] src, int rank, int dim, boolean isW,
			int n1, int n3) {
		// Row-oriented format (new, since 2026-06-03):
		//   {"0": {"i": [0, 4, 8], "c": [1, -1, 1]}, "1": {...}, ...}
		// One product per outer-map entry, indices and coefficients held as
		// parallel arrays for compact one-line-per-product display.
		StringBuilder sb = new StringBuilder();
		sb.append('{');
		// Per product k, gather non-zero (pos, coef) pairs in pos-order.
		StringBuilder idx = new StringBuilder();
		StringBuilder coef = new StringBuilder();
		for (int k = 0; k < rank; k++) {
			if (k > 0) sb.append(',');
			sb.append('"').append(k).append("\":{\"i\":[");
			idx.setLength(0);
			coef.setLength(0);
			boolean first = true;
			for (int p = 0; p < dim; p++) {
				double v = src[p][k];
				if (v == 0.0) continue;
				int outPos;
				if (isW) {
					int i = p / n3;
					int j = p % n3;
					outPos = j * n1 + i;
				} else {
					outPos = p;
				}
				if (!first) {
					idx.append(',');
					coef.append(',');
				}
				first = false;
				idx.append(outPos);
				coef.append(formatCoef(v));
			}
			sb.append(idx).append("],\"c\":[").append(coef).append("]}");
		}
		sb.append('}');
		return sb.toString();
	}

	private static String writeTransposed(double[][] src, int rank, int dim) {
		// src is [dim][rank]; output rank rows of dim ints each.
		StringBuilder sb = new StringBuilder();
		sb.append("[\n");
		for (int k = 0; k < rank; k++) {
			sb.append("    [");
			for (int j = 0; j < dim; j++) {
				if (j > 0) sb.append(", ");
				sb.append(formatCoef(src[j][k]));
			}
			sb.append(']');
			if (k < rank - 1) sb.append(',');
			sb.append('\n');
		}
		sb.append("  ]");
		return sb.toString();
	}

	/**
	 * Write W in dronperminov's col-major C-flatten convention. Reverse of
	 * {@link #readWColMajor}.
	 */
	private static String writeWColMajor(double[][] src, int rank, int n1, int n3) {
		int dim = n1 * n3;
		StringBuilder sb = new StringBuilder();
		sb.append("[\n");
		for (int k = 0; k < rank; k++) {
			sb.append("    [");
			for (int pDron = 0; pDron < dim; pDron++) {
				int i = pDron % n1;
				int j = pDron / n1;
				int pUs = i * n3 + j;
				if (pDron > 0) sb.append(", ");
				sb.append(formatCoef(src[pUs][k]));
			}
			sb.append(']');
			if (k < rank - 1) sb.append(',');
			sb.append('\n');
		}
		sb.append("  ]");
		return sb.toString();
	}

	/**
	 * Emit a coefficient as an <strong>exact</strong> JSON token — never a rounded
	 * decimal. Integers serialise as bare numbers ({@code 3}); non-integer rationals
	 * as a quoted fraction string ({@code "1/17"}), recovered from the {@code double}
	 * by continued-fraction rationalisation (the engine is {@code double}, but a
	 * scheme born of e.g. {@code 1/(n/2+1)} stores the exact ratio, not
	 * {@code 0.058823529411764705}). The reader round-trips this via
	 * {@link #parseCoefString}. A value that is not a small rational within
	 * {@link #RATIONALIZE_DEN_CAP} (a genuine irrational, or one needing a huge
	 * denominator) falls back to {@link Double#toString} — still lossless for the
	 * {@code double}, just not symbolic. The string form is forward-compatible with
	 * symbolic coefficients (e.g. {@code "sqrt(2)"}) once a constructor emits them.
	 */
	private static String formatCoef(double v) {
		long rounded = Math.round(v);
		if (Math.abs(v - rounded) < INTEGER_TOL) return Long.toString(rounded);
		String frac = rationalize(v);
		if (frac != null) return '"' + frac + '"';
		return Double.toString(v);  // irrational / huge-denominator: lossless decimal
	}

	private static final double INTEGER_TOL = 1e-9;
	/** Largest denominator we will recover; catalog rationals are simple (≤ a few
	 *  hundred), this cap just bounds the search and rejects spurious matches. */
	private static final long RATIONALIZE_DEN_CAP = 1_000_000L;
	/** A convergent must reproduce the double this tightly (relative) to be accepted
	 *  as its exact value — loose enough for the double round-trip, tight enough that
	 *  a true irrational never spuriously matches a small rational. */
	private static final double RATIONALIZE_TOL = 1e-12;

	/**
	 * Recover the exact {@code "p/q"} a {@code double} represents via continued-fraction
	 * convergents (denominator-bounded). Returns {@code null} when no rational within
	 * {@link #RATIONALIZE_DEN_CAP} reproduces {@code v} to {@link #RATIONALIZE_TOL} —
	 * the caller then keeps the lossless decimal. Sign is carried on the numerator.
	 */
	private static String rationalize(double v) {
		boolean neg = v < 0;
		double x = Math.abs(v);
		long hPrev = 0, h = 1;   // numerator convergents (h_{-1}=0, h_0 forming)
		long kPrev = 1, k = 0;   // denominator convergents
		double frac = x;
		for (int iter = 0; iter < 64; iter++) {
			long a = (long) Math.floor(frac);
			long hNext = a * h + hPrev;
			long kNext = a * k + kPrev;
			if (kNext > RATIONALIZE_DEN_CAP || kNext <= 0) {
				break;  // denominator blew the cap (or overflowed) — give up
			}
			hPrev = h; h = hNext;
			kPrev = k; k = kNext;
			if (k > 0 && Math.abs((double) h / (double) k - x) <= RATIONALIZE_TOL * Math.max(1.0, x)) {
				if (k == 1) {
					return null;  // integer — caller handles via the integer branch
				}
				return (neg ? "-" : "") + h + "/" + k;
			}
			double rem = frac - a;
			if (rem <= 1e-15) {
				break;  // terminated exactly but k==1 path already handled above
			}
			frac = 1.0 / rem;
		}
		return null;
	}

	/**
	 * Rewrite a scheme file's coefficients <strong>in place</strong> from rounded
	 * decimals to exact rational strings (e.g. {@code 0.058823529411764705} →
	 * {@code "1/17"}), preserving every other field (lineage, fields, source,
	 * additions, buds, …) byte-for-byte modulo the touched values. Walks only the
	 * factor blocks ({@code u/v/w} and {@code u_sparse/v_sparse/w_sparse}, dense or
	 * sparse, real or complex) and replaces each <em>floating-point</em> node with
	 * its exact rational (integers and index arrays are integral nodes, untouched).
	 *
	 * <p>Unlike {@link #write}, which re-serialises from a {@link
	 * NonCubicBilinearAlgorithm} and drops curated metadata, this is a surgical
	 * value-level rewrite. Idempotent: a file already exact reports no change.
	 *
	 * @return true if any coefficient was rewritten (and, when {@code apply}, the
	 *         file was overwritten); false if nothing changed.
	 */
	public static boolean exactifyCoefficients(File f, boolean apply) throws IOException {
		String existing = java.nio.file.Files.readString(f.toPath());
		ObjectNode root = (ObjectNode) MAPPER.readTree(existing);
		boolean[] changed = { false };
		for (String key : new String[] { "u", "v", "w", "u_sparse", "v_sparse", "w_sparse" }) {
			JsonNode sub = root.get(key);
			if (sub != null) {
				exactifyTree(sub, changed);
			}
		}
		if (!changed[0]) {
			return false;
		}
		String formatted = MatrixJsonFormatter.format(root);
		if (formatted.equals(existing)) {
			return false;
		}
		if (apply) {
			java.nio.file.Files.writeString(f.toPath(), formatted);
			FieldAwareLookup.onSchemeWritten(f);
		}
		return true;
	}

	/** Recursively replace floating-point coefficient nodes with exact rationals. */
	private static void exactifyTree(JsonNode node, boolean[] changed) {
		if (node instanceof ArrayNode arr) {
			for (int i = 0; i < arr.size(); i++) {
				JsonNode c = arr.get(i);
				JsonNode repl = exactifyLeaf(c, changed);
				if (repl != null) {
					arr.set(i, repl);
				} else {
					exactifyTree(c, changed);
				}
			}
		} else if (node instanceof ObjectNode obj) {
			// Collect keys first — we mutate the node while iterating.
			java.util.List<String> keys = obj.properties().stream()
					.map(java.util.Map.Entry::getKey).toList();
			for (String k : keys) {
				JsonNode c = obj.get(k);
				JsonNode repl = exactifyLeaf(c, changed);
				if (repl != null) {
					obj.set(k, repl);
				} else {
					exactifyTree(c, changed);
				}
			}
		}
	}

	/**
	 * If {@code c} is a non-integer floating-point coefficient, return its exact
	 * replacement node ({@code "p/q"} text, or an integer node when the float is a
	 * whole number); otherwise {@code null} (recurse / leave as-is — integers, index
	 * arrays, and irrationals untouched).
	 */
	private static JsonNode exactifyLeaf(JsonNode c, boolean[] changed) {
		if (!c.isFloatingPointNumber()) {
			return null;
		}
		double v = c.doubleValue();
		String frac = rationalize(v);
		if (frac != null) {
			changed[0] = true;
			return JsonNodeFactory.instance.textNode(frac);
		}
		long rounded = Math.round(v);
		if (Math.abs(v - rounded) < INTEGER_TOL) {
			changed[0] = true;
			return JsonNodeFactory.instance.numberNode(rounded);
		}
		return null;  // genuine irrational / huge denominator — keep the lossless double
	}

	/**
	 * Enforce the {@link #SPARSE_DIM_THRESHOLD} on an existing file: if it is a
	 * <strong>dense</strong> bilinear scheme with {@code maxDim ≥ threshold}, rewrite
	 * its {@code u/v/w} arrays into the {@code u_sparse/v_sparse/w_sparse} format
	 * (the same writers {@link #write} uses, so coefficients come out exact per
	 * {@link #formatCoef}), preserving every other field. No-op for stubs, complex
	 * schemes, already-sparse files, or dense files below the threshold.
	 *
	 * @return true if the file was (or would be) converted.
	 */
	public static boolean convertDenseToSparse(File f, boolean apply) throws IOException {
		String existing = java.nio.file.Files.readString(f.toPath());
		ObjectNode root = (ObjectNode) MAPPER.readTree(existing);
		if (isStub(root) || (root.has("complex") && root.get("complex").asBoolean())) {
			return false;
		}
		if (root.has("u_sparse") || !root.has("u")) {
			return false;  // already sparse, or not a dense bilinear scheme
		}
		JsonNode nArr = root.get("n");
		if (nArr == null || !nArr.isArray() || nArr.size() != 3) {
			return false;
		}
		int maxDim = Math.max(nArr.get(0).asInt(), Math.max(nArr.get(1).asInt(), nArr.get(2).asInt()));
		if (maxDim < SPARSE_DIM_THRESHOLD) {
			return false;  // correctly dense
		}
		NonCubicBilinearAlgorithm alg = readBilinear(f);
		root.remove("u");
		root.remove("v");
		root.remove("w");
		root.set("u_sparse", MAPPER.readTree(
				writeSparseFactor(alg.denseU(), alg.r, alg.dimU(), false, alg.n, alg.p)));
		root.set("v_sparse", MAPPER.readTree(
				writeSparseFactor(alg.denseV(), alg.r, alg.dimV(), false, alg.n, alg.p)));
		root.set("w_sparse", MAPPER.readTree(
				writeSparseFactor(alg.denseW(), alg.r, alg.dimW(), true, alg.n, alg.p)));
		String formatted = MatrixJsonFormatter.format(root);
		if (formatted.equals(existing)) {
			return false;
		}
		if (apply) {
			java.nio.file.Files.writeString(f.toPath(), formatted);
			FieldAwareLookup.onSchemeWritten(f);
		}
		return true;
	}

	private static String jsonStringArray(List<String> xs) {
		StringBuilder sb = new StringBuilder("[");
		for (int i = 0; i < xs.size(); i++) {
			if (i > 0) sb.append(", ");
			sb.append('"').append(xs.get(i).replace("\\", "\\\\")
					.replace("\"", "\\\"")).append('"');
		}
		sb.append("]");
		return sb.toString();
	}
}
