package eu.solven.matmul.catalog;

import java.io.File;
import java.io.IOException;
import java.nio.file.Files;
import java.nio.file.Path;
import java.util.Locale;
import java.util.Map;

import tools.jackson.core.JacksonException;
import tools.jackson.databind.JsonNode;
import tools.jackson.databind.ObjectMapper;
import tools.jackson.databind.json.JsonMapper;

/**
 * Pretty-prints JSON with a matrix-friendly rule:
 *
 * <ul>
 *   <li><b>All-numeric arrays</b> (numbers / booleans / nulls) → <i>inline</i>
 *       as {@code [a, b, c, d]}. Keeps each U/V/W matrix row, and index lists
 *       like {@code forms}, on one line.</li>
 *   <li><b>SLP op-tuples</b> — one leading string opcode then only numbers,
 *       e.g. {@code ["+", 4, 0, 1]} → <i>inline</i>, so straight-line programs
 *       stay readable.</li>
 *   <li><b>Pure-string arrays</b> ({@code fields}, {@code multiplications}, …)
 *       and <b>arrays holding any nested array/object</b> (the matrices
 *       themselves, the {@code ops} list) → <i>vertical</i>, one element per
 *       line. See {@link #shouldInline}.</li>
 *   <li><b>Objects with only primitive values</b> (e.g. dronperminov's
 *       {@code {"index": 3, "value": 1}}) → inlined as a single-line object.</li>
 *   <li><b>Other objects</b> → vertical, one {@code "key": value} per line.</li>
 * </ul>
 *
 * <p>Indent width is 2 spaces. File ends with a single newline. Idempotent.</p>
 *
 * <p>Backed by Jackson for parsing — {@link ObjectMapper#readTree(String)}
 * produces a {@link JsonNode} tree we then walk recursively.</p>
 */
public final class MatrixJsonFormatter {

	private static final int INDENT_STEP = 2;
	private static final ObjectMapper MAPPER = JsonMapper.builder().build();

	private MatrixJsonFormatter() {}

	/**
	 * THE canonical write entry point. <b>Every procedure that writes a scheme
	 * JSON to disk — creating it from scratch or editing an existing one — MUST
	 * go through {@code write(...)} or {@link #format} here</b>, never through
	 * {@code ObjectNode.toPrettyString()} / {@code writerWithDefaultPrettyPrinter}
	 * / hand-built strings. That is the only way on-disk files stay in one
	 * canonical style and no writer can silently re-introduce a divergent one
	 * (user 2026-06-05). Writes UTF-8 with a single trailing newline.
	 */
	public static void write(File f, JsonNode node) throws IOException {
		Files.writeString(f.toPath(), format(node));
	}

	/** {@link #write(File, JsonNode)} for a {@link Path}. */
	public static void write(Path p, JsonNode node) throws IOException {
		Files.writeString(p, format(node));
	}

	/** Canonical-format an already-parsed JSON tree. */
	public static String format(JsonNode node) {
		StringBuilder sb = new StringBuilder();
		writeNode(node, sb, 0);
		sb.append('\n');
		return sb.toString();
	}

	/** Canonical-format a JSON string (parse, then {@link #format(JsonNode)}). */
	public static String format(String json) throws IOException {
		try {
			return format(MAPPER.readTree(json));
		} catch (JacksonException e) {
			throw new IOException("parse error: " + e.getMessage(), e);
		}
	}

	private static void writeNode(JsonNode node, StringBuilder out, int indent) {
		if (node == null || node.isNull()) {
			out.append("null");
		} else if (node.isBoolean()) {
			out.append(node.booleanValue());
		} else if (node.isNumber()) {
			out.append(formatNumber(node));
		} else if (node.isTextual()) {
			out.append('"').append(escape(node.textValue())).append('"');
		} else if (node.isArray()) {
			writeArray(node, out, indent);
		} else if (node.isObject()) {
			writeObject(node, out, indent);
		} else {
			throw new IllegalArgumentException("unsupported node kind: " + node.getNodeType());
		}
	}

	private static void writeArray(JsonNode arr, StringBuilder out, int indent) {
		if (arr.isEmpty()) {
			out.append("[]");
			return;
		}
		if (shouldInline(arr)) {
			out.append('[');
			for (int i = 0; i < arr.size(); i++) {
				if (i > 0) out.append(", ");
				writeNode(arr.get(i), out, 0);
			}
			out.append(']');
			return;
		}
		out.append("[\n");
		String childIndent = " ".repeat(indent + INDENT_STEP);
		for (int i = 0; i < arr.size(); i++) {
			out.append(childIndent);
			writeNode(arr.get(i), out, indent + INDENT_STEP);
			if (i < arr.size() - 1) out.append(',');
			out.append('\n');
		}
		out.append(" ".repeat(indent)).append(']');
	}

	private static void writeObject(JsonNode obj, StringBuilder out, int indent) {
		if (obj.isEmpty()) {
			out.append("{}");
			return;
		}
		if (allPrimitiveValues(obj)) {
			out.append('{');
			boolean first = true;
			for (Map.Entry<String, JsonNode> e : obj.properties()) {
				if (!first) out.append(", ");
				first = false;
				out.append('"').append(escape(e.getKey())).append("\": ");
				writeNode(e.getValue(), out, 0);
			}
			out.append('}');
			return;
		}
		out.append("{\n");
		String childIndent = " ".repeat(indent + INDENT_STEP);
		int i = 0, last = obj.size() - 1;
		for (Map.Entry<String, JsonNode> e : obj.properties()) {
			out.append(childIndent)
					.append('"').append(escape(e.getKey())).append("\": ");
			writeNode(e.getValue(), out, indent + INDENT_STEP);
			if (i++ < last) out.append(',');
			out.append('\n');
		}
		out.append(" ".repeat(indent)).append('}');
	}

	/**
	 * Decide whether an array is emitted inline (one row) or vertically.
	 * Inline in exactly two cases — both of which are short, scalar, and far
	 * more readable on one line:
	 *
	 * <ol>
	 *   <li><b>All-numeric</b> (numbers / booleans / nulls, no strings) — a
	 *       matrix row {@code [0, 1, -1, 0]} or an index list like
	 *       {@code "forms": [6, 2, 3]} (user 2026-06-03: "true only if we spot
	 *       only numbers").</li>
	 *   <li><b>SLP op-tuple</b>: a single leading string opcode followed by at
	 *       least one number, and nothing else — e.g. {@code ["+", 4, 0, 1]}
	 *       or {@code ["*", 6, 5, -1]} (user 2026-06-05: "if array is one
	 *       string then at least one number and only numbers: put in 1 row,
	 *       hence SLP will be easier to read").</li>
	 *   <li><b>Coefficient row</b>: every element is a number OR a numeric/
	 *       fraction <i>string</i> ({@code "-1/17"}, {@code "2/3"}) — e.g. a
	 *       sparse {@code "c"} list now holding exact rationals, or a U/V/W row.
	 *       Without this, exact-rational coefficient arrays would split one value
	 *       per line (user 2026-06: keep coefficient rows on one line).</li>
	 * </ol>
	 *
	 * <p>Everything else is vertical: any array holding a nested array/object
	 * (e.g. the U/V/W factor matrices, or the {@code ops} list of tuples), and
	 * pure-string arrays like {@code fields} / {@code multiplications} /
	 * {@code elements} (whose strings are NOT numeric, so they stay vertical).</p>
	 */
	private static boolean shouldInline(JsonNode arr) {
		boolean sawNumber = false;
		boolean firstIsText = false;
		boolean allCoeffish = true;  // every element a number / bool / null / numeric-string
		int textualCount = 0;
		int i = 0;
		for (JsonNode v : arr) {
			if (v.isArray() || v.isObject()) return false;  // compound → vertical
			if (v.isTextual()) {
				textualCount++;
				if (i == 0) firstIsText = true;
				if (!isNumericString(v.asString())) allCoeffish = false;
			} else if (v.isNumber()) {
				sawNumber = true;
			}
			i++;
		}
		// (1) all-numeric (no strings).
		if (textualCount == 0) return true;
		// (2) coefficient row: numbers and/or numeric/fraction strings only.
		if (allCoeffish) return true;
		// (3) op-tuple: exactly one leading string opcode + at least one number.
		return textualCount == 1 && firstIsText && sawNumber;
	}

	/** Integer, decimal, scientific, or {@code p/q} fraction — a coefficient token. */
	private static final java.util.regex.Pattern NUMERIC_STRING =
			java.util.regex.Pattern.compile("[+-]?(\\d+(\\.\\d+)?([eE][+-]?\\d+)?|\\d+/-?\\d+)");

	private static boolean isNumericString(String s) {
		return NUMERIC_STRING.matcher(s).matches();
	}

	private static boolean allPrimitiveValues(JsonNode obj) {
		for (JsonNode v : obj.values()) {
			if (v.isArray() || v.isObject()) return false;
		}
		return true;
	}

	private static String formatNumber(JsonNode n) {
		if (n.isIntegralNumber()) {
			return n.asText();
		}
		// Avoid spurious "1.0" / "0.0"; prefer integral when value is integral.
		double d = n.doubleValue();
		if (d == Math.rint(d) && !Double.isInfinite(d) && Math.abs(d) < 1e15) {
			return Long.toString((long) d);
		}
		return String.format(Locale.ROOT, "%s", n.asText());
	}

	private static String escape(String s) {
		StringBuilder sb = new StringBuilder(s.length() + 2);
		for (int i = 0; i < s.length(); i++) {
			char c = s.charAt(i);
			switch (c) {
			case '"': sb.append("\\\""); break;
			case '\\': sb.append("\\\\"); break;
			case '\n': sb.append("\\n"); break;
			case '\r': sb.append("\\r"); break;
			case '\t': sb.append("\\t"); break;
			default:
				if (c < 0x20) {
					sb.append(String.format("\\u%04x", (int) c));
				} else {
					sb.append(c);
				}
			}
		}
		return sb.toString();
	}
}
