package eu.solven.matmul.docs.generate;

import lombok.extern.slf4j.Slf4j;

import java.io.IOException;
import java.io.PrintWriter;
import java.nio.file.Files;
import java.nio.file.Path;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.LinkedHashMap;
import java.util.List;
import java.util.Locale;
import java.util.Map;
import java.util.TreeMap;

import tools.jackson.databind.JsonNode;
import tools.jackson.databind.ObjectMapper;
import tools.jackson.databind.json.JsonMapper;

/**
 * Regenerates {@code paper/tables/*.tex} from the catalog + cited-bounds +
 * derived-bounds JSON inputs.
 *
 * <p><b>Layered inputs.</b>
 * <ul>
 *   <li>{@code docs/catalog.json} — explicit on-disk schemes with verified
 *       {@code (U, V, W)}.</li>
 *   <li>{@code docs/cited-bounds.json} — published rank claims without
 *       on-disk factor matrices (dagger in the table).</li>
 *   <li>{@code docs/derived-from-cited-bounds.json} — formula-derived rank claims our
 *       constructors can re-derive.</li>
 * </ul>
 * Each layer carries a per-entry {@code field} and {@code commutative} tag;
 * cells are populated according to the rules of
 * {@code paper/sota-conventions.md} (transitive field promotion within
 * Z→Q→R→C; no auto-promotion to F₂; commutative-vs-non-commutative kept
 * strictly separate).
 *
 * <p><b>Output.</b> One {@code .tex} file per (field, commutativity)
 * combination, under {@code paper/tables/}. Each file is a self-contained
 * LaTeX {@code \begin{table}...\end{table}} block ready to be
 * {@code \input}'d from {@code paper/sections/tables.tex}.
 *
 * <p><b>Status.</b> Skeleton today: the file structure and JSON parsing
 * are in place but the per-cell SOTA-merge logic is left as TODO until
 * the conventions in {@code paper/sota-conventions.md} are finalised.
 * Run with {@code --dry-run} to print the cells it would write without
 * touching disk.
 *
 * <p>Invoke via:
 * <pre>
 *   mvn -q exec:java -Dexec.mainClass=eu.solven.matmul.docs.generate.GeneratePaperTables
 * </pre>
 */
@Slf4j
public final class GeneratePaperTables {

	private static final Path OUT_DIR = Path.of("paper/tables");

	/** All shapes we want to surface in the per-field tables. Initial list
	 *  mirrors DIS09 Tables 3 + 4 and the post-2009 additions covered by
	 *  the catalog. Extend as needed; the regen tooling pads missing cells
	 *  with "TBD" rather than failing. */
	private static final List<int[]> SHAPES = List.of(
			new int[] {2, 2, 2}, new int[] {2, 2, 3}, new int[] {2, 3, 3},
			new int[] {2, 3, 4}, new int[] {3, 3, 3}, new int[] {3, 3, 4},
			new int[] {3, 3, 6}, new int[] {3, 4, 4}, new int[] {4, 4, 4},
			new int[] {4, 4, 5}, new int[] {5, 5, 5}, new int[] {6, 6, 6},
			new int[] {7, 7, 7}, new int[] {8, 8, 8}, new int[] {9, 9, 9},
			new int[] {10, 10, 10}, new int[] {12, 12, 12},
			new int[] {16, 16, 16}, new int[] {17, 17, 17},
			new int[] {20, 20, 20}, new int[] {24, 24, 24},
			new int[] {32, 32, 32});

	private GeneratePaperTables() {}

	public static void main(String[] args) throws IOException {
		boolean dryRun = false;
		for (String a : args) {
			if ("--dry-run".equals(a)) dryRun = true;
			else if ("--help".equals(a) || "-h".equals(a)) {
				System.out.println("Usage: GeneratePaperTables [--dry-run]");
				System.out.println("Regenerates paper/tables/*.tex from docs/catalog.json + cited + derived.");
				return;
			} else {
				throw new IllegalArgumentException("Unknown arg: " + a);
			}
		}

		// ── Read the three layers. ────────────────────────────────────
		ObjectMapper mapper = JsonMapper.builder().build();
		JsonNode catalog = mapper.readTree(Path.of("docs/catalog.json").toFile());
		JsonNode cited = mapper.readTree(Path.of("docs/cited-bounds.json").toFile());
		JsonNode derived = mapper.readTree(Path.of("docs/derived-from-cited-bounds.json").toFile());

		// ── Merge into per-cell SOTA. ─────────────────────────────────
		// TODO: implement the SOTA-merge rules per paper/sota-conventions.md.
		// For now, this stub iterates SHAPES × {Q, R, C, F2} × {NC, C}
		// and reports "TBD" cells unless an exact match is found.
		Map<String, List<Cell>> tables = new LinkedHashMap<>();
		for (String field : List.of("Q", "R", "C", "F2")) {
			for (boolean commutative : List.of(false, true)) {
				String tableKey = (commutative ? "c-" : "nc-") + field;
				List<Cell> cells = new ArrayList<>();
				for (int[] shape : SHAPES) {
					cells.add(lookupCell(shape, field, commutative,
							catalog, cited, derived));
				}
				tables.put(tableKey, cells);
			}
		}

		// ── Emit one tex file per table. ──────────────────────────────
		if (!dryRun) {
			Files.createDirectories(OUT_DIR);
		}
		for (Map.Entry<String, List<Cell>> e : tables.entrySet()) {
			String tableKey = e.getKey();
			List<Cell> cells = e.getValue();
			String tex = renderTable(tableKey, cells);
			Path out = OUT_DIR.resolve(tableKey + ".tex");
			if (dryRun) {
				log.info("[dry-run] would write {} ({} cells)", out, cells.size());
				for (Cell c : cells) log.info("    {}", c);
			} else {
				try (PrintWriter pw = new PrintWriter(out.toFile())) {
					pw.print(tex);
				}
				log.info("wrote {} ({} cells)", out, cells.size());
			}
		}
	}

	/** Per-cell record: shape, rank, attribution, dagger flag. */
	record Cell(int n, int m, int p, Integer rank, String attribution,
			Integer year, boolean dagger) {
		@Override public String toString() {
			return String.format("⟨%d,%d,%d⟩ %s%s by %s %s",
					n, m, p,
					rank == null ? "TBD" : rank.toString(),
					dagger ? "†" : "",
					attribution == null ? "—" : attribution,
					year == null ? "" : ("(" + year + ")"));
		}
	}

	/** Stub. TODO: implement SOTA merge per paper/sota-conventions.md.
	 *  Currently returns a TBD cell unless an exact match is found in the
	 *  catalog at the requested (shape, field, commutativity). */
	private static Cell lookupCell(int[] shape, String field, boolean commutative,
			JsonNode catalog, JsonNode cited, JsonNode derived) {
		// Walk catalog.json entries — current schema lists schemes under a
		// "schemes" or top-level array; the exact key depends on the
		// manifest version. Accept either.
		JsonNode entries = catalog.has("schemes") ? catalog.get("schemes") : catalog;
		if (entries.isArray()) {
			Integer bestRank = null;
			String bestAttr = null;
			Integer bestYear = null;
			for (JsonNode e : entries) {
				if (!matchShape(e, shape)) continue;
				if (!matchField(e, field)) continue;
				if (!matchCommutative(e, commutative)) continue;
				int r = e.has("rank") ? e.get("rank").asInt()
						: e.has("r") ? e.get("r").asInt() : Integer.MAX_VALUE;
				if (bestRank == null || r < bestRank) {
					bestRank = r;
					// Attribution is the scheme's own source (the importer/discoverer);
					// the rank-only attribution_for_rank property has been retired.
					bestAttr = e.has("source") ? e.get("source").asText() : "?";
					bestYear = e.has("year") ? e.get("year").asInt() : extractYear(bestAttr);
				}
			}
			if (bestRank != null) {
				return new Cell(shape[0], shape[1], shape[2], bestRank,
						bestAttr, bestYear, false);
			}
		}
		// TODO: cited + derived layers; for now emit a dagger TBD.
		return new Cell(shape[0], shape[1], shape[2], null,
				null, null, true);
	}

	private static boolean matchShape(JsonNode entry, int[] shape) {
		// catalog.json: shape is in "format" array [n, m, p]. Compare
		// sorted-shape since the catalog canonicalises (n ≤ m ≤ p).
		JsonNode fmt = entry.path("format");
		if (!fmt.isArray() || fmt.size() != 3) {
			// Older / alternate schemas use n/m/p fields.
			int n = entry.path("n").asInt(-1);
			int m = entry.path("m").asInt(-1);
			int p = entry.path("p").asInt(-1);
			if (n == -1 || m == -1 || p == -1) return false;
			int[] q = sortedCopy(shape);
			int[] r = sortedCopy(new int[] {n, m, p});
			return q[0] == r[0] && q[1] == r[1] && q[2] == r[2];
		}
		int[] q = sortedCopy(shape);
		int[] r = sortedCopy(new int[] {
				fmt.get(0).asInt(), fmt.get(1).asInt(), fmt.get(2).asInt() });
		return q[0] == r[0] && q[1] == r[1] && q[2] == r[2];
	}

	/** Extract a 4-digit year from a free-text attribution (e.g.
	 *  "Strassen 1969" → 1969). Returns null if no year is present. */
	private static Integer extractYear(String text) {
		if (text == null) return null;
		java.util.regex.Matcher m =
				java.util.regex.Pattern.compile("(19|20)\\d{2}").matcher(text);
		return m.find() ? Integer.parseInt(m.group()) : null;
	}

	private static int[] sortedCopy(int[] s) {
		int[] c = s.clone();
		java.util.Arrays.sort(c);
		return c;
	}

	private static boolean matchField(JsonNode entry, String field) {
		// catalog.json field can be a slash-separated list ("R/Q/Z",
		// "C", "F2") representing the wider field plus all narrower
		// ones it's valid over. Tokenise and apply field-promotion
		// logic per paper/sota-conventions.md.
		String f = entry.path("field").asText("");
		if (f.isEmpty()) return false;
		java.util.Set<String> tokens = new java.util.HashSet<>(
				java.util.Arrays.asList(f.split("/")));
		switch (field) {
			case "Q": return tokens.contains("Q") || tokens.contains("Z") || tokens.contains("ZT");
			case "R": return tokens.contains("R") || tokens.contains("Q")
					|| tokens.contains("Z") || tokens.contains("ZT");
			case "C": return tokens.contains("C") || tokens.contains("R")
					|| tokens.contains("Q") || tokens.contains("Z")
					|| tokens.contains("ZT");
			case "F2":
				for (String t : tokens) {
					if (t.equalsIgnoreCase("F2") || t.equalsIgnoreCase("Z2")
							|| t.equalsIgnoreCase("Z/2")) return true;
				}
				return false;
			default: return tokens.contains(field);
		}
	}

	private static boolean matchCommutative(JsonNode entry, boolean commutative) {
		boolean entryComm = entry.path("commutative").asBoolean(false);
		return entryComm == commutative;
	}

	private static String renderTable(String key, List<Cell> cells) {
		String fieldName = key.replaceFirst("^(nc|c)-", "");
		String fieldLatex = switch (fieldName) {
			case "Q" -> "\\Rationals";
			case "R" -> "\\Reals";
			case "C" -> "\\Complex";
			case "F2" -> "\\Ftwo";
			default -> fieldName;
		};
		boolean commutative = key.startsWith("c-");
		// Each table is per-field (the field is named in the caption), so cells
		// use the field-less 3-arg shape macros, not the 4-arg \nmpfield variants.
		String shapeMacro = commutative ? "\\nmpshapec" : "\\nmpshape";
		StringBuilder sb = new StringBuilder();
		sb.append("% AUTO-GENERATED by eu.solven.matmul.docs.generate.GeneratePaperTables.\n");
		sb.append("% Source: docs/catalog.json + cited + derived as of regeneration time.\n");
		sb.append("% Edit the generator (or the source JSON), not this file.\n");
		sb.append("\n");
		sb.append("\\begin{table}[h]\n");
		sb.append("\\centering\n");
		sb.append(String.format(Locale.ROOT,
				"\\caption{%s rank over %s for selected shapes. \\dag\\ marks a "
						+ "cited bound whose explicit scheme is not in the catalog.}%n",
				commutative ? "Commutative" : "Non-commutative", fieldLatex));
		sb.append(String.format("\\label{tab:%s}%n", key));
		sb.append("\\begin{tabular}{l r l l}\n");
		sb.append("\\toprule\n");
		sb.append("Shape & SOTA rank & Attribution & Year \\\\\n");
		sb.append("\\midrule\n");
		for (Cell c : cells) {
			sb.append(String.format(Locale.ROOT,
					"%s{%d}{%d}{%d} & %s%s & %s & %s \\\\%n",
					shapeMacro, c.n, c.m, c.p,
					c.rank == null ? "TBD" : c.rank.toString(),
					c.dagger ? "\\dag" : "",
					c.attribution == null ? "---" : escapeLatex(c.attribution),
					c.year == null ? "---" : c.year.toString()));
		}
		sb.append("\\bottomrule\n");
		sb.append("\\end{tabular}\n");
		sb.append("\\end{table}\n");
		return sb.toString();
	}

	private static String escapeLatex(String s) {
		return s.replace("_", "\\_").replace("&", "\\&");
	}
}
