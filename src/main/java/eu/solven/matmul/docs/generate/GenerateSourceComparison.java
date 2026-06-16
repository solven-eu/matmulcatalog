package eu.solven.matmul.docs.generate;

import java.nio.file.Files;
import java.nio.file.Path;
import java.util.ArrayList;
import java.util.LinkedHashMap;
import java.util.List;
import java.util.Locale;
import java.util.Map;
import java.util.TreeMap;

import lombok.extern.slf4j.Slf4j;

import tools.jackson.databind.JsonNode;
import tools.jackson.databind.json.JsonMapper;

/**
 * DATA layer of the paper-table pipeline (LaTeX-unaware by design): per-shape
 * best non-commutative ranks split by provenance — THIS WORK (further split into
 * {@code us_derived} = our own constructors / closure, and {@code us_imported} =
 * known schemes we merely carry) vs FMM-Lille vs Perminov — emitted
 * <b>once per field</b> as {@code docs/comparison/us-vs-fmm-vs-perminov-{K}.json}
 * (machine) and {@code .md} (human), for each {@code K} in {@link #FIELDS}.
 *
 * <h3>Why per-field (field discipline)</h3>
 * <p>A bare {@code ⟨n,m,p⟩} rank is meaningless without its field — the
 * canonical ⟨4,4,4⟩ shows why: 47/F₂ (AlphaTensor), 48/C (AlphaEvolve) = 48/Q/R
 * (DPS-2025), 49/Z (Strassen²). The previous single, <i>field-blind</i> table
 * took the global min per shape, so ⟨4,4,4⟩'s comparator column became the
 * F₂ 47 and our Q/Z 49 read as a "loss" against a scheme in a different
 * algebra — while the genuine Q result (48) was masked. Each table here is
 * scoped to one field: a scheme contributes to field {@code K}'s table only if
 * it is valid over {@code K}.</p>
 *
 * <h3>Two "this work" columns (honesty)</h3>
 * <ul>
 *   <li>{@code us_derived} — the best rank our own constructors / closure
 *       produced (sources {@code constructed/}, {@code derived/}, {@code curated/},
 *       {@code solven*}). This is the only column we ever claim as a discovery.</li>
 *   <li>{@code us_imported} — the best rank among schemes we <i>carry</i> as
 *       imports (AlphaTensor, Strassen, Laderman, DPS-2025, …). It answers
 *       "do we properly include the known schemes?" without conflating a
 *       carried import with our own work.</li>
 * </ul>
 *
 * <h3>External-catalog field validity</h3>
 * <ul>
 *   <li><b>Our catalog</b>: field membership read directly from {@code fields[]}.</li>
 *   <li><b>Perminov</b> ({@code references/catalogs/perminov-catalog.json}): each entry
 *       carries a single {@code field} (Q / Z / ZT / …); expanded by the
 *       characteristic-0 inclusion {@code Z ⇒ {F2,F3,Z,Q,R,C}}, {@code Q ⇒ {Q,R,C}},
 *       {@code R ⇒ {R,C}}, {@code C ⇒ {C}} — the same policy {@code StampFields}
 *       uses (an integer scheme reduces mod 2/3; a rational one does not
 *       unconditionally, so Q is NOT widened to F₂/F₃).</li>
 *   <li><b>FMM-Lille</b> ({@code references/catalogs/fmm-lille-catalog.json}): entries
 *       carry <i>no</i> field tag. FMM is a characteristic-0 (rational) catalog
 *       by nature — it does not publish F₂/F₃ schemes — so its entries are
 *       treated as valid over {@code {Q,R,C}} only. They are deliberately
 *       <b>excluded from the Z / F₂ / F₃ tables</b>: we cannot tell from rank
 *       alone whether an FMM scheme is integer or reduces. This is an
 *       assumption, documented here and in each table's header.</li>
 * </ul>
 *
 * <pre>mvn -q -ntp exec:java -Dexec.mainClass=eu.solven.matmul.docs.generate.GenerateSourceComparison</pre>
 */
@Slf4j
public final class GenerateSourceComparison {

	private GenerateSourceComparison() {}

	/** Fields that get their own comparison table, narrowest-characteristic first. */
	static final List<String> FIELDS = List.of("F2", "F3", "Z", "Q", "R", "C");

	public static void main(String[] args) throws Exception {
		JsonMapper mapper = JsonMapper.builder().build();
		JsonNode catalog = mapper.readTree(Files.readString(Path.of("docs/catalog.json")));
		JsonNode schemes = catalog.isArray() ? catalog : catalog.get("schemes");

		JsonNode fmmEntries = entries(mapper, "references/catalogs/fmm-lille-catalog.json");
		JsonNode perEntries = entries(mapper, "references/catalogs/perminov-catalog.json");
		// Perminov's serendipitous 17–32 band — a SEPARATE additional catalog
		// (status.json stops at 16), folded into the same "perminov" column by min.
		JsonNode perSerendEntries = entries(mapper, "references/catalogs/perminov-serendipitous-catalog.json");

		Path dir = Path.of("docs/comparison");
		Files.createDirectories(dir);

		for (String field : FIELDS) {
			emitField(mapper, schemes, fmmEntries, perEntries, perSerendEntries, field, dir);
		}
	}

	/** Build + write the {@code .json} and {@code .md} for a single field. */
	private static void emitField(JsonMapper mapper, JsonNode schemes, JsonNode fmmEntries,
			JsonNode perEntries, JsonNode perSerendEntries, String field, Path dir) throws Exception {
		// shapeKey → column → best rank.  Columns: us_derived, us_imported, fmm, perminov.
		Map<String, Map<String, Object>> byShape = new TreeMap<>();

		// --- THIS WORK (our catalog), split derived vs imported, filtered to `field` ---
		for (JsonNode s : schemes) {
			JsonNode fmt = s.get("format");
			if (fmt == null || !fmt.isArray() || fmt.size() != 3) continue;
			if (s.get("rank") == null) continue;
			if (s.get("commutative") != null && s.get("commutative").asBoolean(false)) continue; // NC only
			if (!fieldsContain(s.get("fields"), field)) continue; // not valid over this field
			int rank = s.get("rank").asInt();
			String column = "us".equals(classify(s)) ? "us_derived" : "us_imported";
			put(byShape, fmt, column, rank);
		}

		// --- External references, field-filtered by their own validity rules ---
		overlay(perEntries, byShape, "perminov",
				e -> perminovValidOver(e.get("field"), field));
		overlay(perSerendEntries, byShape, "perminov",
				e -> perminovValidOver(e.get("field"), field));
		overlay(fmmEntries, byShape, "fmm",
				e -> fmmValidOver(field));

		// --- Derive per-row status + headline counts ---
		// Status reflects whether WE hold the best rank, and HOW:
		//   win          — our DERIVED work beats both external catalogs.
		//   tie_derived  — our DERIVED work matches the external best.
		//   tie_imported — we hold the external best, but only via an IMPORT (our
		//                  own derivation is worse or absent). NOT a loss: we are
		//                  at the frontier, we just didn't re-derive it.
		//   loss         — an external catalog strictly beats EVERYTHING we hold
		//                  (both our derivation and our imports).
		//   carried      — import-only coverage with NO external comparator (we
		//                  are the sole holder of the shape, via an import).
		int wins = 0, tieDerived = 0, tieImported = 0, lossDerived = 0, lossImport = 0, carried = 0, usOnly = 0, extOnly = 0;
		List<Map<String, Object>> comparable = new ArrayList<>();
		for (Map<String, Object> row : byShape.values()) {
			Integer der = (Integer) row.get("us_derived");
			Integer imp = (Integer) row.get("us_imported");
			Integer fmm = (Integer) row.get("fmm");
			Integer pr = (Integer) row.get("perminov");
			Integer ext = min(fmm, pr);
			boolean hasUs = der != null || imp != null;
			if (!hasUs) { if (ext != null) extOnly++; continue; }
			int status;
			if (ext == null) {
				// No external comparator → coverage-only, not a head-to-head row.
				if (der != null) { usOnly++; continue; }
				status = 3; carried++; // import-only unique coverage
			} else if (der != null && der < ext) {
				status = 0; wins++;
			} else if (der != null && der.intValue() == ext.intValue()) {
				status = 2; tieDerived++;
			} else if (imp != null && imp <= ext) {
				// Derived absent or worse than best, but an import holds the best.
				status = 4; tieImported++;
			} else {
				// External strictly beats all we hold. Split the gap by the WINNER's
				// provenance: a PUBLISHED scheme (Perminov, or an FMM entry carrying a
				// cited reference) is an IMPORT gap (we should carry it); an FMM-DERIVED
				// composition (no reference) is a DERIVATION gap (we should re-derive it).
				boolean perminovWins = pr != null && (fmm == null || pr <= fmm);
				boolean importGap = perminovWins || Boolean.TRUE.equals(row.get("fmm_published"));
				if (importGap) { status = 5; lossImport++; }
				else { status = 1; lossDerived++; }
			}
			row.put("_status", status);
			comparable.add(row);
		}
		comparable.sort(GenerateSourceComparison::byGrowingShape);

		writeJson(mapper, dir, field, comparable, wins, tieDerived, tieImported, lossDerived, lossImport, carried, usOnly, extOnly);
		writeMarkdown(dir, field, comparable, wins, tieDerived, tieImported, lossDerived, lossImport, carried, usOnly, extOnly);
		log.info("field {}: {} comparable ({} wins / {} tie_derived / {} tie_imported / {} loss_derived "
				+ "/ {} loss_import / {} carried; {} derived-only, {} external-only)",
				field, comparable.size(), wins, tieDerived, tieImported, lossDerived, lossImport, carried, usOnly, extOnly);
	}

	// ──────────────────────────────────────────────────────────────────────
	//  Field-validity rules
	// ──────────────────────────────────────────────────────────────────────

	/** Our catalog rows carry the full inclusion set in {@code fields[]}. */
	private static boolean fieldsContain(JsonNode fields, String field) {
		if (fields == null || !fields.isArray()) return false;
		for (JsonNode f : fields) {
			if (f.isTextual() && f.asString().equals(field)) return true;
		}
		return false;
	}

	/** Perminov entry's single {@code field}, expanded by char-0 inclusion. */
	private static boolean perminovValidOver(JsonNode fieldNode, String target) {
		if (fieldNode == null || !fieldNode.isTextual()) return false;
		return expand(fieldNode.asString()).contains(target);
	}

	/** FMM carries no field; treat as characteristic-0 rational → {Q,R,C} only. */
	private static boolean fmmValidOver(String target) {
		return target.equals("Q") || target.equals("R") || target.equals("C");
	}

	/** Inclusion set of a narrowest field tag (mirrors {@code StampFields#expand}). */
	private static List<String> expand(String field) {
		switch (field) {
			case "Z":
			case "ZT": return FIELDS;                       // integer ⇒ all
			case "Q":  return List.of("Q", "R", "C");
			case "R":  return List.of("R", "C");
			case "C":  return List.of("C");
			case "F2": return List.of("F2");
			case "F3": return List.of("F3");
			default:   return List.of();                    // unknown → contributes nowhere
		}
	}

	// ──────────────────────────────────────────────────────────────────────
	//  Catalog plumbing
	// ──────────────────────────────────────────────────────────────────────

	private static JsonNode entries(JsonMapper mapper, String path) throws Exception {
		JsonNode root = mapper.readTree(Files.readString(Path.of(path)));
		JsonNode e = root.get("entries");
		if (e == null || !e.isArray()) {
			log.warn("reference catalog {} has no entries[] — column will be empty", path);
		}
		return e;
	}

	private interface EntryFilter { boolean keep(JsonNode entry); }

	/** Overlay an external reference catalog onto {@code column}, keeping the best
	 *  (lowest) rank among entries that pass {@code filter}, keyed by canonical
	 *  sorted shape. Creates the shape row if only the reference covers it. */
	private static void overlay(JsonNode entries, Map<String, Map<String, Object>> byShape,
			String column, EntryFilter filter) {
		if (entries == null || !entries.isArray()) return;
		for (JsonNode e : entries) {
			JsonNode fmt = e.get("format");
			JsonNode rk = e.get("rank");
			if (fmt == null || !fmt.isArray() || fmt.size() != 3 || rk == null) continue;
			if (!filter.keep(e)) continue;
			int[] d = { fmt.get(0).asInt(), fmt.get(1).asInt(), fmt.get(2).asInt() };
			java.util.Arrays.sort(d); // canonical n≤m≤p to match the manifest's sorted keys
			Map<String, Object> row = byShape.computeIfAbsent(d[0] + "x" + d[1] + "x" + d[2], k -> {
				Map<String, Object> m = new LinkedHashMap<>();
				m.put("format", List.of(d[0], d[1], d[2]));
				return m;
			});
			Integer prev = (Integer) row.get(column);
			int rank = rk.asInt();
			if (prev == null || rank < prev) {
				row.put(column, rank);
				// Provenance of the winning entry: an FMM entry with a cited
				// `references` is a PUBLISHED result (import target); without, it is
				// an FMM-DERIVED composition. (Perminov entries carry no references —
				// they are his own published schemes, handled as import targets.)
				row.put(column + "_published", hasReferences(e));
			}
		}
	}

	/** True iff an external entry cites a published prior result (non-empty `references`). */
	private static boolean hasReferences(JsonNode e) {
		JsonNode r = e.get("references");
		return r != null && r.isArray() && r.size() > 0;
	}

	private static void put(Map<String, Map<String, Object>> byShape, JsonNode fmt,
			String column, int rank) {
		String key = fmt.get(0).asInt() + "x" + fmt.get(1).asInt() + "x" + fmt.get(2).asInt();
		Map<String, Object> row = byShape.computeIfAbsent(key, k -> {
			Map<String, Object> m = new LinkedHashMap<>();
			m.put("format", List.of(fmt.get(0).asInt(), fmt.get(1).asInt(), fmt.get(2).asInt()));
			return m;
		});
		Integer prev = (Integer) row.get(column);
		if (prev == null || rank < prev) row.put(column, rank);
	}

	/** Order by GROWING shape: largest dim first, then smallest, then middle. */
	private static int byGrowingShape(Map<String, Object> a, Map<String, Object> b) {
		@SuppressWarnings("unchecked") List<Integer> fa = (List<Integer>) a.get("format");
		@SuppressWarnings("unchecked") List<Integer> fb = (List<Integer>) b.get("format");
		int c = Integer.compare(fa.get(2), fb.get(2));
		if (c != 0) return c;
		c = Integer.compare(fa.get(0), fb.get(0));
		if (c != 0) return c;
		return Integer.compare(fa.get(1), fb.get(1));
	}

	// ──────────────────────────────────────────────────────────────────────
	//  Output
	// ──────────────────────────────────────────────────────────────────────

	private static void writeJson(JsonMapper mapper, Path dir, String field,
			List<Map<String, Object>> rows, int wins, int tieDerived, int tieImported, int lossDerived,
			int lossImport, int carried, int usOnly, int extOnly) throws Exception {
		Map<String, Object> out = new LinkedHashMap<>();
		out.put("_what", "Per-shape best NON-COMMUTATIVE rank over " + field
				+ ", by provenance: us_derived (our constructors) and us_imported (schemes we "
				+ "carry) vs FMM-Lille vs Perminov. Field-filtered: a scheme appears only if "
				+ "valid over " + field + ". FMM has no field tag and is treated as char-0 "
				+ "(valid Q/R/C only).");
		out.put("_field", field);
		out.put("_generator", GenerateSourceComparison.class.getName());
		out.put("_inputs", List.of("docs/catalog.json",
				"references/catalogs/fmm-lille-catalog.json", "references/catalogs/perminov-catalog.json"));
		Map<String, Object> counts = new LinkedHashMap<>();
		counts.put("wins", wins); counts.put("tie_derived", tieDerived);
		counts.put("tie_imported", tieImported);
		counts.put("loss_derived", lossDerived); counts.put("loss_import", lossImport);
		counts.put("carried", carried); counts.put("derived_only", usOnly);
		counts.put("external_only", extOnly);
		out.put("_counts", counts);
		out.put("rows", rows);
		Path json = dir.resolve("us-vs-fmm-vs-perminov-" + field + ".json");
		Files.writeString(json, mapper.writerWithDefaultPrettyPrinter().writeValueAsString(out));
	}

	private static void writeMarkdown(Path dir, String field, List<Map<String, Object>> rows,
			int wins, int tieDerived, int tieImported, int lossDerived, int lossImport,
			int carried, int usOnly, int extOnly) throws Exception {
		StringBuilder md = new StringBuilder();
		md.append("# This work vs FMM-Lille vs Perminov over ").append(field)
				.append(" (non-commutative best ranks)\n\n");
		md.append("Auto-generated by `").append(GenerateSourceComparison.class.getSimpleName())
				.append("` from `docs/catalog.json` + `references/catalogs/fmm-lille-catalog.json` + ")
				.append("`references/catalogs/perminov-catalog.json` — do not edit. ")
				.append("Machine twin: `us-vs-fmm-vs-perminov-").append(field).append(".json`.\n\n");
		md.append("**Field-scoped to ").append(field).append("**: a scheme appears only if valid over ")
				.append(field).append(". `us_derived` = best rank our own constructors/closure produced ")
				.append("(the only column we claim as a discovery); `us_imported` = best rank among ")
				.append("schemes we carry as imports. FMM-Lille carries no per-entry field tag and is ")
				.append("treated as characteristic-0 (valid over Q/R/C only), so it is absent from the ")
				.append("Z/F₂/F₃ tables.\n\n");
		md.append("`vs best` legend: **win** = our derived work beats both external catalogs; ")
				.append("**tie_derived** = our derived work matches the external best; ")
				.append("**tie_imported** = we hold the external best but only via an import (our derivation ")
				.append("is worse or absent — at the frontier, just not re-derived, so NOT a loss); ")
				.append("**loss_derived** = an external catalog beats everything we hold AND the winner is an ")
				.append("FMM-*derived* composition (no cited reference) — a gap in our DERIVATION, re-derivable; ")
				.append("**loss_import** = the winner is a *published* scheme (Perminov, or FMM with a cited ")
				.append("reference) — a gap in our IMPORTS, we should carry it; ")
				.append("**carried** = import-only coverage with no external comparator. The `FMM src` column ")
				.append("marks whether the FMM-Lille entry is `pub` (cites a published result) or `deriv` ")
				.append("(FMM's own composition).\n\n");
		md.append(String.format(Locale.ROOT,
				"Head-to-head over %s: **%d wins, %d tie_derived, %d tie_imported, %d loss_derived, "
				+ "%d loss_import**, plus **%d carried** (we are the sole holder, via an import). "
				+ "(Separately, %d shapes only our derived work covers and %d only an external catalog "
				+ "covers.)%n%n",
				field, wins, tieDerived, tieImported, lossDerived, lossImport, carried, usOnly, extOnly));

		// Indexed by _status: 0 win, 1 loss_derived, 2 tie_derived, 3 carried,
		// 4 tie_imported, 5 loss_import.
		String[] label = { "win", "loss_derived", "tie_derived", "carried", "tie_imported", "loss_import" };
		List<String[]> table = new ArrayList<>();
		table.add(new String[] { "shape", "us_derived", "us_imported", "FMM-Lille", "FMM src", "Perminov", "vs best" });
		for (Map<String, Object> row : rows) {
			@SuppressWarnings("unchecked") List<Integer> f = (List<Integer>) row.get("format");
			int st = (int) row.get("_status");
			Object der = row.get("us_derived");
			// Only bold us_derived when it is a strict win.
			String derCell = der == null ? "—" : (st == 0 ? "**" + der + "**" : String.valueOf(der));
			String fmmSrc = row.get("fmm") == null ? "—"
					: (Boolean.TRUE.equals(row.get("fmm_published")) ? "pub" : "deriv");
			table.add(new String[] {
					"⟨" + f.get(0) + "," + f.get(1) + "," + f.get(2) + "⟩",
					derCell, String.valueOf(dash(row.get("us_imported"))),
					String.valueOf(dash(row.get("fmm"))), fmmSrc, String.valueOf(dash(row.get("perminov"))),
					label[st] });
		}
		int cols = 7;
		int[] w = new int[cols];
		for (String[] r : table) {
			for (int i = 0; i < cols; i++) w[i] = Math.max(w[i], r[i].codePointCount(0, r[i].length()));
		}
		appendAlignedRow(md, table.get(0), w);
		md.append('|');
		for (int i = 0; i < cols; i++) md.append(' ').append("-".repeat(w[i])).append(" |");
		md.append('\n');
		for (int r = 1; r < table.size(); r++) appendAlignedRow(md, table.get(r), w);
		Files.writeString(dir.resolve("us-vs-fmm-vs-perminov-" + field + ".md"), md.toString());
	}

	private static void appendAlignedRow(StringBuilder sb, String[] cells, int[] w) {
		sb.append('|');
		for (int i = 0; i < cells.length; i++) {
			int pad = w[i] - cells[i].codePointCount(0, cells[i].length());
			sb.append(' ').append(cells[i]);
			for (int k = 0; k < pad; k++) sb.append(' ');
			sb.append(" |");
		}
		sb.append('\n');
	}

	private static Object dash(Object v) {
		return v == null ? "—" : v;
	}

	private static Integer min(Integer a, Integer b) {
		if (a == null) return b;
		if (b == null) return a;
		return Math.min(a, b);
	}

	/** Provenance bucket of a manifest row, or null when un-classifiable. */
	static String classify(JsonNode s) {
		String src = lower(s, "source");
		String file = lower(s, "file");
		String prior = lower(s, "prior_art");
		String all = src + "|" + file + "|" + prior;
		if (all.contains("fmm") || all.contains("lille")) return "fmm";
		if (all.contains("perminov")) return "perminov";
		if (src.contains("this work") || src.contains("solven")
				|| file.contains("/constructed/") || file.contains("/derived/") || file.contains("/curated/")
				|| file.contains("derived") || file.contains("constructive")) {
			return "us";
		}
		// Direct historical imports (strassen, smirnov, laderman, alphatensor…).
		return "other";
	}

	private static String lower(JsonNode s, String field) {
		JsonNode v = s.get(field);
		return v != null && v.isTextual() ? v.asString().toLowerCase(Locale.ROOT) : "";
	}
}
