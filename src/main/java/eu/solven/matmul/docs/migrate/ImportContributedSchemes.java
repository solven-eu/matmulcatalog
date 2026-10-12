package eu.solven.matmul.docs.migrate;

import java.io.File;
import java.io.IOException;
import java.nio.file.Files;
import java.nio.file.Path;
import java.util.ArrayList;
import java.util.LinkedHashMap;
import java.util.List;
import java.util.Map;
import java.util.stream.Stream;

import eu.solven.matmul.NonCubicBilinearAlgorithm;
import eu.solven.matmul.algebra.Field;
import eu.solven.matmul.catalog.FieldAwareLookup;
import eu.solven.matmul.catalog.SchemeAnalysis;
import eu.solven.matmul.catalog.SchemeIO;
import eu.solven.matmul.verifiers.Verifier;

import lombok.extern.slf4j.Slf4j;

/**
 * Import <b>externally contributed explicit schemes</b> — factor matrices handed to
 * us through a GitHub issue / mail / zip (first use: issue #7, Marcos Adriano's six
 * integer {@code ⟨2,p,n⟩} schemes at the exact Hopcroft–Kerr formula, task #9) —
 * into {@code schemes/known/}, behind the same fail-loud gate the HK emitter uses.
 *
 * <p>Per input file (any of the three on-disk encodings {@link SchemeIO#readBilinear}
 * accepts — {@code u_sparse}, dense, Perminov-reduced; shape/rank/coefficients are
 * read from CONTENT, the contributor's own metadata is ignored):</p>
 * <ol>
 *   <li><b>Exact symbolic proof</b> — {@link Verifier#isExactNonCubic}: every
 *       coefficient scaled to a {@code BigInteger} by the common denominator, the
 *       full bilinear identity checked term-wise. No floating point in the
 *       decision.</li>
 *   <li><b>Random matmul spot-check</b> — {@link Verifier#passesRandomMatmulSpotCheck}:
 *       an independent numeric path (runs the scheme on random {@code A,B}).</li>
 *   <li><b>Field stamping from coefficients</b> — {@code Z} iff all-integer;
 *       {@code F2}/{@code F3} iff the mod-2 / mod-3 reductions verify exactly
 *       ({@link Verifier#isExactNonCubicFp}); {@code Q,R,C} always (char-0 exact
 *       rationals). Never trusts the contributor's {@code fields[]}.</li>
 *   <li><b>Disk round-trip</b> — the WRITTEN canonical file is re-read and must
 *       re-verify exactly and hash identically (guards write-path lossiness).</li>
 * </ol>
 *
 * <p>Then stamps provenance ({@code source, year, source_url, discovery,
 * attribution_for_rank, discovery_note, notes}) and the standard Phase-2 metrics
 * ({@link SchemeAnalysis#defaults()}: {@code verified, additions, has_buds/buds,
 * projection_margin, zt}) through {@link SchemeIO#updateFields} — the single
 * canonical writer. The catalog's current best rank for the shape (over the
 * {@code --field}, default {@code Q}) is reported as BEATS / ties / above; a
 * worse-than-catalog scheme is still importable (historical value — see
 * CLAUDE.md "don't be selective about improvements only") but is logged loudly.</p>
 *
 * <p>Dry-run by default (full gate, nothing written); {@code --execute} writes to
 * {@code schemes/{category}/section{maxDim}/{n}x{m}x{p}-r{r}-{note}-{hash7}.json}.
 * Any gate failure aborts the whole run with a non-zero exit — a contributed batch
 * is imported all-or-nothing.</p>
 *
 * <pre>{@code
 * mvn -q -ntp exec:java -Dexec.mainClass=eu.solven.matmul.docs.migrate.ImportContributedSchemes \
 *     -Dexec.args="path/to/dir-or-files --note=adriano_2026 --source='Adriano 2026' --year=2026 \
 *        --source-url=https://github.com/solven-eu/matmulcatalog/issues/7 \
 *        --discovery=false --attribution='Hopcroft & Kerr 1971 …' [--discovery-note=…] [--notes=…] \
 *        [--category=known] [--field=Q] [--execute]"
 * # or, with the provenance in a UTF-8 JSON file (same keys as stamped; flags override):
 *     -Dexec.args="path/to/dir --meta=references/contributions/issue7-adriano-2026.meta.json --execute"
 * }</pre>
 */
@Slf4j
public final class ImportContributedSchemes {
	private ImportContributedSchemes() {}

	private static final Path SCHEMES_ROOT = Path.of("src/main/resources/schemes");

	/** Parsed CLI options: provenance to stamp + run switches. */
	public record Options(List<Path> inputs, String note, String source, Integer year, String sourceUrl,
			String discovery, String attribution, String discoveryNote, String notes, String category,
			String fieldTag, boolean execute) {
	}

	/** Outcome of the verification gate for one input file. */
	public record Gate(NonCubicBilinearAlgorithm alg, boolean exactSymbolic, boolean spotCheck,
			List<String> fields, boolean integer, boolean ternary, Integer catalogBest) {
		public boolean passed() {
			return exactSymbolic && spotCheck;
		}

		public String vsCatalog() {
			if (catalogBest == null) return "ABSENT";
			if (alg.r < catalogBest) return "BEATS " + catalogBest + " (−" + (catalogBest - alg.r) + ")";
			if (alg.r == catalogBest) return "ties " + catalogBest;
			return "ABOVE " + catalogBest + " (+" + (alg.r - catalogBest) + ")";
		}
	}

	public static void main(String[] args) throws Exception {
		Options opts = parseArgs(args);
		if (opts.inputs().isEmpty()) {
			throw new IllegalArgumentException("usage: ImportContributedSchemes <file|dir>... --note=<token>"
					+ " --source=<label> [--year=YYYY] [--source-url=…] [--discovery=true|false|TBD]"
					+ " [--attribution=…] [--discovery-note=…] [--notes=…] [--category=known] [--field=Q] [--execute]");
		}
		List<File> files = expand(opts.inputs());
		if (files.isEmpty()) {
			throw new IllegalArgumentException("no *.json input under " + opts.inputs());
		}
		FieldAwareLookup lookup = new FieldAwareLookup(Field.fromTag(opts.fieldTag()), SCHEMES_ROOT);

		long t0 = System.nanoTime();
		List<String> failures = new ArrayList<>();
		List<String> report = new ArrayList<>();
		// Phase 1 — gate EVERY input before writing ANY (all-or-nothing: a batch
		// with one bad file is not trusted, and no partial import lands on disk).
		List<Gate> passed = new ArrayList<>();
		for (File in : files) {
			long s = System.nanoTime();
			Gate g;
			try {
				g = gate(in, lookup);
			} catch (Exception e) {
				failures.add(in.getName() + ": unreadable — " + e);
				continue;
			}
			String line = String.format("%s → ⟨%d,%d,%d⟩ r=%d  exact=%b spot=%b  fields=%s%s  catalog(%s): %s  (%d ms)",
					in.getName(), g.alg().n, g.alg().m, g.alg().p, g.alg().r, g.exactSymbolic(), g.spotCheck(),
					g.fields(), g.ternary() ? " ZT" : "", opts.fieldTag(), g.vsCatalog(),
					(System.nanoTime() - s) / 1_000_000L);
			if (!g.passed()) {
				failures.add(in.getName() + ": GATE FAILED — " + line);
				log.error("{}", line);
				continue;
			}
			log.info("{}", line);
			report.add(line);
			if (g.catalogBest() != null && g.alg().r > g.catalogBest()) {
				log.warn("{}: rank {} is ABOVE the catalog's {} — importing anyway (historical record); check the source",
						in.getName(), g.alg().r, g.catalogBest());
			}
			passed.add(g);
		}

		// Phase 2 — write, only when the whole batch passed and --execute was given.
		int written = 0;
		if (opts.execute() && failures.isEmpty()) {
			for (Gate g : passed) {
				File out = write(g, opts);
				written++;
				log.info("wrote {}", SCHEMES_ROOT.relativize(out.toPath()));
			}
		}

		log.info("==================================================================");
		for (String r : report) log.info("  OK  {}", r);
		for (String f : failures) log.error("  FAIL {}", f);
		log.info("ImportContributedSchemes: {} inputs, {} passed, {} failed, {} written ({}) — {} ms",
				files.size(), report.size(), failures.size(), written,
				opts.execute() ? (failures.isEmpty() ? "--execute" : "--execute REFUSED: batch not clean")
						: "DRY-RUN, pass --execute to write",
				(System.nanoTime() - t0) / 1_000_000L);
		if (!failures.isEmpty()) {
			throw new IllegalStateException(failures.size() + " input(s) failed the import gate — nothing written");
		}
	}

	/**
	 * The verification gate for one contributed file: exact symbolic proof + random
	 * spot-check + coefficient-driven field stamping + catalog comparison. Pure
	 * (writes nothing); the tests exercise it on a committed scheme and on a
	 * sabotaged copy.
	 */
	public static Gate gate(File in, FieldAwareLookup lookup) throws IOException {
		return gate(SchemeIO.readBilinear(in), lookup);
	}

	/** {@link #gate(File, FieldAwareLookup)} for a scheme already in memory (a parsed
	 *  upstream artifact that is not one of the JSON encodings). */
	public static Gate gate(NonCubicBilinearAlgorithm alg, FieldAwareLookup lookup) {
		boolean exact = Verifier.isExactNonCubic(alg);
		boolean spot = Verifier.passesRandomMatmulSpotCheck(alg);
		boolean integer = allIntegers(alg);
		List<String> fields = new ArrayList<>();
		if (Verifier.isExactNonCubicFp(alg, 2)) fields.add("F2");
		if (Verifier.isExactNonCubicFp(alg, 3)) fields.add("F3");
		if (integer) fields.add("Z");
		fields.add("Q");
		fields.add("R");
		fields.add("C");
		boolean ternary = integer && SchemeIO.isTernary(alg);
		Integer best = null;
		if (lookup != null) {
			// findRank falls back to the naïve n·m·p (a valid formula bound) when the
			// catalog holds nothing for the shape — report that as ABSENT, not as a tie.
			int v = lookup.findRank(alg.n, alg.m, alg.p);
			long naive = (long) alg.n * alg.m * alg.p;
			if (v > 0 && v < naive) best = v;
		}
		return new Gate(alg, exact, spot, fields, integer, ternary, best);
	}

	/**
	 * Write the canonical file + provenance + Phase-2 metrics, then prove the DISK
	 * artifact (re-read, exact re-verify, identical content hash). Fail-loud.
	 */
	static File write(Gate g, Options opts) throws IOException {
		NonCubicBilinearAlgorithm alg = g.alg();
		int maxDim = Math.max(alg.n, Math.max(alg.m, alg.p));
		Path dir = SCHEMES_ROOT.resolve(opts.category()).resolve("section" + maxDim);
		Files.createDirectories(dir);
		File out = dir.resolve(SchemeIO.canonicalName(alg, opts.note())).toFile();
		if (out.isFile()) {
			log.warn("{} already exists (identical content hash) — re-stamping metadata only", out.getName());
		}
		SchemeIO.write(alg, out);

		Map<String, Object> meta = new LinkedHashMap<>();
		meta.put("fields", g.fields());
		meta.put("fields_not", List.of());
		meta.put("commutative", false);
		meta.put("source", opts.source());
		if (opts.year() != null) meta.put("year", opts.year());
		if (opts.sourceUrl() != null) meta.put("source_url", opts.sourceUrl());
		if (opts.discovery() != null) {
			Object d = "true".equals(opts.discovery()) ? Boolean.TRUE
					: "false".equals(opts.discovery()) ? Boolean.FALSE : opts.discovery();
			meta.put("discovery", d);
		}
		if (opts.attribution() != null) meta.put("attribution_for_rank", opts.attribution());
		if (opts.discoveryNote() != null) meta.put("discovery_note", opts.discoveryNote());
		if (opts.notes() != null) meta.put("notes", opts.notes());
		meta.put("verified", true);
		// Standard Phase-2 metrics (additions, buds, projection_margin, zt, verified):
		// the same analyses EnrichSchemeMetrics stamps catalog-wide, run on the
		// in-memory scheme we just proved — so the file is complete on import and
		// the manifest never has to re-expand it.
		Field field = Field.fromTag(g.fields().contains("Z") ? "Z" : "Q");
		for (SchemeAnalysis an : SchemeAnalysis.defaults()) {
			for (Map.Entry<String, Object> e : an.analyse(alg, field).entrySet()) {
				meta.putIfAbsent(e.getKey(), e.getValue());
			}
		}
		SchemeIO.updateFields(out, meta, List.of(), true);

		// Disk round-trip: the published bits, not the in-memory object, are what we vouch for.
		NonCubicBilinearAlgorithm back = SchemeIO.readBilinear(out);
		if (!Verifier.isExactNonCubic(back)) {
			throw new IllegalStateException("DISK ROUND-TRIP not exact: " + out);
		}
		if (!SchemeIO.contentHash(back).equals(SchemeIO.contentHash(alg))) {
			throw new IllegalStateException("DISK ROUND-TRIP content hash drifted: " + out);
		}
		return out;
	}

	private static boolean allIntegers(NonCubicBilinearAlgorithm alg) {
		for (double[][] f : new double[][][] { alg.denseU(), alg.denseV(), alg.denseW() }) {
			for (double[] row : f) {
				for (double v : row) {
					if (Math.abs(v - Math.round(v)) > 1e-9) return false;
				}
			}
		}
		return true;
	}

	private static List<File> expand(List<Path> inputs) throws IOException {
		List<File> files = new ArrayList<>();
		for (Path p : inputs) {
			if (Files.isDirectory(p)) {
				try (Stream<Path> ls = Files.list(p)) {
					ls.filter(f -> f.toString().endsWith(".json")).sorted().forEach(f -> files.add(f.toFile()));
				}
			} else if (Files.isRegularFile(p)) {
				files.add(p.toFile());
			} else {
				throw new IllegalArgumentException("no such file or directory: " + p);
			}
		}
		return files;
	}

	static Options parseArgs(String[] args) throws IOException {
		List<Path> inputs = new ArrayList<>();
		String note = null, source = null, sourceUrl = null, discovery = "TBD", attribution = null;
		String discoveryNote = null, notes = null, category = "known", fieldTag = "Q";
		Integer year = null;
		boolean execute = false;
		// --meta=<file.json> — the same keys as the JSON they get stamped as
		// (note, source, year, source_url, discovery, attribution_for_rank,
		// discovery_note, notes, category). Read as UTF-8, so prose with ⌈…⌉ / ⟨…⟩
		// survives a Windows command line (sun.jnu.encoding=Cp1252 mangles CLI
		// args). Explicit --flags below override the file.
		for (String a : args) {
			if (!a.startsWith("--meta=")) continue;
			var meta = SchemeIO.parseJson(new File(a.substring("--meta=".length())));
			if (meta.has("note")) note = meta.get("note").asString();
			if (meta.has("source")) source = meta.get("source").asString();
			if (meta.has("year")) year = meta.get("year").asInt();
			if (meta.has("source_url")) sourceUrl = meta.get("source_url").asString();
			if (meta.has("discovery")) discovery = meta.get("discovery").asString();
			if (meta.has("attribution_for_rank")) attribution = meta.get("attribution_for_rank").asString();
			if (meta.has("discovery_note")) discoveryNote = meta.get("discovery_note").asString();
			if (meta.has("notes")) notes = meta.get("notes").asString();
			if (meta.has("category")) category = meta.get("category").asString();
		}
		for (String a : args) {
			if ("--execute".equals(a)) execute = true;
			else if (a.startsWith("--meta=")) continue;
			else if (a.startsWith("--note=")) note = a.substring("--note=".length());
			else if (a.startsWith("--source=")) source = a.substring("--source=".length());
			else if (a.startsWith("--year=")) year = Integer.parseInt(a.substring("--year=".length()));
			else if (a.startsWith("--source-url=")) sourceUrl = a.substring("--source-url=".length());
			else if (a.startsWith("--discovery=")) discovery = a.substring("--discovery=".length());
			else if (a.startsWith("--attribution=")) attribution = a.substring("--attribution=".length());
			else if (a.startsWith("--discovery-note=")) discoveryNote = a.substring("--discovery-note=".length());
			else if (a.startsWith("--notes=")) notes = a.substring("--notes=".length());
			else if (a.startsWith("--category=")) category = a.substring("--category=".length());
			else if (a.startsWith("--field=")) fieldTag = a.substring("--field=".length());
			else if (a.startsWith("--")) throw new IllegalArgumentException("unknown option " + a);
			else inputs.add(Path.of(a));
		}
		if (execute && (note == null || source == null)) {
			throw new IllegalArgumentException("--note and --source are required with --execute");
		}
		if (note != null && !note.matches("[a-z0-9_]+")) {
			throw new IllegalArgumentException("--note must be a filename token [a-z0-9_]+, got " + note);
		}
		if (!List.of("true", "false", "TBD").contains(discovery)) {
			throw new IllegalArgumentException("--discovery must be true|false|TBD, got " + discovery);
		}
		if ("false".equals(discovery) && execute && attribution == null) {
			throw new IllegalArgumentException("--discovery=false requires --attribution=<earliest source of the bound>");
		}
		return new Options(inputs, note, source, year, sourceUrl, discovery, attribution, discoveryNote, notes,
				category, fieldTag, execute);
	}
}
