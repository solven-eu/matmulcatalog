package eu.solven.matmul.docs.verify;

import java.io.IOException;
import java.nio.file.Files;
import java.nio.file.Path;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.LinkedHashMap;
import java.util.LinkedHashSet;
import java.util.List;
import java.util.Map;
import java.util.Set;
import java.util.concurrent.atomic.AtomicInteger;
import java.util.regex.Matcher;
import java.util.regex.Pattern;
import java.util.stream.Stream;

import eu.solven.matmul.NonCubicBilinearAlgorithm;
import eu.solven.matmul.SymbolicVerifier;
import eu.solven.matmul.SymbolicVerifier.Algebra;
import eu.solven.matmul.Verifier;
import eu.solven.matmul.catalog.MatrixJsonFormatter;
import eu.solven.matmul.catalog.SchemeIO;
import lombok.extern.slf4j.Slf4j;
import tools.jackson.databind.JsonNode;
import tools.jackson.databind.json.JsonMapper;
import tools.jackson.databind.node.ArrayNode;
import tools.jackson.databind.node.ObjectNode;

/**
 * Combined per-JSON sanity migration (task #174).
 *
 * <p>For each bilinear scheme JSON under {@code src/main/resources/schemes/}:
 * <ol>
 *   <li>Run {@link SymbolicVerifier} across {F2, F3, Z, Q, R, C} and collect
 *       the fields where the trilinear identity verifies vs fails.</li>
 *   <li>Write {@code fields: [...]} (verifies) and {@code fields_not: [...]}
 *       (tested, failed). Drop the fragmented legacy {@code field} +
 *       {@code z2} + {@code complex} tags.</li>
 *   <li>Set {@code commutative: false} for bilinear (U, V, W) schemes — the
 *       encoding structurally forbids commutativity dependencies (each m_r
 *       is a pure A-side × B-side product; no A·A or B·B possible).</li>
 *   <li>Re-emit {@code u_sparse}/{@code v_sparse}/{@code w_sparse} in the
 *       new row-oriented map format ({@code {"k": {"i": [...], "c": [...]}}})
 *       introduced 2026-06-03.</li>
 * </ol>
 *
 * <p>Scope of this driver: BILINEAR u_sparse files only. Stubs and
 * non-bilinear (Waksman-style) files are handled by follow-up drivers —
 * stubs inherit {@code fields[]} from lineage atoms; non-bilinear schemes
 * need a matrix-scalar Verifier variant to confirm commutativity tags.</p>
 *
 * <h2>Usage</h2>
 * <pre>
 *   # Pilot on a single file:
 *   mvn -q exec:java -Dexec.mainClass=eu.solven.matmul.docs.verify.SanityCatalogMigration \
 *       -Dexec.args="--pilot src/main/resources/schemes/section2/strassen-2x2x2_m7_a18.json"
 *
 *   # Pilot multiple files (no in-place writes; emits target/migration-preview/):
 *   mvn -q exec:java -Dexec.mainClass=eu.solven.matmul.docs.verify.SanityCatalogMigration \
 *       -Dexec.args="--pilot path1.json path2.json ..."
 *
 *   # Full sweep with in-place writes:
 *   mvn -q exec:java -Dexec.mainClass=eu.solven.matmul.docs.verify.SanityCatalogMigration \
 *       -Dexec.args="--apply"
 *
 *   # Dry-run report (no writes, prints summary):
 *   mvn -q exec:java -Dexec.mainClass=eu.solven.matmul.docs.verify.SanityCatalogMigration
 * </pre>
 */
@Slf4j
public final class SanityCatalogMigration {

	private static final Path SCHEMES_ROOT = Path.of("src/main/resources/schemes");
	private static final Path PILOT_OUT = Path.of("target/migration-preview");

	/**
	 * Candidate fields, tested in this order. Order in {@code fields[]} /
	 * {@code fields_not[]} follows this iteration so output is stable.
	 */
	private static final List<Algebra> CANDIDATES =
			List.of(Algebra.F2, Algebra.F3, Algebra.Z, Algebra.Q, Algebra.R, Algebra.C);

	private static final JsonMapper MAPPER = JsonMapper.builder().build();

	private SanityCatalogMigration() {}

	public static void main(String[] args) throws Exception {
		List<String> argv = Arrays.asList(args);
		boolean apply = argv.contains("--apply");
		boolean verify = argv.contains("--verify");
		int pilotIdx = argv.indexOf("--pilot");
		int rootIdx = argv.indexOf("--root");
		List<Path> files;
		if (pilotIdx >= 0) {
			files = new ArrayList<>();
			for (int i = pilotIdx + 1; i < argv.size(); i++) {
				if (argv.get(i).startsWith("--")) break;
				files.add(Path.of(argv.get(i)));
			}
			if (files.isEmpty()) {
				throw new IllegalArgumentException("--pilot requires at least one path");
			}
			Files.createDirectories(PILOT_OUT);
		} else {
			Path root = rootIdx >= 0 && rootIdx + 1 < argv.size()
					? Path.of(argv.get(rootIdx + 1))
					: SCHEMES_ROOT;
			try (Stream<Path> walk = Files.walk(root)) {
				files = walk.filter(p -> p.toString().endsWith(".json")).toList();
			}
		}

		AtomicInteger bilinear = new AtomicInteger();
		AtomicInteger skippedStub = new AtomicInteger();
		AtomicInteger skippedNB = new AtomicInteger();
		AtomicInteger skippedComplex = new AtomicInteger();
		AtomicInteger skippedReduced = new AtomicInteger();
		AtomicInteger errors = new AtomicInteger();
		AtomicInteger rewritten = new AtomicInteger();
		AtomicInteger seen = new AtomicInteger();
		long startNs = System.nanoTime();
		int total = files.size();

		boolean writeInPlace = apply;
		boolean writeToPreview = pilotIdx >= 0;
		Path previewDir = PILOT_OUT;

		files.parallelStream().forEach(file -> {
			int s = seen.incrementAndGet();
			if (s % 200 == 0) {
				long elapsedMs = (System.nanoTime() - startNs) / 1_000_000L;
				double rate = s * 1000.0 / Math.max(1, elapsedMs);
				long etaMs = (long) ((total - s) / rate * 1000.0);
				log.info("[progress] {}/{} processed ({}b rewritten, {}err, {}stub, {}nb, {}red), {}ms elapsed, ETA {}s",
						s, total, rewritten.get(), errors.get(), skippedStub.get(),
						skippedNB.get(), skippedReduced.get(), elapsedMs, etaMs / 1000);
			}
			try {
				JsonNode root = SchemeIO.parseJson(file.toFile());
				if (SchemeIO.isStub(root)) { skippedStub.incrementAndGet(); return; }
				if (SchemeIO.isNonBilinear(root)) { skippedNB.incrementAndGet(); return; }
				if (SchemeIO.isComplex(root)) { skippedComplex.incrementAndGet(); return; }
				if (SchemeIO.isReduced(root)) { skippedReduced.incrementAndGet(); return; }
				if (!root.has("u_sparse") && !root.has("u")) { skippedStub.incrementAndGet(); return; }

				bilinear.incrementAndGet();
				String migrated = migrateBilinear(file, root, verify);
				if (migrated == null) { errors.incrementAndGet(); return; }

				Path outFile = writeToPreview ? previewDir.resolve(file.getFileName()) : file;
				if (writeInPlace || writeToPreview) {
					Files.writeString(outFile, migrated);
					rewritten.incrementAndGet();
				}
			} catch (Exception e) {
				log.warn("failed on {}: {}", file, e.getMessage());
				errors.incrementAndGet();
			}
		});

		log.info("Sanity migration summary: bilinear={} rewritten={} stubs={} nb={} complex={} reduced={} errors={}",
				bilinear.get(), rewritten.get(), skippedStub.get(), skippedNB.get(),
				skippedComplex.get(), skippedReduced.get(), errors.get());
	}

	/**
	 * Migrates a single bilinear JSON: runs the per-field Verifier sweep,
	 * writes {@code fields[]}, {@code fields_not[]}, {@code commutative:
	 * false}, and re-emits {@code u_sparse}/{@code v_sparse}/{@code w_sparse}
	 * in the new row-oriented map format. Returns the rewritten JSON, or
	 * {@code null} if the scheme failed to verify in any candidate algebra
	 * (in which case the file should be left alone for manual triage).
	 */
	public static String migrateBilinear(Path file, JsonNode root, boolean verify) throws IOException {
		NonCubicBilinearAlgorithm alg = SchemeIO.read(root);
		FieldSweepResult sweep = verify
				? sweepFields(file.getFileName().toString(), root, alg)
				: fieldsFromTags(root, file.getFileName().toString());
		if (sweep.fields.isEmpty()) {
			log.warn("no field signal: {}", file.getFileName());
			return null;
		}
		return rewriteJson(root, alg, sweep);
	}

	/**
	 * Structural-only path: read existing {@code field} / {@code z2} /
	 * {@code complex} tags, project to a single-element {@code fields[]}.
	 * Use when the full Verifier-based narrowing is delegated to a separate
	 * job (e.g. {@code FieldWideningSweep}). Cheap — no per-field identity check.
	 */
	public static FieldSweepResult fieldsFromTags(JsonNode root, String filename) {
		List<String> fields = new ArrayList<>();
		// F2 / complex status is often encoded ONLY in the filename token (the
		// AlphaTensor-F2 corpus has no `z2` flag in its JSON body), so consult
		// the filename before falling back to the JSON `field`/`z2`/`complex`
		// tags — otherwise an F2 scheme leaks into the char-0 lookup. See the
		// 2026-06-03 ⟨8,8,8⟩=329 / ⟨2,3,3⟩=15 mis-tag incident.
		String lname = filename == null ? "" : filename.toLowerCase();
		boolean f2Name = lname.contains("_f2-") || lname.contains("_f2_")
				|| lname.contains("-f2-") || lname.contains("-f2_")
				|| lname.contains("atf2") || lname.contains("_z2-") || lname.contains("_z2_");
		boolean complexName = lname.contains("0.5xc") || lname.contains("xc-") || lname.contains("_c-");
		if (root.has("z2") && root.get("z2").asBoolean(false) || f2Name) {
			fields.add("F2");
		} else if (root.has("complex") && root.get("complex").asBoolean(false) || complexName) {
			fields.add("C");
		} else if (root.has("field") && root.get("field").isTextual()) {
			String tag = root.get("field").asString();
			// Canonicalise common shorthands; preserve ZT/Z distinction by collapsing to Z
			// (Z and ZT differ in catalog provenance, not in algebra).
			String canon = switch (tag) {
			case "Z", "ZT" -> "Z";
			case "Q", "R", "C", "F2", "F3" -> tag;
			default -> {
				if (tag.contains("F2")) yield "F2";
				if (tag.contains("0.5") && tag.contains("C")) yield "C";
				if (tag.contains("0.5") && tag.contains("Z")) yield "Q";
				if (tag.contains("C")) yield "C";
				if (tag.contains("R")) yield "R";
				if (tag.contains("Q")) yield "Q";
				if (tag.contains("Z")) yield "Z";
				yield "Z";  // default: characteristic-0 catch-all
			}
			};
			fields.add(canon);
		} else {
			fields.add("Z");  // default
		}
		return new FieldSweepResult(fields, List.of(), Map.of());
	}

	public record FieldSweepResult(List<String> fields, List<String> fieldsNot, Map<String, String> reasons) {}

	public static FieldSweepResult sweepFields(String filenameHint, JsonNode root,
			NonCubicBilinearAlgorithm alg) {
		List<String> fields = new ArrayList<>();
		List<String> fieldsNot = new ArrayList<>();
		Map<String, String> reasons = new LinkedHashMap<>();
		for (Algebra cand : CANDIDATES) {
			SymbolicVerifier.Result r = verifyOne(alg, cand);
			if (r.verified()) {
				fields.add(cand.name());
			} else {
				fieldsNot.add(cand.name());
				reasons.put(cand.name(), r.reason());
			}
		}
		return new FieldSweepResult(fields, fieldsNot, reasons);
	}

	private static SymbolicVerifier.Result verifyOne(NonCubicBilinearAlgorithm alg, Algebra cand) {
		try {
			return switch (cand) {
			case F2 -> {
				if (!allIntegers(alg)) yield SymbolicVerifier.Result.fail(cand, "non-integer coefficient");
				yield Verifier.isExactNonCubicF2(alg)
						? SymbolicVerifier.Result.ok(cand)
						: SymbolicVerifier.Result.fail(cand, "F2 identity fails");
			}
			case F3 -> {
				if (!allIntegers(alg)) yield SymbolicVerifier.Result.fail(cand, "non-integer coefficient");
				yield Verifier.isExactNonCubicF3(alg)
						? SymbolicVerifier.Result.ok(cand)
						: SymbolicVerifier.Result.fail(cand, "F3 identity fails");
			}
			default -> SymbolicVerifier.verifyBilinear(alg, cand);
			};
		} catch (Exception e) {
			return SymbolicVerifier.Result.fail(cand, "exception: " + e.getMessage());
		}
	}

	private static boolean allIntegers(NonCubicBilinearAlgorithm alg) {
		for (double[][] mat : new double[][][] { alg.denseU(), alg.denseV(), alg.denseW() }) {
			for (double[] row : mat) {
				for (double v : row) {
					if (v != Math.rint(v)) return false;
				}
			}
		}
		return true;
	}

	/**
	 * Build the rewritten JSON: preserve all original metadata, replace the
	 * fragmented field tags with the unified pair, set {@code commutative:
	 * false} for bilinear schemes, and re-emit sparse factor matrices.
	 */
	public static String rewriteJson(JsonNode root, NonCubicBilinearAlgorithm alg,
			FieldSweepResult sweep) {
		ObjectNode out = MAPPER.createObjectNode();
		// Drop legacy field tags; keep everything else.
		Set<String> drop = Set.of("field", "z2", "complex", "fields", "fields_not");
		for (Map.Entry<String, JsonNode> e : ((ObjectNode) root).properties()) {
			if (drop.contains(e.getKey())) continue;
			if (e.getKey().equals("u_sparse") || e.getKey().equals("v_sparse")
					|| e.getKey().equals("w_sparse")) {
				continue;  // re-emitted below
			}
			if (e.getKey().equals("u") || e.getKey().equals("v") || e.getKey().equals("w")) {
				continue;  // dense factors re-emitted via SchemeIO
			}
			out.set(e.getKey(), e.getValue());
		}
		// Insert the new field tags + commutative flag at the top alongside `n` / `m`.
		ArrayNode fieldsArr = MAPPER.createArrayNode();
		sweep.fields.forEach(fieldsArr::add);
		ArrayNode fieldsNotArr = MAPPER.createArrayNode();
		sweep.fieldsNot.forEach(fieldsNotArr::add);
		out.set("fields", fieldsArr);
		out.set("fields_not", fieldsNotArr);
		if (!out.has("commutative")) {
			out.put("commutative", false);  // bilinear (U,V,W) is structurally NC
		}
		// Emit the bilinear factors in the new sparse format. SchemeIO
		// already writes the new shape post-2026-06-03.
		String factorJson = SchemeIO.toJsonSparse(alg);
		JsonNode factorNode = MAPPER.readTree(factorJson);
		out.set("u_sparse", factorNode.get("u_sparse"));
		out.set("v_sparse", factorNode.get("v_sparse"));
		out.set("w_sparse", factorNode.get("w_sparse"));

		// Minified → MatrixJsonFormatter (matches the rest of the catalog's
		// matrix-friendly indentation, so diffs stay focused on the actual
		// schema change rather than whitespace churn).
		String minified;
		try {
			minified = MAPPER.writeValueAsString(out);
		} catch (Exception e) {
			throw new IllegalStateException(e);
		}
		try {
			// format() already appends the canonical trailing newline.
			return MatrixJsonFormatter.format(minified);
		} catch (IOException e) {
			throw new IllegalStateException("MatrixJsonFormatter failed", e);
		}
	}

	// Reserved: source-name → cmt flag heuristic (unused so far — bilinear
	// is uniformly NC). Will be wired when this driver gains non-bilinear support.
	private static final Pattern SOURCE_PREFIX =
			Pattern.compile("^([a-z][a-z0-9_\\-]*?)[-_]\\d");

	@SuppressWarnings("unused")
	private static boolean inferCommutativeFromSource(String filename) {
		Matcher m = SOURCE_PREFIX.matcher(filename);
		if (!m.find()) return false;
		String source = m.group(1).toLowerCase();
		return source.contains("waksman") || source.contains("rosowski")
				|| source.contains("makarov") || source.contains("islam")
				|| source.contains("probert") || source.contains("smith")
				|| source.contains("schachtel");
	}
}
