package eu.solven.matmul.docs.migrate;

import java.io.File;
import java.io.IOException;
import java.nio.file.Files;
import java.util.Arrays;
import java.util.List;
import java.util.Map;

import eu.solven.matmul.NonCubicBilinearAlgorithm;
import eu.solven.matmul.catalog.MatrixJsonFormatter;
import eu.solven.matmul.catalog.SchemeIO;

import tools.jackson.databind.JsonNode;
import tools.jackson.databind.ObjectMapper;
import tools.jackson.databind.json.JsonMapper;
import tools.jackson.databind.node.ObjectNode;

/**
 * Sanitization procedure that stamps the {@code "zt"} flag into every scheme
 * JSON, so the catalog (and anything else) can simply READ it instead of
 * recomputing.
 *
 * <p><b>What ZT is.</b> ZT is NOT a field — it was historically, wrongly,
 * conflated with F₂/Z₂ (characteristic 2). ZT is the <i>sub-class of {@code Z}</i>
 * (integer) schemes whose every U/V/W coefficient is in {@code {-1,0,1}} —
 * "ternary integer" in Perminov's terminology. (Beware: {@code F3}/{@code Z3}
 * is also "ternary" but ternary <i>modular</i> over GF(3) — a different algebra.)</p>
 *
 * <p><b>What this does.</b> For each scheme:</p>
 * <ul>
 *   <li>If {@code Z} is among the scheme's verified {@code fields[]} AND the
 *       factor matrices are readable → set {@code "zt": true|false} from
 *       {@link SchemeIO#isTernary(NonCubicBilinearAlgorithm)}.</li>
 *   <li>If {@code Z} is NOT among the fields → remove any stale {@code "zt"}
 *       (the flag is only meaningful for integer schemes).</li>
 *   <li>If {@code Z} is present but the matrices can't be read (stub scheme
 *       with stripped U/V/W) → leave the file untouched and count it skipped.</li>
 * </ul>
 *
 * <p>The flag is inserted just before the first matrix (the first compound
 * array), keeping it in the metadata block. The file is rewritten through
 * {@link MatrixJsonFormatter} so it is also canonically formatted (primitives
 * inline, compound vertical) as a side effect. Idempotent.</p>
 *
 * <p>Run via the {@code scripts/mmc.sh sanitize zt} operation, or directly:
 * <pre>
 *   mvn -q test-compile
 *   java -cp target/classes:target/test-classes:$CLASSPATH \
 *        eu.solven.matmul.docs.migrate.MaterialiseZT [path]
 * </pre>
 * {@code path} defaults to {@code src/main/resources/schemes}; pass a single
 * file or sub-directory to scope the run.</p>
 */
public final class MaterialiseZT {

	private static final ObjectMapper MAPPER = JsonMapper.builder().build();

	private MaterialiseZT() {}

	private enum Outcome { TRUE, FALSE, CLEARED, SKIPPED, UNCHANGED }

	public static void main(String[] args) throws IOException {
		File target = args.length > 0 ? new File(args[0]) : new File("src/main/resources/schemes");

		File[] files;
		if (target.isFile()) {
			files = new File[] { target };
		} else if (target.isDirectory()) {
			try (var s = Files.walk(target.toPath())) {
				files = s.filter(p -> p.toString().endsWith(".json"))
						.map(java.nio.file.Path::toFile)
						.toArray(File[]::new);
			}
			Arrays.sort(files);
		} else {
			System.err.println("not a file or directory: " + target);
			System.exit(1);
			return;
		}

		long t0 = System.currentTimeMillis();
		int stampedTrue = 0, stampedFalse = 0, cleared = 0, skipped = 0, errors = 0, unchanged = 0;
		for (int i = 0; i < files.length; i++) {
			File f = files[i];
			try {
				switch (process(f)) {
				case TRUE -> stampedTrue++;
				case FALSE -> stampedFalse++;
				case CLEARED -> cleared++;
				case SKIPPED -> skipped++;
				case UNCHANGED -> unchanged++;
				}
			} catch (Exception e) {
				errors++;
				System.err.printf("ERROR %s: %s%n", f.getName(), e.getMessage());
			}
			if (files.length > 200 && (i + 1) % 1000 == 0) {
				long ms = System.currentTimeMillis() - t0;
				System.out.printf("[progress] %d/%d processed (zt+ %d, zt- %d, skip %d) %dms elapsed%n",
						i + 1, files.length, stampedTrue, stampedFalse, skipped, ms);
			}
		}

		System.out.printf(
				"MaterialiseZT: zt:true=%d, zt:false=%d, cleared=%d, skipped(stub/unreadable)=%d, "
						+ "already-canonical=%d, errors=%d (total %d files, %dms)%n",
				stampedTrue, stampedFalse, cleared, skipped, unchanged, errors, files.length,
				System.currentTimeMillis() - t0);
		if (errors > 0) System.exit(2);
	}

	private static Outcome process(File f) throws IOException {
		String before = Files.readString(f.toPath());
		JsonNode root = SchemeIO.parseJson(f);
		if (!(root instanceof ObjectNode obj)) {
			return Outcome.SKIPPED;  // not a scheme object
		}

		List<String> tags = SchemeIO.fieldTags(root);
		boolean hasZ = tags.contains("Z");

		Boolean zt;
		Outcome intent;
		if (!hasZ) {
			// Not an integer scheme: the flag is meaningless — drop any stale one.
			zt = null;
			intent = obj.has("zt") ? Outcome.CLEARED : Outcome.UNCHANGED;
		} else {
			NonCubicBilinearAlgorithm alg;
			try {
				alg = SchemeIO.isReduced(root) ? SchemeIO.readReduced(root) : SchemeIO.read(root);
			} catch (Exception readFailed) {
				// Stub / stripped matrices: can't compute — leave the file alone.
				return Outcome.SKIPPED;
			}
			zt = SchemeIO.isTernary(alg);
			intent = zt ? Outcome.TRUE : Outcome.FALSE;
		}

		String after = MatrixJsonFormatter.format(withZt(obj, zt));
		if (after.equals(before)) {
			return Outcome.UNCHANGED;
		}
		Files.writeString(f.toPath(), after);
		// A pure reformat (zt unchanged) still rewrites the file but isn't a
		// zt change — report it as UNCHANGED for zt-accounting purposes only if
		// the flag itself didn't move; otherwise report the zt intent.
		return intent == Outcome.UNCHANGED ? Outcome.UNCHANGED : intent;
	}

	/**
	 * Return a copy of {@code obj} with {@code "zt"} set to {@code zt} (when
	 * non-null) inserted just before the first compound array (the first
	 * matrix), or with any existing {@code "zt"} removed (when {@code zt} is
	 * null). Property order is otherwise preserved.
	 */
	private static ObjectNode withZt(ObjectNode obj, Boolean zt) {
		ObjectNode out = MAPPER.createObjectNode();
		boolean inserted = false;
		for (Map.Entry<String, JsonNode> e : obj.properties()) {
			String k = e.getKey();
			if (k.equals("zt")) {
				continue;  // drop the old one; we re-insert in the canonical slot
			}
			if (!inserted && zt != null && isCompoundArray(e.getValue())) {
				out.put("zt", zt.booleanValue());
				inserted = true;
			}
			out.set(k, e.getValue());
		}
		if (!inserted && zt != null) {
			out.put("zt", zt.booleanValue());  // no matrix found — append at end
		}
		return out;
	}

	/** An array that holds at least one array/object element (i.e. a matrix). */
	private static boolean isCompoundArray(JsonNode node) {
		if (node == null || !node.isArray()) return false;
		for (JsonNode e : node) {
			if (e.isArray() || e.isObject()) return true;
		}
		return false;
	}
}
