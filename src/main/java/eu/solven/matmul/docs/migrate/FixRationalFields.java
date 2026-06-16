package eu.solven.matmul.docs.migrate;

import eu.solven.matmul.docs.verify.SanityCatalogMigration;

import java.nio.file.Files;
import java.nio.file.Path;
import java.util.List;
import java.util.concurrent.atomic.AtomicInteger;
import java.util.stream.Stream;

import eu.solven.matmul.NonCubicBilinearAlgorithm;
import eu.solven.matmul.catalog.MatrixJsonFormatter;
import eu.solven.matmul.catalog.SchemeIO;
import lombok.extern.slf4j.Slf4j;
import tools.jackson.databind.JsonNode;
import tools.jackson.databind.json.JsonMapper;
import tools.jackson.databind.node.ArrayNode;
import tools.jackson.databind.node.ObjectNode;

/**
 * Repairs a field-tag bug introduced by the structural #174 migration: schemes
 * whose field lived only in the FILENAME (e.g. {@code dis09_Q-…}, {@code
 * perminov_Q-…}, {@code fmm_lille-…} with ½/⅛ coefficients, AlphaEvolve, DPS)
 * were tagged {@code fields=["Z"]} because the migration read only the JSON
 * body (which had no {@code field} tag) and defaulted to Z.
 *
 * <p>A scheme with a non-integer (rational) coefficient <strong>cannot</strong>
 * be in Z / F₂ / F₃. This driver inspects the actual U/V/W coefficients and,
 * when a non-integer is present, re-tags the scheme as {@code fields=["Q","R","C"]}
 * / {@code fields_not=["F2","F3","Z"]}. Integer-coefficient schemes (correctly
 * tagged Z or F2 already) are left untouched.</p>
 *
 * <p><strong>Complex schemes</strong> (authoritative {@code complex:true} flag)
 * are valid <em>only</em> over C. Some AlphaEvolve imports — including the
 * canonical C-only ⟨4,4,4⟩=48 — were mis-tagged {@code fields=["Z"]}, and
 * {@code SchemeIO.isComplex} reads {@code fields[]} (not the flag) so it failed to
 * catch them. We re-tag those to {@code fields=["C"]} /
 * {@code fields_not=["F2","F3","Z","Q","R"]}, rewriting the field tags ONLY and
 * preserving the {@code [re,im]} complex matrices verbatim.</p>
 *
 * <p>Identity correctness is NOT re-checked here — these are already-verified
 * catalog schemes; only their field <em>membership</em> tag was wrong. Files are
 * rewritten only when {@code fields[]} actually changes, keeping the diff minimal.</p>
 *
 * <pre>mvn -q exec:java -Dexec.mainClass=eu.solven.matmul.docs.migrate.FixRationalFields -Dexec.args="--apply"</pre>
 */
@Slf4j
public final class FixRationalFields {

	private static final Path SCHEMES_ROOT = Path.of("src/main/resources/schemes");

	private FixRationalFields() {}

	public static void main(String[] args) throws Exception {
		boolean apply = List.of(args).contains("--apply");
		List<Path> files;
		try (Stream<Path> walk = Files.walk(SCHEMES_ROOT)) {
			files = walk.filter(p -> p.toString().endsWith(".json")).toList();
		}
		AtomicInteger scanned = new AtomicInteger();
		AtomicInteger retagged = new AtomicInteger();
		AtomicInteger retaggedComplex = new AtomicInteger();
		AtomicInteger errors = new AtomicInteger();

		files.parallelStream().forEach(file -> {
			scanned.incrementAndGet();
			try {
				JsonNode root = SchemeIO.parseJson(file.toFile());
				if (SchemeIO.isStub(root) || SchemeIO.isNonBilinear(root)) {
					return;
				}
				if (!root.has("u_sparse") && !root.has("u")) return;
				List<String> tags = SchemeIO.fieldTags(root);

				// Complex schemes (authoritative `complex` flag): a complex factor
				// has a non-integer/non-real coefficient, so it is valid ONLY over C.
				// Some AlphaEvolve imports (incl. the canonical C-only ⟨4,4,4⟩=48)
				// were mis-tagged fields=["Z"]; SchemeIO.isComplex reads fields[],
				// not the flag, so it failed to catch them. Rewrite field-tags ONLY,
				// preserving the [re,im] complex matrices verbatim.
				if (root.path("complex").asBoolean(false)) {
					boolean overclaims = tags.stream().anyMatch(t -> !t.equals("C"));
					if (!overclaims) return;  // already [C]
					ObjectNode obj = (ObjectNode) root;
					ArrayNode f = MAPPER.createArrayNode();
					f.add("C");
					ArrayNode fn = MAPPER.createArrayNode();
					List.of("F2", "F3", "Z", "Q", "R").forEach(fn::add);
					obj.set("fields", f);
					obj.set("fields_not", fn);
					if (apply) {
						// format() already appends the canonical trailing newline.
						Files.writeString(file, MatrixJsonFormatter.format(obj));
					}
					retaggedComplex.incrementAndGet();
					log.info("retag {} (complex) -> [C]", file.getFileName());
					return;
				}

				// Only touch schemes currently claiming a characteristic-strict tag.
				boolean claimsStrict = tags.contains("Z") || tags.contains("F2") || tags.contains("F3");
				if (!claimsStrict) return;

				NonCubicBilinearAlgorithm alg =
						SchemeIO.isReduced(root) ? SchemeIO.readReduced(root) : SchemeIO.read(root);
				if (allIntegers(alg)) return;  // genuinely integer → leave as-is

				// Non-integer rational coefficients → cannot be Z/F2/F3.
				String migrated = SanityCatalogMigration.rewriteJson(
						root, alg,
						new SanityCatalogMigration.FieldSweepResult(
								List.of("Q", "R", "C"),
								List.of("F2", "F3", "Z"),
								java.util.Map.of()));
				if (apply) {
					Files.writeString(file, migrated);
				}
				retagged.incrementAndGet();
				log.info("retag {} -> [Q,R,C]", file.getFileName());
			} catch (Exception e) {
				errors.incrementAndGet();
				log.warn("failed on {}: {}", file.getFileName(), e.getMessage());
			}
		});

		log.info("FixRationalFields: scanned={} retagged(rational)={} retagged(complex)={} errors={} (apply={})",
				scanned.get(), retagged.get(), retaggedComplex.get(), errors.get(), apply);
	}

	private static final JsonMapper MAPPER = JsonMapper.builder().build();

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
}
