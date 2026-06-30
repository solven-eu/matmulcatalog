package eu.solven.matmul.docs.migrate;

import java.nio.file.Files;
import java.nio.file.Path;
import java.util.ArrayList;
import java.util.List;
import java.util.Optional;
import java.util.concurrent.atomic.AtomicInteger;
import java.util.stream.Collectors;

import eu.solven.matmul.NonCubicBilinearAlgorithm;
import eu.solven.matmul.verifiers.Verifier;
import eu.solven.matmul.algebra.Field;
import eu.solven.matmul.catalog.FieldAwareLookup;
import eu.solven.matmul.catalog.Lineage;
import eu.solven.matmul.catalog.MatrixJsonFormatter;
import eu.solven.matmul.catalog.SchemeIO;

import tools.jackson.databind.JsonNode;
import tools.jackson.databind.node.ArrayNode;
import tools.jackson.databind.node.ObjectNode;

/**
 * Backfill a {@code fields[]} (and, when content-derived, {@code fields_not[]})
 * onto every scheme JSON that lacks one — closing the data hole that let the
 * manifest invent a bogus narrowest field.
 *
 * <p>The canonical example is the derived ⟨3,3,8⟩=55 stub: it has no
 * {@code fields[]}, so {@link eu.solven.matmul.docs.generate.GenerateCatalogManifest}
 * fell back to {@code expandFieldCluster("R/Q/Z")} → {@code [Z,Q,R]}, claiming
 * <b>Z</b> even though the scheme holds {@code 1/8} coefficients and depends on
 * the Q-only {@code smirnov13} ⟨3,3,6⟩ leaf. Its correct field set is
 * {@code [F3,Q,R,C]} (matching the curated twin).</p>
 *
 * <h3>How fields are derived</h3>
 * <ul>
 *   <li><b>Matrix-bearing, non-complex, bilinear</b> → from the actual U/V/W
 *       <i>content</i>, using the same logic as
 *       {@code GenerateHk2npConstructed}: F₂ / F₃ via
 *       {@link Verifier#isExactNonCubicFp} (handles rational coefficients whose
 *       denominator is coprime to the prime — so {@code 1/8} grants F₃ but not
 *       F₂), Z iff every coefficient is an integer, and Q/R/C by inclusion for
 *       any verified characteristic-0 scheme. Content is GROUND TRUTH for the
 *       Z question — a {@code 1/8} coefficient definitively excludes Z no matter
 *       what the lineage suggests.</li>
 *   <li><b>Matrix-less stub / non-bilinear / complex</b> → fall back to LINEAGE
 *       inference ({@link FieldAwareLookup#inferFieldFromLineage}, no dense
 *       replay) and expand the narrowest field to its inclusion set — the same
 *       conservative path as {@link StampFields}.</li>
 * </ul>
 *
 * <p>Default DRY-RUN. {@code --execute} writes. Re-run the manifest afterwards.</p>
 *
 * <pre>mvn -q -ntp exec:java -Dexec.mainClass=eu.solven.matmul.docs.migrate.BackfillMissingFields [-Dexec.args=--execute]</pre>
 */
public final class BackfillMissingFields {
	private BackfillMissingFields() {}

	private static final List<String> ORDER = List.of("F2", "F3", "Z", "Q", "R", "C");
	private static final List<String> ALL = ORDER;

	public static void main(String[] args) throws Exception {
		boolean execute = List.of(args).contains("--execute");
		Path root = Path.of("src/main/resources/schemes");
		FieldAwareLookup lookup = new FieldAwareLookup("Q");

		List<Path> files;
		try (var s = Files.walk(root)) {
			files = s.filter(p -> p.toString().endsWith(".json")).collect(Collectors.toList());
		}
		System.out.println("Scanning " + files.size() + " scheme files (mode=" + (execute ? "EXECUTE" : "DRY-RUN") + ")…");

		AtomicInteger already = new AtomicInteger(), content = new AtomicInteger();
		AtomicInteger skipped = new AtomicInteger(), errors = new AtomicInteger();
		for (Path f : files) {
			try {
				JsonNode parsed = SchemeIO.parseJson(f.toFile());
				if (!SchemeIO.fieldTags(parsed).isEmpty()) { already.incrementAndGet(); continue; }
				if (!(parsed instanceof ObjectNode obj)) { errors.incrementAndGet(); continue; }

				boolean hasMatrices = parsed.has("u") || parsed.has("u_sparse");
				boolean contentEligible = hasMatrices && !SchemeIO.isComplex(parsed)
						&& !SchemeIO.isNonBilinear(parsed) && !SchemeIO.isStub(parsed);
				if (!contentEligible) {
					// NO FABRICATION (user 2026-06-12, "no fields → no fields"): only
					// stamp fields we can derive from ACTUAL coefficients. A matrix-less
					// stub / complex / non-bilinear scheme is left WITHOUT fields[] — we
					// NEVER write a lineage char-0 floor (it silently drops F₂/F₃) nor a
					// blind [Z] default (e.g. the smirnov13-backed ⟨3,3,13⟩=89 stub is
					// rational — 1/8 — and is NOT valid over Z). Materialise the stub
					// first if you want its real, content-derived fields.
					skipped.incrementAndGet();
					System.out.println("[skip]  " + f.getFileName() + " -> not content-eligible; left WITHOUT fields[]");
					continue;
				}
				NonCubicBilinearAlgorithm alg = SchemeIO.read(parsed);
				List<String> fields = fieldsFromContent(alg);
				List<String> fieldsNot = ALL.stream().filter(t -> !fields.contains(t)).collect(Collectors.toList());
				content.incrementAndGet();

				System.out.println((execute ? "[write] " : "[plan]  ") + f.getFileName()
						+ " -> fields=" + order(fields) + (fieldsNot != null ? " not=" + order(fieldsNot) : " (lineage)"));

				if (execute) {
					ArrayNode fa = obj.arrayNode();
					order(fields).forEach(fa::add);
					obj.set("fields", fa);
					if (fieldsNot != null) {
						ArrayNode fn = obj.arrayNode();
						order(fieldsNot).forEach(fn::add);
						obj.set("fields_not", fn);
					}
					Files.writeString(f, MatrixJsonFormatter.format(obj));
				}
			} catch (Exception e) {
				errors.incrementAndGet();
				System.out.println("[error] " + f.getFileName() + ": " + e);
			}
		}

		System.out.println("\n=== " + (execute ? "BACKFILLED" : "PLAN") + " ===");
		System.out.println("already had fields[]:                 " + already.get());
		System.out.println("stamped from CONTENT (matrices):      " + content.get());
		System.out.println("left WITHOUT fields[] (no content):   " + skipped.get());
		System.out.println("errors (skipped):                     " + errors.get());
		if (!execute) System.out.println("\n(DRY-RUN — pass --execute to write)");
	}

	/**
	 * Content-derived field set for a verified, non-complex bilinear scheme —
	 * the same policy as {@code GenerateHk2npConstructed}: a coefficient like
	 * {@code 1/8} excludes Z and F₂ but (being coprime to 3) still reduces into
	 * F₃, so {@code [F3,Q,R,C]}. F₂/F₃ are decided by the modular identity
	 * actually holding, NOT by an all-integer pre-gate.
	 */
	static List<String> fieldsFromContent(NonCubicBilinearAlgorithm alg) {
		List<String> fields = new ArrayList<>();
		if (Verifier.isExactNonCubicFp(alg, 2)) fields.add("F2");
		if (Verifier.isExactNonCubicFp(alg, 3)) fields.add("F3");
		if (allIntegers(alg)) fields.add("Z");
		// Any verified scheme with these (rational/real) coefficients is valid
		// over Q ⊆ R ⊆ C. Complex-only schemes are routed to the lineage branch.
		fields.add("Q");
		fields.add("R");
		fields.add("C");
		return fields;
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

	/** Narrowest field → inclusion-correct set (mirrors {@link StampFields#expand}). */
	private static List<String> expand(Field f) {
		switch (f.name()) {
			case "F2": return List.of("F2");
			case "F3": return List.of("F3");
			case "Z": return List.of("F2", "F3", "Z", "Q", "R", "C");
			case "Q": return List.of("Q", "R", "C");
			case "R": return List.of("R", "C");
			case "C": return List.of("C");
			default: throw new IllegalArgumentException("Unknown field, refusing to silently widen: " + f);
		}
	}

	private static List<String> order(List<String> fields) {
		return fields.stream().sorted((a, b) -> Integer.compare(
				idx(ORDER.indexOf(a)), idx(ORDER.indexOf(b)))).collect(Collectors.toList());
	}

	private static int idx(int i) { return i < 0 ? 99 : i; }
}
