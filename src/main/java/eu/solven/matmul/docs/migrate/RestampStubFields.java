package eu.solven.matmul.docs.migrate;

import java.nio.file.Files;
import java.nio.file.Path;
import java.util.ArrayList;
import java.util.List;
import java.util.Locale;
import java.util.Map;
import java.util.Set;

import eu.solven.matmul.NonCubicBilinearAlgorithm;
import eu.solven.matmul.verifiers.Verifier;
import eu.solven.matmul.catalog.FieldAwareLookup;
import eu.solven.matmul.catalog.SchemeIO;
import eu.solven.matmul.search.LineageReplayer;
import lombok.extern.slf4j.Slf4j;
import tools.jackson.databind.JsonNode;

/**
 * Recompute a stub's {@code fields[]} from its REPLAYED coefficients — repairs the
 * field drift left by field-scoped sweeps. Canonical case: {@code 17x19x20-r3780}
 * was stamped {@code [R, C]} by an R-scoped run although its replay is a
 * composition of Q/Z atoms — invisible to Q lookups until restamped (see the
 * "Sweeps default to --field=Q" rule in CLAUDE.md).
 *
 * <p>For each stub file at {@code --shape=NxMxP[,…]}: replay the stored lineage →
 * dense U/V/W → classify by the actual coefficients:</p>
 * <ul>
 *   <li><b>integer</b> entries + char-0 spot-check → {@code Z, Q, R, C}, plus
 *       {@code F2}/{@code F3} when the corresponding mod-p spot-check passes
 *       (each modular claim is spot-checked on its own, not granted by theorem);</li>
 *   <li><b>non-integer</b> entries + char-0 spot-check → {@code Q, R, C}, plus
 *       {@code F2}/{@code F3} per mod-p spot-check (denominators permitting —
 *       an unrepresentable scheme simply fails the check).</li>
 * </ul>
 *
 * <p>Spot-check-based, matching the bar {@code VerifySchemes --check=full} applies
 * to maxDim&gt;16 stubs (an exact symbolic verify at these dims is prohibitive) —
 * so the stamp is "verified (randomised)", the same policy as the enrich pipeline.
 * Doubles ARE rationals, so char-0-valid always implies at least {@code [Q,R,C]}:
 * a replayable stub can never legitimately be R-without-Q.</p>
 *
 * <p>Dry-run by default (reports would-change); {@code --apply} rewrites
 * {@code fields[]} through the canonical formatter. Re-runnable whenever drift
 * recurs.</p>
 *
 * <pre>mvn -q -ntp exec:java -Dexec.mainClass=eu.solven.matmul.docs.migrate.RestampStubFields \
 *     -Dexec.args="--shape=17x19x20 --apply"</pre>
 */
@Slf4j
public final class RestampStubFields {

	private static final Path SCHEMES_ROOT = Path.of("src/main/resources/schemes");
	/** Canonical fields[] order (mirrors StampFields.ORDER). */
	private static final List<String> ORDER = List.of("F2", "F3", "Z", "Q", "R", "C");

	private RestampStubFields() {}

	public static void main(String[] args) throws Exception {
		boolean apply = List.of(args).contains("--apply");
		Set<String> shapes = new java.util.LinkedHashSet<>();
		for (String a : args) {
			if (a.startsWith("--shape=")) {
				for (String tok : a.substring(8).split(",")) {
					if (!tok.isBlank()) shapes.add(canon(tok.trim()));
				}
			}
		}
		if (shapes.isEmpty()) {
			throw new IllegalArgumentException("--shape=NxMxP[,NxMxP…] is mandatory "
					+ "(this is a targeted repair tool, not a tree-wide stamper).");
		}

		// Widest char-0 lookup (C admits Z⊂Q⊂R⊂C): the repair tool must resolve
		// leaves of R/C-drifted stubs, so this is a DOCUMENTED exception to the
		// Q-by-default sweep rule — nothing is searched here, and the resulting
		// stamp comes from the coefficients, not from this lookup's field.
		FieldAwareLookup lookup = new FieldAwareLookup("C");
		LineageReplayer replayer = LineageReplayer.withDefaultPool(lookup);

		List<Path> files;
		try (var s = Files.walk(SCHEMES_ROOT)) {
			files = s.filter(p -> p.toString().endsWith(".json")).sorted().toList();
		}
		int scanned = 0, changed = 0, errors = 0;
		for (Path f : files) {
			JsonNode root = SchemeIO.parseJson(f.toFile());
			JsonNode nArr = root.get("n");
			if (nArr == null || nArr.size() != 3) continue;
			String shape = nArr.get(0).asInt() + "x" + nArr.get(1).asInt() + "x" + nArr.get(2).asInt();
			if (!shapes.contains(canon(shape))) continue;
			if (!SchemeIO.isStub(root)) continue;   // explicit-matrix files: NarrowFields' turf
			scanned++;
			try {
				NonCubicBilinearAlgorithm alg = replayer.replayFromFile(f.toFile());
				if (alg.r != root.get("m").asInt()) {
					log.warn("{}: replayed rank {} ≠ stamped m={} — rank drift, NOT restamping",
							f.getFileName(), alg.r, root.get("m").asInt());
					errors++;
					continue;
				}
				List<String> fields = coefficientFields(alg);
				if (fields.isEmpty()) {
					log.warn("{}: replay does not compute matmul over char-0 — NOT restamping",
							f.getFileName());
					errors++;
					continue;
				}
				List<String> cur = SchemeIO.fieldTags(root);
				if (new java.util.HashSet<>(cur).equals(new java.util.HashSet<>(fields))) {
					log.info("{}: fields already correct {}", f.getFileName(), cur);
					continue;
				}
				changed++;
				log.info("{}: fields {} → {}{}", f.getFileName(), cur, fields,
						apply ? " (APPLIED)" : " (dry-run)");
				if (apply) {
					SchemeIO.updateFields(f.toFile(), Map.of("fields", fields),
							List.of("fields"), true);
				}
			} catch (RuntimeException e) {
				log.warn("{}: replay failed ({}) — NOT restamping", f.getFileName(), e.getMessage());
				errors++;
			}
		}
		log.info("RestampStubFields: scanned={} changed={} errors={} apply={}",
				scanned, changed, errors, apply);
	}

	/** Field set from the replayed coefficients — every claim spot-checked. */
	static List<String> coefficientFields(NonCubicBilinearAlgorithm alg) {
		if (!Verifier.passesRandomMatmulSpotCheck(alg)) {
			return List.of();
		}
		List<String> out = new ArrayList<>();
		if (Verifier.passesRandomMatmulSpotCheckFp(alg, 2)) out.add("F2");
		if (Verifier.passesRandomMatmulSpotCheckFp(alg, 3)) out.add("F3");
		if (allIntegers(alg)) out.add("Z");
		out.add("Q");
		out.add("R");
		out.add("C");
		// Already emitted in ORDER; assert the invariant cheaply in case of edits.
		assert isOrdered(out) : "fields not in canonical order: " + out;
		return out;
	}

	private static boolean isOrdered(List<String> fields) {
		int last = -1;
		for (String f : fields) {
			int i = ORDER.indexOf(f);
			if (i < last) return false;
			last = i;
		}
		return true;
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

	/** Sorted canonical "NxMxP" (orientation-agnostic match). */
	private static String canon(String shape) {
		int[] d = java.util.stream.Stream.of(shape.toLowerCase(Locale.ROOT).split("x"))
				.mapToInt(Integer::parseInt).toArray();
		java.util.Arrays.sort(d);
		return d[0] + "x" + d[1] + "x" + d[2];
	}
}
