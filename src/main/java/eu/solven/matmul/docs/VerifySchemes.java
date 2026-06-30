package eu.solven.matmul.docs;

import java.nio.file.Files;
import java.nio.file.Path;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.LinkedHashSet;
import java.util.List;
import java.util.Locale;
import java.util.Optional;
import java.util.Set;
import java.util.concurrent.atomic.AtomicInteger;
import java.util.stream.Collectors;

import eu.solven.matmul.NonCubicBilinearAlgorithm;
import eu.solven.matmul.verifiers.Verifier;
import eu.solven.matmul.algebra.Field;
import eu.solven.matmul.catalog.FieldAwareLookup;
import eu.solven.matmul.catalog.SchemeIO;
import eu.solven.matmul.docs.verify.VerifyLineageFieldCompat;

import lombok.extern.slf4j.Slf4j;
import tools.jackson.databind.JsonNode;

/**
 * Single CLI entry-point for VERIFYING catalog schemes — the read-side
 * counterpart of {@link SchemeSweep} (which GENERATES schemes). One driver,
 * three independent checks selectable via {@code --check=}, cheapest first:
 *
 * <ol>
 *   <li><b>{@code lineage}</b> (cheap, NO dense) — a derived scheme must not
 *       declare a characteristic-0 field <i>narrower</i> than its lineage atoms
 *       support ({@code Z < Q < R < C}). This is the check that catches the
 *       derived ⟨3,3,8⟩=55 "declares Z but its smirnov13 ⟨3,3,6⟩ leaf is Q-only"
 *       over-claim. Pure metadata: reads each atom's own {@code fields[]}, never
 *       materializes a tensor. Delegates to {@link VerifyLineageFieldCompat#check}.</li>
 *   <li><b>{@code coefficients}</b> (cheap, content) — for a scheme that carries
 *       explicit matrices (an atom / materialized scheme), the declared fields
 *       must be consistent with the EXACT coefficient denominators: a declared
 *       {@code Z} demands all-integer entries (a {@code 1/8} coefficient
 *       disproves it); declared {@code F2} demands odd denominators and
 *       {@code F3} demands denominators coprime to 3 ({@code 1/8} ∈ F₃ but ∉ F₂).
 *       These are NECESSARY conditions read straight off the raw tokens via
 *       {@link SchemeIO#fieldsContradictedByCoefficients} (no dense
 *       materialization); the exact modular identity is the {@code full} check.</li>
 *   <li><b>{@code full}</b> (EXPENSIVE, dense) — the scheme must actually COMPUTE
 *       matmul over each declared field: an exact symbolic residual for
 *       F₂/F₃ ({@link Verifier#isExactNonCubicF2}/{@code F3}) and a random
 *       matmul spot-check for characteristic-0 fields
 *       ({@link Verifier#passesRandomMatmulSpotCheck}). This subsumes what
 *       {@code VerifyAllSchemes} does.</li>
 * </ol>
 *
 * <h3>CLI</h3>
 * <pre>
 *   --check=lineage,coefficients       checks to run (default: this pair — the cheap ones)
 *   --check=full                       add the dense identity verification
 *   --check=all                        lineage,coefficients,full
 *   --field=Z|Q|R|C|F2|F3              restrict to schemes declaring this field
 *   --shape=NxMxP[,NxMxP…]             restrict to these shapes (sorted-agnostic)
 *   --schemes-root=PATH                default src/main/resources/schemes
 *   --max-report=N                     cap printed failures (default 50)
 * </pre>
 *
 * <pre>mvn -q -ntp exec:java -Dexec.mainClass=eu.solven.matmul.docs.VerifySchemes -Dexec.args="--check=lineage,coefficients"</pre>
 *
 * Exit code 1 if any selected check fails on any scheme (CI merge gate).
 */
@Slf4j
public final class VerifySchemes {
	private VerifySchemes() {}

	enum Check { LINEAGE, COEFFICIENTS, FULL }

	public static void main(String[] args) throws Exception {
		List<String> argv = Arrays.asList(args);
		Set<Check> checks = parseChecks(arg(argv, "--check", "lineage,coefficients"));
		Optional<Field> fieldFilter = parseField(arg(argv, "--field", null));
		Set<String> shapeFilter = parseShapes(arg(argv, "--shape", null));
		Path root = Path.of(arg(argv, "--schemes-root", "src/main/resources/schemes"));
		int maxReport = Integer.parseInt(arg(argv, "--max-report", "50"));

		FieldAwareLookup lookup = new FieldAwareLookup("Q");
		List<Path> files;
		try (var s = Files.walk(root)) {
			files = s.filter(p -> p.toString().endsWith(".json")).collect(Collectors.toList());
		}
		log.info("VerifySchemes: {} files under {} | checks={} field={} shapes={}",
				files.size(), root, checks, fieldFilter.map(Enum::name).orElse("*"),
				shapeFilter.isEmpty() ? "*" : shapeFilter);

		List<String> failures = new ArrayList<>();
		AtomicInteger scanned = new AtomicInteger(), skipped = new AtomicInteger();
		long t0 = System.nanoTime();
		int processed = 0;
		for (Path f : files) {
			try {
				JsonNode root2 = SchemeIO.parseJson(f.toFile());
				if (!matchesShape(root2, shapeFilter)) { skipped.incrementAndGet(); continue; }
				List<String> declared = SchemeIO.fieldTags(root2);
				if (fieldFilter.isPresent() && !declared.contains(fieldFilter.get().name())) {
					skipped.incrementAndGet();
					continue;
				}
				scanned.incrementAndGet();
				failures.addAll(checkOne(f, root2, declared, checks, lookup));
			} catch (java.io.FileNotFoundException | java.nio.file.NoSuchFileException e) {
				// File vanished between the directory walk and the read (e.g. a
				// regenerator deleting an old emission). Not a verification failure.
				skipped.incrementAndGet();
			} catch (Exception e) {
				failures.add(f.getFileName() + " [UNREADABLE] " + e);
			}
			if (++processed % 2000 == 0) {
				log.info("[progress] {}/{} processed, {} failure(s), {}ms",
						processed, files.size(), failures.size(), (System.nanoTime() - t0) / 1_000_000);
			}
		}

		log.info("=== VerifySchemes done: scanned={} skipped={} failures={} ===",
				scanned.get(), skipped.get(), failures.size());
		failures.stream().limit(maxReport).forEach(m -> log.error("  [FAIL] {}", m));
		if (failures.size() > maxReport) {
			log.error("  … and {} more (raise --max-report to see all)", failures.size() - maxReport);
		}
		if (!failures.isEmpty()) {
			System.exit(1);
		}
		log.info("OK: all selected checks passed.");
	}

	/** Run the selected checks on one scheme; returns the (possibly empty) failure messages. */
	static List<String> checkOne(Path f, JsonNode root, List<String> declared,
			Set<Check> checks, FieldAwareLookup lookup) {
		List<String> out = new ArrayList<>();
		String name = f.getFileName().toString();

		// (1) lineage field-compat — derived scheme vs its atoms (no dense).
		if (checks.contains(Check.LINEAGE)) {
			Optional<VerifyLineageFieldCompat.Violation> v =
					VerifyLineageFieldCompat.check(f, root, lookup);
			v.ifPresent(viol -> out.add("lineage: " + viol));
		}

		// (2) coefficients — declared field vs actual U/V/W content (matrix-bearing only).
		boolean hasMatrices = root.has("u") || root.has("u_sparse");
		boolean bilinear = hasMatrices && !SchemeIO.isComplex(root)
				&& !SchemeIO.isNonBilinear(root) && !SchemeIO.isReduced(root);
		if (checks.contains(Check.COEFFICIENTS) && bilinear) {
			// EXACT, cheap necessary-condition gate read straight off the raw
			// coefficient tokens (no dense materialization): a declared Z must be
			// all-integer, F2 must have odd denominators, F3 denominators coprime
			// to 3. This is the over-claim that has bitten the catalog repeatedly
			// (a 1/8 coefficient tagged Z/F2). See
			// SchemeIO.fieldsContradictedByCoefficients.
			for (String why : SchemeIO.fieldsContradictedByCoefficients(root)) {
				out.add(name + " coefficients: " + why);
			}
		}
		// (3) full identity verification per declared field (dense, expensive).
		if (checks.contains(Check.FULL) && bilinear) {
			NonCubicBilinearAlgorithm alg;
			try {
				alg = SchemeIO.read(root);
			} catch (Exception e) {
				out.add(name + " full: cannot read matrices: " + e);
				return out;
			}
			String why = fullVerifyViolation(alg, declared);
			if (why != null) out.add(name + " full: " + why);
		}
		return out;
	}

	/** Dense identity verification over each declared field. Returns a reason or null. */
	private static String fullVerifyViolation(NonCubicBilinearAlgorithm alg, List<String> declared) {
		boolean char0Checked = false;
		for (String tag : declared) {
			switch (tag) {
				case "F2":
					// Fp (not the integer-only isExactNonCubicF2) so a rational coprime
					// to 2 is judged on the actual modular identity, not its denominator.
					if (!Verifier.isExactNonCubicFp(alg, 2)) return "F2 identity fails";
					break;
				case "F3":
					if (!Verifier.isExactNonCubicFp(alg, 3)) return "F3 identity fails";
					break;
				case "Z": case "Q": case "R": case "C":
					// One char-0 spot-check covers the whole Z⊂Q⊂R⊂C chain.
					if (!char0Checked) {
						if (!Verifier.passesRandomMatmulSpotCheck(alg)) return "char-0 matmul spot-check fails";
						char0Checked = true;
					}
					break;
				default:
					break;  // ZT and friends: not a verification field
			}
		}
		return null;
	}

	// ── arg parsing ─────────────────────────────────────────────────────────

	private static Set<Check> parseChecks(String spec) {
		Set<Check> out = new LinkedHashSet<>();
		for (String tok : spec.split(",")) {
			switch (tok.trim().toLowerCase(Locale.ROOT)) {
				case "lineage" -> out.add(Check.LINEAGE);
				case "coefficients", "coeff", "coeffs" -> out.add(Check.COEFFICIENTS);
				case "full", "identity" -> out.add(Check.FULL);
				case "all" -> { out.add(Check.LINEAGE); out.add(Check.COEFFICIENTS); out.add(Check.FULL); }
				case "" -> { }
				default -> throw new IllegalArgumentException("Unknown --check token: " + tok);
			}
		}
		if (out.isEmpty()) throw new IllegalArgumentException("--check selected no checks");
		return out;
	}

	private static Optional<Field> parseField(String s) {
		return s == null ? Optional.empty() : Optional.of(Field.valueOf(s.trim()));
	}

	private static Set<String> parseShapes(String s) {
		if (s == null) return Set.of();
		Set<String> out = new LinkedHashSet<>();
		for (String tok : s.split(",")) {
			int[] d = Arrays.stream(tok.trim().split("x")).mapToInt(Integer::parseInt).toArray();
			Arrays.sort(d);
			out.add(d[0] + "x" + d[1] + "x" + d[2]);
		}
		return out;
	}

	private static boolean matchesShape(JsonNode root, Set<String> shapeFilter) {
		if (shapeFilter.isEmpty()) return true;
		JsonNode n = root.get("n");
		if (n == null || !n.isArray() || n.size() != 3) return false;
		int[] d = { n.get(0).asInt(), n.get(1).asInt(), n.get(2).asInt() };
		Arrays.sort(d);
		return shapeFilter.contains(d[0] + "x" + d[1] + "x" + d[2]);
	}

	private static String arg(List<String> argv, String key, String dflt) {
		for (String a : argv) {
			if (a.startsWith(key + "=")) return a.substring(key.length() + 1);
		}
		return dflt;
	}
}
