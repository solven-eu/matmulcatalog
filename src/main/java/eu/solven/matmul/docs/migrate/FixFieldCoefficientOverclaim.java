package eu.solven.matmul.docs.migrate;

import java.nio.file.Files;
import java.nio.file.Path;
import java.util.LinkedHashSet;
import java.util.List;
import java.util.Map;
import java.util.Set;
import java.util.concurrent.atomic.AtomicInteger;
import java.util.stream.Stream;

import eu.solven.matmul.catalog.SchemeIO;
import lombok.extern.slf4j.Slf4j;
import tools.jackson.databind.JsonNode;

/**
 * Strips field tags that a scheme's coefficient denominators make impossible —
 * the recurring "{@code 1/8} coefficient tagged {@code Z}/{@code F2}" over-claim
 * (see {@link SchemeIO#fieldsContradictedByCoefficients}). The materialiser's
 * "integer base ⇒ all fields" stamping mechanically emitted
 * {@code fields=[F2,F3,Z,Q,R,C]} on derived schemes that picked up rational
 * (½-polarization / recombination) coefficients; the narrowest correct field is
 * {@code Q}, not {@code Z}, and {@code F2} is invalid whenever a denominator is
 * even.
 *
 * <p>The fix is deterministic and verification-free: remove exactly the tags the
 * <em>coefficients</em> contradict ({@code Z} when any coefficient is rational,
 * {@code F2} when any denominator is even, {@code F3} when any denominator is
 * divisible by 3) and keep the rest ({@code F3}/{@code Q}/{@code R}/{@code C}
 * survive {@code 1/8}). Only {@code fields[]} is rewritten — a minimal textual
 * patch that preserves the canonical formatting and the U/V/W matrices verbatim.</p>
 *
 * <p>This is the sibling repair to {@link FixRationalFields} / {@link NarrowFields}
 * but driven purely by the coefficient denominators (not a lineage base tag, which
 * is what mis-expanded in the first place). Re-runnable; rewrites a file only when
 * {@code fields[]} actually changes.</p>
 *
 * <pre>mvn -q -ntp exec:java -Dexec.mainClass=eu.solven.matmul.docs.migrate.FixFieldCoefficientOverclaim -Dexec.args="--apply"</pre>
 */
@Slf4j
public final class FixFieldCoefficientOverclaim {

	private static final Path SCHEMES_ROOT = Path.of("src/main/resources/schemes");

	private FixFieldCoefficientOverclaim() {}

	public static void main(String[] args) throws Exception {
		boolean apply = List.of(args).contains("--apply");
		List<Path> files;
		try (Stream<Path> walk = Files.walk(SCHEMES_ROOT)) {
			files = walk.filter(p -> p.toString().endsWith(".json")).toList();
		}
		AtomicInteger scanned = new AtomicInteger();
		AtomicInteger fixed = new AtomicInteger();
		AtomicInteger errors = new AtomicInteger();

		for (Path file : files) {
			scanned.incrementAndGet();
			try {
				JsonNode root = SchemeIO.parseJson(file.toFile());
				List<String> contradicted = contradictedTags(root);
				if (contradicted.isEmpty()) {
					continue;
				}
				List<String> tags = SchemeIO.fieldTags(root);
				Set<String> kept = new LinkedHashSet<>(tags);
				kept.removeAll(contradicted);
				if (kept.equals(new LinkedHashSet<>(tags))) {
					continue;  // nothing actually removed
				}

				// Canonical, drift-free re-emit: remove the old fields[] then re-add
				// the narrowed list (updateFields never clobbers, so the remove must
				// run first). Goes through MatrixJsonFormatter — preserves formatting
				// and the U/V/W matrices verbatim, and notifies FieldAwareLookup.
				log.info("{} : {} → {} (drop {})",
						file.getFileName(), tags, kept, contradicted);
				boolean changed = SchemeIO.updateFields(file.toFile(),
						Map.of("fields", List.copyOf(kept)), List.of("fields"), apply);
				if (changed) {
					fixed.incrementAndGet();
				}
			} catch (Exception e) {
				errors.incrementAndGet();
				log.warn("error on {}: {}", file, e.toString());
			}
		}

		log.info("FixFieldCoefficientOverclaim: scanned={} {}={} errors={} apply={}",
				scanned.get(), apply ? "fixed" : "would-fix", fixed.get(), errors.get(), apply);
		if (!apply && fixed.get() > 0) {
			log.info("Dry run — re-run with --apply to rewrite.");
		}
	}

	/** The field tags contradicted by this scheme's coefficient denominators, de-duplicated. */
	private static List<String> contradictedTags(JsonNode root) {
		Set<String> bad = new LinkedHashSet<>();
		for (String reason : SchemeIO.fieldsContradictedByCoefficients(root)) {
			// reasons read "declares <TAG> but …" — recover the tag token.
			int i = reason.indexOf("declares ");
			if (i < 0) continue;
			String rest = reason.substring(i + "declares ".length());
			int sp = rest.indexOf(' ');
			bad.add(sp < 0 ? rest : rest.substring(0, sp));
		}
		return List.copyOf(bad);
	}
}
