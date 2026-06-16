package eu.solven.matmul.docs.verify;

import java.nio.file.Path;
import java.util.List;
import java.util.Optional;

import eu.solven.matmul.algebra.Field;
import eu.solven.matmul.catalog.FieldAwareLookup;
import eu.solven.matmul.catalog.Lineage;
import eu.solven.matmul.catalog.SchemeIO;

import tools.jackson.databind.JsonNode;

/**
 * CI guard: a derived scheme must not declare a characteristic-0 field
 * <i>narrower</i> than the field its lineage atoms can actually support — a
 * <b>lightweight, no-dense-materialization</b> check.
 *
 * <h3>The bug this catches</h3>
 * <p>The derived ⟨3,3,8⟩=55 stub declared {@code [Z,Q,R]} (claiming <b>Z</b>)
 * while it is built by concatenating an alphatensor ⟨3,3,2⟩ (Z) with the
 * <b>Q-only</b> {@code smirnov13} ⟨3,3,6⟩ (which holds {@code 1/8}
 * coefficients). The composed scheme lives in the smallest field containing
 * both leaves — Q, not Z — so the Z claim is false. Nothing checked this:
 * {@code VerifyAllSchemes} only verifies the bilinear <i>identity</i>, which
 * holds over every characteristic-0 field regardless of whether the
 * coefficients are integers, so it can never disprove a bogus Z.</p>
 *
 * <h3>What it checks (and what it deliberately does NOT)</h3>
 * <p>For every scheme carrying a {@code lineage}, we infer the char-0
 * <i>floor</i> from the atoms via
 * {@link FieldAwareLookup#inferFieldFromLineage} (pure metadata — it reads the
 * leaves' own {@code fields[]}, never materializes a tensor) and assert that
 * the scheme declares no char-0 field strictly below that floor
 * ({@code Z < Q < R < C}).</p>
 *
 * <p>F₂/F₃ membership is <b>intentionally excluded</b>: it is a coefficient
 * mod-p representability fact, not derivable from the lineage (a Q leaf with
 * denominators coprime to 3 still reduces into F₃). Those are validated from
 * content by {@code BackfillMissingFields}/{@code SanityCatalogMigration}, not
 * here. We also skip schemes whose lineage field is "unknown" (an unresolved
 * leaf) or finite-only — this check is purely about the char-0 over-claim.</p>
 *
 * <p>This class is a LIBRARY — the CLI lives in {@code VerifySchemes}
 * ({@code --check=lineage}). {@link #check} is the reusable per-scheme predicate.</p>
 */
public final class VerifyLineageFieldCompat {
	private VerifyLineageFieldCompat() {}

	/** A char-0 over-claim: {@code declared} names a field narrower than the lineage {@code floor}. */
	public record Violation(Path file, String declaredNarrowest, String floor, String lineage) {
		@Override public String toString() {
			return file.getFileName() + ": declares " + declaredNarrowest
					+ " but lineage floor is " + floor + "  [" + lineage + "]";
		}
	}

	/**
	 * Returns a {@link Violation} iff {@code root}'s declared {@code fields[]}
	 * names a characteristic-0 field strictly narrower than what its lineage
	 * atoms support. Empty for atoms (no lineage), unknown-leaf lineages, and
	 * finite-field floors. Pure metadata — no dense replay.
	 */
	public static Optional<Violation> check(Path file, JsonNode root, FieldAwareLookup lookup) {
		Optional<Lineage.Node> lin = SchemeIO.readLineage(root);
		if (lin.isEmpty()) return Optional.empty();
		Optional<Field> inferred = lookup.inferFieldFromLineage(lin.get()).field();
		if (inferred.isEmpty() || !isChar0(inferred.get())) return Optional.empty();
		int floor = char0Rank(inferred.get());

		// Narrowest char-0 field the scheme actually declares.
		List<String> declared = SchemeIO.fieldTags(root);
		Field narrowest = null;
		for (String tag : declared) {
			Field d = parseChar0(tag);
			if (d != null && (narrowest == null || char0Rank(d) < char0Rank(narrowest))) {
				narrowest = d;
			}
		}
		if (narrowest == null) return Optional.empty();  // declares no char-0 field at all
		if (char0Rank(narrowest) < floor) {
			String lstr = root.has("lineage_compact") ? root.get("lineage_compact").asString()
					: root.has("lineage_str") ? root.get("lineage_str").asString() : "lineage";
			return Optional.of(new Violation(file, narrowest.name(), inferred.get().name(), lstr));
		}
		return Optional.empty();
	}

	private static boolean isChar0(Field f) {
		return f == Field.Z || f == Field.Q || f == Field.R || f == Field.C;
	}

	private static Field parseChar0(String tag) {
		return switch (tag) {
			case "Z" -> Field.Z;
			case "Q" -> Field.Q;
			case "R" -> Field.R;
			case "C" -> Field.C;
			default -> null;  // F2 / F3 are not part of the char-0 floor check
		};
	}

	private static int char0Rank(Field f) {
		return switch (f) {
			case Z -> 0;
			case Q -> 1;
			case R -> 2;
			case C -> 3;
			default -> -1;
		};
	}
}
