package eu.solven.matmul.docs.migrate;

import java.nio.file.Files;
import java.nio.file.Path;
import java.util.ArrayList;
import java.util.List;
import java.util.Optional;
import java.util.concurrent.atomic.AtomicInteger;
import java.util.stream.Collectors;

import eu.solven.matmul.NonCubicBilinearAlgorithm;
import eu.solven.matmul.algebra.Field;
import eu.solven.matmul.catalog.FieldAwareLookup;
import eu.solven.matmul.catalog.Lineage;
import eu.solven.matmul.catalog.MatrixJsonFormatter;
import eu.solven.matmul.catalog.SchemeIO;

import tools.jackson.databind.JsonNode;
import tools.jackson.databind.node.ArrayNode;
import tools.jackson.databind.node.ObjectNode;

/**
 * Normalize the {@code fields[]} / {@code fields_not[]} of schemes against their
 * LINEAGE — the systemic fix for the derived ⟨3,3,8⟩-style char-0 over-claim
 * (~1.7k matrix-less stubs declared <b>Z</b> while a Q-only leaf forces a Q
 * floor, so their rational coefficients propagate and Z is false, not merely
 * unprovable).
 *
 * <h3>Policy</h3>
 * <ul>
 *   <li><b>matrix-bearing &amp; over-claiming</b> (declared char-0 narrowest below
 *       the lineage floor) → recompute {@code fields[]} from CONTENT (ground
 *       truth) via {@link BackfillMissingFields#fieldsFromContent}.</li>
 *   <li><b>matrix-less stub</b> → claim exactly the lineage char-0 floor,
 *       expanded ({@code Q→[Q,R,C]}); everything else is {@code fields_not[]}
 *       ({@code [F2,F3,Z]}). The below-floor char-0 field (Z) is a hard
 *       exclusion (a non-integer leaf coefficient propagates); F₂/F₃ are
 *       excluded too because a stub cannot materialize to certify them, so it
 *       claims only what its atoms prove. Idempotent: only rewrites a stub whose
 *       fields/fields_not already differ from this honest form.</li>
 * </ul>
 *
 * <pre>mvn -q -ntp exec:java -Dexec.mainClass=eu.solven.matmul.docs.migrate.RestampLineageOverclaims [-Dexec.args=--execute]</pre>
 */
public final class RestampLineageOverclaims {
	private RestampLineageOverclaims() {}

	private static final List<String> ORDER = List.of("F2", "F3", "Z", "Q", "R", "C");

	public static void main(String[] args) throws Exception {
		boolean execute = List.of(args).contains("--execute");
		Path root = Path.of("src/main/resources/schemes");
		FieldAwareLookup lookup = new FieldAwareLookup("Q");

		List<Path> files;
		try (var s = Files.walk(root)) {
			files = s.filter(p -> p.toString().endsWith(".json")).collect(Collectors.toList());
		}
		System.out.println("Scanning " + files.size() + " scheme files (mode=" + (execute ? "EXECUTE" : "DRY-RUN") + ")…");

		AtomicInteger fromContent = new AtomicInteger(), fromFloor = new AtomicInteger();
		AtomicInteger errors = new AtomicInteger();
		int processed = 0;
		for (Path f : files) {
			try {
				JsonNode parsed = SchemeIO.parseJson(f.toFile());
				Optional<Lineage.Node> lin = SchemeIO.readLineage(parsed);
				if (lin.isEmpty()) continue;
				Optional<Field> floorOpt = lookup.inferFieldFromLineage(lin.get()).field();
				if (floorOpt.isEmpty() || !isChar0(floorOpt.get())) continue;
				if (!(parsed instanceof ObjectNode obj)) { errors.incrementAndGet(); continue; }
				Field floor = floorOpt.get();

				List<String> declared = SchemeIO.fieldTags(parsed);
				Field narrowest = narrowestChar0(declared);
				boolean overClaim = narrowest != null && char0Rank(narrowest) < char0Rank(floor);

				boolean hasMatrices = parsed.has("u") || parsed.has("u_sparse");
				boolean contentEligible = hasMatrices && !SchemeIO.isComplex(parsed)
						&& !SchemeIO.isNonBilinear(parsed) && !SchemeIO.isReduced(parsed);

				List<String> newFields, newFieldsNot;
				if (contentEligible) {
					if (!overClaim) continue;  // content already consistent or wider — leave it
					NonCubicBilinearAlgorithm alg = SchemeIO.read(parsed);
					newFields = order(BackfillMissingFields.fieldsFromContent(alg));
					final List<String> nf = newFields;
					newFieldsNot = ORDER.stream().filter(t -> !nf.contains(t)).collect(Collectors.toList());
					if (newFields.equals(order(declared)) && newFieldsNot.equals(readFieldsNot(parsed))) continue;
					fromContent.incrementAndGet();
				} else {
					// Stub: claim the lineage char-0 floor (Q→[Q,R,C]); everything else
					// is fields_not ([F2,F3,Z]). A stub cannot materialize, so it
					// certifies neither F₂ nor F₃ — excluded alongside the below-floor
					// char-0 field (Z, a hard exclusion since the rational propagates).
					//
					// SAFETY: only RESTRICT (floor Q/R/C → never claims F2/F3/Z), never
					// WIDEN. A floor-Z inference is NOT trustworthy for a stub: it is
					// read from LEAF fields, but a recombination can introduce rational
					// multipliers of its own (½-symmetrization) that leaf inference can't
					// see — so claiming Z/F2/F3 could over-claim. Leave floor-Z stubs
					// untouched (their integrality needs materialization to certify).
					if (floor == Field.Z) continue;
					newFields = order(expand(floor));
					final List<String> nf = newFields;
					newFieldsNot = ORDER.stream().filter(t -> !nf.contains(t)).collect(Collectors.toList());
					if (newFields.equals(order(declared)) && newFieldsNot.equals(readFieldsNot(parsed))) continue;
					fromFloor.incrementAndGet();
				}

				System.out.println((execute ? "[write] " : "[plan]  ") + f.getFileName()
						+ "  floor " + floor.name() + (overClaim ? " (over-claim " + narrowest + ")" : "")
						+ " -> fields=" + newFields + " not=" + newFieldsNot);

				if (execute) {
					ArrayNode fa = obj.arrayNode();
					newFields.forEach(fa::add);
					obj.set("fields", fa);
					if (newFieldsNot.isEmpty()) {
						obj.remove("fields_not");
					} else {
						ArrayNode fn = obj.arrayNode();
						newFieldsNot.forEach(fn::add);
						obj.set("fields_not", fn);
					}
					Files.writeString(f, MatrixJsonFormatter.format(obj));
				}
			} catch (Exception e) {
				errors.incrementAndGet();
				System.out.println("[error] " + f.getFileName() + ": " + e);
			}
			if (++processed % 4000 == 0) System.out.println("[progress] " + processed + "/" + files.size());
		}

		System.out.println("\n=== " + (execute ? "RE-STAMPED" : "PLAN") + " ===");
		System.out.println("re-stamped from CONTENT (matrices):  " + fromContent.get());
		System.out.println("re-stamped to LINEAGE FLOOR (stub):  " + fromFloor.get());
		System.out.println("errors (skipped):                    " + errors.get());
		if (!execute) System.out.println("\n(DRY-RUN — pass --execute to write)");
	}

	private static List<String> readFieldsNot(JsonNode root) {
		List<String> out = new ArrayList<>();
		JsonNode fn = root.get("fields_not");
		if (fn != null && fn.isArray()) fn.forEach(t -> out.add(t.asString()));
		return order(out);
	}

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

	private static boolean isChar0(Field f) {
		return f == Field.Z || f == Field.Q || f == Field.R || f == Field.C;
	}

	private static Field narrowestChar0(List<String> declared) {
		Field n = null;
		for (String tag : declared) {
			Field d = switch (tag) {
				case "Z" -> Field.Z; case "Q" -> Field.Q; case "R" -> Field.R; case "C" -> Field.C;
				default -> null;
			};
			if (d != null && (n == null || char0Rank(d) < char0Rank(n))) n = d;
		}
		return n;
	}

	private static int char0Rank(Field f) {
		return switch (f) { case Z -> 0; case Q -> 1; case R -> 2; case C -> 3; default -> -1; };
	}

	private static List<String> order(List<String> fields) {
		return fields.stream().sorted((a, b) -> Integer.compare(
				idx(ORDER.indexOf(a)), idx(ORDER.indexOf(b)))).collect(Collectors.toList());
	}

	private static int idx(int i) { return i < 0 ? 99 : i; }
}
