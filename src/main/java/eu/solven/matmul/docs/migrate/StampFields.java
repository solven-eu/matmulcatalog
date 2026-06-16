package eu.solven.matmul.docs.migrate;

import java.nio.file.Files;
import java.nio.file.Path;
import java.util.List;
import java.util.Optional;
import java.util.concurrent.atomic.AtomicInteger;
import java.util.stream.Collectors;

import eu.solven.matmul.algebra.Field;
import eu.solven.matmul.catalog.FieldAwareLookup;
import eu.solven.matmul.catalog.Lineage;
import eu.solven.matmul.catalog.MatrixJsonFormatter;
import eu.solven.matmul.catalog.SchemeIO;

import tools.jackson.databind.JsonNode;
import tools.jackson.databind.node.ArrayNode;
import tools.jackson.databind.node.ObjectNode;

/**
 * Stamp a content-derived {@code fields[]} property into every scheme JSON that
 * lacks one — the enabler that let us drop filename field-parsing entirely
 * (the former {@code FieldAwareLookup.classifyFilenameField}) and rely on JSON
 * content as the single source of truth.
 *
 * <p>Imports already carry verified {@code fields[]}. Derived stubs (matrix-less)
 * don't, so we infer the narrowest field from the LINEAGE (no replay) via
 * {@link FieldAwareLookup#inferFieldFromLineage} and expand it to the inclusion-correct
 * set ({@code Z→[Z,Q,R,C]}, {@code Q→[Q,R,C]}, …; F₂/F₃ are separate prime fields).
 * Inference unknown → conservative {@code [Z]}.</p>
 *
 * <p>Default DRY-RUN (reports counts). {@code --execute} writes. Re-run the manifest
 * afterwards to confirm the catalog is unchanged.</p>
 *
 * <pre>mvn -q -ntp exec:java -Dexec.mainClass=eu.solven.matmul.docs.migrate.StampFields [-Dexec.args=--execute]</pre>
 */
public final class StampFields {
	private StampFields() {}

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

		AtomicInteger already = new AtomicInteger(), inferred = new AtomicInteger();
		AtomicInteger unknownInfer = new AtomicInteger(), noLineage = new AtomicInteger(), errors = new AtomicInteger();
		int processed = 0;
		for (Path f : files) {
			try {
				JsonNode parsed = SchemeIO.parseJson(f.toFile());
				if (!SchemeIO.fieldTags(parsed).isEmpty()) { already.incrementAndGet(); continue; }
				if (!(parsed instanceof ObjectNode obj)) { errors.incrementAndGet(); continue; }

				// Field set = INTERSECTION of the lineage leaves' stamped fields[]
				// (content-only, FieldAwareLookup#fieldNamesFromLineage). NO [Z]
				// fallback — that over-claimed Z on rational-leaf schemes (e.g.
				// ⟨5,32,32⟩ via the rational ⟨3,8,8⟩=145). Inference unknown /
				// no lineage → leave UNSTAMPED rather than guess.
				List<String> fields = List.of();
				Optional<Lineage.Node> lin = SchemeIO.readLineage(parsed);
				if (lin.isPresent()) {
					fields = lookup.fieldNamesFromLineage(lin.get());
					if (!fields.isEmpty()) { inferred.incrementAndGet(); }
					else { unknownInfer.incrementAndGet(); }
				} else {
					noLineage.incrementAndGet();
				}

				if (execute && !fields.isEmpty()) {
					ArrayNode fa = obj.arrayNode();
					fields.forEach(fa::add);  // fieldNamesFromLineage already returns canonical order
					obj.set("fields", fa);
					Files.writeString(f, MatrixJsonFormatter.format(obj));
				}
			} catch (Exception e) {
				errors.incrementAndGet();
			}
			if (++processed % 2000 == 0) System.out.println("[progress] " + processed + "/" + files.size());
		}

		System.out.println("\n=== " + (execute ? "STAMPED" : "PLAN") + " ===");
		System.out.println("already had fields[] (imports + newer): " + already.get());
		System.out.println("stamped via leaf-fields intersection:    " + inferred.get());
		System.out.println("inference unknown → left UNSTAMPED:      " + unknownInfer.get());
		System.out.println("no lineage → left UNSTAMPED:             " + noLineage.get());
		System.out.println("errors (skipped):                       " + errors.get());
		if (!execute) System.out.println("\n(DRY-RUN — pass --execute to write)");
	}

	/**
	 * Narrowest field → inclusion-correct verified set. F₂/F₃ are separate prime
	 * fields and never lift to characteristic 0; but a Z (integer-exact) scheme
	 * reduces mod ANY prime, so Z grants F₂ and F₃ as a theorem — the same policy
	 * as {@link NarrowFields#expandFromBase}. (These two stampers historically
	 * disagreed: this one emitted [Z,Q,R,C], leaving ~8k integer schemes invisible
	 * under the SPA's F2/F3 selectors. The expansion is exactly as trustworthy as
	 * the input Z claim itself — schemes that over-claim Z while holding rational
	 * coefficients are a separate, tracked data bug.)
	 */
	private static List<String> expand(Field f) {
		// Canonical expander (single source of truth) — see FieldAwareLookup.
		return FieldAwareLookup.inclusionFieldNames(f);
	}

	private static List<String> order(List<String> fields) {
		return fields.stream().sorted((a, b) -> Integer.compare(
				idx(ORDER.indexOf(a)), idx(ORDER.indexOf(b)))).collect(Collectors.toList());
	}

	private static int idx(int i) { return i < 0 ? 99 : i; }
}
