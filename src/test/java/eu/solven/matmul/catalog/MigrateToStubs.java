package eu.solven.matmul.catalog;

import java.io.IOException;
import java.nio.file.Files;
import java.nio.file.Path;
import java.util.ArrayList;
import java.util.Comparator;
import java.util.List;
import java.util.regex.Matcher;
import java.util.regex.Pattern;
import java.util.stream.Stream;

import tools.jackson.databind.JsonNode;

/**
 * One-shot migration: strip the {@code u}/{@code v}/{@code w}
 * (and {@code u_sparse}/{@code u_fresh}/…) arrays from large scheme
 * files whose lineage is sufficient to reconstruct them. The remaining
 * "stub" carries only:
 *
 * <ul>
 *   <li>{@code n}: shape</li>
 *   <li>{@code rank} (alias {@code m}), {@code additions} (alias {@code s})
 *       — the load-bearing invariants verifiers must match on replay</li>
 *   <li>{@code field}, attribution, {@code lineage}</li>
 *   <li>{@code scheme_type: "stub"}</li>
 * </ul>
 *
 * <p>Reproducibility: a CI workflow replays the lineage through
 * {@code RecursiveMaterialiser} and verifies the materialised algorithm
 * matches {@code (rank, additions)} — drift is caught immediately. Bit
 * identity is NOT required; only the invariants.</p>
 *
 * <h3>Selection criteria</h3>
 *
 * <p>A file is eligible iff ALL of:</p>
 * <ol>
 *   <li>{@code maxDim > MAX_FLAT_DIM} (default 16, matching FMM-Lille's
 *       per-page coverage).</li>
 *   <li>The JSON has a {@code lineage} field.</li>
 *   <li>The {@code u}/{@code v}/{@code w} (or sparse equivalents) are
 *       still present (i.e. not already a stub).</li>
 * </ol>
 *
 * <h3>Modes</h3>
 *
 * <p>Default: <strong>dry-run</strong> — reports candidates + estimated
 * size savings, modifies nothing. Pass {@code --apply} as the first arg
 * to perform the rewrite.</p>
 *
 * <p>Run:</p>
 * <pre>
 *   # Dry run
 *   mvn -q -o exec:java -Dexec.classpathScope=test \
 *       -Dexec.mainClass=eu.solven.matmul.catalog.MigrateToStubs
 *
 *   # Apply
 *   mvn -q -o exec:java -Dexec.classpathScope=test \
 *       -Dexec.mainClass=eu.solven.matmul.catalog.MigrateToStubs \
 *       -Dexec.args="--apply"
 * </pre>
 */
public final class MigrateToStubs {

	private MigrateToStubs() {}

	private static final Path SCHEMES_ROOT = Path.of("src/main/resources/schemes");
	private static final int MAX_FLAT_DIM = 16;
	private static final Pattern SHAPE = Pattern.compile("_(\\d+)x(\\d+)x(\\d+)_");

	/** Fields preserved on a stub. Everything else is stripped. */
	private static final java.util.Set<String> STUB_FIELDS = java.util.Set.of(
			"n", "rank", "m", "additions", "z2", "field", "complex", "commutative",
			"scheme_type", "source", "attribution_for_rank", "discovery", "year",
			"reference", "notes", "lineage"
	);

	public static void main(String[] args) throws IOException {
		boolean apply = args.length > 0 && "--apply".equals(args[0]);
		System.out.printf("MigrateToStubs %s (threshold maxDim > %d)%n",
				apply ? "APPLY" : "DRY-RUN", MAX_FLAT_DIM);

		List<Path> files;
		try (Stream<Path> walk = Files.walk(SCHEMES_ROOT)) {
			files = walk.filter(p -> p.toString().endsWith(".json"))
					.sorted()
					.toList();
		}

		List<String> skipped = new ArrayList<>();
		long bytesBefore = 0, bytesAfter = 0;
		int stubbed = 0, alreadyStub = 0, noLineage = 0, tooSmall = 0;

		for (Path p : files) {
			long sizeBefore = Files.size(p);
			Matcher m = SHAPE.matcher(p.getFileName().toString());
			if (!m.find()) { tooSmall++; continue; }
			int n = Integer.parseInt(m.group(1));
			int mm = Integer.parseInt(m.group(2));
			int pp = Integer.parseInt(m.group(3));
			int maxDim = Math.max(n, Math.max(mm, pp));
			if (maxDim <= MAX_FLAT_DIM) { tooSmall++; continue; }

			JsonNode root;
			try {
				root = SchemeIO.parseJson(p.toFile());
			} catch (Exception e) {
				skipped.add(p + " (parse failed: " + e.getMessage() + ")");
				continue;
			}
			if (SchemeIO.isStub(root)) { alreadyStub++; continue; }
			if (root.get("lineage") == null || root.get("lineage").isNull()) {
				noLineage++;
				if (noLineage <= 10) {
					skipped.add(p.toString() + " (no lineage)");
				}
				continue;
			}

			String stubJson = renderStub(root);
			bytesBefore += sizeBefore;
			bytesAfter += stubJson.length();
			stubbed++;
			if (apply) {
				Files.writeString(p, stubJson);
			}
		}

		System.out.println();
		System.out.printf("Files walked          : %d%n", files.size());
		System.out.printf("Below threshold       : %d%n", tooSmall);
		System.out.printf("Already stub          : %d%n", alreadyStub);
		System.out.printf("No lineage (skipped)  : %d%n", noLineage);
		System.out.printf("Stubbable             : %d%n", stubbed);
		System.out.printf("Bytes before          : %.1f MB%n", bytesBefore / 1024.0 / 1024);
		System.out.printf("Bytes after (stubs)   : %.1f MB%n", bytesAfter / 1024.0 / 1024);
		System.out.printf("Bytes saved           : %.1f MB%n",
				(bytesBefore - bytesAfter) / 1024.0 / 1024);
		if (!apply) {
			System.out.println();
			System.out.println("DRY-RUN — no files modified. Re-run with --apply to perform the migration.");
		}
		if (!skipped.isEmpty()) {
			System.out.println();
			System.out.println("First 10 skipped (no lineage):");
			skipped.stream().limit(10).sorted(Comparator.naturalOrder())
					.forEach(s -> System.out.println("  " + s));
		}
	}

	/**
	 * Build the stub JSON: keep {@link #STUB_FIELDS}, drop everything else
	 * (most importantly {@code u}/{@code v}/{@code w} and their sparse /
	 * fresh variants). Tag {@code scheme_type: "stub"}.
	 */
	private static String renderStub(JsonNode root) {
		// Tools-jackson ObjectNode is mutable; build via JsonMapper for safety.
		tools.jackson.databind.json.JsonMapper mapper =
				tools.jackson.databind.json.JsonMapper.builder()
						.enable(tools.jackson.databind.SerializationFeature.INDENT_OUTPUT)
						.build();
		tools.jackson.databind.node.ObjectNode stub = mapper.createObjectNode();
		// Preserve fields in stable order.
		for (String key : new String[]{
				"n", "scheme_type", "field", "z2", "complex", "commutative",
				"rank", "m", "additions",
				"source", "attribution_for_rank", "discovery", "year", "reference",
				"notes", "lineage"}) {
			JsonNode v = root.get(key);
			if (v != null) stub.set(key, v);
		}
		// Tag as stub (overwrite if the input already had a scheme_type).
		stub.put("scheme_type", "stub");
		try {
			// Canonical formatting via the single shared formatter (NOT Jackson's
			// default), so stubs match every other scheme writer's style.
			return MatrixJsonFormatter.format(mapper.writeValueAsString(stub));
		} catch (Exception e) {
			throw new RuntimeException("stub serialisation failed", e);
		}
	}
}
