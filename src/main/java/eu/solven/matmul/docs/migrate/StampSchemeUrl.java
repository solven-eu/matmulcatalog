package eu.solven.matmul.docs.migrate;

import java.nio.file.Files;
import java.nio.file.Path;
import java.util.LinkedHashMap;
import java.util.List;
import java.util.Map;
import java.util.concurrent.atomic.AtomicInteger;
import java.util.stream.Collectors;

import eu.solven.matmul.catalog.SchemeIO;

import tools.jackson.databind.JsonNode;

/**
 * Backfill {@code source_scheme_url} — a pointer to the scheme's file in the
 * originating author's own repository — onto every imported scheme that records
 * an {@code original_source_path} but predates the field.
 *
 * <p>Today this only fires for Perminov imports: every Perminov scheme carries
 * {@code original_source_path} (the relative path inside
 * {@code dronperminov/FastMatrixMultiplication}); the URL is just
 * {@link ImportPerminovSchemes#BLOB_BASE} + that path. The field is distinct from
 * {@code source_paper_url} (the author's paper) — it points at the JSON file, not
 * the publication — and survives round-trips via {@link SchemeIO#addFields}
 * (which never clobbers an existing value). {@link ImportPerminovSchemes} sets it
 * inline going forward, so re-running this after the next import is a no-op.</p>
 *
 * <pre>mvn -q -ntp exec:java -Dexec.mainClass=eu.solven.matmul.docs.migrate.StampSchemeUrl [-Dexec.args=--execute]</pre>
 */
public final class StampSchemeUrl {
	private StampSchemeUrl() {}

	public static void main(String[] args) throws Exception {
		boolean execute = List.of(args).contains("--execute");
		Path root = Path.of("src/main/resources/schemes");
		List<Path> files;
		try (var s = Files.walk(root)) {
			files = s.filter(p -> p.toString().endsWith(".json")).collect(Collectors.toList());
		}
		System.out.println("Scanning " + files.size() + " scheme files (mode="
				+ (execute ? "EXECUTE" : "DRY-RUN") + ")…");

		AtomicInteger already = new AtomicInteger(), stamped = new AtomicInteger(), noPath = new AtomicInteger();
		int processed = 0;
		for (Path f : files) {
			try {
				JsonNode root2 = SchemeIO.parseJson(f.toFile());
				if (root2.has("source_scheme_url") && !root2.get("source_scheme_url").asString().isBlank()) {
					already.incrementAndGet();
					continue;
				}
				JsonNode osp = root2.get("original_source_path");
				if (osp == null || osp.asString().isBlank()) {
					noPath.incrementAndGet();
					continue;
				}
				String url = ImportPerminovSchemes.BLOB_BASE + osp.asString();
				if (execute) {
					Map<String, Object> add = new LinkedHashMap<>();
					add.put("source_scheme_url", url);
					// addFields preserves every other field and won't clobber.
					SchemeIO.addFields(f.toFile(), add, /* apply */ true);
				}
				stamped.incrementAndGet();
			} catch (Exception e) {
				System.out.println("[ERR] " + f + ": " + e);
			}
			if (++processed % 2000 == 0) {
				System.out.println("[progress] " + processed + "/" + files.size()
						+ " (stamped=" + stamped.get() + ")");
			}
		}
		System.out.println("\n=== " + (execute ? "STAMPED" : "PLAN") + " ===");
		System.out.println("already had source_scheme_url: " + already.get());
		System.out.println("stamped from original_source_path: " + stamped.get());
		System.out.println("no original_source_path (skipped): " + noPath.get());
		if (!execute) {
			System.out.println("\n(DRY-RUN — pass --execute to write)");
		}
	}
}
