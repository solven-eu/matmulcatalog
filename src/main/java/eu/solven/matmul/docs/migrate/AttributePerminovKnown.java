package eu.solven.matmul.docs.migrate;

import java.nio.file.Files;
import java.nio.file.Path;
import java.util.List;
import java.util.Map;
import java.util.TreeMap;
import java.util.concurrent.atomic.AtomicInteger;
import java.util.stream.Collectors;

import eu.solven.matmul.catalog.MatrixJsonFormatter;
import eu.solven.matmul.catalog.PerminovKnownAttribution;
import eu.solven.matmul.catalog.PerminovKnownAttribution.Attribution;
import eu.solven.matmul.catalog.SchemeIO;

import tools.jackson.databind.JsonNode;
import tools.jackson.databind.node.ObjectNode;

/**
 * One-shot re-attribution of schemes that Perminov files under his
 * {@code schemes/known/<sub-source>/} subtree: those originate elsewhere (he flags
 * them "known"), so crediting {@code source="Perminov 2023"} misattributes the
 * historical record. We rewrite {@code source} to the true originator (via
 * {@link PerminovKnownAttribution}) and record the importer hop in
 * {@code imported_via}, while LEAVING {@code source_scheme_url} /
 * {@code original_source_path} intact (the Perminov-file provenance link).
 *
 * <p>Canonical example: {@code ⟨2,7,7⟩=76} under {@code known/meta_flip_graph/} is
 * Kauers &amp; Wood 2025 (arXiv:2510.19787), not Perminov.</p>
 *
 * <p>Only the {@code source}/{@code imported_via} metadata changes — the U/V/W
 * matrices and content hash are untouched, so filenames and lineage refs stay
 * valid. Schemes under {@code known/tensor} (FMM mirror) and
 * {@code known/matmulcatalog} (our own) are re-attributed too, NOT deleted: they
 * are load-bearing recombination bases for hundreds of {@code derived/} schemes.</p>
 *
 * <pre>mvn -q -ntp exec:java \
 *   -Dexec.mainClass=eu.solven.matmul.docs.migrate.AttributePerminovKnown \
 *   [-Dexec.args=--execute]</pre>
 */
public final class AttributePerminovKnown {
	private AttributePerminovKnown() {}

	private static final String IMPORTED_VIA = "Perminov FastMatrixMultiplication";

	public static void main(String[] args) throws Exception {
		boolean execute = List.of(args).contains("--execute");
		Path root = Path.of("src/main/resources/schemes");
		List<Path> files;
		try (var s = Files.walk(root)) {
			files = s.filter(p -> p.toString().endsWith(".json")).collect(Collectors.toList());
		}
		System.out.println("Scanning " + files.size() + " scheme files (mode="
				+ (execute ? "EXECUTE" : "DRY-RUN") + ")…");

		AtomicInteger reattributed = new AtomicInteger(), keptPerminov = new AtomicInteger(),
				notPerminov = new AtomicInteger(), noPath = new AtomicInteger();
		// source label → count, for a readable summary.
		Map<String, Integer> bySource = new TreeMap<>();
		int processed = 0;
		for (Path f : files) {
			try {
				JsonNode parsed = SchemeIO.parseJson(f.toFile());
				if (!(parsed instanceof ObjectNode obj)) continue;
				String source = obj.has("source") ? obj.get("source").asString() : "";
				// Only touch schemes currently credited to Perminov.
				if (!source.toLowerCase().contains("perminov")) { notPerminov.incrementAndGet(); continue; }
				if (!obj.has("original_source_path")) { noPath.incrementAndGet(); continue; }
				String osp = obj.get("original_source_path").asString();

				Attribution attr = PerminovKnownAttribution.forPath(osp).orElse(null);
				if (attr == null || attr.isPerminovOwn()) {
					// schemes/results/* — genuinely Perminov's own; leave as-is.
					keptPerminov.incrementAndGet();
					continue;
				}
				bySource.merge(attr.source(), 1, Integer::sum);
				if (execute) {
					obj.put("source", attr.source());
					obj.put("imported_via", IMPORTED_VIA);
					Files.writeString(f, MatrixJsonFormatter.format(obj));
				}
				reattributed.incrementAndGet();
			} catch (Exception e) {
				System.out.println("[skip] " + f + ": " + e);
			}
			if (++processed % 2000 == 0) {
				System.out.println("[progress] " + processed + "/" + files.size());
			}
		}

		System.out.println("\n=== " + (execute ? "RE-ATTRIBUTED" : "PLAN") + " ===");
		System.out.println("re-attributed (known/<sub> → true origin): " + reattributed.get());
		bySource.forEach((src, n) -> System.out.printf("    %5d  → %s%n", n, src));
		System.out.println("kept Perminov (schemes/results/* — his own): " + keptPerminov.get());
		System.out.println("not Perminov-sourced (skipped): " + notPerminov.get());
		System.out.println("Perminov but no original_source_path (skipped): " + noPath.get());
		if (!execute) System.out.println("\n(DRY-RUN — pass --execute to write)");
		else System.out.println("\nNext: GenerateCatalogManifest to refresh docs/catalog.json.");
	}
}
