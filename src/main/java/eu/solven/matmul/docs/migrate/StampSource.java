package eu.solven.matmul.docs.migrate;

import eu.solven.matmul.docs.generate.GenerateCatalogManifest;

import java.nio.file.Files;
import java.nio.file.Path;
import java.util.List;
import java.util.concurrent.atomic.AtomicInteger;
import java.util.regex.Matcher;
import java.util.regex.Pattern;
import java.util.stream.Collectors;

import eu.solven.matmul.catalog.MatrixJsonFormatter;
import eu.solven.matmul.catalog.SchemeIO;

import tools.jackson.databind.JsonNode;
import tools.jackson.databind.node.ObjectNode;

/**
 * Stamp a content {@code source} (display attribution) into every scheme JSON that
 * lacks one — the last filename-derived property the manifest reads, so that after
 * this the manifest can be fully content-driven and the catalog can be renamed.
 *
 * <p>The source is derived from the CURRENT filename prefix (everything before the
 * {@code NxMxP} token) via {@link GenerateCatalogManifest#normalizeSource}
 * (e.g. {@code alphatensor_Z} → "AlphaTensor 2022", {@code derived_recursive} →
 * "Derived_recursive"). Run BEFORE any rename, while filenames still carry the prefix.</p>
 *
 * <pre>mvn -q -ntp exec:java -Dexec.mainClass=eu.solven.matmul.docs.migrate.StampSource [-Dexec.args=--execute]</pre>
 */
public final class StampSource {
	private StampSource() {}

	private static final Pattern SHAPE = Pattern.compile("[-_](\\d+)x(\\d+)x(\\d+)");

	public static void main(String[] args) throws Exception {
		boolean execute = List.of(args).contains("--execute");
		Path root = Path.of("src/main/resources/schemes");
		List<Path> files;
		try (var s = Files.walk(root)) {
			files = s.filter(p -> p.toString().endsWith(".json")).collect(Collectors.toList());
		}
		System.out.println("Scanning " + files.size() + " scheme files (mode=" + (execute ? "EXECUTE" : "DRY-RUN") + ")…");

		AtomicInteger already = new AtomicInteger(), stamped = new AtomicInteger(), noPrefix = new AtomicInteger();
		int processed = 0;
		for (Path f : files) {
			try {
				JsonNode parsed = SchemeIO.parseJson(f.toFile());
				if (parsed.has("source") && !parsed.get("source").asString().isBlank()) { already.incrementAndGet(); continue; }
				if (!(parsed instanceof ObjectNode obj)) continue;
				String stem = f.getFileName().toString().replaceFirst("\\.json$", "");
				Matcher m = SHAPE.matcher(stem);
				String prefix = m.find() ? stem.substring(0, m.start()) : stem;
				if (prefix.isBlank()) { noPrefix.incrementAndGet(); prefix = "unknown"; }
				String source = GenerateCatalogManifest.normalizeSource(prefix);
				if (execute) {
					obj.put("source", source);
					Files.writeString(f, MatrixJsonFormatter.format(obj));
				}
				stamped.incrementAndGet();
			} catch (Exception e) {
				// skip
			}
			if (++processed % 2000 == 0) System.out.println("[progress] " + processed + "/" + files.size());
		}
		System.out.println("\n=== " + (execute ? "STAMPED" : "PLAN") + " ===");
		System.out.println("already had source: " + already.get());
		System.out.println("stamped from filename prefix: " + stamped.get() + "  (" + noPrefix.get() + " had no prefix → 'unknown')");
		if (!execute) System.out.println("\n(DRY-RUN — pass --execute to write)");
	}
}
