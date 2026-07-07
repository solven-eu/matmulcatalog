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
 * <p>New-convention filenames ({@code {n}x{m}x{p}-r{rank}-{note}-{hash7}}) carry no
 * source prefix — the shape comes first — so they are SKIPPED, never stamped. (The
 * 2026-07-07 pre-fix version fell back to the whole stem and stamped 7.6k filename
 * echoes like {@code "10x3x3-R69-Derived-79abb2c"}; {@link FixFilenameEchoSources}
 * reverted them.)</p>
 *
 * <pre>mvn -q -ntp exec:java -Dexec.mainClass=eu.solven.matmul.docs.migrate.StampSource [-Dexec.args=--execute]</pre>
 */
public final class StampSource {
	private StampSource() {}

	private static final Pattern SHAPE = Pattern.compile("[-_](\\d+)x(\\d+)x(\\d+)");
	private static final Pattern NEW_CONVENTION = Pattern.compile("^\\d+x\\d+x\\d+-");

	/**
	 * Source attribution for an old-convention filename stem, or empty for a
	 * new-convention stem (shape-first: filenames are pure labels, nothing to stamp).
	 */
	static java.util.Optional<String> sourceForStem(String stem) {
		if (NEW_CONVENTION.matcher(stem).find()) {
			return java.util.Optional.empty();
		}
		Matcher m = SHAPE.matcher(stem);
		String prefix = m.find() ? stem.substring(0, m.start()) : stem;
		if (prefix.isBlank()) {
			prefix = "unknown";
		}
		return java.util.Optional.of(GenerateCatalogManifest.normalizeSource(prefix));
	}

	public static void main(String[] args) throws Exception {
		boolean execute = List.of(args).contains("--execute");
		Path root = Path.of("src/main/resources/schemes");
		List<Path> files;
		try (var s = Files.walk(root)) {
			files = s.filter(p -> p.toString().endsWith(".json")).collect(Collectors.toList());
		}
		System.out.println("Scanning " + files.size() + " scheme files (mode=" + (execute ? "EXECUTE" : "DRY-RUN") + ")…");

		AtomicInteger already = new AtomicInteger(), stamped = new AtomicInteger(), newConvention = new AtomicInteger();
		int processed = 0;
		for (Path f : files) {
			try {
				JsonNode parsed = SchemeIO.parseJson(f.toFile());
				if (parsed.has("source") && !parsed.get("source").asString().isBlank()) { already.incrementAndGet(); continue; }
				if (!(parsed instanceof ObjectNode obj)) continue;
				String stem = f.getFileName().toString().replaceFirst("\\.json$", "");
				java.util.Optional<String> source = sourceForStem(stem);
				if (source.isEmpty()) { newConvention.incrementAndGet(); continue; }
				if (execute) {
					obj.put("source", source.get());
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
		System.out.println("new-convention (shape-first) filenames skipped: " + newConvention.get());
		System.out.println("stamped from filename prefix: " + stamped.get());
		if (!execute) System.out.println("\n(DRY-RUN — pass --execute to write)");
	}
}
