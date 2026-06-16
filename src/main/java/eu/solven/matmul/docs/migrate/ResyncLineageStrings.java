package eu.solven.matmul.docs.migrate;

import java.nio.file.Files;
import java.nio.file.Path;
import java.util.List;
import java.util.concurrent.atomic.AtomicInteger;
import java.util.regex.Pattern;
import java.util.stream.Collectors;

import eu.solven.matmul.catalog.MatrixJsonFormatter;
import eu.solven.matmul.catalog.SchemeIO;

import tools.jackson.databind.JsonNode;
import tools.jackson.databind.node.ObjectNode;

/**
 * Strip content hashes from the human-readable {@code lineage_str} /
 * {@code lineage_compact} strings, leaving the BARE shape (user 2026-06-13: "we
 * should not have the hash in the human lineage_str"). The hash lives only in the
 * structured {@code lineage} DAG, where it pins the exact base for replay; the
 * human form is for reading the COMPOSITION (shapes + operators), matching the
 * dominant convention already in the catalog (e.g. {@code R[2x4x4; 6,6 | …]}).
 *
 * <p>Two hashed forms are reduced to {@code {shape}}:
 * <ul>
 *   <li>{@code {shape}@{hash}}            → {@code {shape}}   (pinned-ref display)</li>
 *   <li>{@code {shape}-r{rank}-{note}-{hash}} → {@code {shape}}   (legacy filename stem)</li>
 * </ul>
 * Text-surgery on the two string values only — NOT a full {@code prettyCompact}
 * regeneration (that has format-drifted and would churn thousands of untouched
 * files). Idempotent: a string with no hash is a no-op.
 *
 * <pre>mvn -q -ntp exec:java -Dexec.mainClass=eu.solven.matmul.docs.migrate.ResyncLineageStrings [-Dexec.args=--execute]</pre>
 */
public final class ResyncLineageStrings {
	private ResyncLineageStrings() {}

	/** {@code {shape}@{hash}} → {@code {shape}}. */
	private static final Pattern AT_HASH = Pattern.compile("(\\d+x\\d+x\\d+)@[0-9a-f]{4,}");
	/** {@code {shape}-r{rank}-{note}-{hash}} → {@code {shape}} (legacy filename stem). */
	private static final Pattern STEM = Pattern.compile("(\\d+x\\d+x\\d+)-r\\d+-[A-Za-z0-9_]+-[0-9a-f]{4,}");

	private static String strip(String s) {
		String out = AT_HASH.matcher(s).replaceAll("$1");
		out = STEM.matcher(out).replaceAll("$1");
		return out;
	}

	public static void main(String[] args) throws Exception {
		boolean execute = List.of(args).contains("--execute");
		Path root = Path.of("src/main/resources/schemes");
		List<Path> files;
		try (var s = Files.walk(root)) {
			files = s.filter(p -> p.toString().endsWith(".json")).collect(Collectors.toList());
		}
		System.out.println("Scanning " + files.size() + " files (mode=" + (execute ? "EXECUTE" : "DRY-RUN") + ")…");

		AtomicInteger changed = new AtomicInteger(), errors = new AtomicInteger();
		int processed = 0;
		for (Path f : files) {
			try {
				JsonNode parsed = SchemeIO.parseJson(f.toFile());
				if (!(parsed instanceof ObjectNode obj)) { continue; }
				boolean fileChanged = false;
				for (String key : new String[] { "lineage_str", "lineage_compact" }) {
					if (!obj.has(key)) { continue; }
					String before = obj.get(key).asString();
					String after = strip(before);
					if (!after.equals(before)) {
						obj.put(key, after);
						fileChanged = true;
					}
				}
				if (fileChanged) {
					changed.incrementAndGet();
					if (changed.get() <= 6) {
						System.out.println("  " + f.getFileName() + " -> " + obj.path("lineage_compact").asString(""));
					}
					if (execute) {
						Files.writeString(f, MatrixJsonFormatter.format(obj));
					}
				}
			} catch (Exception e) {
				errors.incrementAndGet();
				System.out.println("[err] " + f.getFileName() + ": " + e);
			}
			if (++processed % 4000 == 0) System.out.println("[progress] " + processed + "/" + files.size());
		}
		System.out.println("\n=== " + (execute ? "STRIPPED" : "PLAN") + " ===");
		System.out.println("human strings de-hashed: " + changed.get());
		System.out.println("errors:                  " + errors.get());
		if (!execute) System.out.println("\n(DRY-RUN — pass --execute to write)");
	}
}
