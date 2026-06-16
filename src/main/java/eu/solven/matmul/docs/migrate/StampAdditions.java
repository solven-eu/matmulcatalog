package eu.solven.matmul.docs.migrate;

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
 * Stamp the {@code additions} count into every scheme JSON that lacks it, reading it
 * from the current filename's {@code _a{N}} token (the historical authoritative home).
 * The last filename-derived property the manifest reads — after this + StampFields +
 * StampSource, the manifest is fully content-driven and the catalog can be renamed.
 * Run BEFORE any rename.
 *
 * <pre>mvn -q -ntp exec:java -Dexec.mainClass=eu.solven.matmul.docs.migrate.StampAdditions [-Dexec.args=--execute]</pre>
 */
public final class StampAdditions {
	private StampAdditions() {}

	private static final Pattern A = Pattern.compile("_a(\\d+)");

	public static void main(String[] args) throws Exception {
		boolean execute = List.of(args).contains("--execute");
		Path root = Path.of("src/main/resources/schemes");
		List<Path> files;
		try (var s = Files.walk(root)) {
			files = s.filter(p -> p.toString().endsWith(".json")).collect(Collectors.toList());
		}
		System.out.println("Scanning " + files.size() + " (mode=" + (execute ? "EXECUTE" : "DRY-RUN") + ")…");
		AtomicInteger already = new AtomicInteger(), stamped = new AtomicInteger(), noToken = new AtomicInteger();
		int processed = 0;
		for (Path f : files) {
			try {
				JsonNode parsed = SchemeIO.parseJson(f.toFile());
				if (parsed.has("additions") && parsed.get("additions").isInt()) { already.incrementAndGet(); continue; }
				if (!(parsed instanceof ObjectNode obj)) continue;
				Matcher m = A.matcher(f.getFileName().toString());
				if (!m.find()) { noToken.incrementAndGet(); continue; }
				if (execute) {
					obj.put("additions", Integer.parseInt(m.group(1)));
					Files.writeString(f, MatrixJsonFormatter.format(obj));
				}
				stamped.incrementAndGet();
			} catch (Exception e) {
				// skip
			}
			if (++processed % 2000 == 0) System.out.println("[progress] " + processed + "/" + files.size());
		}
		System.out.println("\n=== " + (execute ? "STAMPED" : "PLAN") + " ===");
		System.out.println("already had additions: " + already.get());
		System.out.println("stamped from _a token:  " + stamped.get());
		System.out.println("no _a token (left null): " + noToken.get());
		if (!execute) System.out.println("\n(DRY-RUN — pass --execute to write)");
	}
}
