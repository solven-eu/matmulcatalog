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
 * Remove the bogus filename-echo {@code source} values stamped by the pre-fix
 * {@link StampSource} (2026-07-07): on new-convention filenames
 * ({@code {n}x{m}x{p}-r{rank}-{note}-{hash7}}) its old-convention prefix
 * extraction fell back to the whole stem, title-cased it, and wrote e.g.
 * {@code "source": "10x3x3-R69-Derived-79abb2c"} — a label echo, not an
 * attribution. No legitimate source starts with a shape token, so any
 * {@code source} matching {@code ^\d+x\d+x\d+-} is stripped.
 *
 * <pre>mvn -q -ntp exec:java -Dexec.mainClass=eu.solven.matmul.docs.migrate.FixFilenameEchoSources [-Dexec.args=--execute]</pre>
 */
public final class FixFilenameEchoSources {
	private FixFilenameEchoSources() {}

	static final Pattern ECHO = Pattern.compile("^\\d+x\\d+x\\d+-");

	static boolean isFilenameEcho(String source) {
		return source != null && ECHO.matcher(source).find();
	}

	public static void main(String[] args) throws Exception {
		boolean execute = List.of(args).contains("--execute");
		Path root = Path.of("src/main/resources/schemes");
		List<Path> files;
		try (var s = Files.walk(root)) {
			files = s.filter(p -> p.toString().endsWith(".json")).collect(Collectors.toList());
		}
		System.out.println("Scanning " + files.size() + " scheme files (mode=" + (execute ? "EXECUTE" : "DRY-RUN") + ")…");

		AtomicInteger stripped = new AtomicInteger(), kept = new AtomicInteger();
		int processed = 0;
		for (Path f : files) {
			try {
				JsonNode parsed = SchemeIO.parseJson(f.toFile());
				if (!(parsed instanceof ObjectNode obj) || !obj.has("source")) continue;
				String source = obj.get("source").asString();
				if (!isFilenameEcho(source)) { kept.incrementAndGet(); continue; }
				if (execute) {
					obj.remove("source");
					Files.writeString(f, MatrixJsonFormatter.format(obj));
				}
				stripped.incrementAndGet();
			} catch (Exception e) {
				System.err.println("SKIP (parse error) " + f + " : " + e);
			}
			if (++processed % 2000 == 0) System.out.println("[progress] " + processed + "/" + files.size());
		}
		System.out.println("\n=== " + (execute ? "STRIPPED" : "PLAN") + " ===");
		System.out.println("filename-echo sources " + (execute ? "removed" : "to remove") + ": " + stripped.get());
		System.out.println("legitimate sources kept: " + kept.get());
		if (!execute) System.out.println("\n(DRY-RUN — pass --execute to write)");
	}
}
