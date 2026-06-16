package eu.solven.matmul.docs.migrate;

import java.nio.file.Files;
import java.nio.file.Path;
import java.util.List;
import java.util.concurrent.atomic.AtomicInteger;
import java.util.stream.Collectors;

import eu.solven.matmul.catalog.MatrixJsonFormatter;
import eu.solven.matmul.catalog.SchemeIO;

import tools.jackson.databind.JsonNode;

/**
 * Re-serialise every scheme JSON in the canonical {@link MatrixJsonFormatter} format
 * (the same one {@code SchemeIO.write}/{@code updateFields} emit). The stamping +
 * rename passes wrote files with Jackson's {@code toPrettyString} (spaces around
 * colons, inline matrix rows), which drifts from the canonical style. This normalises
 * all files to one format so diffs are clean and future writes don't churn. Content is
 * preserved exactly (parse → reformat), including stamped {@code fields}/{@code source}/
 * {@code additions}.
 *
 * <pre>mvn -q -ntp exec:java -Dexec.mainClass=eu.solven.matmul.docs.migrate.ReformatSchemes [-Dexec.args=--execute]</pre>
 */
public final class ReformatSchemes {
	private ReformatSchemes() {}

	public static void main(String[] args) throws Exception {
		boolean execute = List.of(args).contains("--execute");
		Path root = Path.of("src/main/resources/schemes");
		List<Path> files;
		try (var s = Files.walk(root)) {
			files = s.filter(p -> p.toString().endsWith(".json")).collect(Collectors.toList());
		}
		System.out.println("Scanning " + files.size() + " scheme files (mode=" + (execute ? "EXECUTE" : "DRY-RUN") + ")…");

		AtomicInteger changed = new AtomicInteger(), same = new AtomicInteger(), errors = new AtomicInteger();
		int processed = 0;
		for (Path f : files) {
			try {
				JsonNode node = SchemeIO.parseJson(f.toFile());
				String formatted = MatrixJsonFormatter.format(node);
				String current = Files.readString(f);
				if (formatted.equals(current)) {
					same.incrementAndGet();
				} else {
					if (execute) Files.writeString(f, formatted);
					changed.incrementAndGet();
				}
			} catch (Exception e) {
				errors.incrementAndGet();
			}
			if (++processed % 2000 == 0) System.out.println("[progress] " + processed + "/" + files.size());
		}
		System.out.println("\n=== " + (execute ? "REFORMATTED" : "PLAN") + " ===");
		System.out.println("already canonical: " + same.get());
		System.out.println("reformatted:       " + changed.get());
		System.out.println("errors (skipped):  " + errors.get());
		if (!execute) System.out.println("\n(DRY-RUN — pass --execute to write)");
	}
}
