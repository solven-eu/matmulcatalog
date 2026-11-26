package eu.solven.matmul.docs.migrate;

import java.nio.file.Files;
import java.nio.file.Path;
import java.util.ArrayList;
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
 * {@code additions}; only whitespace and the canonical key order
 * ({@link MatrixJsonFormatter#canonicalOrder}) change.
 *
 * <pre>
 * mvn -q -ntp exec:java -Dexec.mainClass=eu.solven.matmul.docs.migrate.ReformatSchemes [-Dexec.args="--execute"]
 * mvn -q -ntp exec:java -Dexec.mainClass=eu.solven.matmul.docs.migrate.ReformatSchemes -Dexec.args="--check a.json b.json"
 * </pre>
 *
 * <ul>
 *   <li>default: dry-run over {@code src/main/resources/schemes}, report counts;</li>
 *   <li>{@code --execute}: rewrite the non-canonical files;</li>
 *   <li>{@code --check}: list each non-canonical (or unparsable) file and exit 1 — the CI
 *       gate ({@code test.yml}) runs it on the scheme files a PR changes;</li>
 *   <li>any other argument is a {@code .json} file — or {@code @list.txt}, one path per
 *       line — to process instead of the whole tree (missing files — deleted by the PR —
 *       are skipped).</li>
 * </ul>
 */
public final class ReformatSchemes {
	private ReformatSchemes() {}

	public static void main(String[] args) throws Exception {
		boolean execute = List.of(args).contains("--execute");
		boolean check = List.of(args).contains("--check");
		List<Path> files = new ArrayList<>();
		for (String a : args) {
			if (a.startsWith("--")) continue;
			// "@list.txt": one path per line (a PR's changed files can exceed the arg limit).
			List<String> paths = a.startsWith("@") ? Files.readAllLines(Path.of(a.substring(1))) : List.of(a);
			for (String s : paths) {
				Path p = Path.of(s.strip());
				if (s.strip().endsWith(".json") && Files.isRegularFile(p)) files.add(p);
			}
		}
		boolean explicit = List.of(args).stream().anyMatch(a -> !a.startsWith("--"));
		if (!explicit) {
			try (var s = Files.walk(Path.of("src/main/resources/schemes"))) {
				files = s.filter(p -> p.toString().endsWith(".json")).collect(Collectors.toList());
			}
		}
		String mode = execute ? "EXECUTE" : check ? "CHECK" : "DRY-RUN";
		System.out.println("Scanning " + files.size() + " scheme files (mode=" + mode + ")…");

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
					if (check) System.out.println("NOT CANONICAL: " + f);
					changed.incrementAndGet();
				}
			} catch (Exception e) {
				if (check) System.out.println("UNPARSABLE:    " + f + " (" + e.getMessage() + ")");
				errors.incrementAndGet();
			}
			if (++processed % 2000 == 0) System.out.println("[progress] " + processed + "/" + files.size());
		}
		System.out.println("\n=== " + (execute ? "REFORMATTED" : "PLAN") + " ===");
		System.out.println("already canonical: " + same.get());
		System.out.println("reformatted:       " + changed.get());
		System.out.println("errors (skipped):  " + errors.get());
		if (check && (changed.get() > 0 || errors.get() > 0)) {
			System.out.println("\nRun ReformatSchemes --execute (or write through MatrixJsonFormatter) and commit.");
			System.exit(1);
		}
		if (!execute && !check) System.out.println("\n(DRY-RUN — pass --execute to write)");
	}
}
