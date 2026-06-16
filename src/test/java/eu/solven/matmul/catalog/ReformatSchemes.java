package eu.solven.matmul.catalog;

import java.io.File;
import java.io.IOException;
import java.nio.file.Files;
import java.util.Arrays;

import eu.solven.matmul.catalog.MatrixJsonFormatter;

/**
 * Canonicalize every {@code src/**\/*.json} file using
 * {@link MatrixJsonFormatter} — primitives inline, compound vertical.
 *
 * <p>Run via:
 * <pre>
 *   mvn -q test-compile
 *   java -cp target/classes:target/test-classes:&dollar;CLASSPATH \
 *        eu.solven.matmul.catalog.ReformatSchemes
 * </pre>
 *
 * <p>Idempotent (a second run after a clean one is a no-op). Intended for use
 * before commit; you can also wire it into your editor's save hook.</p>
 */
public class ReformatSchemes {

	public static void main(String[] args) throws IOException {
		File target = args.length > 0 ? new File(args[0]) : new File("src/main/resources/schemes");

		// Single-file fast path — useful while iterating on the formatter rules.
		if (target.isFile()) {
			String before = Files.readString(target.toPath());
			String after = MatrixJsonFormatter.format(before);
			if (!before.equals(after)) {
				Files.writeString(target.toPath(), after);
				System.out.printf("ReformatSchemes: rewrote %s%n", target);
			} else {
				System.out.printf("ReformatSchemes: %s already canonical%n", target);
			}
			return;
		}

		if (!target.isDirectory()) {
			System.err.println("not a file or directory: " + target);
			System.exit(1);
		}

		int changed = 0, unchanged = 0, errors = 0;
		// Recurse into section{N}/ subdirectories.
		File[] files;
		try (var s = java.nio.file.Files.walk(target.toPath())) {
			files = s.filter(p -> p.toString().endsWith(".json"))
					.map(java.nio.file.Path::toFile)
					.toArray(File[]::new);
		}
		if (files == null) {
			System.err.println("no files in " + target);
			System.exit(1);
		}
		Arrays.sort(files);

		for (File f : files) {
			try {
				String before = Files.readString(f.toPath());
				String after = MatrixJsonFormatter.format(before);
				if (!before.equals(after)) {
					Files.writeString(f.toPath(), after);
					changed++;
				} else {
					unchanged++;
				}
			} catch (Exception e) {
				errors++;
				System.err.printf("ERROR %s: %s%n", f.getName(), e.getMessage());
			}
		}

		System.out.printf("ReformatSchemes: %d changed, %d unchanged, %d errors (total %d files)%n",
				changed, unchanged, errors, files.length);
		if (errors > 0) System.exit(2);
	}
}
