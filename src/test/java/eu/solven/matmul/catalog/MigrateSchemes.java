package eu.solven.matmul.catalog;

import java.io.File;
import java.io.IOException;
import java.nio.file.Files;
import java.nio.file.Path;
import java.nio.file.StandardCopyOption;
import java.util.Arrays;
import java.util.regex.Matcher;
import java.util.regex.Pattern;

import eu.solven.matmul.ComplexNonCubicBilinearAlgorithm;
import eu.solven.matmul.NonCubicBilinearAlgorithm;
import eu.solven.matmul.catalog.SchemeIO;

/**
 * One-shot migration: walks {@code src/main/resources/schemes/} (flat layout)
 * and reorganizes into per-max-dim subdirectories
 * ({@code section2/}, {@code section3/}, …, {@code section32/}). For
 * dense schemes with {@code max(n,m,p) ≥ SchemeIO.SPARSE_DIM_THRESHOLD} the
 * file is also converted to the sparse {@code u_sparse}/{@code v_sparse}/
 * {@code w_sparse} representation (≈5× smaller for typical schemes).
 *
 * <p>Idempotent: re-running after a clean migration is a no-op. Run via:</p>
 * <pre>
 *   mvn -q test-compile
 *   java -cp target/classes:target/test-classes:&dollar;CLASSPATH \
 *        eu.solven.matmul.catalog.MigrateSchemes
 * </pre>
 */
public class MigrateSchemes {

	private static final File ROOT = new File("src/main/resources/schemes");
	private static final Pattern NAME = Pattern.compile(
			"(?<source>.*)_(?<n>\\d+)x(?<m>\\d+)x(?<p>\\d+)_(?:r|m)(?<rank>\\d+)[^/]*\\.json");

	public static void main(String[] args) throws IOException {
		File[] flat = ROOT.listFiles((d, n) -> n.endsWith(".json"));
		if (flat == null) {
			System.err.println("no JSON files at " + ROOT);
			System.exit(1);
		}
		Arrays.sort(flat);

		int moved = 0, converted = 0, errors = 0;
		for (File src : flat) {
			Matcher m = NAME.matcher(src.getName());
			if (!m.matches()) {
				System.err.println("skip (name mismatch): " + src.getName());
				continue;
			}
			int n = Integer.parseInt(m.group("n"));
			int mm = Integer.parseInt(m.group("m"));
			int p = Integer.parseInt(m.group("p"));
			int maxDim = Math.max(n, Math.max(mm, p));

			Path targetDir = ROOT.toPath().resolve("section" + maxDim);
			Files.createDirectories(targetDir);
			File dst = targetDir.resolve(src.getName()).toFile();

			try {
				if (shouldConvertToSparse(src, maxDim)) {
					// Load + re-write through SchemeIO (auto-selects sparse).
					tools.jackson.databind.JsonNode root = SchemeIO.parseJson(src);
					NonCubicBilinearAlgorithm alg = SchemeIO.read(root);
					SchemeIO.write(alg, dst);
					Files.delete(src.toPath());
					converted++;
				} else {
					Files.move(src.toPath(), dst.toPath(), StandardCopyOption.REPLACE_EXISTING);
					moved++;
				}
			} catch (Exception e) {
				System.err.printf("ERROR %s → %s: %s%n", src.getName(), dst, e.getMessage());
				errors++;
			}
		}

		System.out.printf("Migration: moved=%d converted=%d errors=%d%n",
				moved, converted, errors);
		if (errors > 0) System.exit(2);
	}

	/**
	 * Convert iff max-dim ≥ threshold AND not already sparse AND not
	 * {@code _reduced} (those have their own sparse-list format) AND
	 * not complex (complex schemes use the {@code [re, im]} pair encoding).
	 */
	private static boolean shouldConvertToSparse(File src, int maxDim) throws IOException {
		if (maxDim < SchemeIO.SPARSE_DIM_THRESHOLD) return false;
		if (src.getName().contains("_reduced")) return false;
		try {
			if (SchemeIO.isComplex(SchemeIO.parseJson(src))) return false;
		} catch (IOException ignored) {
			return false;
		}
		String body = Files.readString(src.toPath());
		return !body.contains("\"u_sparse\"");
	}

	/**
	 * Walk all JSON files under {@link #ROOT}, recursing into subdirectories.
	 * Useful for code that previously did a flat {@code listFiles} — switch to
	 * this once migration has run.
	 */
	public static File[] allSchemes() throws IOException {
		try (var s = Files.walk(ROOT.toPath())) {
			return s.filter(p -> p.toString().endsWith(".json"))
					.map(Path::toFile)
					.sorted()
					.toArray(File[]::new);
		}
	}

	// Touch ComplexNonCubicBilinearAlgorithm/NonCubicBilinearAlgorithm so static
	// init doesn't bite us downstream when the SchemeIO factory paths run.
	@SuppressWarnings("unused")
	private static final Class<?>[] KEEP = {
			ComplexNonCubicBilinearAlgorithm.class, NonCubicBilinearAlgorithm.class };
}
