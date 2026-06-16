package eu.solven.matmul.docs.migrate;

import java.io.File;
import java.nio.file.Files;
import java.nio.file.Path;
import java.util.List;

import eu.solven.matmul.NonCubicBilinearAlgorithm;
import eu.solven.matmul.Verifier;
import eu.solven.matmul.catalog.SchemeIO;
import tools.jackson.databind.JsonNode;

/**
 * Re-emit scheme files with <strong>exact rational</strong> coefficients: rounded
 * decimals like {@code 0.058823529411764705} become the string {@code "1/17"}.
 * The JSON is more expressive (and forward-compatible with a future
 * rational/symbolic compute backend) while parsing to the same value today.
 *
 * <p>Surgical + metadata-preserving: delegates to
 * {@link SchemeIO#exactifyCoefficients(File, boolean)}, which rewrites only the
 * floating-point coefficient nodes and leaves every other field intact. Integer /
 * F₂ / ternary schemes are untouched (no floating-point coeffs → no change).</p>
 *
 * <p>Optional first arg = a filename substring filter (e.g. {@code schwartz_zwecher});
 * default processes the whole catalog. Pass {@code --apply} to write (default: dry
 * run). On {@code --apply} each rewritten file is re-read and spot-checked so a bad
 * rewrite can never land silently.</p>
 *
 * <pre>
 *   mvn -q -o exec:java -Dexec.mainClass=eu.solven.matmul.docs.migrate.ExactifyRationalCoeffs \
 *       -Dexec.args="schwartz_zwecher --apply"
 * </pre>
 */
public final class ExactifyRationalCoeffs {
	private ExactifyRationalCoeffs() {}

	public static void main(String[] args) throws Exception {
		List<String> a = List.of(args);
		boolean apply = a.contains("--apply");
		String filter = a.stream().filter(s -> !s.startsWith("--")).findFirst().orElse("");
		Path root = Path.of("src/main/resources/schemes");

		List<Path> files;
		try (var s = Files.walk(root)) {
			files = s.filter(p -> p.getFileName().toString().endsWith(".json"))
					.filter(p -> filter.isEmpty() || p.getFileName().toString().contains(filter))
					.sorted().toList();
		}

		int changed = 0, skipped = 0, verifyFail = 0, scanned = 0;
		for (Path p : files) {
			File f = p.toFile();
			JsonNode rn = SchemeIO.parseJson(f);
			// Only standard real bilinear files carry rational coeffs we can re-verify.
			boolean complex = rn.has("complex") && rn.get("complex").asBoolean();
			boolean hasFactors = rn.has("u_sparse") || rn.has("u");
			if (SchemeIO.isStub(rn) || complex || !hasFactors) {
				skipped++;
				continue;
			}
			scanned++;
			boolean wouldChange = SchemeIO.exactifyCoefficients(f, apply);
			if (!wouldChange) {
				continue;
			}
			changed++;
			String tag = apply ? "EXACT" : "(dry)";
			if (apply) {
				try {
					NonCubicBilinearAlgorithm back = SchemeIO.readBilinear(f);
					boolean ok = Verifier.passesRandomMatmulSpotCheck(back);
					if (!ok) { verifyFail++; tag = "VERIFY-FAIL"; }
				} catch (RuntimeException e) {
					verifyFail++;
					tag = "READ-FAIL(" + e.getMessage() + ")";
				}
			}
			System.out.printf("%-6s %s%n", tag, f.getName());
		}
		System.out.printf("%n%s: %d scanned, %d rewritten, %d skipped (stub/complex/integer-only), %d verify-fail%n",
				apply ? "APPLIED" : "DRY RUN", scanned, changed, skipped, verifyFail);
		if (!apply) {
			System.out.println("Re-run with --apply to write.");
		}
	}
}
