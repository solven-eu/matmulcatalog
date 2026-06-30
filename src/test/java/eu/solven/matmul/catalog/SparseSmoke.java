package eu.solven.matmul.catalog;

import java.io.File;

import eu.solven.matmul.NonCubicBilinearAlgorithm;
import eu.solven.matmul.verifiers.Verifier;
import eu.solven.matmul.catalog.SchemeIO;

/**
 * Smoke test: pick one max-dim-9 scheme, round-trip via the new sparse format,
 * confirm read+verify still works. For diagnostic use; not a unit test.
 */
public class SparseSmoke {
	public static void main(String[] args) throws Exception {
		File in = eu.solven.matmul.catalog.SchemeResolver.byHint("src/main/resources/schemes/derived_laderman2-9x9x9_r529.json");
		NonCubicBilinearAlgorithm orig = SchemeIO.read(in);
		System.out.printf("loaded %s: ⟨%d,%d,%d⟩ r=%d residual=%.2e%n",
				in.getName(), orig.n, orig.m, orig.p, orig.r, Verifier.residualNonCubic(orig));

		File tmp = File.createTempFile("sparse-smoke-", ".json");
		tmp.deleteOnExit();
		SchemeIO.write(orig, tmp);
		long bytesIn = in.length(), bytesOut = tmp.length();
		System.out.printf("wrote sparse to %s: %d → %d bytes (%.1f%%)%n",
				tmp.getName(), bytesIn, bytesOut, 100.0 * bytesOut / bytesIn);

		NonCubicBilinearAlgorithm round = SchemeIO.read(tmp);
		System.out.printf("re-loaded: ⟨%d,%d,%d⟩ r=%d residual=%.2e%n",
				round.n, round.m, round.p, round.r, Verifier.residualNonCubic(round));

		System.out.println("\n--- head of sparse output ---");
		try (var lines = java.nio.file.Files.lines(tmp.toPath())) {
			lines.limit(15).forEach(System.out::println);
		}
	}
}
