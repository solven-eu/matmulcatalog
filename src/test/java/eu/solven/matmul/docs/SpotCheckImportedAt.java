package eu.solven.matmul.docs;

import java.io.File;
import java.nio.file.Files;
import java.nio.file.Path;
import java.util.stream.Stream;

import eu.solven.matmul.NonCubicBilinearAlgorithm;
import eu.solven.matmul.Verifier;
import eu.solven.matmul.catalog.SchemeIO;

/**
 * Sanity: spot-check all 27 newly-imported AT files from
 * factorizations_r.npz via our Java verifier. Numpy verified them
 * already; this confirms our Java {@link SchemeIO} reads the format
 * correctly + the matmul holds end-to-end.
 */
public final class SpotCheckImportedAt {

	public static void main(String[] args) throws Exception {
		Path root = Path.of("src/main/resources/schemes");
		int total = 0, pass = 0, fail = 0;
		try (Stream<Path> walk = Files.walk(root)) {
			for (Path p : (Iterable<Path>) walk::iterator) {
				String name = p.getFileName().toString();
				if (!name.startsWith("alphatensor-Q_") || !name.endsWith("_aN.json")) continue;
				total++;
				try {
					NonCubicBilinearAlgorithm alg = SchemeIO.read(p.toFile());
					boolean ok = Verifier.passesRandomMatmulSpotCheck(alg);
					if (ok) {
						pass++;
					} else {
						fail++;
						System.out.println("FAIL " + name + " rank=" + alg.r);
					}
				} catch (Exception e) {
					fail++;
					System.out.println("FAIL-READ " + name + ": " + e.getMessage());
				}
			}
		}
		System.out.printf("Total: %d, Pass: %d, Fail: %d%n", total, pass, fail);
	}
}
