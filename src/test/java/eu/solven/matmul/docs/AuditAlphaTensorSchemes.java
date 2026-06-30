package eu.solven.matmul.docs;
import eu.solven.matmul.catalog.FieldAwareLookup;

import java.io.File;
import java.nio.file.Files;
import java.nio.file.Path;
import java.util.ArrayList;
import java.util.List;
import java.util.stream.Stream;

import eu.solven.matmul.NonCubicBilinearAlgorithm;
import eu.solven.matmul.verifiers.Verifier;
import eu.solven.matmul.catalog.SchemeIO;

/**
 * Audit every {@code alphatensor-Z_*.json} and {@code alphatensor-F2_*.json}
 * file in the catalog: spot-check each. Any scheme that FAILS is corrupt
 * (wrong convention, bad import, or simply broken in the source) and
 * should be quarantined so {@code FieldAwareLookup} doesn't return it.
 */
public final class AuditAlphaTensorSchemes {

	public static void main(String[] args) throws Exception {
		Path root = Path.of("src/main/resources/schemes");
		List<Path> bad = new ArrayList<>();
		int total = 0, passed = 0;
		try (Stream<Path> walk = Files.walk(root)) {
			for (Path p : (Iterable<Path>) walk::iterator) {
				String n = p.getFileName().toString();
				if (!n.endsWith(".json")) continue;
				if (!n.toLowerCase().startsWith("alphatensor")) continue;
				total++;
				NonCubicBilinearAlgorithm alg;
				try {
					alg = SchemeIO.read(p.toFile());
				} catch (Exception e) {
					System.out.println("READ-ERROR " + n + " : " + e.getMessage());
					bad.add(p);
					continue;
				}
				boolean isF2 = n.contains("F2");
				boolean ok;
				if (isF2) {
					// F2 schemes use modular arithmetic — real-number spot-check
					// would give false negatives. Use the F2-specific verifier.
					ok = Verifier.isExactNonCubicF2(alg);
				} else {
					ok = Verifier.passesRandomMatmulSpotCheck(alg);
				}
				if (ok) {
					passed++;
				} else {
					System.out.println("BROKEN " + n + " rank=" + alg.r);
					bad.add(p);
				}
			}
		}
		System.out.println();
		System.out.printf("Audited %d AT files: %d PASS, %d BROKEN%n", total, passed, bad.size());
	}
}
