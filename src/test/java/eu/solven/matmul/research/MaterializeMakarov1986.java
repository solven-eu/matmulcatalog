package eu.solven.matmul.research;

import java.nio.file.Files;
import java.nio.file.Path;

import eu.solven.matmul.NonBilinearAlgorithm;
import eu.solven.matmul.Verifier;
import eu.solven.matmul.catalog.Lineage;
import eu.solven.matmul.catalog.SparseNonBilinearWriter;
import eu.solven.matmul.catalog.SparseNonBilinearWriter.Metadata;
import eu.solven.matmul.papers.makarov1986.Makarov22;

/**
 * Materialise Makarov 1986's explicit {@code ⟨3,3,3⟩ = 22} non-bilinear
 * commutative scheme (one better than Laderman 1976's bilinear
 * ⟨3,3,3⟩=23).
 *
 * <p>Source: Makarov 1986, "An algorithm for multiplication of 3×3
 * matrices", <em>USSR Comput. Math. Math. Phys.</em> 26(2):293–294.
 * Cross-checked against the Russian original (Zh. Vychisl. Mat. i
 * Mat. Fiz. 26(2)) — Islam 2009 §3.3.1's transcription has a
 * single-index typo in γ18 (see {@link Makarov22}).</p>
 *
 * <p>Distinct from Rosowski Corollary 1's {@code ⟨3,3,3⟩=21} (also on
 * disk) — Makarov 22 sits between Laderman 23 and Rosowski 21
 * chronologically, preserving the historical record.</p>
 */
public final class MaterializeMakarov1986 {

	private MaterializeMakarov1986() {}

	private static final String SCHEMES_DIR = "src/main/resources/schemes";

	public static void main(String[] args) throws Exception {
		NonBilinearAlgorithm alg = Makarov22.buildDefault();
		if (!Verifier.isExactNonBilinear(alg)) {
			throw new IllegalStateException("Makarov22.build() did not produce an exact scheme — "
					+ "residual: " + Verifier.residualNonBilinear(alg));
		}
		Path dir = Path.of(SCHEMES_DIR, "section3");
		Files.createDirectories(dir);
		Path file = dir.resolve("makarov_1986-3x3x3_m22_a105.json");
		Metadata meta = new Metadata(
				"Q",
				/*commutative*/ true,
				"Makarov 1986",
				"Makarov 1986",
				/*discovery*/ true,
				1986,
				"Makarov 1986 — explicit ⟨3,3,3⟩=22 non-bilinear, commutative. "
						+ "USSR Comput. Math. Math. Phys. 26(2):293–294, "
						+ "DOI:10.1016/0041-5553(86)90203-X. "
						+ "Cross-checked against Russian original (mathnet.ru "
						+ "zvmmf4056); Islam 2009 §3.3.1's transcription has a "
						+ "single-index typo in γ18 (b32 → b23).");
		Lineage.Node lineage = new Lineage.Atom("Makarov22()");
		SparseNonBilinearWriter.write(alg, file.toFile(), meta, lineage);
		System.out.printf("Wrote %s (rank=%d, additions=%d)%n",
				file, alg.r, SparseNonBilinearWriter.additionCount(alg));
	}
}
