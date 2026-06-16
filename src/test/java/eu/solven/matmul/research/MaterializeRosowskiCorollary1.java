package eu.solven.matmul.research;

import java.nio.file.Files;
import java.nio.file.Path;

import eu.solven.matmul.NonBilinearAlgorithm;
import eu.solven.matmul.Verifier;
import eu.solven.matmul.catalog.Lineage;
import eu.solven.matmul.catalog.SparseNonBilinearWriter;
import eu.solven.matmul.catalog.SparseNonBilinearWriter.Metadata;
import eu.solven.matmul.papers.rosowski2019.Rosowski21;

/**
 * One-shot: materialise the explicit Rosowski 2019 <strong>Corollary 1</strong>
 * non-bilinear commutative scheme {@code ⟨3,3,3⟩ = 21}. Writes to
 * {@code section3/rosowski_2019_cor1-3x3x3_r21.json}.
 *
 * <p>The bound {@code R_c(⟨3,3,3⟩) ≤ 21} is also derivable from Rosowski
 * Theorem 3 (odd-{@code n} formula) — see
 * {@link eu.solven.matmul.papers.rosowski2019.RosowskiBound}. The Cor 1
 * construction is the explicit set of 21 products published in the paper
 * (Section 5). It coincides in rank with the {@code n = 3} specialisation
 * of {@link eu.solven.matmul.papers.rosowski2019.RosowskiAlgorithm1} (which
 * is the transpose-involution construction for {@code ⟨n,3,3⟩ = 6n+3}),
 * but the two schemes are distinct: different product indexing, different
 * algebraic structure inside the rank-1 atoms.</p>
 *
 * <p>Keeping BOTH schemes in the catalog documents the alternative
 * construction (per CLAUDE.md's "research-grade catalog of fast matrix
 * multiplication algorithms" goal — historical chronology + alternative
 * derivations are valuable even at the same rank).</p>
 *
 * <p>Run:</p>
 * <pre>
 *   mvn -q -o exec:java -Dexec.classpathScope=test \
 *       -Dexec.mainClass=eu.solven.matmul.research.MaterializeRosowskiCorollary1
 * </pre>
 */
public final class MaterializeRosowskiCorollary1 {

	private MaterializeRosowskiCorollary1() {}

	// Formula constructors live under schemes/constructed/ (alongside HK71).
	private static final String SCHEMES_DIR = "src/main/resources/schemes/constructed";

	public static void main(String[] args) throws Exception {
		NonBilinearAlgorithm alg = Rosowski21.build();
		if (!Verifier.isExactNonBilinear(alg)) {
			throw new IllegalStateException("Rosowski21.build() did not produce an exact scheme");
		}
		Path dir = Path.of(SCHEMES_DIR, "section3");
		Files.createDirectories(dir);
		Path file = dir.resolve("rosowski_2019_cor1-3x3x3_m21.json");
		Metadata meta = new Metadata(
				"Q",
				/*commutative*/ true,
				"Rosowski 2019",
				"Rosowski 2019",
				/*discovery*/ true,
				2019,
				"arXiv:1904.07683 Corollary 1 — explicit ⟨3,3,3⟩=21 commutative");
		Lineage.Node lineage = new Lineage.Atom("Rosowski21()");
		SparseNonBilinearWriter.write(alg, file.toFile(), meta, lineage);
		System.out.printf("Wrote %s (rank=%d, additions=%d)%n",
				file, alg.r, SparseNonBilinearWriter.additionCount(alg));
	}
}
