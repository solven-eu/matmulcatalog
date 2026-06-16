package eu.solven.matmul.research;

import java.nio.file.Files;
import java.nio.file.Path;

import eu.solven.matmul.NonBilinearAlgorithm;
import eu.solven.matmul.Verifier;
import eu.solven.matmul.catalog.CatalogPolicy;
import eu.solven.matmul.catalog.Lineage;
import eu.solven.matmul.catalog.SparseNonBilinearWriter;
import eu.solven.matmul.catalog.SparseNonBilinearWriter.Metadata;
import eu.solven.matmul.papers.rosowski2019.RosowskiTheorem2;

/**
 * One-shot: materialise Rosowski 2019 Theorem 2 ({@code ⟨2,2,p⟩ = 3p+1},
 * divisions-free non-bilinear commutative) for {@code p = 3 ..
 * MATERIALISE_MAX_DIM}, one JSON scheme file per format. Each file is verified
 * before being written.
 *
 * <p>Why {@code p ≥ 3}: {@code ⟨2,2,2⟩} gives rank 7 — tied with (and worse in
 * field coverage than) the bilinear all-field Strassen ⟨2,2,2⟩=7 already in the
 * catalog, so it adds no information. From {@code p = 3} the family
 * (⟨2,2,3⟩=10, ⟨2,2,4⟩=13, …) is the commutative non-bilinear record this
 * project tracks.</p>
 *
 * <p><strong>Field = Z</strong>: Theorem 2's whole point over Waksman 1970 [19]
 * / Islam 2009 [10] (which establish the same {@code 3p+1} rank but with
 * divisions by 2, needing 2 invertible) is that it is DIVISIONS-FREE — all
 * coefficients are {@code ±1}, so it is valid over any commutative ring
 * including {@code Z} (and {@code F₂}, {@code F₃}). The rank itself is prior
 * art, so {@code discovery=false} and {@code attribution_for_rank=Waksman
 * 1970}; the divisions-free construction is Rosowski's.</p>
 *
 * <p>The closed-form bound {@code 3p+1} is also exposed via
 * {@link eu.solven.matmul.papers.rosowski2019.RosowskiBound#commutativeBoundBilinear}
 * and emitted to {@code docs/derived-bounds.json}.</p>
 *
 * <p>Run:</p>
 * <pre>
 *   mvn -q -o exec:java -Dexec.classpathScope=test \
 *       -Dexec.mainClass=eu.solven.matmul.research.MaterializeRosowskiTheorem2
 * </pre>
 */
public final class MaterializeRosowskiTheorem2 {

	private MaterializeRosowskiTheorem2() {}

	private static final int MIN_P = 3;
	private static final int MAX_P = CatalogPolicy.MATERIALISE_MAX_DIM;
	// Formula constructors live under schemes/constructed/ (alongside HK71).
	private static final String SCHEMES_DIR = "src/main/resources/schemes/constructed";

	public static void main(String[] args) throws Exception {
		int wrote = 0, skipped = 0;
		for (int p = MIN_P; p <= MAX_P; p++) {
			int rank = 3 * p + 1;
			NonBilinearAlgorithm alg = RosowskiTheorem2.build22p(p);
			if (!Verifier.passesRandomMatmulSpotCheckNB(alg)) {
				System.err.printf("SKIP ⟨2,2,%d⟩: spot-check failed%n", p);
				skipped++;
				continue;
			}
			int maxDim = Math.max(2, p);
			Path dir = Path.of(SCHEMES_DIR, "section" + maxDim);
			Files.createDirectories(dir);
			Path file = dir.resolve("rosowski_2019_thm2-2x2x" + p + "_m" + rank + ".json");
			Metadata meta = new Metadata(
					"Z",
					/*commutative*/ true,
					"Rosowski 2019",
					/*attributionForRank*/ "Waksman 1970",
					/*discovery*/ false,
					2019,
					"arXiv:1904.07683 Theorem 2 (⟨2,2," + p + "⟩=" + rank
							+ ", divisions-free; rank due to Waksman 1970 [19] / Islam 2009 [10])");
			Lineage.Node lineage = new Lineage.Atom("RosowskiTheorem2(l=2,n=2,m=" + p + ")");
			SparseNonBilinearWriter.write(alg, file.toFile(), meta, lineage);
			wrote++;
		}
		System.out.printf("Materialised %d Rosowski Theorem 2 schemes (⟨2,2,p⟩, p=%d..%d), skipped %d%n",
				wrote, MIN_P, MAX_P, skipped);
	}
}
