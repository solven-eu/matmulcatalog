package eu.solven.matmul.research;

import java.nio.file.Files;
import java.nio.file.Path;

import eu.solven.matmul.NonBilinearAlgorithm;
import eu.solven.matmul.Verifier;
import eu.solven.matmul.catalog.Lineage;
import eu.solven.matmul.catalog.SparseNonBilinearWriter;
import eu.solven.matmul.catalog.SparseNonBilinearWriter.Metadata;
import eu.solven.matmul.papers.rosowski2019.RosowskiAlgorithm1;

/**
 * One-shot: materialise Rosowski 2019 Algorithm 1 ({@code ⟨n,3,3⟩ = 6n+3},
 * commutative non-bilinear) for {@code n = 2..32}, writing one JSON
 * scheme file per format. Each file is verified before being written.
 *
 * <p>Writing goes through {@link SparseNonBilinearWriter} so the JSON
 * layout + coefficient formatting + {@code additions} count are
 * consistent with other non-bilinear materialisers
 * ({@link MaterializeWaksman1970}, future {@code MaterializeRosowskiThm2}, …).</p>
 *
 * <p>For {@code n > 32} the formula bound {@code 6n+3} is already exposed
 * via {@link eu.solven.matmul.papers.rosowski2019.RosowskiBound} and
 * emitted to {@code docs/derived-bounds.json}.</p>
 *
 * <p>Run:</p>
 * <pre>
 *   mvn -q -o exec:java -Dexec.classpathScope=test \
 *       -Dexec.mainClass=eu.solven.matmul.research.MaterializeRosowskiAlgorithm1
 * </pre>
 */
public final class MaterializeRosowskiAlgorithm1 {

	private MaterializeRosowskiAlgorithm1() {}

	private static final int MAX_N = 32;
	// Formula constructors live under schemes/constructed/ (alongside HK71).
	private static final String SCHEMES_DIR = "src/main/resources/schemes/constructed";

	public static void main(String[] args) throws Exception {
		int wrote = 0, skipped = 0;
		for (int n = 2; n <= MAX_N; n++) {
			NonBilinearAlgorithm alg = RosowskiAlgorithm1.build(n);
			if (!Verifier.passesRandomMatmulSpotCheckNB(alg)) {
				System.err.printf("SKIP n=%d: spot-check failed%n", n);
				skipped++;
				continue;
			}
			int maxDim = Math.max(n, 3);
			Path dir = Path.of(SCHEMES_DIR, "section" + maxDim);
			Files.createDirectories(dir);
			Path file = dir.resolve("rosowski_2019_alg1_" + n + "x3x3_r" + (6 * n + 3) + ".json");
			Metadata meta = new Metadata(
					"Q",
					/*commutative*/ true,
					"Rosowski 2019",
					"Rosowski 2019",
					/*discovery*/ true,
					2019,
					"arXiv:1904.07683 Algorithm 1 (⟨" + n + ",3,3⟩=" + (6 * n + 3) + ")");
			Lineage.Node lineage = new Lineage.Atom("RosowskiAlgorithm1(n=" + n + ")");
			SparseNonBilinearWriter.write(alg, file.toFile(), meta, lineage);
			wrote++;
		}
		System.out.printf("Materialised %d Rosowski Algorithm 1 schemes (n=2..%d), skipped %d%n",
				wrote, MAX_N, skipped);
	}
}
