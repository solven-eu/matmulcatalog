package eu.solven.matmul.research;

import java.nio.file.Files;
import java.nio.file.Path;

import eu.solven.matmul.NonCubicBilinearAlgorithm;
import eu.solven.matmul.Verifier;
import eu.solven.matmul.catalog.Lineage;
import eu.solven.matmul.catalog.SchemeIO;
import eu.solven.matmul.papers.dis2009.PanTrilinearAggregation;

/**
 * Materialise the Pan/Islam trilinear-aggregation cubic ⟨n,n,n⟩ scheme
 * (DIS09 Appendix Lemma 4, even + odd cases) for every {@code n} where
 * the formula could conceivably win against our current SOTA — meaning,
 * for now, the full range {@code n = 4..32} (we filter post-hoc against
 * existing catalog ranks in a downstream comparison pass).
 *
 * <p>Each scheme is verified via
 * {@link Verifier#passesRandomMatmulSpotCheck} (sufficient to
 * catch any indexing/coefficient bug; the construction itself is unit-
 * tested for symbolic exactness in {@code TestPanTrilinearAggregationBuild}
 * for {@code n ∈ {3, 4, 5, 6}}).</p>
 *
 * <p>Files are tagged {@code "field": "Q"}, {@code "commutative": false}
 * (the construction is genuinely non-commutative — it lifts to recursive
 * matmul). Naming convention:
 * {@code dis09-Q_NxNxN_rRANK.json}.</p>
 *
 * <p>Run:</p>
 * <pre>
 *   mvn -q -o exec:java -Dexec.classpathScope=test \
 *       -Dexec.mainClass=eu.solven.matmul.research.MaterializeDIS09TA
 * </pre>
 */
public final class MaterializeDIS09TA {

	private MaterializeDIS09TA() {}

	private static final int MIN_N = 4;
	private static final int MAX_N = 32;
	private static final String SCHEMES_DIR = "src/main/resources/schemes";

	public static void main(String[] args) throws Exception {
		int wrote = 0, skipped = 0, failed = 0;
		long t0 = System.nanoTime();
		for (int n = MIN_N; n <= MAX_N; n++) {
			long tn = System.nanoTime();
			NonCubicBilinearAlgorithm alg;
			try {
				alg = PanTrilinearAggregation.build(n);
			} catch (RuntimeException e) {
				System.err.printf("BUILD FAILED n=%d: %s%n", n, e.getMessage());
				failed++;
				continue;
			}
			boolean ok = Verifier.passesRandomMatmulSpotCheck(alg);
			if (!ok) {
				System.err.printf("SKIP n=%d (rank %d): spot-check failed%n", n, alg.r);
				skipped++;
				continue;
			}
			Path dir = Path.of(SCHEMES_DIR, "section" + n);
			Files.createDirectories(dir);
			String fname = String.format("dis09-Q_%dx%dx%d_r%d.json", n, n, n, alg.r);
			Path file = dir.resolve(fname);
			Lineage.Node lineage = new Lineage.Atom("DIS09Lemma4(n=" + n + ")");
			SchemeIO.write(alg, file.toFile(), lineage);
			long stepMs = (System.nanoTime() - tn) / 1_000_000L;
			System.out.printf("  n=%d  r=%d  (%s case)  %dms%n",
					n, alg.r,
					n % 2 == 0 ? "even" : "odd",
					stepMs);
			wrote++;
		}
		long ms = (System.nanoTime() - t0) / 1_000_000L;
		System.out.printf("Materialised %d DIS09 TA schemes (skipped %d, failed %d) in %dms%n",
				wrote, skipped, failed, ms);
	}
}
