package eu.solven.matmul.research;

import java.nio.file.Files;
import java.nio.file.Path;

import eu.solven.matmul.NonBilinearAlgorithm;
import eu.solven.matmul.verifiers.Verifier;
import eu.solven.matmul.catalog.Lineage;
import eu.solven.matmul.catalog.SparseNonBilinearWriter;
import eu.solven.matmul.catalog.SparseNonBilinearWriter.Metadata;
import eu.solven.matmul.papers.waksman1970.Waksman1970;

/**
 * One-shot: materialise the Waksman 1970 generic commutative algorithm for
 * {@code ⟨n,n,n⟩} for {@code n = 2..32}, writing one JSON scheme file per
 * format. Each file is verified before being written.
 *
 * <p>Rank formula: {@code R_c(⟨n,n,n⟩) = (n²+2n−1)·⌊n/2⌋ + (n mod 2)·n²}.
 * Concrete values: 2→7, 3→23, 4→46, 5→93, 6→138, 8→312, 16→2576, 32→21472.</p>
 *
 * <p>Files are tagged {@code "commutative": true} — Waksman's algorithm
 * exploits scalar commutativity inside each rank-1 product, so it does
 * NOT lift to recursive matmul over non-commutative rings (matrices). It
 * DOES give a valid scheme for scalar matmul over any commutative ring,
 * plus an upper bound on the <em>commutative</em> asymptotic exponent
 * {@code ω_c ≤ log_n(R_c(n))}.</p>
 *
 * <p>Run:</p>
 * <pre>
 *   mvn -q -o exec:java -Dexec.classpathScope=test \
 *       -Dexec.mainClass=eu.solven.matmul.research.MaterializeWaksman1970
 * </pre>
 */
public final class MaterializeWaksman1970 {

	private MaterializeWaksman1970() {}

	private static final int MAX_N = 32;
	// Formula constructors live under schemes/constructed/ (alongside HK71).
	private static final String SCHEMES_DIR = "src/main/resources/schemes/constructed";

	public static void main(String[] args) throws Exception {
		int wrote = 0, skipped = 0;
		long t0 = System.nanoTime();
		for (int n = 2; n <= MAX_N; n++) {
			long tn = System.nanoTime();
			NonBilinearAlgorithm alg = Waksman1970.build(n);
			// Full residualNonBilinear is O(n⁴·r) which becomes ~minutes at
			// n≥19. The fast spot-check at default samples is O(n²·r) and is
			// enough to catch any indexing/coefficient bug — the generator
			// is unit-tested for exactness up to n≤8 in TestWaksman1970, this
			// is a regression guard during materialisation.
			if (!Verifier.passesRandomMatmulSpotCheckNB(alg)) {
				System.err.printf("SKIP n=%d: spot-check failed%n", n);
				skipped++;
				continue;
			}
			int adds = SparseNonBilinearWriter.additionCount(alg);
			Path dir = Path.of(SCHEMES_DIR, "section" + n);
			Files.createDirectories(dir);
			String fname = String.format("waksman-1970_%dx%dx%d_r%d_a%d_commutative.json",
					n, n, n, alg.r, adds);
			Path file = dir.resolve(fname);
			Metadata meta = new Metadata(
					"Q",
					/*commutative*/ true,
					"Waksman 1970",
					"Waksman 1970",
					/*discovery*/ true,
					1970,
					String.format("A. Waksman, On Winograd's algorithm for inner products, "
							+ "IEEE Trans. Comp. C-19(4), 1970, ⟨%d,%d,%d⟩=%d commutative",
							n, n, n, alg.r));
			// Generator-style lineage: this scheme is the deterministic output
			// of `Waksman1970.build(n)`. The ref string conventionally encodes
			// the constructor + parameter so a future LineageReplayer with
			// Generator-leaf support can call the build() method directly.
			Lineage.Node lineage = new Lineage.Atom("Waksman1970(n=" + n + ")");
			SparseNonBilinearWriter.write(alg, file.toFile(), meta, lineage);
			long stepMs = (System.nanoTime() - tn) / 1_000_000L;
			System.out.printf("  n=%d  r=%d  a=%d  %dms%n", n, alg.r, adds, stepMs);
			wrote++;
		}
		long ms = (System.nanoTime() - t0) / 1_000_000L;
		System.out.printf("Materialised %d Waksman 1970 schemes, skipped %d, %dms%n",
				wrote, skipped, ms);
	}
}
