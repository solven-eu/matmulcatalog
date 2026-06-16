package eu.solven.matmul.f2.sat;

import java.util.Optional;

import eu.solven.matmul.f2.sat.BreakIdBridge;
import eu.solven.matmul.f2.sat.KissatSolver;
import eu.solven.matmul.f2.sat.SatMatmulPipeline;

/**
 * Phase 1.6 standalone runner: search for a Z/2 rank-r decomposition of
 * the non-cubic ⟨2, 3, 3⟩ matmul tensor.
 *
 * <p>The research target is {@code R_{Z/2}(⟨2,3,3⟩) = 14} (Hopcroft–Kerr
 * 1971 gives ≥ 15 over fields of characteristic ≠ 2; over Z/2 the bound
 * is 14 ≤ r ≤ 15 — a decisive UNSAT at r=14 would prove the bound is 15).</p>
 *
 * <p>Run via:
 * <pre>
 *   mvn -q test-compile
 *   java -cp target/classes:target/test-classes \
 *        eu.solven.matmul.f2.sat.Phase16Runner [r=15]
 * </pre>
 * Default {@code r = 15} is SAT (verifies the encoder). Pass {@code r = 14}
 * for the research target (likely UNSAT, slow).</p>
 */
public final class Phase16Runner {

	public static void main(String[] args) {
		int r = args.length > 0 ? Integer.parseInt(args[0]) : 15;
		int n = 2, m = 3, p = 3;

		System.out.printf("Phase 1.6: searching Z/2 decomposition of ⟨%d,%d,%d⟩ at r=%d%n", n, m, p, r);
		System.out.printf("  dimU = %d, dimV = %d, dimW = %d%n", n * m, m * p, n * p);
		System.out.printf("  solver: kissat=%b, BreakID=%b%n",
				KissatSolver.isAvailable(), BreakIdBridge.isAvailable());

		int[][][] target = SatMatmulPipeline.z2DenseMatmulTensor(n, m, p);
		long t0 = System.currentTimeMillis();
		Optional<double[][][]> result =
				SatMatmulPipeline.findZ2NonCubicDecomposition(n, m, p, r, target);
		long dt = System.currentTimeMillis() - t0;

		if (result.isPresent()) {
			boolean ok = SatMatmulPipeline.verifyZ2NonCubic(result.get(), target);
			System.out.printf("SAT in %d ms (verified=%b)%n", dt, ok);
			System.out.printf("U shape [%d][%d], V shape [%d][%d], W shape [%d][%d]%n",
					result.get()[0].length, result.get()[0][0].length,
					result.get()[1].length, result.get()[1][0].length,
					result.get()[2].length, result.get()[2][0].length);
		} else {
			System.out.printf("UNSAT in %d ms — R_{Z/2}(⟨%d,%d,%d⟩) > %d%n", dt, n, m, p, r);
		}
	}
}
