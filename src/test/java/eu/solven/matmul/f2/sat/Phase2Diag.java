package eu.solven.matmul.f2.sat;

import java.util.List;
import java.util.Optional;

import eu.solven.matmul.BilinearAlgorithm;
import eu.solven.matmul.f2.sat.KissatSolver;
import eu.solven.matmul.f2.sat.SatMatmulPipeline;
import eu.solven.matmul.f2.sat.Z2CnfEncoder;

/**
 * Phase 2 calibration: dense ⟨3,3,3⟩ over Z/2, varying r. Reports CNF size
 * and SAT4J wall-clock per rank so we know where the boundary of "feasible
 * without BreakID/kissat" lies.
 *
 * Hard wall-clock cap per rank (configurable below). Skips and prints a
 * TIMEOUT marker rather than hanging indefinitely.
 */
public class Phase2Diag {

	public static void main(String[] args) {
		System.out.println("--- Dense ⟨3,3,3⟩ over Z/2 (solver=kissat) ---");
		int[][][] target = SatMatmulPipeline.z2DenseMatmulTensor(3);
		for (int r : new int[] { 30, 27, 25, 23, 22 }) {
			run(target, r);
		}
	}

	private static void run(int[][][] target, int r) {
		Z2CnfEncoder encoder = new Z2CnfEncoder(3, r, target);
		List<int[]> clauses = encoder.getClauses();
		System.out.printf("r=%2d: CNF %d vars, %d clauses — ", r, encoder.getVarCount(), clauses.size());
		System.out.flush();

		long start = System.nanoTime();
		Optional<boolean[]> model = KissatSolver.solve(encoder.getVarCount(), clauses);
		long ms = (System.nanoTime() - start) / 1_000_000;
		System.out.printf("%s in %d ms", model.isPresent() ? "SAT" : "UNSAT", ms);
		if (model.isPresent()) {
			BilinearAlgorithm alg = encoder.decode(model.get());
			System.out.printf("  verify=%s%n", SatMatmulPipeline.verifyZ2(alg, target) ? "OK" : "FAIL");
		} else {
			System.out.println();
		}
		System.out.flush();
	}
}
