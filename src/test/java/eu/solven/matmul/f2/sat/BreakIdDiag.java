package eu.solven.matmul.f2.sat;

import java.util.List;
import java.util.Optional;

import eu.solven.matmul.BilinearAlgorithm;
import eu.solven.matmul.f2.sat.BreakIdBridge;
import eu.solven.matmul.f2.sat.Cnf;
import eu.solven.matmul.f2.sat.CryptoMiniSatSolver;
import eu.solven.matmul.f2.sat.KissatSolver;
import eu.solven.matmul.f2.sat.SatMatmulPipeline;
import eu.solven.matmul.f2.sat.Z2CnfEncoder;

/**
 * Diagnostic: dump CNF size and solver wall-clock for the BreakID pipeline
 * variants on dense ⟨2,2,2⟩ Z/2 at r=7.
 *
 * Compares:
 *   1) bare CNF (no hand-lex)            — for reference (expected to hang)
 *   2) hand-lex only                     — Phase 1 baseline
 *   3) BreakID alone (no hand-lex)       — BreakID does all symmetry breaking
 *   4) BreakID + hand-lex                — redundant but checks for harm
 *
 * Uses kissat as the underlying solver so we can compare fairly across all
 * four variants (SAT4J chokes on variants 1 and 3 without symmetry breaking).
 */
public class BreakIdDiag {

	public static void main(String[] args) throws Exception {
		System.out.println("BreakID available: " + BreakIdBridge.isAvailable());
		System.out.println("Kissat available: " + KissatSolver.isAvailable());
		System.out.println("CryptoMiniSat available: " + CryptoMiniSatSolver.isAvailable());
		System.out.println();

		// Phase 1 validation: dense ⟨2,2,2⟩ Z/2 at r=7
		System.out.println("### Dense ⟨2,2,2⟩ Z/2 at r=7 ###");
		int[][][] t2 = SatMatmulPipeline.z2DenseMatmulTensor(2);
		runAll("baseline (hand-lex only)", new Z2CnfEncoder(2, 7, t2, true), t2);
		runAll("BreakID only (no hand-lex)", encodeAndBreak(2, 7, t2, false), t2);
		runAll("hand-lex + BreakID", encodeAndBreak(2, 7, t2, true), t2);

		// Phase 2 spot: dense ⟨3,3,3⟩ Z/2 at r=27 (naive rank, was 3.9 min with kissat alone)
		System.out.println("### Dense ⟨3,3,3⟩ Z/2 at r=27 ###");
		int[][][] t3 = SatMatmulPipeline.z2DenseMatmulTensor(3);
		runAll("baseline (hand-lex only)", new Z2CnfEncoder(3, 27, t3, true), t3);
		runAll("hand-lex + BreakID", encodeAndBreak(3, 27, t3, true), t3);
	}

	private static void runAll(String label, Object cnfHolder, int[][][] target) {
		run(label + " / kissat", cnfHolder, target, "kissat");
		run(label + " / CMS", cnfHolder, target, "cms");
		System.out.println();
	}

	private static void run(String label, Object cnfHolder, int[][][] target, String solver) {
		int varCount;
		List<int[]> clauses;
		Z2CnfEncoder encoder;
		if (cnfHolder instanceof Z2CnfEncoder e) {
			encoder = e;
			varCount = e.getVarCount();
			clauses = e.getClauses();
		} else {
			BreakIdRun run = (BreakIdRun) cnfHolder;
			encoder = run.encoder;
			varCount = run.varCount;
			clauses = run.clauses;
		}
		long start = System.nanoTime();
		Optional<boolean[]> model = solver.equals("kissat")
				? KissatSolver.solve(varCount, clauses)
				: CryptoMiniSatSolver.solve(varCount, clauses);
		long ms = (System.nanoTime() - start) / 1_000_000;
		System.out.printf("[%-50s] %d vars %5d clauses — %s in %d ms",
				label, varCount, clauses.size(),
				model.isPresent() ? "SAT" : "UNSAT", ms);
		if (model.isPresent()) {
			BilinearAlgorithm alg = encoder.decode(model.get());
			System.out.printf("  verify=%s%n", SatMatmulPipeline.verifyZ2(alg, target) ? "OK" : "FAIL");
		} else {
			System.out.println();
		}
		System.out.flush();
	}

	private static BreakIdRun encodeAndBreak(int n, int r, int[][][] target, boolean handLex) throws Exception {
		Z2CnfEncoder encoder = new Z2CnfEncoder(n, r, target, handLex);
		Cnf.ReadResult augmented = BreakIdBridge.preprocess(encoder.getVarCount(), encoder.getClauses());
		return new BreakIdRun(encoder, augmented.varCount, augmented.clauses);
	}

	private record BreakIdRun(Z2CnfEncoder encoder, int varCount, List<int[]> clauses) {}
}
