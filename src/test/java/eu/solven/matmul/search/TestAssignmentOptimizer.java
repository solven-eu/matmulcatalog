package eu.solven.matmul.search;

import static org.assertj.core.api.Assertions.assertThat;

import java.io.File;

import org.junit.jupiter.api.Test;

import eu.solven.matmul.NonCubicBilinearAlgorithm;
import eu.solven.matmul.catalog.FieldAwareLookup;
import eu.solven.matmul.catalog.Recombination.SotaResolver;
import eu.solven.matmul.catalog.SchemeIO;

/**
 * The partition+assignment optimiser must return the SAME global optimum as the
 * flat-composition {@link AllocationOptimizer} branch-and-bound — it is a
 * re-organisation of the exact search, not an approximation.
 */
public class TestAssignmentOptimizer {

	private static SotaResolver sota() {
		FieldAwareLookup lk = new FieldAwareLookup("R");
		return (p, q, r) -> {
			if (p == 0 || q == 0 || r == 0) return 0;
			if (p == 1) return q * r;
			if (q == 1) return p * r;
			if (r == 1) return p * q;
			int v = lk.findRank(p, q, r);
			return v >= Integer.MAX_VALUE / 100 ? p * q * r : v;
		};
	}

	private static NonCubicBilinearAlgorithm strassen() throws Exception {
		return SchemeIO.read(eu.solven.matmul.catalog.SchemeResolver.byHint("src/main/resources/schemes/known/section2/strassen-2x2x2_m7_a18.json"));
	}

	@Test
	public void matchesAllocationOptimizer_cubic() throws Exception {
		NonCubicBilinearAlgorithm S = strassen();
		SotaResolver sota = sota();

		for (int[] tgt : new int[][] { { 8, 8, 8 }, { 13, 13, 13 }, { 9, 9, 9 } }) {
			long flat = AllocationOptimizer.optimize(S, sota, tgt[0], tgt[1], tgt[2]).rank();
			AssignmentOptimizer.Result asg = AssignmentOptimizer.optimize(S, sota, tgt[0], tgt[1], tgt[2]);
			assertThat(asg.rank()).as("⟨%d,%d,%d⟩", tgt[0], tgt[1], tgt[2]).isEqualTo(flat);
		}
	}

	@Test
	public void knownOptima() throws Exception {
		NonCubicBilinearAlgorithm S = strassen();
		SotaResolver sota = sota();
		assertThat(AssignmentOptimizer.optimize(S, sota, 8, 8, 8).rank()).isEqualTo(336);
		assertThat(AssignmentOptimizer.optimize(S, sota, 13, 13, 13).rank()).isEqualTo(1432);
	}
}
