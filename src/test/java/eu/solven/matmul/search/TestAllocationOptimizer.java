package eu.solven.matmul.search;

import eu.solven.matmul.recombination.AllocationOptimizer;
import eu.solven.matmul.recombination.AnalyticalMaskSearch;

import static org.assertj.core.api.Assertions.assertThat;

import java.io.File;

import org.junit.jupiter.api.Test;

import eu.solven.matmul.NonCubicBilinearAlgorithm;
import eu.solven.matmul.catalog.FieldAwareLookup;
import eu.solven.matmul.recombination.Recombination.SotaResolver;
import eu.solven.matmul.catalog.SchemeIO;
import eu.solven.matmul.recombination.AnalyticalMaskSearch.SchemeSupports;

/**
 * Guards the single-base, non-recursive allocation optimiser (the 2026-06-04
 * balance discussion): branch-and-bound must (1) return the exact global
 * minimum over allocations, including the unbalanced cubic win
 * \langle13,13,13\rangle; and (2) honour the Kronecker upper-bound drop.
 */
public class TestAllocationOptimizer {

	private static NonCubicBilinearAlgorithm strassen() throws Exception {
		return SchemeIO.read(eu.solven.matmul.catalog.SchemeResolver.byHint("src/main/resources/schemes/known/section2/strassen-2x2x2_m7_a18.json"));
	}

	/** Catalog-backed SOTA over R, naive fallback for un-catalogued sub-shapes. */
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

	private static long brute(SchemeSupports sup, SotaResolver sota, int T) {
		long best = Long.MAX_VALUE;
		for (int a = 1; a < T; a++) for (int b = 1; b < T; b++) for (int c = 1; c < T; c++) {
			int[][] sh = AnalyticalMaskSearch.shapesAt(sup, new int[] { a, T - a }, new int[] { b, T - b }, new int[] { c, T - c });
			long tot = 0;
			for (int[] s : sh) tot += sota.getRank(s[0], s[1], s[2]);
			best = Math.min(best, tot);
		}
		return best;
	}

	@Test
	public void bnb_matches_exhaustive_including_unbalanced_cubic_win() throws Exception {
		NonCubicBilinearAlgorithm S = strassen();
		SotaResolver sota = sota();
		SchemeSupports sup = SchemeSupports.extract(S);

		// ⟨8,8,8⟩: balanced optimum (336 = Strassen²⊗DPS).
		AllocationOptimizer.Result r8 = AllocationOptimizer.optimize(S, sota, 8, 8, 8);
		assertThat(r8.rank()).isEqualTo(brute(sup, sota, 8)).isEqualTo(336);

		// ⟨13,13,13⟩: the optimum is UNBALANCED (1432 < balanced 1434).
		AllocationOptimizer.Result r13 = AllocationOptimizer.optimize(S, sota, 13, 13, 13);
		assertThat(r13.rank()).isEqualTo(brute(sup, sota, 13)).isEqualTo(1432);
		boolean balanced = java.util.Arrays.equals(r13.allocA(), new int[] { 7, 6 })
				&& java.util.Arrays.equals(r13.allocB(), new int[] { 7, 6 })
				&& java.util.Arrays.equals(r13.allocC(), new int[] { 7, 6 });
		assertThat(balanced).as("⟨13,13,13⟩ optimum is unbalanced").isFalse();
	}

	@Test
	public void kronecker_bound_drops_when_unbeatable() throws Exception {
		NonCubicBilinearAlgorithm S = strassen();
		SotaResolver sota = sota();

		// A loose bound (above the optimum) → full search, same answer.
		AllocationOptimizer.Result loose = AllocationOptimizer.optimize(S, sota, 8, 8, 8, 10_000);
		assertThat(loose.rank()).isEqualTo(336);
		assertThat(loose.improvedOnBound()).isTrue();

		// A bound at/below the root LB (R(4,4,4)=49 + 6·1 = 55) → early drop,
		// no allocations visited, base declared unable to improve.
		AllocationOptimizer.Result dropped = AllocationOptimizer.optimize(S, sota, 8, 8, 8, 50);
		assertThat(dropped.nodes()).as("root LB ≥ bound → sweep skipped").isZero();
		assertThat(dropped.improvedOnBound()).isFalse();
	}
}
