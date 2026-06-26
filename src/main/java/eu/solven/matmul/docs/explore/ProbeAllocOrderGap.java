package eu.solven.matmul.docs.explore;

import java.io.File;

import eu.solven.matmul.NonCubicBilinearAlgorithm;
import eu.solven.matmul.catalog.FieldAwareLookup;
import eu.solven.matmul.catalog.SchemeIO;
import eu.solven.matmul.recombination.AllocationOptimizer;
import eu.solven.matmul.recombination.AnalyticalMaskSearch;
import eu.solven.matmul.recombination.AnalyticalMaskSearch.SchemeSupports;
import eu.solven.matmul.recombination.Recombination;

/**
 * Throwaway probe for the ⟨21,21,21⟩/⟨23,23,23⟩ allocation-ORDER gap: master holds 5202
 * via Strassen ⟨2,2,2⟩ at allocC=[12,9], the optimizer returns 5240 at [9,12]. Is
 * [12,9]³ actually cheaper (→ the uniqueAxis dim-vector dedup wrongly drops it), or do
 * both orders cost the same (→ master's 5202 came from a different base variant / path)?
 */
public class ProbeAllocOrderGap {

	public static void main(String[] args) throws java.io.IOException {
		FieldAwareLookup lookup = new FieldAwareLookup("Q");
		Recombination.SotaResolver sota = (a, b, c) -> lookup.findRank(a, b, c);

		NonCubicBilinearAlgorithm strassen = SchemeIO.read(
				new File("src/main/resources/schemes/known/section2/strassen-2x2x2_m7_a18.json"));
		SchemeSupports sup = SchemeSupports.extract(strassen);

		for (int[] c : new int[][] { { 21, 21, 21 }, { 23, 23, 23 } }) {
			int N = c[0], M = c[1], P = c[2];
			System.out.printf("%n=== ⟨%d,%d,%d⟩ Strassen<2,2,2> ===%n", N, M, P);
			int hi = (N + 1) / 2, lo = N / 2;   // e.g. 21 -> 11,10 ; we want the master split
			// master split for 21 is [12,9] (imbalance 3); enumerate a few orders explicitly.
			for (int big : new int[] { N - N / 2 + 1, N - N / 3, (N * 4 + 3) / 7 }) {
				int small = N - big;
				if (small <= 0 || big <= 0) continue;
				long bigFirst = cost(sup, sota, new int[] { big, small });
				long smallFirst = cost(sup, sota, new int[] { small, big });
				System.out.printf("  split [%d,%d]: bigFirst=%d  smallFirst=%d  %s%n",
						big, small, bigFirst, smallFirst,
						bigFirst == smallFirst ? "(order-invariant)" : "ORDER MATTERS");
			}
			AllocationOptimizer.Result r = AllocationOptimizer.optimize(strassen, sota, N, M, P);
			System.out.printf("  optimizer best: rank=%d allocC=%s (exhaustive=%s, nodes=%d)%n",
					r.rank(), java.util.Arrays.toString(r.allocC()), r.exhaustive(), r.nodes());
		}
	}

	/** Exact recombination cost with the SAME split on all three axes. */
	private static long cost(SchemeSupports sup, Recombination.SotaResolver sota, int[] split) {
		int[][] shapes = AnalyticalMaskSearch.shapesAt(sup, split, split, split);
		long tot = 0;
		for (int[] s : shapes) tot += sota.getRank(s[0], s[1], s[2]);
		return tot;
	}
}
