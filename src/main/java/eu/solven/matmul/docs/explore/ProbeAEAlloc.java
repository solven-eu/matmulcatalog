package eu.solven.matmul.docs.explore;

import java.io.File;

import eu.solven.matmul.NonCubicBilinearAlgorithm;
import eu.solven.matmul.catalog.FieldAwareLookup;
import eu.solven.matmul.catalog.SchemeIO;
import eu.solven.matmul.recombination.AllocationOptimizer;
import eu.solven.matmul.recombination.Recombination;

/**
 * Throwaway: run the UNBOUNDED {@link AllocationOptimizer} for the AE⟨5,5,5⟩=93 base onto
 * ⟨14,15,31⟩ (master reaches 3839 via alloc [3,2,3,3,3]/[6,6,6,7,6]; the branch search
 * returned recomb=-1). If unbounded optimize reaches 3839 → the bounded search prunes it
 * (cheapUB/stagnation bug); if it reaches only ≥3840 → the optimizer itself misses the
 * asymmetric 5-part allocation (uniqueAxis/axisDims dedup).
 */
public class ProbeAEAlloc {

	public static void main(String[] args) throws java.io.IOException {
		FieldAwareLookup lookup = new FieldAwareLookup("Q");
		Recombination.SotaResolver sota = (a, b, c) -> lookup.findRank(a, b, c);
		NonCubicBilinearAlgorithm ae = SchemeIO.read(
				new File("src/main/resources/schemes/known/section5/alphaevolve-5x5x5_m93_a846.json"));
		System.out.printf("AE base = ⟨%d,%d,%d⟩=%d%n", ae.n, ae.m, ae.p, ae.r);
		for (int[] c : new int[][] { { 14, 15, 31 }, { 14, 15, 29 } }) {
			long t0 = System.currentTimeMillis();
			AllocationOptimizer.Result r = AllocationOptimizer.optimize(ae, sota, c[0], c[1], c[2]);
			System.out.printf("⟨%d,%d,%d⟩ optimize: rank=%d  allocA=%s allocB=%s allocC=%s  exhaustive=%s nodes=%d (%dms)%n",
					c[0], c[1], c[2], r.rank(), java.util.Arrays.toString(r.allocA()),
					java.util.Arrays.toString(r.allocB()), java.util.Arrays.toString(r.allocC()),
					r.exhaustive(), r.nodes(), System.currentTimeMillis() - t0);
		}
	}
}
