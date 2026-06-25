package eu.solven.matmul.docs.explore;

import java.io.File;

import eu.solven.matmul.NonCubicBilinearAlgorithm;
import eu.solven.matmul.algebra.Algebra;
import eu.solven.matmul.algebra.Field;
import eu.solven.matmul.catalog.SchemeIO;
import eu.solven.matmul.recombination.AllocationOptimizer;
import eu.solven.matmul.recombination.AnalyticalMaskSearch.SchemeSupports;
import eu.solven.matmul.recombination.Recombination;
import eu.solven.matmul.recombination.Recombination.SotaResolver;
import eu.solven.matmul.recombination.RecombinationMultisetOrbit;
import eu.solven.matmul.recombination.RecombinationMultisetOrbit.Result;
import eu.solven.matmul.search.SearchBudget;

/**
 * FAIR-COSTED frontier test (no descending-partition handicap): for one base, cost its NATIVE
 * support exactly as the pool does ({@link AllocationOptimizer} free composition), then cost EVERY
 * GL-frontier member the SAME way — {@code materialise(transform)} → real support → AllocationOptimizer.
 * Native ≡ pool by construction; any frontier member that beats native (same costing, same base,
 * same native orientation) is a REAL win, settling whether the orbit's alternative supports help.
 */
public final class ProbeFairFrontier {
	private ProbeFairFrontier() {}

	public static void main(String[] args) throws Exception {
		String path = args.length > 0 ? args[0]
				: "src/main/resources/schemes/known/section3/2x3x3-r15-alphatensor_Z-497eea7.json";
		int bound = args.length > 1 ? Integer.parseInt(args[1]) : 2;
		NonCubicBilinearAlgorithm base = SchemeIO.read(new File(path));
		SotaResolver sota = Recombination.catalogResolver(Algebra.nonCommutative(Field.R));

		int maxDim = Math.max(base.n, Math.max(base.m, base.p));
		System.out.printf("base %s ⟨%d,%d,%d⟩=%d — enumerating frontier (%s)…%n",
				new File(path).getName(), base.n, base.m, base.p, base.r, maxDim <= 3 ? "GL-exact" : "sampled");
		long t0 = System.nanoTime();
		Result orbit = maxDim <= 3 ? RecombinationMultisetOrbit.enumerate(base, bound)
				: RecombinationMultisetOrbit.enumerateSampled(base, 100_000, 2);
		var frontier = orbit.dominanceFrontier();
		System.out.printf("frontier=%d (%.0fs)%n%n", frontier.size(), (System.nanoTime() - t0) / 1e9);

		// Pre-materialise each frontier member and pre-build its ≤6 axis orientations (s3Orbit),
		// each as a real SchemeSupports — the fair, free-comp, orientation-covered frontier.
		java.util.List<SchemeSupports> frontierSupports = new java.util.ArrayList<>();
		for (String key : frontier) {
			int[][][] xyz = orbit.representativeTransforms.get(key);
			if (xyz == null) continue;
			NonCubicBilinearAlgorithm orbited = RecombinationMultisetOrbit.materialise(base, xyz[0], xyz[1], xyz[2]);
			for (NonCubicBilinearAlgorithm ob : eu.solven.matmul.SymmetryTransforms.s3Orbit(orbited))
				frontierSupports.add(SchemeSupports.extract(ob));
		}
		System.out.printf("frontier members × orientations → %d supports%n%n", frontierSupports.size());

		var pool = eu.solven.matmul.recombination.BlockSplitSearch.defaultPool();
		int[][] T = {{8,8,9},{7,8,9},{6,6,7},{3,7,8},{5,7,11},{3,3,13},{6,7,7},{3,5,12},{3,8,10},{3,4,14},{8,9,9},{6,8,11},{9,9,10},{7,7,11}};
		int beatsPool = 0, tiesPool = 0;
		for (int[] t : T) {
			long poolRank = eu.solven.matmul.recombination.BlockSplitSearch
					.findBestStrategy(t[0], t[1], t[2], pool, sota, false).map(s -> s.rank()).orElse(Long.MAX_VALUE);
			long fr = Long.MAX_VALUE;
			for (SchemeSupports sup : frontierSupports)
				fr = Math.min(fr, AllocationOptimizer.optimize(sup, sota, t[0], t[1], t[2], SearchBudget.EXACT, null).rank());
			String mark = fr < poolRank ? "<<< FRONTIER beats POOL by " + (poolRank - fr)
					: (fr == poolRank ? "(tie)" : "(pool better by " + (fr - poolRank) + ")");
			if (fr < poolRank) beatsPool++; else if (fr == poolRank) tiesPool++;
			System.out.printf("⟨%d,%d,%d⟩  pool=%d  frontier(fair,oriented)=%d  %s%n", t[0], t[1], t[2], poolRank, fr, mark);
		}
		System.out.printf("%nFAIR+ORIENTED frontier-of-⟨2,3,3⟩ vs FULL POOL:  beats=%d  ties=%d  of %d%n",
				beatsPool, tiesPool, T.length);
	}
}
