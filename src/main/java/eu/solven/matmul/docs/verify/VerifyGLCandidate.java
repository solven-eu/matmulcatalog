package eu.solven.matmul.docs.verify;

import java.io.File;
import java.nio.file.Path;
import java.util.ArrayList;
import java.util.List;
import java.util.Optional;

import eu.solven.matmul.NonCubicBilinearAlgorithm;
import eu.solven.matmul.SymmetryTransforms;
import eu.solven.matmul.Verifier;
import eu.solven.matmul.catalog.FieldAwareLookup;
import eu.solven.matmul.catalog.SchemeIO;
import eu.solven.matmul.recombination.BlockSplitSearch;
import eu.solven.matmul.recombination.BlockSplitSearch.NamedBase;
import eu.solven.matmul.recombination.BlockSplitSearch.NonCubicStrategy;
import eu.solven.matmul.recombination.Recombination;
import eu.solven.matmul.recombination.Recombination.SotaResolver;
import eu.solven.matmul.recombination.RecombinationMultisetOrbit;
import eu.solven.matmul.recombination.RecombinationMultisetOrbit.Result;
import eu.solven.matmul.search.RecombinationPoolConfig;
import eu.solven.matmul.search.RecursiveClosureSota;
import eu.solven.matmul.search.RecursiveMaterialiser;

/**
 * Build-verify the GLLargeSweep candidates by injecting the full GL-orbit frontier members of the
 * ⟨2,2,2⟩ base (which AXIS_FLIP cannot reach) into the materialiser pool, then materialising the
 * ACTUAL recombined scheme (no disk write) and running {@link Verifier#isExactNonCubic}. A candidate
 * is a real win iff Verifier PASSES and the built rank is strictly below the catalog rank.
 */
public final class VerifyGLCandidate {
	private VerifyGLCandidate() {}

	private static SotaResolver catalogSota(FieldAwareLookup lk) {
		return (p, q, r) -> {
			if (p == 0 || q == 0 || r == 0) return 0;
			if (p == 1) return q * r;
			if (q == 1) return p * r;
			if (r == 1) return p * q;
			int v = lk.findRank(p, q, r);
			return v >= SotaResolver.UNKNOWN_RANK ? p * q * r : v;
		};
	}

	public static void main(String[] args) throws Exception {
		FieldAwareLookup lk = new FieldAwareLookup("R");
		SotaResolver sota = catalogSota(lk);

		// All GL-orbit frontier members of ⟨2,2,2⟩ (× s3 orientations) as candidate bases.
		NonCubicBilinearAlgorithm str = SchemeIO.read(
				new File("src/main/resources/schemes/known/section2/2x2x2-r7-strassen-db11bcc.json"));
		Result orbit = RecombinationMultisetOrbit.enumerate(str, 2);
		List<NonCubicBilinearAlgorithm> glMembers = new ArrayList<>();
		for (String key : orbit.dominanceFrontier()) {
			int[][][] xyz = orbit.representativeTransforms.get(key);
			if (xyz == null) continue;
			NonCubicBilinearAlgorithm mem = RecombinationMultisetOrbit.materialise(str, xyz[0], xyz[1], xyz[2]);
			for (NonCubicBilinearAlgorithm o : SymmetryTransforms.s3Orbit(mem)) glMembers.add(o);
		}
		System.out.printf("%d GL ⟨2,2,2⟩ candidate members%n%n", glMembers.size());

		int[][] cands = args.length >= 3
				? new int[][] { { Integer.parseInt(args[0]), Integer.parseInt(args[1]), Integer.parseInt(args[2]) } }
				: new int[][] { { 17, 19, 20 }, { 18, 30, 31 } };
		System.out.printf("%-12s %8s %8s %8s %9s  %s%n", "shape", "catalog", "eval", "matRank", "verified", "verdict");
		System.out.println("-".repeat(78));
		for (int[] c : cands) {
			int n = c[0], m = c[1], p = c[2];
			int catalog = lk.findRank(n, m, p);

			// Pick ONLY the single best GL member for THIS shape (injecting all 36 rational bases
			// blows up the recursive search at every node). Score via the same B&B costing.
			NonCubicBilinearAlgorithm bestMem = null; long bestR = Long.MAX_VALUE;
			for (NonCubicBilinearAlgorithm mem : glMembers) {
				long r = eu.solven.matmul.recombination.AllocationOptimizer.optimize(
						eu.solven.matmul.recombination.AnalyticalMaskSearch.SchemeSupports.extract(mem),
						sota, n, m, p, eu.solven.matmul.search.SearchBudget.EXACT, null).rank();
				if (r < bestR) { bestR = r; bestMem = mem; }
			}
			List<NamedBase> pool = new ArrayList<>(BlockSplitSearch.buildPool(RecombinationPoolConfig.simple()));
			pool.add(new NamedBase("GL222_win", bestMem));
			RecursiveClosureSota recSota = new RecursiveClosureSota(lk, pool, false, true);
			RecursiveMaterialiser mat = new RecursiveMaterialiser(lk, pool, recSota, Path.of("target/verify-tmp"), false, false);

			Optional<NonCubicStrategy> strat = BlockSplitSearch.findBestStrategy(n, m, p, pool, sota, false);
			long eval = strat.map(NonCubicStrategy::rank).orElse(-1L);

			long matRank = -1; boolean verified = false; String note = "glBnB=" + bestR + " ";
			try {
				Optional<RecursiveMaterialiser.Result> r = mat.materialise(n, m, p);
				if (r.isPresent()) {
					NonCubicBilinearAlgorithm alg = r.get().alg();
					matRank = alg.r;
					verified = Verifier.isExactNonCubic(alg);
				} else note = "materialise-empty";
			} catch (Throwable t) {
				note = "EXC:" + t.getClass().getSimpleName() + ":" + t.getMessage();
			}

			boolean haveCat = catalog < SotaResolver.UNKNOWN_RANK;
			String verdict = !verified ? "NOT VERIFIED " + note
					: !haveCat ? "verified, no catalog"
					: matRank < catalog ? "✔ REAL WIN (−" + (catalog - matRank) + ")"
					: matRank == catalog ? "tie (eval not reproduced)"
					: "materialiser worse";
			System.out.printf("⟨%d,%d,%d⟩ %8s %8d %8d %9s  %s%n",
					n, m, p, haveCat ? String.valueOf(catalog) : "none", eval, matRank, verified, verdict);
		}
	}
}
