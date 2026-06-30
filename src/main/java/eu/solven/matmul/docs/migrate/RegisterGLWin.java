package eu.solven.matmul.docs.migrate;

import java.io.File;
import java.nio.file.Path;
import java.util.ArrayList;
import java.util.List;
import java.util.Optional;

import eu.solven.matmul.NonCubicBilinearAlgorithm;
import eu.solven.matmul.SymmetryTransforms;
import eu.solven.matmul.verifiers.Verifier;
import eu.solven.matmul.catalog.FieldAwareLookup;
import eu.solven.matmul.catalog.SchemeIO;
import eu.solven.matmul.recombination.AllocationOptimizer;
import eu.solven.matmul.recombination.AnalyticalMaskSearch.SchemeSupports;
import eu.solven.matmul.recombination.BlockSplitSearch;
import eu.solven.matmul.recombination.BlockSplitSearch.NamedBase;
import eu.solven.matmul.recombination.Recombination.SotaResolver;
import eu.solven.matmul.isotropy.RecombinationMultisetOrbit;
import eu.solven.matmul.isotropy.RecombinationMultisetOrbit.Result;
import eu.solven.matmul.search.RecombinationPoolConfig;
import eu.solven.matmul.search.RecursiveMaterialiser;
import eu.solven.matmul.search.SearchBudget;

/**
 * Register a GL-frontier recombination win to the catalog. Injects the single best GL-orbit member
 * of ⟨2,2,2⟩ for the target into the pool and materialises with {@code writeNewSchemes=true} into
 * {@code src/main/resources/schemes/derived/section{maxDim}/}.
 *
 * <p>Costing uses a DISK-READING sota ({@link FieldAwareLookup#findRank}) — NOT
 * {@code RecursiveClosureSota}, which re-runs the full {@code findBestStrategy} for every shape
 * recursively (ignoring the lineaged disk rank) and so takes minutes. With the disk sota the
 * decomposition is chosen in O(1) lookups and only the composition itself runs.</p>
 *
 * <p>Args: {@code n m p} (default 17 19 20). Re-verifies the written file before returning.</p>
 */
public final class RegisterGLWin {
	private RegisterGLWin() {}

	private static SotaResolver diskSota(FieldAwareLookup lk) {
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
		int n = args.length >= 3 ? Integer.parseInt(args[0]) : 17;
		int m = args.length >= 3 ? Integer.parseInt(args[1]) : 19;
		int p = args.length >= 3 ? Integer.parseInt(args[2]) : 20;

		FieldAwareLookup lk = new FieldAwareLookup("Q"); // rational pipeline → stamp over Q (avoids R-floor dropping Q)
		SotaResolver sota = diskSota(lk);
		int catalog = lk.findRank(n, m, p);

		// Find the single best GL ⟨2,2,2⟩ frontier member for this shape.
		NonCubicBilinearAlgorithm str = SchemeIO.read(
				new File("src/main/resources/schemes/known/section2/2x2x2-r7-strassen-db11bcc.json"));
		Result orbit = RecombinationMultisetOrbit.enumerate(str, 2);
		NonCubicBilinearAlgorithm bestMem = null; long bestR = Long.MAX_VALUE;
		for (String key : orbit.dominanceFrontier()) {
			int[][][] xyz = orbit.representativeTransforms.get(key);
			if (xyz == null) continue;
			NonCubicBilinearAlgorithm mem = RecombinationMultisetOrbit.materialise(str, xyz[0], xyz[1], xyz[2]);
			for (NonCubicBilinearAlgorithm o : SymmetryTransforms.s3Orbit(mem)) {
				long r = AllocationOptimizer.optimize(SchemeSupports.extract(o), sota, n, m, p, SearchBudget.EXACT, null).rank();
				if (r < bestR) { bestR = r; bestMem = o; }
			}
		}
		System.out.printf("⟨%d,%d,%d⟩ catalog=%d  best GL-member recomb=%d%n", n, m, p, catalog, bestR);
		if (bestR >= catalog) { System.out.println("no strict win — aborting registration"); return; }

		List<NamedBase> pool = new ArrayList<>(BlockSplitSearch.buildPool(RecombinationPoolConfig.simple()));
		pool.add(new NamedBase("GL222_win", bestMem));

		// writeNewSchemes=true, improveExisting=true → persist the improvement to the real catalog.
		Path schemesRoot = Path.of("src/main/resources/schemes");
		RecursiveMaterialiser mat = new RecursiveMaterialiser(lk, pool, sota, schemesRoot, true, false, true);

		long t0 = System.nanoTime();
		Optional<RecursiveMaterialiser.Result> r = mat.materialise(n, m, p);
		System.out.printf("materialise+write took %.1fs%n", (System.nanoTime() - t0) / 1e9);
		if (r.isEmpty()) { System.out.println("materialise EMPTY — nothing written"); return; }

		NonCubicBilinearAlgorithm alg = r.get().alg();
		boolean verified = Verifier.isExactNonCubic(alg);
		System.out.printf("built ⟨%d,%d,%d⟩ rank=%d verified=%s fromDisk=%s%n",
				alg.n, alg.m, alg.p, alg.r, verified, r.get().fromDisk());

		// Locate + reload the written file to confirm it persisted correctly.
		int maxDim = Math.max(n, Math.max(m, p));
		File dir = schemesRoot.resolve("derived").resolve("section" + maxDim).toFile();
		File[] hits = dir.listFiles((d, name) -> name.startsWith(n + "x" + m + "x" + p + "-") && name.contains("r" + alg.r));
		if (hits != null && hits.length > 0) {
			for (File f : hits) {
				NonCubicBilinearAlgorithm re = SchemeIO.read(f);
				System.out.printf("WROTE %s — reload rank=%d verified=%s%n",
						f.getName(), re.r, Verifier.isExactNonCubic(re));
			}
		} else {
			System.out.println("NO written file found under " + dir + " (rank " + alg.r + ") — check write gate");
		}
	}
}
