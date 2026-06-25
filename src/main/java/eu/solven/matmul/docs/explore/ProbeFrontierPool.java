package eu.solven.matmul.docs.explore;

import java.io.File;
import java.nio.file.Path;
import java.util.ArrayList;
import java.util.List;

import eu.solven.matmul.algebra.Algebra;
import eu.solven.matmul.algebra.Field;
import eu.solven.matmul.recombination.BlockSplitSearch;
import eu.solven.matmul.recombination.FrontierRecombination;
import eu.solven.matmul.recombination.Recombination;
import eu.solven.matmul.recombination.Recombination.SotaResolver;
import eu.solven.matmul.recombination.RecombFrontierIO;

/**
 * FAIR, INDEX-BASED comparison: load precomputed frontier sidecars (one per base, canonical), score
 * each over all 6 axis orientations via column permutation, and compare the pool minimum to the
 * current {@code findBestStrategy}. Coverage is ≥ by construction (frontier ⊇ native); the question
 * is whether the exhaustive orientation scan ever STRICTLY beats the live search.
 */
public final class ProbeFrontierPool {
	private ProbeFrontierPool() {}

	public static void main(String[] args) throws Exception {
		SotaResolver sota = Recombination.catalogResolver(Algebra.nonCommutative(Field.R));
		var pool = BlockSplitSearch.defaultPool();

		File dir = new File("src/main/resources/frontiers");
		File[] files = dir.listFiles((d, n) -> n.endsWith(".json"));
		List<RecombFrontierIO.Loaded> loaded = new ArrayList<>();
		for (File f : files) loaded.add(RecombFrontierIO.read(Path.of(f.getPath())));
		System.out.println("loaded " + loaded.size() + " frontier sidecars: "
				+ java.util.Arrays.stream(files).map(File::getName).sorted().toList());

		int[][] shapes = {{2,9,10},{3,10,10},{8,8,9},{4,5,11},{5,7,9},{6,6,7},{3,7,8},{4,4,9},{2,7,11},{5,5,6},{7,8,9},{4,9,9},{3,8,10},{4,7,7},{2,5,13},{3,4,14},{2,3,17},{5,8,8},{6,7,7},{3,3,13},{4,6,9},{2,6,13},{5,6,10},{7,7,8},{3,5,12},{4,8,8},{2,4,15},{6,8,9},{3,9,9},{5,7,11}};
		int wins = 0, ties = 0;
		long tSearch = 0, tFront = 0;
		for (int[] t : shapes) {
			long s0 = System.nanoTime();
			var best = BlockSplitSearch.findBestStrategy(t[0], t[1], t[2], pool, sota, false);
			long cur = best.map(s -> s.rank()).orElse(Long.MAX_VALUE);
			String poolLabel = best.map(s -> s.label()).orElse("none");
			tSearch += System.nanoTime() - s0;
			long f0 = System.nanoTime();
			long fr = Long.MAX_VALUE;
			for (RecombFrontierIO.Loaded c : loaded)
				fr = Math.min(fr, FrontierRecombination.bestRankOverOrientations(t[0], t[1], t[2], c.dims(), c.frontier(), sota));
			tFront += System.nanoTime() - f0;
			String mark = fr < cur ? "  <<< FRONTIER WINS by " + (cur - fr)
					: (fr == cur ? "  (tie)" : "  (search better by " + (fr - cur) + ")");
			if (fr < cur) wins++; else if (fr == cur) ties++;
			System.out.printf("⟨%d,%d,%d⟩  pool=%d  frontier=%d%s   poolBase=%s%n", t[0], t[1], t[2], cur, fr, mark, poolLabel);
		}
		System.out.printf("%nFRONTIER-WINS=%d  TIES=%d  of %d  (only %d sidecar bases loaded; not the full pool)%n",
				wins, ties, shapes.length, loaded.size());
		double perBaseFront = tFront / 1e6 / (shapes.length * (double) loaded.size());
		System.out.printf("TIMING: findBestStrategy total=%.0fms (%.1fms/shape) | frontier-scan total=%.0fms over %d bases"
				+ " (%.3fms per base·shape) | extrapolated to 78 bases ≈ %.1fms/shape%n",
				tSearch / 1e6, tSearch / 1e6 / shapes.length, tFront / 1e6, loaded.size(),
				perBaseFront, perBaseFront * 78);
	}
}
