package eu.solven.matmul.docs.generate;

import java.io.File;
import java.nio.file.Path;

import eu.solven.matmul.NonCubicBilinearAlgorithm;
import eu.solven.matmul.catalog.SchemeIO;
import eu.solven.matmul.recombination.RecombFrontierIO;
import eu.solven.matmul.isotropy.RecombinationMultisetOrbit;
import eu.solven.matmul.isotropy.RecombinationMultisetOrbit.Result;

/**
 * Emit the recombination-frontier sidecar JSON ({@link RecombFrontierIO}) for one or more base
 * schemes. For a base with all axes ≤3 it uses the exact GL enumeration and certifies
 * exhaustiveness via direction-bound stability; for a dim-≥4 axis (GL intractable) it falls back to
 * the partial structural frontier and stamps {@code exhaustive=false}.
 *
 * <p>Run: {@code mvn -q -ntp exec:java
 * -Dexec.mainClass=eu.solven.matmul.docs.generate.GenerateSchemeFrontiers
 * -Dexec.args="path/to/base1.json path/to/base2.json …"}</p>
 */
public final class GenerateSchemeFrontiers {
	private GenerateSchemeFrontiers() {}

	public static void main(String[] args) throws Exception {
		if (args.length == 0) {
			System.out.println("usage: GenerateSchemeFrontiers <base.json> [<base.json> …]");
			return;
		}
		int bound = 2;
		for (String path : args) {
			File f = new File(path);
			if (!f.exists()) { System.out.println("MISSING " + path); continue; }
			NonCubicBilinearAlgorithm base;
			try {
				base = SchemeIO.read(f);
			} catch (Exception e) {
				System.out.println("UNREADABLE " + path + " : " + e.getMessage());
				continue;
			}
			int maxDim = Math.max(base.n, Math.max(base.m, base.p));
			Result orbit;
			boolean exhaustive;
			String method;
			long t0 = System.nanoTime();
			if (maxDim <= 3) {
				orbit = RecombinationMultisetOrbit.enumerate(base, bound);
				// Exact requires BOTH a non-inflated antichain (orbit.frontierExact — not sampled) AND a
				// dirBound-saturated canonical set (isStable). dim≤2 trivially stable; dim3 certified.
				exhaustive = orbit.frontierExact && RecombinationMultisetOrbit.isStable(base, bound);
				method = "GL-exact";
			} else {
				// d≥4: exact GL odometer is intractable → saturation-sampled anytime frontier
				// (random GL samples until 100k consecutive yield nothing new). Partial menu — the
				// sampled dominance frontier is INFLATED (orbit.frontierExact == false), see Result doc.
				orbit = RecombinationMultisetOrbit.enumerateSampled(base, 100_000, 2);
				exhaustive = false;
				method = "sampled";
			}
			Path out = RecombFrontierIO.write(base, orbit, exhaustive, method, bound);
			System.out.printf("%s ⟨%d,%d,%d⟩ r=%d : frontier=%d %s exhaustive=%s -> %s (%.1fs)%n",
					f.getName(), base.n, base.m, base.p, base.r, orbit.dominanceFrontier().size(),
					method, exhaustive, out, (System.nanoTime() - t0) / 1e9);
		}
	}
}
