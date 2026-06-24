package eu.solven.matmul.docs.explore;

import java.io.File;
import java.util.TreeSet;

import eu.solven.matmul.NonCubicBilinearAlgorithm;
import eu.solven.matmul.catalog.SchemeIO;
import eu.solven.matmul.recombination.RecombinationMultisetOrbit;
import eu.solven.matmul.recombination.RecombinationMultisetOrbit.Result;

/**
 * Validate the GL-free structural enumerators against the GL-odometer oracle where it is
 * tractable (dim ≤ 3), then run the d=4-capable frontier path where it is not (⟨2,4,4⟩).
 *
 * <p>Two structural entry points are checked:
 * <ul>
 * <li>{@code enumerateStructural} — full canonical set; must EQUAL {@code enumerate}.</li>
 * <li>{@code enumerateStructuralFrontier} — per-axis-pruned; its frontier must EQUAL
 *     {@code enumerate(...).dominanceFrontier()}.</li>
 * </ul>
 * Generics are now Vandermonde (node &gt; max|coef|), so {@code genericBound=1} suffices.
 */
public final class ProbeStructuralFrontier {
	private ProbeStructuralFrontier() {}

	public static void main(String[] args) throws Exception {
		String[][] cases = {
				{ "src/main/resources/schemes/known/section2/2x2x2-r7-strassen-db11bcc.json", "2" },
				{ "src/main/resources/schemes/known/section2/2x2x2-r7-winograd_1971-511df05.json", "2" },
				{ "src/main/resources/schemes/known/section3/2x2x3-r11-alphatensor_Z-682e003.json", "2" },
				{ "src/main/resources/schemes/known/section3/2x2x3-r11-kauers_2026-72562fb.json", "2" },
		};
		for (String[] c : cases) {
			File f = new File(c[0]);
			if (!f.exists()) { System.out.println("MISSING " + c[0]); continue; }
			NonCubicBilinearAlgorithm base = SchemeIO.read(f);
			int glBound = Integer.parseInt(c[1]);

			Result gl = RecombinationMultisetOrbit.enumerate(base, glBound);
			Result stFull = RecombinationMultisetOrbit.enumerateStructural(base, 1);
			Result stFr = RecombinationMultisetOrbit.enumerateStructuralFrontier(base, 1);

			boolean fullSame = gl.canonicalMultisets.equals(stFull.canonicalMultisets);
			var glFr = new TreeSet<>(gl.dominanceFrontier());
			var stFrFr = new TreeSet<>(stFr.dominanceFrontier());
			boolean frSame = glFr.equals(stFrFr);

			System.out.printf("%n%s  ⟨%d,%d,%d⟩ r=%d%n", f.getName(), base.n, base.m, base.p, base.r);
			System.out.printf("  full      : GL canonical=%d  STRUCT canonical=%d  -> %s%n",
					gl.canonicalMultisets.size(), stFull.canonicalMultisets.size(), fullSame ? "MATCH ✓" : "DIFF ✗");
			System.out.printf("  frontier  : GL frontier=%d  STRUCT-FR frontier=%d (from %d cand)  -> %s%n",
					glFr.size(), stFrFr.size(), stFr.canonicalMultisets.size(), frSame ? "MATCH ✓" : "DIFF ✗");
			if (!frSame) {
				var onlyGl = new TreeSet<>(glFr); onlyGl.removeAll(stFrFr);
				var onlySt = new TreeSet<>(stFrFr); onlySt.removeAll(glFr);
				System.out.println("    only-in-GL=" + onlyGl.size() + " only-in-STRUCT=" + onlySt.size());
				for (String k : onlyGl) System.out.println("      GL-only: " + RecombinationMultisetOrbit.prettySymbolic(k, "n", "m", "p"));
			}
			System.out.flush();
		}

		// ⟨2,4,4⟩ (dim-4): GL = 9^16 intractable; the frontier path prunes per-axis, so it runs.
		File f244 = new File("src/main/resources/schemes/known/section4/2x4x4-r26-alphatensor_Z-8e5986a.json");
		if (f244.exists()) {
			NonCubicBilinearAlgorithm base = SchemeIO.read(f244);
			System.out.printf("%n%s  ⟨2,4,4⟩ r=%d  [FRONTIER-ONLY, GL=9^16 intractable]%n", f244.getName(), base.r);
			for (int gen : new int[] { 1, 2 }) {
				long t0 = System.nanoTime();
				Result st = RecombinationMultisetOrbit.enumerateStructuralFrontier(base, gen);
				System.out.printf("  gen=%d : per-axis antichains=%s  combined-cand=%d  frontier=%d  (%.1fs)%n",
						gen, java.util.Arrays.toString(st.perAxisPatternCounts), st.canonicalMultisets.size(),
						st.dominanceFrontier().size(), (System.nanoTime() - t0) / 1e9);
				System.out.flush();
			}
			Result st = RecombinationMultisetOrbit.enumerateStructuralFrontier(base, 1);
			System.out.println("  frontier multisets (gen=1):");
			int i = 1;
			for (String k : st.dominanceFrontier())
				System.out.printf("    [%2d] %s%n", i++, RecombinationMultisetOrbit.prettySymbolic(k, "n", "m", "p"));
		}
	}
}
