package eu.solven.matmul.docs.explore;

import java.io.File;
import java.util.List;

import eu.solven.matmul.NonCubicBilinearAlgorithm;
import eu.solven.matmul.catalog.SchemeIO;
import eu.solven.matmul.recombination.RecombinationMultisetOrbit;

/**
 * Per-base recombination-multiset frontier for ⟨2,2,3⟩: for each distinct rank-11
 * scheme, enumerate its GL-orbit, print canonical/frontier counts, and list its
 * frontier multisets symbolically (n₁=largest n-block, …). Different bases sweep
 * different GL-orbits → different (often disjoint) frontiers; the exhaustive menu is
 * the union. Prints per-base so a slow representative doesn't block earlier output.
 */
public final class ProbeMultisetReach {
	private ProbeMultisetReach() {}

	public static void main(String[] args) throws Exception {
		int bound = args.length > 0 ? Integer.parseInt(args[0]) : 2;
		String[] files = {
				"src/main/resources/schemes/known/section3/2x2x3-r11-alphatensor_Z-682e003.json",
				"src/main/resources/schemes/known/section3/2x2x3-r11-kauers_2026-72562fb.json",
				"src/main/resources/schemes/derived/section3/2x2x3-r11-derived-eaa4b53.json",
				"src/main/resources/schemes/known/section3/2x2x3-r11-alphatensor_F2-279f789.json",
				"src/main/resources/schemes/known/section3/2x2x3-r11-perminov_cr20_cn20_ZT_reduced-b7f24e6.json",
		};
		for (String fp : files) {
			File f = new File(fp);
			if (!f.exists()) { System.out.println("MISSING " + fp + "\n"); continue; }
			NonCubicBilinearAlgorithm base;
			try { base = SchemeIO.read(f); } catch (Exception e) { System.out.println("UNREADABLE " + fp + "\n"); continue; }
			RecombinationMultisetOrbit.Result r = RecombinationMultisetOrbit.enumerate(base, bound);
			List<String> fr = r.dominanceFrontier();
			System.out.printf("%n########## %s  (canonical=%d, frontier=%d) ##########%n",
					f.getName(), r.canonicalMultisets.size(), fr.size());
			int i = 1;
			for (String key : fr) {
				System.out.printf("  [%2d] %s%n", i++,
						RecombinationMultisetOrbit.prettySymbolic(key, "n", "m", "p"));
			}
			System.out.flush();
		}
	}
}
