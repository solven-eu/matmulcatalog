package eu.solven.matmul.search;

import java.io.File;
import java.util.ArrayList;
import java.util.LinkedHashMap;
import java.util.List;
import java.util.Map;

import org.junit.jupiter.api.Test;

import eu.solven.matmul.NonCubicBilinearAlgorithm;
import eu.solven.matmul.SymmetryTransforms;
import eu.solven.matmul.catalog.SchemeIO;

/**
 * Audit the catalog's rank-7 ⟨2,2,2⟩ schemes by their axis-flip canonical
 * signature. Schemes with the same canonical signature are in the same
 * discrete orbit and are cost-redundant (one of them is enough for
 * orbit-expanding search).
 *
 * <p>Diagnostic-only — dumps to stdout. Doesn't assert specific counts; the
 * catalog may grow.
 */
class TestCatalog2x2x2DiscreteOrbits {

	private static final String SECTION2 = "src/main/resources/schemes/known/section2/";

	@Test
	void audit_section2_rank7_discrete_orbits() throws Exception {
		File dir = new File(SECTION2);
		File[] files = dir.listFiles((d, n) -> n.endsWith(".json") && n.contains("_2x2x2_m7_"));
		if (files == null || files.length == 0) {
			System.out.println("No rank-7 ⟨2,2,2⟩ schemes found");
			return;
		}
		java.util.Arrays.sort(files, java.util.Comparator.comparing(File::getName));

		// Group by axis-flip canonical signature.
		Map<String, List<String>> orbits = new LinkedHashMap<>();
		List<String> unloadable = new ArrayList<>();
		for (File f : files) {
			try {
				NonCubicBilinearAlgorithm alg = SchemeIO.readBilinear(f);
				if (alg.n != 2 || alg.m != 2 || alg.p != 2 || alg.r != 7) continue;
				String canonical = SymmetryTransforms.axisFlipCanonicalSignature(alg);
				orbits.computeIfAbsent(canonical, k -> new ArrayList<>()).add(f.getName());
			} catch (Exception e) {
				unloadable.add(f.getName() + " (" + e.getClass().getSimpleName() + ")");
			}
		}

		System.out.println();
		System.out.println("==============================================================");
		System.out.println(" Section2 rank-7 ⟨2,2,2⟩ discrete-orbit audit");
		System.out.println("==============================================================");
		System.out.println();
		System.out.printf("Total loadable: %d schemes in %d distinct discrete orbits%n",
				orbits.values().stream().mapToInt(List::size).sum(), orbits.size());
		if (!unloadable.isEmpty()) {
			System.out.printf("Unloadable (different format / non-bilinear): %d%n", unloadable.size());
			for (String n : unloadable) System.out.println("    - " + n);
		}
		System.out.println();

		int orbitIdx = 0;
		for (var entry : orbits.entrySet()) {
			orbitIdx++;
			List<String> members = entry.getValue();
			System.out.printf("Orbit #%d (%d member%s):%n",
					orbitIdx, members.size(), members.size() == 1 ? "" : "s");
			for (String name : members) {
				System.out.println("    " + name);
			}
			if (members.size() > 1) {
				System.out.println("    → all members are axis-flip-equivalent (cost-redundant); "
						+ "keep one canonical, treat others as orbit-derived.");
			}
			System.out.println();
		}

		System.out.println("==============================================================");
		System.out.println(" Interpretation: each distinct orbit = potential outer-base for");
		System.out.println(" recombination. ⟨17,17,17⟩=2930 at (9,8)³ requires being in");
		System.out.println(" Winograd's orbit; Strassen's orbit only reaches 2940 there.");
		System.out.println(" The more distinct orbits, the more unbalanced-allocation paths");
		System.out.println(" the search can explore.");
		System.out.println("==============================================================");
	}
}
