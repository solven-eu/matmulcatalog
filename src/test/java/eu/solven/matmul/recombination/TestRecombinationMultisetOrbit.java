package eu.solven.matmul.recombination;


import static org.assertj.core.api.Assertions.assertThat;

import java.io.File;
import java.util.ArrayList;
import java.util.List;

import org.junit.jupiter.api.Tag;
import org.junit.jupiter.api.Test;

import eu.solven.matmul.NonCubicBilinearAlgorithm;
import eu.solven.matmul.catalog.SchemeIO;

/**
 * Drives {@link RecombinationMultisetOrbit} — the base-agnostic recombination
 * multiset enumerator — on ⟨2,2,2⟩ (cross-checked to 40) and a non-cubic base.
 */
class TestRecombinationMultisetOrbit {

	private static final String KNOWN = "src/main/resources/schemes/known/";

	private static NonCubicBilinearAlgorithm load(String relPath) throws Exception {
		return SchemeIO.readBilinear(new File(KNOWN + relPath));
	}

	@Test
	@Tag("slow")
	void strassen_2x2x2_reproduces_40() throws Exception {
		NonCubicBilinearAlgorithm strassen = load("section2/2x2x2-r7-strassen-db11bcc.json");
		RecombinationMultisetOrbit.Result res = RecombinationMultisetOrbit.enumerate(strassen, 2);
		System.out.printf("%n⟨2,2,2⟩ Strassen: per-axis patterns=%s, %d combinations → %d canonical multisets%n",
				java.util.Arrays.toString(res.perAxisPatternCounts), res.combinations, res.canonicalMultisets.size());
		assertThat(res.canonicalMultisets).hasSize(40);
		assertThat(RecombinationMultisetOrbit.isStable(strassen, 2)).as("stable at bound 2").isTrue();

		// Symbolic rendering (n₁≥n₂ etc.) — the result is generic in the block sizes.
		List<String> keys = new ArrayList<>(res.canonicalMultisets);
		java.util.Collections.sort(keys);
		System.out.println("All 40 over ℚ (symbolic in block sizes n₁≥n₂):");
		int i = 0;
		for (String k : keys) {
			i++;
			System.out.printf("  #%-2d %s%n", i, RecombinationMultisetOrbit.prettySymbolic(k, "n", "n", "n"));
		}
	}

	@Test
	@Tag("slow")
	void noncubic_base_2x3x3() throws Exception {
		// A non-cubic base: only the m↔p axis-swap is a symmetry, so the
		// canonicalisation group is Z₂ (not S₃).
		NonCubicBilinearAlgorithm base = load("section3/2x3x3-r15-alphatensor_Z-497eea7.json");
		System.out.printf("%n⟨%d,%d,%d⟩ base, rank %d%n", base.n, base.m, base.p, base.r);
		RecombinationMultisetOrbit.Result res = RecombinationMultisetOrbit.enumerate(base, 2);
		System.out.printf("per-axis patterns=%s, %,d combinations → %d canonical multisets, stabilizer=%d perms%n",
				java.util.Arrays.toString(res.perAxisPatternCounts), res.combinations,
				res.canonicalMultisets.size(),
				RecombinationMultisetOrbit.shapeStabilizer(base.n, base.m, base.p).length);
		assertThat(res.canonicalMultisets).isNotEmpty();
	}
}
