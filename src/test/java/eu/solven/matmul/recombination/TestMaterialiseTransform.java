package eu.solven.matmul.recombination;

import static org.assertj.core.api.Assertions.assertThat;

import java.io.File;
import java.util.List;

import org.junit.jupiter.api.Test;

import eu.solven.matmul.NonCubicBilinearAlgorithm;
import eu.solven.matmul.verifiers.Verifier;
import eu.solven.matmul.catalog.SchemeIO;
import eu.solven.matmul.recombination.AnalyticalMaskSearch.SchemeSupports;
import eu.solven.matmul.isotropy.RecombinationMultisetOrbit.Result;
import eu.solven.matmul.isotropy.RecombinationMultisetOrbit;

/**
 * Step (a) guard: the GL transform captured per frontier multiset, applied via
 * {@link RecombinationMultisetOrbit#materialise}, must (1) still compute the same matmul and
 * (2) actually realise the target multiset (its native resolved support == that frontier key).
 * This certifies the "1 scheme + transform config" foundation for the frontier index.
 */
public class TestMaterialiseTransform {

	@Test
	public void transforms_realise_frontier_and_compute_matmul_2x2x2() {
		check("src/main/resources/schemes/known/section2/2x2x2-r7-strassen-db11bcc.json", 2);
	}

	@Test
	public void transforms_realise_frontier_and_compute_matmul_2x2x3() {
		check("src/main/resources/schemes/known/section3/2x2x3-r11-alphatensor_Z-682e003.json", 2);
	}

	private void check(String path, int bound) {
		File f = new File(path);
		assertThat(f).exists();
		NonCubicBilinearAlgorithm base;
		try {
			base = SchemeIO.read(f);
		} catch (Exception e) {
			throw new AssertionError(e);
		}
		Result orbit = RecombinationMultisetOrbit.enumerate(base, bound);
		List<String> frontier = orbit.dominanceFrontier();
		assertThat(frontier).isNotEmpty();

		int checked = 0;
		for (String key : frontier) {
			int[][][] xyz = orbit.representativeTransforms.get(key);
			assertThat(xyz).as("transform present for frontier key " + key).isNotNull();
			NonCubicBilinearAlgorithm orbited =
					RecombinationMultisetOrbit.materialise(base, xyz[0], xyz[1], xyz[2]);

			// (1) still computes the SAME matmul, at the SAME rank
			assertThat(orbited.r).isEqualTo(base.r);
			assertThat(Verifier.isExactNonCubic(orbited))
					.as("orbited scheme for " + key + " computes ⟨" + base.n + "," + base.m + "," + base.p + "⟩")
					.isTrue();

			// (2) its native resolved support == the target frontier multiset (canonical)
			String realised = resolvedCanonical(orbited);
			assertThat(realised).as("orbited scheme realises its target multiset").isEqualTo(key);
			checked++;
		}
		assertThat(checked).isEqualTo(frontier.size());
	}

	/** Canonical multiset key of a scheme's NATIVE support, by the same min/max index rule as the orbit. */
	private static String resolvedCanonical(NonCubicBilinearAlgorithm alg) {
		SchemeSupports s = SchemeSupports.extract(alg);
		int[][] shapes = new int[alg.r][3];
		for (int k = 0; k < alg.r; k++) {
			shapes[k][0] = Math.max(min(s.uRowSupport[k]), min(s.wRowSupport[k]));
			shapes[k][1] = Math.max(min(s.uColSupport[k]), min(s.vRowSupport[k]));
			shapes[k][2] = Math.max(min(s.vColSupport[k]), min(s.wColSupport[k]));
		}
		return RecombinationMultisetOrbit.canonicalKey(shapes,
				RecombinationMultisetOrbit.shapeStabilizer(alg.n, alg.m, alg.p));
	}

	private static int min(int[] a) {
		int mn = Integer.MAX_VALUE;
		for (int v : a) mn = Math.min(mn, v);
		return mn == Integer.MAX_VALUE ? 0 : mn;
	}
}
