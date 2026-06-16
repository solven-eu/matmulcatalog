package eu.solven.matmul;

import static org.assertj.core.api.Assertions.assertThat;

import java.io.File;

import org.junit.jupiter.api.Test;

import eu.solven.matmul.catalog.SchemeIO;

/**
 * Verifies axis-flip canonicalisation:
 * <ul>
 *   <li>An axis-flip variant of a scheme has the same canonical form as the
 *       canonical scheme — so they're recognised as same-orbit.</li>
 *   <li>Two schemes from different GL orbits (e.g. Strassen vs Winograd) have
 *       different canonical forms.</li>
 *   <li>The canonical form is idempotent: canonical(canonical(x)) ==
 *       canonical(x).</li>
 * </ul>
 */
class TestAxisFlipCanonical {

	private static NonCubicBilinearAlgorithm load(String name) throws Exception {
		return SchemeIO.read(eu.solven.matmul.catalog.SchemeResolver.byHint("src/main/resources/schemes/known/section2/" + name));
	}

	@Test
	void winograd_canonical_and_axflip1_cousin_are_orbit_equivalent() throws Exception {
		NonCubicBilinearAlgorithm winograd = load("winograd_1971-2x2x2_m7_a24.json");
		NonCubicBilinearAlgorithm cousinAxflip1 = load("solven_winograd_cousin_axflip1-2x2x2_m7_a24.json");
		assertThat(SymmetryTransforms.axisFlipEquivalent(winograd, cousinAxflip1))
				.as("winograd-1971 and solven-winograd-cousin-axflip1 are in the same axis-flip orbit")
				.isTrue();
		// Pool-dedup implication: keeping both is cost-redundant once
		// AnalyticalMaskSearch enumerates masks at search time.
	}

	@Test
	void strassen_and_winograd_have_distinct_canonical_forms() throws Exception {
		NonCubicBilinearAlgorithm strassen = load("strassen-2x2x2_m7_a18.json");
		NonCubicBilinearAlgorithm winograd = load("winograd_1971-2x2x2_m7_a24.json");
		assertThat(SymmetryTransforms.axisFlipEquivalent(strassen, winograd))
				.as("Strassen and Winograd are NOT axis-flip-equivalent (different discrete orbits)")
				.isFalse();
	}

	@Test
	void canonical_form_is_idempotent() throws Exception {
		NonCubicBilinearAlgorithm strassen = load("strassen-2x2x2_m7_a18.json");
		NonCubicBilinearAlgorithm c1 = SymmetryTransforms.axisFlipCanonical(strassen);
		NonCubicBilinearAlgorithm c2 = SymmetryTransforms.axisFlipCanonical(c1);
		assertThat(SymmetryTransforms.axisFlipCanonicalSignature(c1))
				.isEqualTo(SymmetryTransforms.axisFlipCanonicalSignature(c2));
	}

	@Test
	void canonical_form_invariant_under_full_orbit() throws Exception {
		NonCubicBilinearAlgorithm winograd = load("winograd_1971-2x2x2_m7_a24.json");
		String canonical = SymmetryTransforms.axisFlipCanonicalSignature(winograd);
		for (NonCubicBilinearAlgorithm variant : SymmetryTransforms.axisFlipOrbit(winograd)) {
			assertThat(SymmetryTransforms.axisFlipCanonicalSignature(variant))
					.as("orbit member must canonicalise to the same form")
					.isEqualTo(canonical);
		}
	}
}
