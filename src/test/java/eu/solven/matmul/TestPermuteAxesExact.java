package eu.solven.matmul;

import static org.assertj.core.api.Assertions.assertThat;

import java.util.List;

import org.junit.jupiter.api.Test;

import eu.solven.matmul.catalog.Compose;
import eu.solven.matmul.catalog.SchemeIO;
import eu.solven.matmul.catalog.SchemeResolver;

/**
 * Regression guard for the recombination base-orientation divergence
 * ({@code project_recomb_base_orientation_not_pinned}).
 *
 * <p>A recombination base with two equal-sized axes (⟨4,4,3⟩ from ⟨3,4,4⟩) has
 * MORE THAN ONE axis-permutation landing on the same shape, with DIFFERENT
 * U/V/W. The search scored one orientation but the lineage recorded only the
 * shape, so {@code orientAs}-on-replay picked a different (worse) one → the
 * built scheme's rank exceeded the evaluated rank and the predict/build guard
 * rejected an otherwise-valid win.
 *
 * <p>The fix records the EXACT axis-relabel ({@link SymmetryTransforms#s3OrbitWithPerms})
 * and replays it deterministically ({@link SymmetryTransforms#permuteAxes}, also
 * used by {@code LineageReplayer.applyTranspose}). This test pins the invariant
 * that makes the lineage faithful: {@code permuteAxes(canonical, v.perm())} must
 * reproduce {@code v.alg()} bit-for-bit, including for equal-axis bases.
 */
public class TestPermuteAxesExact {

	private static NonCubicBilinearAlgorithm loadStrassen() throws Exception {
		return SchemeIO.read(SchemeResolver.byHint(
				"src/main/resources/schemes/known/section2/strassen-2x2x2_m7_a18.json"));
	}

	private static void assertSame(NonCubicBilinearAlgorithm actual, NonCubicBilinearAlgorithm expected) {
		assertThat(new int[] { actual.n, actual.m, actual.p })
				.as("shape").isEqualTo(new int[] { expected.n, expected.m, expected.p });
		assertThat(actual.denseU()).as("U").isDeepEqualTo(expected.denseU());
		assertThat(actual.denseV()).as("V").isDeepEqualTo(expected.denseV());
		assertThat(actual.denseW()).as("W").isDeepEqualTo(expected.denseW());
	}

	@Test
	public void permuteAxes_realizes_each_of_six_canonical_perms() throws Exception {
		NonCubicBilinearAlgorithm s = loadStrassen();
		assertSame(SymmetryTransforms.permuteAxes(s, "ABC->ABC"), s);
		assertSame(SymmetryTransforms.permuteAxes(s, "ABC->CBA"), s.transpose());
		assertSame(SymmetryTransforms.permuteAxes(s, "ABC->BCA"), s.cyclicShift());
		assertSame(SymmetryTransforms.permuteAxes(s, "ABC->CAB"), s.cyclicShift().cyclicShift());
		assertSame(SymmetryTransforms.permuteAxes(s, "ABC->ACB"), s.cyclicShift().transpose());
		assertSame(SymmetryTransforms.permuteAxes(s, "ABC->BAC"), s.transpose().cyclicShift());
	}

	@Test
	public void permuteAxes_accepts_NMP_letters_and_rejects_nonbijective() throws Exception {
		NonCubicBilinearAlgorithm s = loadStrassen();
		// NMP letters are the modern alias of ABC; PMN == CBA == transpose.
		assertSame(SymmetryTransforms.permuteAxes(s, "NMP->PMN"), s.transpose());
		// Legacy / malformed perms are not bijections of {A,B,C} → null, caller falls back.
		assertThat(SymmetryTransforms.permuteAxes(s, "ABC->ABA")).isNull();
		assertThat(SymmetryTransforms.permuteAxes(s, "garbage")).isNull();
	}

	/**
	 * THE guard: for an EQUAL-AXIS base, every orbit variant's recorded perm
	 * reproduces that variant bit-exactly — so {@code Transpose(canonical, perm)}
	 * replays to exactly the orientation the search scored. The old shape-based
	 * {@code orientAs} replay could not distinguish the equal-axis orientations
	 * and so failed this for at least one variant.
	 */
	@Test
	public void s3OrbitWithPerms_each_perm_reproduces_its_variant_bit_exact() throws Exception {
		// ⟨2,4,4⟩ = Strassen ⊗ naive(1,2,2): asymmetric AND two equal axes (m=p=4),
		// exactly the ambiguous family that triggered the ⟨13,13,17⟩ divergence.
		NonCubicBilinearAlgorithm base = Compose.kroneckerGeneral(
				loadStrassen(), NonCubicBilinearAlgorithm.naive(1, 2, 2));
		assertThat(base.m).as("equal axes present").isEqualTo(base.p);

		List<SymmetryTransforms.S3Variant> orbit = SymmetryTransforms.s3OrbitWithPerms(base);
		assertThat(orbit).as("non-trivial orbit").hasSizeGreaterThan(1);
		for (SymmetryTransforms.S3Variant v : orbit) {
			NonCubicBilinearAlgorithm rebuilt = SymmetryTransforms.permuteAxes(base, v.perm());
			assertThat(rebuilt).as("perm %s must be a recognised bijection", v.perm()).isNotNull();
			assertSame(rebuilt, v.alg());
		}
	}

	/**
	 * Demonstrates WHY the exact perm is load-bearing: swapping the two equal
	 * axes ("ABC-&gt;ACB") yields a scheme of the SAME shape but DIFFERENT
	 * coefficients. A shape-only record cannot tell the two apart, which is the
	 * bug the exact-perm pinning fixes.
	 */
	@Test
	public void swapping_equal_axes_changes_coefficients_same_shape() throws Exception {
		NonCubicBilinearAlgorithm base = Compose.kroneckerGeneral(
				loadStrassen(), NonCubicBilinearAlgorithm.naive(1, 2, 2));
		NonCubicBilinearAlgorithm swapped = SymmetryTransforms.permuteAxes(base, "ABC->ACB");
		// Same shape (m and p both 4) …
		assertThat(new int[] { swapped.n, swapped.m, swapped.p })
				.isEqualTo(new int[] { base.n, base.m, base.p });
		// … but genuinely a different scheme, so orient-by-shape is ambiguous.
		boolean identical = java.util.Arrays.deepEquals(swapped.denseV(), base.denseV())
				&& java.util.Arrays.deepEquals(swapped.denseW(), base.denseW());
		assertThat(identical).as("equal-axis swap is non-trivial").isFalse();
		// Both still compute the matmul.
		assertThat(Verifier.isExactNonCubic(swapped)).isTrue();
	}
}
