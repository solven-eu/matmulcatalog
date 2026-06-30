package eu.solven.matmul.catalog;

import static org.assertj.core.api.Assertions.assertThat;

import java.io.File;

import org.junit.jupiter.api.Test;

import eu.solven.matmul.NonCubicBilinearAlgorithm;
import eu.solven.matmul.verifiers.Verifier;
import eu.solven.matmul.catalog.SchemeIO;

/**
 * Guards the SPARSE (no-dense-round-trip) rewrite of {@link Compose} operators (task b:
 * stop materialising {@code denseU/V/W} in the build path). The correctness criterion is
 * the strongest one: the composed scheme must be an EXACT matmul of the composed shape —
 * any row-reindexing slip makes {@link Verifier#isExactNonCubic} fail.
 */
public class TestComposeSparse {

	@Test
	public void concatRight_of_naive_blocks_is_exact() {
		// ⟨2,3,2⟩ + ⟨2,3,3⟩ → ⟨2,3,5⟩ (P-axis concat); naive blocks are exact by construction.
		NonCubicBilinearAlgorithm left = NonCubicBilinearAlgorithm.naive(2, 3, 2);
		NonCubicBilinearAlgorithm right = NonCubicBilinearAlgorithm.naive(2, 3, 3);
		NonCubicBilinearAlgorithm cat = Compose.concatRight(left, right);

		assertThat(cat.n).isEqualTo(2);
		assertThat(cat.m).isEqualTo(3);
		assertThat(cat.p).isEqualTo(5);
		assertThat(cat.r).isEqualTo(left.r + right.r);
		assertThat(Verifier.isExactNonCubic(cat)).isTrue();
	}

	@Test
	public void concatRight_of_real_strassen_blocks_is_exact() throws java.io.IOException {
		// Strassen ⟨2,2,2⟩=7 concatenated with itself on P → exact ⟨2,2,4⟩=14.
		NonCubicBilinearAlgorithm s = SchemeIO.read(
				new File("src/main/resources/schemes/known/section2/strassen-2x2x2_m7_a18.json"));
		NonCubicBilinearAlgorithm cat = Compose.concatRight(s, s);

		assertThat(cat.n).isEqualTo(2);
		assertThat(cat.m).isEqualTo(2);
		assertThat(cat.p).isEqualTo(4);
		assertThat(cat.r).isEqualTo(14);
		assertThat(Verifier.isExactNonCubic(cat)).isTrue();
	}

	@Test
	public void concatBelow_of_naive_blocks_is_exact() {
		// ⟨1,3,4⟩ + ⟨2,3,4⟩ → ⟨3,3,4⟩ (N-axis concat, uneven).
		NonCubicBilinearAlgorithm cat = Compose.concatBelow(
				NonCubicBilinearAlgorithm.naive(1, 3, 4), NonCubicBilinearAlgorithm.naive(2, 3, 4));
		assertThat(cat.n).isEqualTo(3);
		assertThat(cat.m).isEqualTo(3);
		assertThat(cat.p).isEqualTo(4);
		assertThat(Verifier.isExactNonCubic(cat)).isTrue();
	}

	@Test
	public void concatInner_of_naive_blocks_is_exact() {
		// ⟨3,2,4⟩ + ⟨3,5,4⟩ → ⟨3,7,4⟩ (M-axis contraction sum, uneven).
		NonCubicBilinearAlgorithm cat = Compose.concatInner(
				NonCubicBilinearAlgorithm.naive(3, 2, 4), NonCubicBilinearAlgorithm.naive(3, 5, 4));
		assertThat(cat.n).isEqualTo(3);
		assertThat(cat.m).isEqualTo(7);
		assertThat(cat.p).isEqualTo(4);
		assertThat(Verifier.isExactNonCubic(cat)).isTrue();
	}

	@Test
	public void project_drops_axis_index_exact() {
		// Drop N-index 2 from naive ⟨3,3,3⟩ → exact ⟨2,3,3⟩ (sparse DCE projection).
		NonCubicBilinearAlgorithm proj = Compose.project(
				NonCubicBilinearAlgorithm.naive(3, 3, 3),
				new int[] { 0, 1 }, new int[] { 0, 1, 2 }, new int[] { 0, 1, 2 });
		assertThat(proj.n).isEqualTo(2);
		assertThat(proj.m).isEqualTo(3);
		assertThat(proj.p).isEqualTo(3);
		assertThat(Verifier.isExactNonCubic(proj)).isTrue();
	}

	@Test
	public void project_drops_each_axis_exact() {
		// Drop one index per axis from naive ⟨4,4,4⟩ → exact ⟨3,3,3⟩.
		NonCubicBilinearAlgorithm proj = Compose.project(
				NonCubicBilinearAlgorithm.naive(4, 4, 4),
				new int[] { 0, 1, 2 }, new int[] { 0, 1, 3 }, new int[] { 1, 2, 3 });
		assertThat(proj.n).isEqualTo(3);
		assertThat(proj.m).isEqualTo(3);
		assertThat(proj.p).isEqualTo(3);
		assertThat(Verifier.isExactNonCubic(proj)).isTrue();
	}

	@Test
	public void concatRight_uneven_split_is_exact() {
		// Asymmetric split exercises the distinct p1/p2 row-reindexing paths.
		NonCubicBilinearAlgorithm cat = Compose.concatRight(
				NonCubicBilinearAlgorithm.naive(3, 4, 1), NonCubicBilinearAlgorithm.naive(3, 4, 6));
		assertThat(cat.p).isEqualTo(7);
		assertThat(Verifier.isExactNonCubic(cat)).isTrue();
	}
}
