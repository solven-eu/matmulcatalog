package eu.solven.matmul.catalog;

import static org.assertj.core.api.Assertions.assertThat;

import java.io.File;

import org.junit.jupiter.api.Test;

import eu.solven.matmul.NonCubicBilinearAlgorithm;
import eu.solven.matmul.Verifier;

public class TestComposeConcat {

	@Test
	public void concatRight_2x10x6_and_2x10x10_gives_2x10x16() throws Exception {
		// FMM ⟨2,10,16⟩=248 decomposes as ⟨2,10,6⟩ + ⟨2,10,10⟩. Verify
		// concatRight applied to our local schemes (Perminov ⟨2,6,10⟩=94
		// + FMM ⟨2,10,10⟩=155) yields a valid ⟨2,10,16⟩ scheme of rank 249.
		File leftFile = eu.solven.matmul.catalog.SchemeResolver.byHint("src/main/resources/schemes/known/section10/perminov_cr325_fv153_cn668_ZT_reduced-2x6x10_m94_a668.json");
		NonCubicBilinearAlgorithm left = SchemeIO.readBilinear(leftFile);
		// ⟨2,6,10⟩ is our convention's ⟨2,6,10⟩ — but we need ⟨2,10,6⟩
		// (m=10, p=6). Use orientAs to permute.
		NonCubicBilinearAlgorithm leftReordered = left.orientAs(2, 10, 6)
				.orElseThrow(() -> new AssertionError("orientAs ⟨2,10,6⟩ failed"));

		// ⟨2,10,10⟩=155 — Hopcroft-Kerr 1971 (was imported via fmm-lille; the
		// 2026-06 rename re-attributed it to hk71 — hint updated accordingly).
		NonCubicBilinearAlgorithm right = SchemeIO.read(eu.solven.matmul.catalog.SchemeResolver.byHint(
				"src/main/resources/schemes/known/section10/2x10x10-r155-hk71-177ef8b.json"));

		NonCubicBilinearAlgorithm combined = Compose.concatRight(leftReordered, right);
		assertThat(combined.n).isEqualTo(2);
		assertThat(combined.m).isEqualTo(10);
		assertThat(combined.p).isEqualTo(16);
		assertThat(combined.r).isEqualTo(94 + 155);
		assertThat(Verifier.passesRandomMatmulSpotCheck(combined)).isTrue();
	}

	@Test
	public void concatBelow_strassen_with_strassen_gives_4x2x2() throws Exception {
		// Sanity: two ⟨2,2,2⟩=7 stacked vertically should give ⟨4,2,2⟩=14.
		NonCubicBilinearAlgorithm s = SchemeIO.read(eu.solven.matmul.catalog.SchemeResolver.byHint(
				"src/main/resources/schemes/known/section2/strassen-2x2x2_m7_a18.json"));
		NonCubicBilinearAlgorithm combined = Compose.concatBelow(s, s);
		assertThat(combined.n).isEqualTo(4);
		assertThat(combined.m).isEqualTo(2);
		assertThat(combined.p).isEqualTo(2);
		assertThat(combined.r).isEqualTo(14);
		assertThat(Verifier.passesRandomMatmulSpotCheck(combined)).isTrue();
	}

	@Test
	public void concatInner_strassen_with_strassen_gives_2x4x2() throws Exception {
		// m-axis (contraction) sum: two ⟨2,2,2⟩=7 along the inner dimension
		// give ⟨2,4,2⟩=14, computing C = A1·B1 + A2·B2. This is the third
		// sibling — it accumulates rather than tiles, but is still exact.
		NonCubicBilinearAlgorithm s = SchemeIO.read(eu.solven.matmul.catalog.SchemeResolver.byHint(
				"src/main/resources/schemes/known/section2/strassen-2x2x2_m7_a18.json"));
		NonCubicBilinearAlgorithm combined = Compose.concatInner(s, s);
		assertThat(combined.n).isEqualTo(2);
		assertThat(combined.m).isEqualTo(4);
		assertThat(combined.p).isEqualTo(2);
		assertThat(combined.r).isEqualTo(14);
		assertThat(Verifier.passesRandomMatmulSpotCheck(combined)).isTrue();
	}

	@Test
	public void concatInner_asymmetric_gives_correct_shape_and_verifies() throws Exception {
		// Mixed inner split: ⟨2,2,2⟩=7 (m1=2) + ⟨2,3,2⟩ naive (m2=3) → ⟨2,5,2⟩,
		// exercising the m1 != m2 path. ⟨2,3,2⟩ naive has rank 12.
		NonCubicBilinearAlgorithm s = SchemeIO.read(eu.solven.matmul.catalog.SchemeResolver.byHint(
				"src/main/resources/schemes/known/section2/strassen-2x2x2_m7_a18.json"));
		NonCubicBilinearAlgorithm naive232 =
				eu.solven.matmul.NaiveMatMul.ofNonCubic(2, 3, 2);
		NonCubicBilinearAlgorithm combined = Compose.concatInner(s, naive232);
		assertThat(combined.n).isEqualTo(2);
		assertThat(combined.m).isEqualTo(5);
		assertThat(combined.p).isEqualTo(2);
		assertThat(combined.r).isEqualTo(7 + 12);
		assertThat(Verifier.passesRandomMatmulSpotCheck(combined)).isTrue();
	}
}
