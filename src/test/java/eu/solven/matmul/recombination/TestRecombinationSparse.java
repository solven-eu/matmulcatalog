package eu.solven.matmul.recombination;

import static org.assertj.core.api.Assertions.assertThat;

import java.io.File;
import java.util.Optional;

import org.junit.jupiter.api.Test;

import eu.solven.matmul.NonCubicBilinearAlgorithm;
import eu.solven.matmul.Verifier;
import eu.solven.matmul.catalog.SchemeIO;

/**
 * Guards the SPARSE rewrite of {@link Recombination#constructWithAllocation} /
 * {@code constructFromResult} (task b: no more dense {@code [dim²][totalRank]} result, the
 * recombination-replay OOM source). Criterion: the recombined scheme must be EXACT matmul
 * of the composed shape — any block-embed / row-reindex slip fails {@code isExactNonCubic}.
 */
public class TestRecombinationSparse {

	private static NonCubicBilinearAlgorithm strassen() throws java.io.IOException {
		return SchemeIO.read(new File(
				"src/main/resources/schemes/known/section2/strassen-2x2x2_m7_a18.json"));
	}

	@Test
	public void strassen_recombined_to_4x4x4_is_exact() throws java.io.IOException {
		NonCubicBilinearAlgorithm s = strassen();
		// [2,2]³ on a ⟨2,2,2⟩ base → ⟨4,4,4⟩; each of the 7 products is a ⟨2,2,2⟩ block matmul.
		Recombination.AlgorithmLookup lookup =
				(n, m, p) -> (n == 2 && m == 2 && p == 2) ? Optional.of(s) : Optional.empty();
		NonCubicBilinearAlgorithm cat = Recombination.constructWithAllocation(
				s, lookup, new int[] { 2, 2 }, new int[] { 2, 2 }, new int[] { 2, 2 });

		assertThat(cat.n).isEqualTo(4);
		assertThat(cat.m).isEqualTo(4);
		assertThat(cat.p).isEqualTo(4);
		assertThat(cat.r).isEqualTo(49); // 7 base products × 7 (Strassen leaf)
		assertThat(Verifier.isExactNonCubic(cat)).isTrue();
	}

	@Test
	public void ta_fusion_naive_grid_is_exact() {
		// naive-⟨1,2,2⟩ grid, allocA=[10] allocB=[3,10] allocC=[3,10] → ⟨10,13,13⟩. The
		// ⟨10,3,10⟩ & ⟨10,10,3⟩ products are a rot² pair whose TA fusion saves (2·300 > the
		// fused 460), so it fuses — exercising the SPARSE TA build (embedTaPair +
		// embedProductSparse, no dense result). Naive leaves keep the whole thing exact.
		NonCubicBilinearAlgorithm base = NonCubicBilinearAlgorithm.naive(1, 2, 2);
		Recombination.SubResolver resolve = sz -> NonCubicBilinearAlgorithm.naive(sz[0], sz[1], sz[2]);
		Recombination.SotaResolver sota = (a, b, c) -> a * b * c;
		Recombination.TaFusedConstruction tc = Recombination.constructWithTaFusion(
				base, resolve, sota, new int[] { 10 }, new int[] { 3, 10 }, new int[] { 3, 10 });

		assertThat(tc.fusedPairs()).isNotEmpty();   // at least one TA pair fused
		NonCubicBilinearAlgorithm cat = tc.alg();
		assertThat(cat.n).isEqualTo(10);
		assertThat(cat.m).isEqualTo(13);
		assertThat(cat.p).isEqualTo(13);
		assertThat(Verifier.isExactNonCubic(cat)).isTrue();
	}

	@Test
	public void strassen_unbalanced_3_2_split_with_padding_is_exact() throws java.io.IOException {
		NonCubicBilinearAlgorithm s = strassen();
		// [3,2]³ → ⟨5,5,5⟩ with UNEQUAL blocks (3 vs 2): Strassen's block-combining products
		// (e.g. A11+A22) stitch a 3×3 and a 2×2 block, exercising the block-size padding skip.
		// Naive leaves of every shape keep it exact while stressing the sparse embed.
		Recombination.AlgorithmLookup naive =
				(n, m, p) -> Optional.of(NonCubicBilinearAlgorithm.naive(n, m, p));
		NonCubicBilinearAlgorithm cat = Recombination.constructWithAllocation(
				s, naive, new int[] { 3, 2 }, new int[] { 3, 2 }, new int[] { 3, 2 });

		assertThat(cat.n).isEqualTo(5);
		assertThat(cat.m).isEqualTo(5);
		assertThat(cat.p).isEqualTo(5);
		assertThat(Verifier.isExactNonCubic(cat)).isTrue();
	}
}
