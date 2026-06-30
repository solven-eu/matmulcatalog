package eu.solven.matmul.catalog;

import static org.assertj.core.api.Assertions.assertThat;

import java.io.File;

import org.junit.jupiter.api.Test;

import eu.solven.matmul.NonCubicBilinearAlgorithm;
import eu.solven.matmul.verifiers.Verifier;

/**
 * Regression guard for the {@code read()} dispatch bug: the common entry point
 * {@link SchemeIO#read(File)} previously only branched on {@code u_sparse}, never on the Perminov
 * "reduced" sparse-list format ({@code u}/{@code v}/{@code w} of {@code {index,value}} objects with
 * optional {@code *_fresh} shared intermediates). So all 168 {@code *_reduced} files — including the
 * best NC integer recombination bases (perminov ⟨2,4,4⟩=26, ⟨4,4,4⟩=49) — threw
 * "row 0 has 1 cols, expected …" when loaded via {@code read}, silently dropping them from any pool
 * that used {@code read} rather than {@code readBilinear}. The fix makes {@code read(JsonNode)}
 * dispatch {@code isReduced} → {@code readReduced} like {@code readBilinear} does.
 */
final class TestSchemeIOReducedDispatch {

	private static final String P244 =
			"src/main/resources/schemes/known/section4/2x4x4-r26-perminov_cr92_cn130_ZT_reduced-8dabbe0.json";
	private static final String P444 =
			"src/main/resources/schemes/known/section4/4x4x4-r49-perminov_cr159_fv100_cn474_ZT_reduced-2c050f4.json";

	@Test
	void read_loads_reduced_perminov_244_and_verifies() throws Exception {
		NonCubicBilinearAlgorithm b = SchemeIO.read(new File(P244));
		assertThat(new int[] { b.n, b.m, b.p, b.r }).containsExactly(2, 4, 4, 26);
		assertThat(Verifier.isExactNonCubic(b)).as("perminov ⟨2,4,4⟩=26 computes matmul").isTrue();
	}

	@Test
	void read_loads_reduced_perminov_444_and_verifies() throws Exception {
		NonCubicBilinearAlgorithm b = SchemeIO.read(new File(P444));
		assertThat(new int[] { b.n, b.m, b.p, b.r }).containsExactly(4, 4, 4, 49);
		assertThat(Verifier.isExactNonCubic(b)).as("perminov ⟨4,4,4⟩=49 computes matmul").isTrue();
	}

	@Test
	void read_and_readBilinear_agree_on_reduced_shape() throws Exception {
		NonCubicBilinearAlgorithm viaRead = SchemeIO.read(new File(P244));
		NonCubicBilinearAlgorithm viaBilinear = SchemeIO.readBilinear(new File(P244));
		assertThat(new int[] { viaRead.n, viaRead.m, viaRead.p, viaRead.r })
				.containsExactly(viaBilinear.n, viaBilinear.m, viaBilinear.p, viaBilinear.r);
	}
}
