package eu.solven.matmul.catalog;

import static org.assertj.core.api.Assertions.assertThat;

import java.io.File;
import java.io.IOException;
import java.nio.file.Files;
import java.nio.file.Path;
import java.util.Optional;
import java.util.stream.Stream;

import org.junit.jupiter.api.Test;

import eu.solven.matmul.NonCubicBilinearAlgorithm;
import eu.solven.matmul.verifiers.Verifier;
import eu.solven.matmul.recombination.Recombination;
import eu.solven.matmul.catalog.SchemeIO;
import eu.solven.matmul.recombination.Recombination.AlgorithmLookup;

/**
 * Probes whether {@link Recombination#construct} on Strassen ⟨2,2,2⟩ base
 * with target ⟨7,7,7⟩ produces Sedoglavic's 250-mult algorithm
 * (the algebraic saving comes from {@code _process_additions}'s
 * cross-block min-reduction). Expected to land somewhere between
 * the naïve block-decomposition rank 273 and Sedoglavic's 250.
 */
public class TestStrassenRecombination777 {

	// Content-driven catalog-best resolver over Q (the char-0 default field; excludes
	// F2-only schemes via fields[], not filenames). Replaces a Files.walk+findFirst
	// FILENAME matcher whose pick was walk-order-dependent: as the catalog grew it
	// started resolving DOMINATED files (⟨3,4,4⟩ → a r=40 derived stub instead of the
	// 38-rank best), silently inflating the built rank (250 → 256).
	private static final FieldAwareLookup LOOKUP = new FieldAwareLookup("Q");

	private static AlgorithmLookup catalogLookup() {
		return LOOKUP::find;
	}

	@Test
	public void strassen_recombination_777_with_forced_split_4_3_yields_sedoglavic_rank() throws IOException {
		NonCubicBilinearAlgorithm strassen = SchemeIO.read(
				eu.solven.matmul.catalog.SchemeResolver.byHint("src/main/resources/schemes/known/section2/strassen-2x2x2_m7_a18.json"));
		assertThat(strassen.r).isEqualTo(7);

		Recombination.SotaResolver sota = (a, b, c) -> {
			if (a == 0 || b == 0 || c == 0) return 0;
			// Skip catalog entries for the full ⟨7,7,7⟩ target itself: we want
			// to force a constructive build, not a direct lookup.
			if (a == 7 && b == 7 && c == 7) return Integer.MAX_VALUE / 100;
			return catalogLookup().find(a, b, c).map(alg -> alg.r).orElse(Integer.MAX_VALUE / 100);
		};

		// Force split [4, 3] each axis — the Sedoglavic recipe.
		Recombination.Result rec = Recombination.recombineWithAllocation(
				strassen, sota, new int[] { 4, 3 }, new int[] { 4, 3 }, new int[] { 4, 3 });
		System.out.println("recombine ⟨7,7,7⟩ via Strassen with FORCED alloc [4,3]³: rank=" + rec.totalRank);
		for (int k = 0; k < strassen.r; k++) {
			System.out.println("  base mult " + k + " → sub-shape "
					+ java.util.Arrays.toString(rec.smallMatrixSizes[k]));
		}
		// Expectation: with split [4,3] and Strassen's cross-block sharing, recombination
		// should find a rank lower than the naïve block decomposition (273) — ideally 250.
		assertThat(rec.totalRank).isLessThanOrEqualTo(273);
	}

	@Test
	public void constructive_sedoglavic_777_via_strassen_4_3_split_yields_250_and_verifies()
			throws IOException {
		NonCubicBilinearAlgorithm strassen = SchemeIO.read(
				eu.solven.matmul.catalog.SchemeResolver.byHint("src/main/resources/schemes/known/section2/strassen-2x2x2_m7_a18.json"));

		// Constructive Sedoglavic: Strassen ⟨2,2,2⟩ outer × mixed-shape inner via
		// forced split [4, 3] on each axis. Bounds-checked embedding correctly
		// handles the padded zero positions where block (1,1) is 3×3 but sub-mults
		// see 4×4 inputs.
		NonCubicBilinearAlgorithm built = Recombination.constructWithAllocation(
				strassen, catalogLookup(),
				new int[] { 4, 3 }, new int[] { 4, 3 }, new int[] { 4, 3 });
		assertThat(built.n).isEqualTo(7);
		// The deterministic forced-[4,3] constructive build — Sedoglavic 2017's
		// recombination, 250 with 2017-era leaves (⟨4,4,4⟩=49, ⟨3,4,4⟩=38, …).
		// SOTA-or-better, not equality (repo guard rule): catalog-best leaves can
		// only improve the total (e.g. the DPS-2025 ⟨4,4,4⟩=48 rationalisation),
		// and an improvement must never break the guard. The global free-allocation
		// SOTA ⟨7,7,7⟩=249 (Perminov) is a different construction.
		assertThat(built.r).isLessThanOrEqualTo(250);
		assertThat(Verifier.isExactNonCubic(built)).isTrue();
	}
}
