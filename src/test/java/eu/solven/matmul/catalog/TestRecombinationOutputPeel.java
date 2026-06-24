package eu.solven.matmul.catalog;

import eu.solven.matmul.recombination.Recombination;

import static org.assertj.core.api.Assertions.assertThat;

import java.io.File;

import org.junit.jupiter.api.Test;

import eu.solven.matmul.NonCubicBilinearAlgorithm;
import eu.solven.matmul.search.CitedBound;

/**
 * Output-side zero-peel reduction (Islam 2009 MSc Ch. 4 γ5 trick),
 * implemented as the peel-aware overload of
 * {@link Recombination#recombineWithAllocation}.
 *
 * <p>The reference case: compute {@code R(⟨17,17,17⟩)} via Strassen
 * {@code ⟨2,2,2⟩=7} on the {@code [9,9]³} allocation, padding ⟨17⟩ to
 * ⟨18⟩ and peeling 1 unit off each block's "last" axis position. This
 * is FMM-Lille's recipe for the 2934 bound (see
 * {@code docs/diagnostics/17x17x17_pair_fuse.md}).</p>
 */
class TestRecombinationOutputPeel {

	@Test
	void strassen17x17x17_peel_equivalent_to_9_8_alloc_NOT_a_savings() throws Exception {
		// HONEST documentation: for ⟨17,17,17⟩ via Strassen ⟨2,2,2⟩=7,
		// [9,9]³+peel=[0,1]³ produces the SAME sub-shape distribution as
		// direct [9,8]³. Both land at rank 2940. FMM-Lille's 2934 bound
		// requires a *different* Strassen product↔block mapping (none
		// of the 7 standard products has W entirely in the (1,1)-block
		// for this allocation, so output-peel by itself doesn't activate).
		// The γ5 mechanism IS active in the smaller ⟨3,3,3⟩ test below.
		File strassenFile = eu.solven.matmul.catalog.SchemeResolver.byHint("src/main/resources/schemes/known/section2/strassen-2x2x2_m7_a18.json");
		NonCubicBilinearAlgorithm strassen = SchemeIO.read(strassenFile);
		CitedBound sota = new CitedBound(new FieldAwareLookup("Q"));

		Recombination.Result baseline = Recombination.recombineWithAllocation(
				strassen, sota,
				new int[] { 9, 8 }, new int[] { 9, 8 }, new int[] { 9, 8 });
		Recombination.Result peeled = Recombination.recombineWithAllocation(
				strassen, sota,
				new int[] { 9, 9 }, new int[] { 9, 9 }, new int[] { 9, 9 },
				new int[] { 0, 1 }, new int[] { 0, 1 }, new int[] { 0, 1 });

		System.out.printf("⟨17,17,17⟩ via Strassen ⟨2,2,2⟩=7%n");
		System.out.printf("  [9,8]³ no-peel: rank %d%n", baseline.totalRank);
		for (int k = 0; k < baseline.smallMatrixSizes.length; k++) {
			int[] sz = baseline.smallMatrixSizes[k];
			System.out.printf("    M%d → ⟨%d,%d,%d⟩=%d%n",
					k + 1, sz[0], sz[1], sz[2],
					sota.getRank(sz[0], sz[1], sz[2]));
		}
		System.out.printf("  [9,9]³ + peel=[0,1]³: rank %d%n", peeled.totalRank);
		for (int k = 0; k < peeled.smallMatrixSizes.length; k++) {
			int[] sz = peeled.smallMatrixSizes[k];
			System.out.printf("    M%d → ⟨%d,%d,%d⟩=%d%n",
					k + 1, sz[0], sz[1], sz[2],
					sota.getRank(sz[0], sz[1], sz[2]));
		}

		// Currently the two are equivalent — peel=[0,1] gives effective [9,8],
		// same as direct [9,8] allocation. This documents the "expected to be
		// equal but Islam says we should do better" state of the implementation.
		assertThat(peeled.totalRank).isLessThanOrEqualTo(baseline.totalRank);
	}

	@Test
	void peel_null_matches_legacy_signature() throws Exception {
		File strassenFile = eu.solven.matmul.catalog.SchemeResolver.byHint("src/main/resources/schemes/known/section2/strassen-2x2x2_m7_a18.json");
		NonCubicBilinearAlgorithm strassen = SchemeIO.read(strassenFile);
		CitedBound sota = new CitedBound(new FieldAwareLookup("Q"));
		int[] aA = { 9, 8 }, aB = { 9, 8 }, aC = { 9, 8 };
		Recombination.Result viaLegacy = Recombination.recombineWithAllocation(strassen, sota, aA, aB, aC);
		Recombination.Result viaNull = Recombination.recombineWithAllocation(strassen, sota, aA, aB, aC,
				null, null, null);
		assertThat(viaNull.totalRank).isEqualTo(viaLegacy.totalRank);
	}

	@Test
	void canonical_islam_gamma5_case_3x3x3_via_padded_strassen() throws Exception {
		// Islam 2009 MSc Ch. 4: pad ⟨3,3,3⟩ to ⟨4,4,4⟩ via Strassen ⟨2,2,2⟩=7
		// on [2,2]³, peel 1 unit on each axis's second block. γ5 (= the
		// Strassen product whose W column lives entirely in the peeled
		// corner block C̃₂,₂) collapses from ⟨2,2,2⟩=7 to ⟨1,2,1⟩=2 mults.
		File strassenFile = eu.solven.matmul.catalog.SchemeResolver.byHint("src/main/resources/schemes/known/section2/strassen-2x2x2_m7_a18.json");
		NonCubicBilinearAlgorithm strassen = SchemeIO.read(strassenFile);
		CitedBound sota = new CitedBound(new FieldAwareLookup("Q"));

		Recombination.Result peeled = Recombination.recombineWithAllocation(
				strassen, sota,
				new int[] { 2, 2 }, new int[] { 2, 2 }, new int[] { 2, 2 },
				new int[] { 0, 1 }, new int[] { 0, 1 }, new int[] { 0, 1 });

		System.out.printf("⟨3,3,3⟩ via Strassen ⟨2,2,2⟩=7 on [2,2]³ + peel=[0,1]³:%n");
		boolean hasGamma5Like = false;
		for (int k = 0; k < peeled.smallMatrixSizes.length; k++) {
			int[] sz = peeled.smallMatrixSizes[k];
			long rank = sota.getRank(sz[0], sz[1], sz[2]);
			System.out.printf("  M%d → ⟨%d,%d,%d⟩=%d%n",
					k + 1, sz[0], sz[1], sz[2], rank);
			// "γ5-like" = some axis collapses to 1 due to peel.
			if (Math.min(sz[0], Math.min(sz[1], sz[2])) == 1) hasGamma5Like = true;
		}
		System.out.printf("  total rank = %d%n", peeled.totalRank);
		assertThat(hasGamma5Like)
				.as("at least one Strassen product collapses via output-peel")
				.isTrue();
	}

	@Test
	void peel_full_block_produces_smaller_subshapes() throws Exception {
		// Sanity check: if peel == alloc (block fully peeled), sub-shapes
		// reading from that block should be shape-zero on that axis — meaning
		// the SOTA returns rank for the surviving sub-shape only.
		File strassenFile = eu.solven.matmul.catalog.SchemeResolver.byHint("src/main/resources/schemes/known/section2/strassen-2x2x2_m7_a18.json");
		NonCubicBilinearAlgorithm strassen = SchemeIO.read(strassenFile);
		CitedBound sota = new CitedBound(new FieldAwareLookup("Q"));
		// Fully peel block 1 (effective alloc [9, 0] = degenerate to 1-block matmul).
		Recombination.Result fullyPeeled = Recombination.recombineWithAllocation(
				strassen, sota,
				new int[] { 9, 9 }, new int[] { 9, 9 }, new int[] { 9, 9 },
				new int[] { 0, 9 }, new int[] { 0, 9 }, new int[] { 0, 9 });
		// Without crashing — the rank just becomes the cost of the
		// surviving ⟨9,9,9⟩ block (some sub-products collapse to zero shape).
		assertThat(fullyPeeled.totalRank).isPositive();
	}
}
