package eu.solven.matmul.search;

import static org.assertj.core.api.Assertions.assertThat;

import java.io.File;

import org.junit.jupiter.api.Test;

import eu.solven.matmul.NonCubicBilinearAlgorithm;
import eu.solven.matmul.catalog.FieldAwareLookup;
import eu.solven.matmul.isotropy.PairedSubProducts;
import eu.solven.matmul.recombination.Recombination;

/**
 * Guards for leaf-level Pan pairing (fmm-gap 2026-07-09, the ⟨20,23,23⟩=5906
 * decode). The silent limitation: the pairing-aware sweep ran only over
 * naive-grid bases on the (false) belief that combining bases cannot carry a
 * fusable pair — a recombination's leaves are formally independent bilinear
 * problems, so any two same-shape cubic leaves fuse via PanPairProduct
 * regardless of which blocks feed them. The device closed the two biggest
 * remaining FMM gaps (⟨20,23,23⟩ −39 and ⟨19,19,22⟩ −2) at exact ties.
 */
public class TestLeafPairFusion {

	/** applyPairingWithMatching pairs same-shape cubic slots and prices them at
	 *  pairCost; the matching covers every slot exactly once. */
	@Test
	public void matching_pairs_cubic_same_shape_slots()  {
		int[][] shapes = { {11,11,11}, {9,12,12}, {11,11,11}, {5,5,5} };
		Recombination.SotaResolver sota = (a, b, c) -> a * b * c; // toy: naive ranks
		PairedSubProducts.Matching m = PairedSubProducts.applyPairingWithMatching(shapes, sota);
		assertThat(m.pairs()).hasDimensions(1, 2);
		assertThat(m.pairs()[0]).containsExactlyInAnyOrder(0, 2);
		assertThat(m.solo()).containsExactlyInAnyOrder(1, 3);
		// paired cubes at pairCost(11)=1694 < 2·1331; solos at naive.
		assertThat(m.cost()).isEqualTo(PairedSubProducts.pairCost(11, 11, 11) + 9*12*12 + 5*5*5);
	}

	/** The persisted leaf-paired stubs replay to exactly their claimed ranks
	 *  through the new RecombinationWithPairN replay path. */
	@Test
	public void leaf_paired_stubs_replay_exactly() {
		FieldAwareLookup lookup = new FieldAwareLookup("Q");
		LineageReplayer replayer = LineageReplayer.withDefaultPool(lookup);
		File s1 = new File("src/main/resources/schemes/derived/section23/20x23x23-r5906-derived-e86f451.json");
		File s2 = new File("src/main/resources/schemes/derived/section22/19x19x22-r4536-derived-e96d27a.json");
		assertThat(s1).exists();
		assertThat(s2).exists();
		NonCubicBilinearAlgorithm a1 = replayer.replayFromFile(s1);
		assertThat(a1.r).isEqualTo(5906);
		NonCubicBilinearAlgorithm a2 = replayer.replayFromFile(s2);
		assertThat(a2.r).isEqualTo(4536);
	}
}
