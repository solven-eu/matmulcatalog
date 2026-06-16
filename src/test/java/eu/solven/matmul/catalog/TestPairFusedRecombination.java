package eu.solven.matmul.catalog;

import eu.solven.matmul.search.PairFusedRecombination;

import eu.solven.matmul.search.BlockSplitSearch;

import static org.assertj.core.api.Assertions.assertThat;

import java.util.Optional;

import org.junit.jupiter.api.Test;

public class TestPairFusedRecombination {

	private static Recombination.SotaResolver fixedSota(int rank777, int rank111111) {
		return (a, b, c) -> {
			if (a == 7 && b == 7 && c == 7) return rank777;
			if (a == 11 && b == 11 && c == 11) return rank111111;
			return -1;
		};
	}

	@Test
	public void predict_14cube_strassen_pair_fused_3pairs_1solo() {
		// ⟨14,14,14⟩ via Strassen-recomb [k=7,k=7]³: 7 outer products → 3 pairs + 1 solo
		// pair-cost  = 7³ + 3·49 = 490
		// solo-cost  = R(⟨7,7,7⟩) = 249  (Sedoglavic 2017 / Perminov import)
		// pair-fused = 3·490 + 249 = 1719   (vs 7·249 = 1743 solo recomb — saves 24)
		Recombination.SotaResolver sota = fixedSota(249, /* unused */ -1);
		Optional<PairFusedRecombination.Prediction> pred =
				PairFusedRecombination.predict(14, 7, 7, sota);
		assertThat(pred).isPresent();
		assertThat(pred.get().pairs()).isEqualTo(3);
		assertThat(pred.get().solos()).isEqualTo(1);
		assertThat(pred.get().pairCost()).isEqualTo(490);
		assertThat(pred.get().soloCost()).isEqualTo(249);
		assertThat(pred.get().totalRank()).isEqualTo(1719);
	}

	@Test
	public void predict_22cube_strassen_pair_fused() {
		// ⟨22,22,22⟩ via Strassen-recomb [k=11,k=11]³: 7 outer products → 3 pairs + 1 solo
		// pair-cost  = 11³ + 3·121 = 1331 + 363 = 1694
		// solo-cost  = R(⟨11,11,11⟩) = 896  (Smirnov 2013)
		// pair-fused = 3·1694 + 896 = 5978  (vs 7·896 = 6272 solo recomb — saves 294)
		Recombination.SotaResolver sota = fixedSota(/* unused */ -1, 896);
		Optional<PairFusedRecombination.Prediction> pred =
				PairFusedRecombination.predict(22, 11, 7, sota);
		assertThat(pred).isPresent();
		assertThat(pred.get().totalRank()).isEqualTo(5978);
	}

	@Test
	public void predict_returns_empty_when_not_profitable() {
		// 2·R = 1024 (=2·512); pair-cost = 8³ + 3·64 = 704 < 1024 — profitable
		// Edge case: when solo-cost is too low, pair-cost ≥ 2·solo and empty.
		Recombination.SotaResolver tinyLeaf = (a, b, c) -> 100;  // wildly low solo
		Optional<PairFusedRecombination.Prediction> pred =
				PairFusedRecombination.predict(16, 8, 7, tinyLeaf);
		assertThat(pred).isEmpty();
	}

	@Test
	public void predictBalancedCubic_rejects_non_cubic_base() {
		// Non-cubic base ⟨2,3,4⟩ — predictBalancedCubic should return empty.
		double[][] U = new double[6][1], V = new double[12][1], W = new double[8][1];
		U[0][0] = V[0][0] = W[0][0] = 1;
		eu.solven.matmul.NonCubicBilinearAlgorithm noncubic =
				new eu.solven.matmul.NonCubicBilinearAlgorithm(2, 3, 4, U, V, W);
		Optional<PairFusedRecombination.Prediction> pred =
				PairFusedRecombination.predictBalancedCubic(24, noncubic, fixedSota(249, 896));
		assertThat(pred).isEmpty();
	}

	@Test
	public void chooseBest_smallLeaf_pairFusedWins() {
		// 7 leaves of ⟨7,7,7⟩ (the canonical small-leaf case):
		//   solo     = 7·249             = 1743
		//   pair     = 3·490 + 1·249     = 1719   ← wins
		//   full-TA  = 7·390 (Islam odd) = 2730   (corrections dominate at k=7)
		Recombination.SotaResolver sota = fixedSota(249, -1);
		Optional<PairFusedRecombination.TaChoice> choice =
				PairFusedRecombination.chooseBest(7, 7, sota);
		assertThat(choice).isPresent();
		assertThat(choice.get().soloTotal()).isEqualTo(1743);
		assertThat(choice.get().pairTotal()).isEqualTo(1719);
		assertThat(choice.get().fullTaTotal()).isEqualTo(2730);
		assertThat(choice.get().best()).isEqualTo(PairFusedRecombination.TaStrategy.PAIR_FUSED);
		assertThat(choice.get().bestTotal()).isEqualTo(1719);
	}

	@Test
	public void chooseBest_largeEvenLeaf_fullTaWins() {
		// 7 leaves of ⟨20,20,20⟩ with an (artificially high) catalog rank so
		// the deep single-product TA dominates:
		//   pair/2  = (20³+3·20²)/2 = 4600 per product
		//   full-TA = bestPanTaBound(20) = 4340 per product (Islam even)  ← wins
		Recombination.SotaResolver sota = (a, b, c) -> (a == 20 && b == 20 && c == 20) ? 9000 : -1;
		Optional<PairFusedRecombination.TaChoice> choice =
				PairFusedRecombination.chooseBest(20, 7, sota);
		assertThat(choice).isPresent();
		assertThat(choice.get().fullTaTotal()).isEqualTo(7L * 4340);
		assertThat(choice.get().best()).isEqualTo(PairFusedRecombination.TaStrategy.FULL_TA);
	}

	@Test
	public void chooseBest_cheapLeaf_soloWins() {
		// When the catalog rank is already very low, neither TA option helps:
		//   solo = 7·100 = 700; pair-cost 490 ≥ 2·100 ⇒ pairing inert (=700);
		//   full-TA = 7·390 = 2730. SOLO wins.
		Recombination.SotaResolver sota = fixedSota(100, -1);
		Optional<PairFusedRecombination.TaChoice> choice =
				PairFusedRecombination.chooseBest(7, 7, sota);
		assertThat(choice).isPresent();
		assertThat(choice.get().best()).isEqualTo(PairFusedRecombination.TaStrategy.SOLO);
		assertThat(choice.get().bestTotal()).isEqualTo(700);
	}

	@Test
	public void predict_strategy_visible_in_findBestStrategy() {
		// End-to-end: a fixed-rank SOTA where pair-fused STRICTLY dominates
		// solo recomb. We force R(⟨7,7,7⟩) = 280 → 7·280 = 1960 solo,
		// pair-fused = 3·490 + 280 = 1750 (saves 210). Pair-fused should win.
		Recombination.SotaResolver biased = fixedSota(280, -1);
		java.util.List<BlockSplitSearch.NamedBase> pool = BlockSplitSearch.defaultPool();
		Optional<BlockSplitSearch.NonCubicStrategy> picked =
				BlockSplitSearch.findBestStrategy(14, 14, 14, pool, biased, true);
		assertThat(picked).isPresent();
		// We expect pair-fused to be the chosen strategy under this biased resolver.
		assertThat(picked.get().rank()).isLessThanOrEqualTo(1750L);
	}
}
