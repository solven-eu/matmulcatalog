package eu.solven.matmul.recombination;

import static org.assertj.core.api.Assertions.assertThat;

import org.junit.jupiter.api.Test;

import eu.solven.matmul.NonCubicBilinearAlgorithm;
import eu.solven.matmul.recombination.Recombination.SotaResolver;
import eu.solven.matmul.recombination.Recombination.TaFusionBreakdown;

/**
 * Guards {@link Recombination#describeTaFusion} — the Pan-TA highlight that makes a
 * naïve-grid recombination's final rank explainable (which products fuse, how much TA
 * saves). The canonical case mirrors master's ⟨26,29,29⟩=11693 scheme: a {@code
 * naive-1x2x2} grid whose ⟨26,3,26⟩ and ⟨26,26,3⟩ products are a disjoint cyclic pair
 * that TA fuses.
 */
public class TestTaFusionBreakdown {

	/** Naïve sota — every sub-shape priced at the cubic a·b·c (always buildable). */
	private static final SotaResolver NAIVE = (a, b, c) -> a * b * c;

	@Test
	public void naive_1x2x2_grid_fuses_one_cyclic_pair() {
		NonCubicBilinearAlgorithm base = NonCubicBilinearAlgorithm.naive(1, 2, 2);
		assertThat(Recombination.isNaiveGrid(base)).isTrue();

		// n=26 (unsplit), m=29 -> [3,26], p=29 -> [3,26]: the ⟨1,2,2⟩ grid's 4 products
		// are ⟨26,3,3⟩, ⟨26,3,26⟩, ⟨26,26,3⟩, ⟨26,26,26⟩; the middle two are a rot² pair.
		TaFusionBreakdown bd = Recombination.describeTaFusion(
				base, NAIVE, new int[] { 26 }, new int[] { 3, 26 }, new int[] { 3, 26 });

		assertThat(bd).isNotNull();
		assertThat(bd.hasFusion()).isTrue();
		assertThat(bd.fusedPairs()).hasSize(1);

		TaFusionBreakdown.FusedPair fp = bd.fusedPairs().get(0);
		// fused cost = abc + ab + bc + ca = 2028 + 78 + 78 + 676 = 2860
		assertThat(fp.fusedCost()).isEqualTo(2860);
		// naive (un-fused) = 2 * 26*3*26 = 4056; saving = 4056 - 2860 = 1196
		assertThat(fp.naiveRank()).isEqualTo(4056);
		assertThat(fp.saving()).isEqualTo(1196);
		assertThat(bd.taSaving()).isEqualTo(1196);

		// unpaired = ⟨26,3,3⟩(234) + ⟨26,26,26⟩(17576) = 17810
		assertThat(bd.unpairedLeafSum()).isEqualTo(17810);
		assertThat(bd.unpairedShapes()).hasSize(2);

		// the load-bearing identity: total = unpaired + fused, and TA is what bought the gap.
		assertThat(bd.totalRank()).isEqualTo(bd.unpairedLeafSum() + bd.fusedCost());
		assertThat(bd.totalRank()).isEqualTo(20670);
		assertThat(bd.summary()).contains("Pan-TA saved 1196").contains("⟨26,3,26⟩");
	}

	@Test
	public void block_combining_base_has_no_ta_structure() {
		// A Strassen 2x2x2 base is NOT a naïve grid (products combine blocks), so there is
		// no disjoint cyclic pair to fuse — describeTaFusion returns null, never a fake "0 saved".
		NonCubicBilinearAlgorithm base = NonCubicBilinearAlgorithm.naive(2, 2, 2); // 8-mult naive cube
		assertThat(Recombination.isNaiveGrid(base)).isTrue(); // naive cube IS a grid
		TaFusionBreakdown bd = Recombination.describeTaFusion(
				base, NAIVE, new int[] { 1, 1 }, new int[] { 1, 1 }, new int[] { 1, 1 });
		// all eight products are ⟨1,1,1⟩ — cyclic pairs exist but every pair shares blocks
		// (not disjoint) OR saves nothing, so no fusion is profitable here.
		assertThat(bd).isNotNull();
		assertThat(bd.hasFusion()).isFalse();
		assertThat(bd.summary()).contains("no Pan-TA pair fused");
	}
}
