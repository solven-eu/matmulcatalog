package eu.solven.matmul.catalog;

import static org.assertj.core.api.Assertions.assertThat;

import java.util.List;
import java.util.Optional;

import org.junit.jupiter.api.Test;

import eu.solven.matmul.recombination.Recombination;
import eu.solven.matmul.recombination.Recombination.TaFusionBreakdown;

/**
 * Guards {@link TaFusionExplainer} — the lineage-node → {@link TaFusionBreakdown} glue
 * the catalog manifest uses to stamp {@code ta_fusion}. Mirrors master's ⟨26,29,29⟩
 * {@code RecombinationTa(base=naive-1x2x2, alloc[26]×[3,26]×[3,26])}.
 */
public class TestTaFusionExplainer {

	private static final Recombination.SotaResolver NAIVE = (a, b, c) -> a * b * c;

	@Test
	public void describes_ta_node_with_naive_grid_base() {
		Lineage.RecombinationTaN node = new Lineage.RecombinationTaN(
				new Lineage.Atom("naive-1x2x2"),
				new int[] { 26 }, new int[] { 3, 26 }, new int[] { 3, 26 },
				List.of());

		Optional<TaFusionBreakdown> bd = TaFusionExplainer.describe(node, NAIVE);

		assertThat(bd).isPresent();
		assertThat(bd.get().hasFusion()).isTrue();
		assertThat(bd.get().fusedPairs()).hasSize(1);
		assertThat(bd.get().taSaving()).isEqualTo(1196);
		assertThat(bd.get().summary()).contains("Pan-TA saved 1196");
	}

	@Test
	public void empty_when_base_is_not_a_naive_grid_atom() {
		// A pinned content-hash base (not "naive-NxMxP") is not a recoverable naïve grid,
		// so the explainer declines rather than guessing.
		Lineage.RecombinationTaN node = new Lineage.RecombinationTaN(
				new Lineage.Atom("2x2x2@deadbeef"),
				new int[] { 1 }, new int[] { 1, 1 }, new int[] { 1, 1 },
				List.of());

		assertThat(TaFusionExplainer.describe(node, NAIVE)).isEmpty();
	}
}
