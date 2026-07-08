package eu.solven.matmul.recombination;

import static org.assertj.core.api.Assertions.assertThat;

import java.util.Arrays;
import java.util.List;
import java.util.stream.Collectors;

import org.junit.jupiter.api.Test;

import eu.solven.matmul.NonCubicBilinearAlgorithm;
import eu.solven.matmul.search.SearchBudget;

/**
 * Human-checkable documentation of HOW the allocation optimiser walks the
 * allocation space (fmm-gap 2026-07-08, the ⟨14,27,27⟩ starvation):
 *
 * <ol>
 * <li>Each axis's composition list is visited <b>balance-first</b>
 *     ({@code imbalanceKey} = Σaᵢ², minimal at balanced) — spelled out
 *     literally in {@link #per_axis_lists_iterate_balance_first()}.</li>
 * <li>The nested A→B→C sweep is LEXICOGRAPHIC over those lists, so on its own
 *     it fully scans B×C (out to extreme imbalance) before advancing A at all —
 *     under a stagnation cap it can die inside the first A-subtree. The
 *     <b>coordinate-descent seeding</b> pass compensates: every ONE-AXIS
 *     deviation from the greedily-updated incumbent is priced BEFORE the sweep,
 *     so single-axis-unbalanced optima are found regardless of the cap —
 *     proven by a budget too starved for the sweep to find anything
 *     ({@link #seeding_finds_single_axis_deviation_under_starved_budget()}).</li>
 * </ol>
 */
public class TestAllocationOrdering {

	/** compositions(6,2) in the exact order the optimiser visits them: balanced
	 *  first, then by growing imbalance (ties keep generation order — stable sort). */
	@Test
	public void per_axis_lists_iterate_balance_first() {
		List<int[]> comps = AllocationOptimizer.compositions(6, 2);
		comps.sort(java.util.Comparator.comparingLong(AllocationOptimizer::imbalanceKey));
		List<String> visited = comps.stream().map(Arrays::toString).collect(Collectors.toList());
		assertThat(visited).containsExactly(
				"[3, 3]",           // Σa² = 18 — balanced, always first
				"[2, 4]", "[4, 2]", // Σa² = 20
				"[1, 5]", "[5, 1]"  // Σa² = 26 — most unbalanced, last
		);
	}

	/**
	 * Synthetic single-axis-optimum case, budget starved to stagnation=1 /
	 * maxNodes=1 so the nested sweep contributes NOTHING: only the seeding pass
	 * can find the answer. Naive ⟨2,3,3⟩ base on target ⟨4,9,9⟩; rank table makes
	 * A=[1,3] (cost 9·1 + 9·5 = 54) beat balanced A=[2,2] (18·10 = 180).
	 */
	@Test
	public void seeding_finds_single_axis_deviation_under_starved_budget() {
		NonCubicBilinearAlgorithm base = NonCubicBilinearAlgorithm.naive(2, 3, 3);
		Recombination.SotaResolver table = (a, b, c) -> {
			if (a == 1) return 1;   // R(1,3,3)-style thin leaves are cheap
			if (a == 2) return 10;  // balanced leaves deliberately expensive
			return 5;               // R(3,3,3)-style
		};
		AllocationOptimizer.Result res = AllocationOptimizer.optimize(base, table, 4, 9, 9,
				new SearchBudget(Long.MAX_VALUE, 1, 1), null);
		assertThat(res.rank()).as("found by seeding only — the sweep had a 1-node budget")
				.isEqualTo(54);
		int[] a = res.allocA().clone();
		Arrays.sort(a);
		assertThat(a).containsExactly(1, 3);
		assertThat(res.allocB()).containsExactly(3, 3, 3);
		assertThat(res.allocC()).containsExactly(3, 3, 3);
	}
}
