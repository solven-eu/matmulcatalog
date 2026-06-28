package eu.solven.matmul.recombination;

import static org.assertj.core.api.Assertions.assertThat;

import org.junit.jupiter.api.Test;

import eu.solven.matmul.NonCubicBilinearAlgorithm;
import eu.solven.matmul.Verifier;
import eu.solven.matmul.catalog.Compose;
import eu.solven.matmul.catalog.FieldAwareLookup;

/**
 * Confirms the engine can REACH the three shapes we had to manually import from master —
 * each via a FAST targeted construction (the right base + allocation), not the slow default
 * sweep. The default full-pool search misses them (the allocation space for high-part bases
 * is huge), but with the proper base+alloc the (now sparse) constructWithAllocation /
 * concat reaches the SOTA fast and passes a random matmul spot-check. (Exactness of the
 * construction operators themselves is guarded separately by TestRecombinationSparse /
 * TestComposeSparse.) Allocations are master's lineage allocations for these shapes.
 */
public class TestManualImportsReachable {

	private final FieldAwareLookup lookup = new FieldAwareLookup("Q");

	private NonCubicBilinearAlgorithm base(int n, int m, int p) {
		return lookup.find(n, m, p).orElseThrow(() -> new AssertionError("no base ⟨" + n + "," + m + "," + p + "⟩"));
	}

	@Test
	public void reach_7x14x28_1769_via_2x3x3_base() {
		// master lineage: Recombination(2x3x3, A=[4,3] B=[5,5,4] C=[10,9,9]).
		NonCubicBilinearAlgorithm s = Recombination.constructWithAllocation(
				base(2, 3, 3), lookup, new int[] { 4, 3 }, new int[] { 5, 5, 4 }, new int[] { 10, 9, 9 });
		assertThat(s.n).isEqualTo(7);
		assertThat(s.m).isEqualTo(14);
		assertThat(s.p).isEqualTo(28);
		assertThat(s.r).isLessThanOrEqualTo(1769);
		assertThat(Verifier.passesRandomMatmulSpotCheck(s)).isTrue();
	}

	@Test
	public void reach_14x15x31_3839_via_5x5x5_base() {
		// master lineage: Recombination(AE<5,5,5>, A=[3,2,3,3,3] B=[3,3,3,3,3] C=[6,6,6,7,6]).
		NonCubicBilinearAlgorithm s = Recombination.constructWithAllocation(
				base(5, 5, 5), lookup,
				new int[] { 3, 2, 3, 3, 3 }, new int[] { 3, 3, 3, 3, 3 }, new int[] { 6, 6, 6, 7, 6 });
		assertThat(s.n).isEqualTo(14);
		assertThat(s.m).isEqualTo(15);
		assertThat(s.p).isEqualTo(31);
		assertThat(s.r).isLessThanOrEqualTo(3839);
		assertThat(Verifier.passesRandomMatmulSpotCheck(s)).isTrue();
	}

	@Test
	public void reach_7x7x32_1069_via_concat_of_strassen_recomb() {
		// master lineage: ConcatCols(⟨7,7,4⟩, Recombination(Strassen<2,2,2>, [3,4]/[3,4]/[14,14])).
		NonCubicBilinearAlgorithm right = Recombination.constructWithAllocation(
				base(2, 2, 2), lookup, new int[] { 3, 4 }, new int[] { 3, 4 }, new int[] { 14, 14 });
		assertThat(right.p).isEqualTo(28);
		NonCubicBilinearAlgorithm s = Compose.concatRight(base(7, 7, 4), right);
		assertThat(s.n).isEqualTo(7);
		assertThat(s.m).isEqualTo(7);
		assertThat(s.p).isEqualTo(32);
		assertThat(s.r).isLessThanOrEqualTo(1069);
		assertThat(Verifier.passesRandomMatmulSpotCheck(s)).isTrue();
	}
}
