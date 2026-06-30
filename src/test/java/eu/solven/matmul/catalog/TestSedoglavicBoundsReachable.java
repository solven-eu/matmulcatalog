package eu.solven.matmul.catalog;

import static org.assertj.core.api.Assertions.assertThat;

import java.util.List;
import java.util.Map;
import java.util.Optional;
import java.util.function.Function;

import org.junit.jupiter.api.Test;

import eu.solven.matmul.recombination.BlockSplitSearch;
import eu.solven.matmul.recombination.Recombination;

/**
 * Regression guard for the removal of {@code SedoglavicProp1Method} (the
 * standalone {@code ConstructiveMethod} enumerator of Sedoglavic 2017 Prop 1).
 *
 * <p>That method was redundant: its bound
 * {@code ⟨u+v,u+v,u+v⟩ ≤ ⟨u,u,u⟩ + 3·⟨u,u,v⟩ + 3·⟨v,v,u⟩} is exactly the
 * recombination multiset of the Strassen {@code ⟨2,2,2⟩} base at a two-part
 * {@code [u,v]} allocation, so the generic recombination path already reaches
 * it (and beats it when a better split exists, e.g. ⟨17,17,17⟩=2930&lt;2940);
 * the {@code u=v} doubling is the {@link eu.solven.matmul.search.PairFusedRecombination}
 * Pan-TA candidate. This test pins that: with the target cubic masked out of the
 * SOTA resolver (so the rank MUST be constructed, not read back), the live
 * {@code findBestStrategy} reaches the published Sedoglavic upper bound OR
 * BETTER — for both the {@code u>v} splits (7,11,17) and the {@code u=v}
 * doubling (14, where Pan-TA fusion strictly beats plain {@code 7·R(⟨7,7,7⟩)}).</p>
 *
 * <p>SOTA-or-better ({@code ≤}) by design: a future improvement must never break
 * this, only a silent regression of the recombination/pair-fused pool would.</p>
 */
public class TestSedoglavicBoundsReachable {

	// Published Sedoglavic Prop 1 cubic bounds (non-commutative, char 0):
	//   (u,v)=(4,3)→⟨7⟩³=250, (6,5)→⟨11⟩³=873, (9,8)→⟨17⟩³=2940,
	//   doubling k=7→⟨14⟩³=1719.
	private static final int[][] SEDOGLAVIC_CUBIC = {
			{ 7, 250 },
			{ 11, 873 },
			{ 14, 1719 }, // u=v=7 doubling — the Pan-TA-vs-plain case
			{ 17, 2940 },
	};

	@Test
	public void findBestStrategy_reaches_sedoglavic_cubic_bounds_without_the_method() {
		List<BlockSplitSearch.NamedBase> pool = BlockSplitSearch.defaultPool();
		Map<String, Integer> ranks = BlockSplitSearch.loadCatalogBestRanksForField("R");
		Function<int[], Optional<Integer>> lookup = BlockSplitSearch.rankLookupFromMap(ranks);

		for (int[] row : SEDOGLAVIC_CUBIC) {
			int n = row[0];
			long bound = row[1];

			// SOTA resolver that returns trivial ranks for degenerate dims and
			// catalog-best otherwise, BUT masks the target ⟨n,n,n⟩ so it cannot be
			// read back directly — forcing a constructive (recombination / pair-fused)
			// derivation, exactly the path the removed method used to shadow.
			int tn = n;
			Recombination.SotaResolver sota = (a, b, c) -> {
				if (a == 0 || b == 0 || c == 0) return 0;
				if (a == 1) return b * c;
				if (b == 1) return a * c;
				if (c == 1) return a * b;
				if (a == tn && b == tn && c == tn) return Integer.MAX_VALUE / 100;
				return lookup.apply(new int[] { a, b, c }).orElse(Integer.MAX_VALUE / 100);
			};

			// balancedOnly=false so the allocation optimizer is free to pick the
			// unbalanced [u,v] split the Sedoglavic identity uses.
			Optional<BlockSplitSearch.NonCubicStrategy> best =
					BlockSplitSearch.findBestStrategy(n, n, n, pool, sota, false);

			assertThat(best)
					.as("⟨%d,%d,%d⟩ must be constructible without SedoglavicProp1Method", n, n, n)
					.isPresent();
			assertThat(best.get().rank())
					.as("⟨%d,%d,%d⟩ reached %s (%s); Sedoglavic published bound is %d — expected SOTA-or-better",
							n, n, n, best.get().rank(), best.get().label(), bound)
					.isLessThanOrEqualTo(bound);
		}
	}
}
