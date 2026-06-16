package eu.solven.matmul.search;

import static org.assertj.core.api.Assertions.assertThat;

import java.util.List;
import java.util.Optional;

import org.junit.jupiter.api.BeforeAll;
import org.junit.jupiter.api.Test;
import org.junit.jupiter.params.ParameterizedTest;
import org.junit.jupiter.params.provider.CsvSource;

import eu.solven.matmul.Verifier;
import eu.solven.matmul.catalog.FieldAwareLookup;

/**
 * Fast, parameter-specific <em>regression guards</em> on the search/materialise
 * pipeline: with the real catalog, the recursive materialiser must still spot the
 * known SOTA (or better) for a spread of shapes that each exercise a distinct
 * mechanism. This is the cheap counterpart to the full {@code SchemeSweep} /
 * {@code VerifyAllSchemes} runs — it builds the lookup once and probes a handful
 * of shapes, so it runs in seconds and fails loudly when an engine silently
 * regresses (e.g. the 2026-06-10 empty-{@code extendedPool} bug, or the
 * bud-ordering bug that hid ⟨8,9,9⟩=430).
 *
 * <p>Each bound is an <em>upper</em> bound (≤): the pipeline may legitimately find
 * something better, but never worse. Keep the shape list small + low-dim so the
 * suite stays fast.</p>
 */
public class TestSweepSpotsSota {

	private static FieldAwareLookup lookup;
	private static RecursiveMaterialiser mat;

	@BeforeAll
	static void setUp() {
		lookup = new FieldAwareLookup("Q");
		List<BlockSplitSearch.NamedBase> pool = BlockSplitSearch.defaultPool();
		RecursiveClosureSota sota = new RecursiveClosureSota(lookup, pool, true, true);
		// dry-run (no disk writes); composes from the real catalog.
		mat = new RecursiveMaterialiser(lookup, pool, sota, null, false, true);
	}

	@ParameterizedTest(name = "⟨{0},{1},{2}⟩ ≤ {3}")
	@CsvSource({
			// shape         SOTA   mechanism exercised
			"2, 2, 2,    7",   // Strassen — direct disk hit
			"3, 3, 3,   23",   // Laderman — direct disk hit
			"4, 4, 4,   49",   // Kronecker ⟨2,2,2⟩⊗⟨2,2,2⟩
			"3, 7, 8,  126",   // concat ⟨3,7,4⟩ +ₚ ⟨3,7,4⟩ (orientation-aware)
			"7, 7, 7,  250",   // Strassen-recursion recombination
			"6, 8, 9,  296",   // serendipitous bud-product (bud-ordering fix)
			"8, 9, 9,  430",   // serendipitous bud-product (size-3 V-bud)
			"4, 8, 12, 272",   // serendipitous
			"8, 8, 12, 504",   // serendipitous
	})
	public void materialise_spots_sota(int n, int m, int p, int sota) {
		Optional<RecursiveMaterialiser.Result> r = mat.materialise(n, m, p);
		assertThat(r).as("⟨%d,%d,%d⟩ should resolve", n, m, p).isPresent();
		assertThat(r.get().alg().r)
				.as("⟨%d,%d,%d⟩ rank must be ≤ SOTA %d (regression if higher)", n, m, p, sota)
				.isLessThanOrEqualTo(sota);
		assertThat(Verifier.passesRandomMatmulSpotCheck(r.get().alg()))
				.as("⟨%d,%d,%d⟩ result must verify", n, m, p).isTrue();
	}

	/**
	 * ⟨17,17,17⟩ has been a long-contested, hard shape: the plain search only
	 * reaches the <b>2940</b> floor, while the catalog holds sub-2940 results
	 * (FMM 2934, LRP/derived 2930) as <b>maxDim&gt;16 lineage-only stubs</b>.
	 * {@code materialise()} deliberately skips those stubs (it returns 2940 here),
	 * so we guard via {@code findRank}, which is stub-inclusive — checking a stub
	 * in a unit test is fine (per the user). The bound 2934 sits below the 2940
	 * search floor, so it fails loudly if the sub-2940 import/derivation is ever
	 * lost (e.g. a folder reorg or over-eager cleanup dropping the stub).
	 */
	@Test
	public void retains_hard_won_17x17x17_below_search_floor() {
		assertThat(lookup.findRank(17, 17, 17))
				.as("⟨17,17,17⟩ must retain a sub-2940 result (2930/2934-class stub); "
						+ "2940 would mean the hard-won import/derivation was lost")
				.isLessThanOrEqualTo(2934);
	}

	/**
	 * The extended template pool must see the whole catalog tree. A
	 * {@code listFiles("section*")} on the schemes root (pre-2026-06-10 bug)
	 * silently returned an empty pool after the known/derived/curated split,
	 * crippling the search; this guards against that regression.
	 */
	@Test
	public void extended_pool_is_not_empty() {
		assertThat(BlockSplitSearch.extendedPool(8))
				.as("extendedPool(8) must load catalog leaves from known/derived/curated, not be seed-only")
				.hasSizeGreaterThan(50);
	}

	/**
	 * ⟨5,32,32⟩ = 3320 (= FMM-Lille) must be reachable through the SchemeSweep
	 * evaluate path — {@code buildPool(includeDerived)} + {@code findBestStrategy}
	 * — via the HK ⟨2,4,4⟩=26 recombination: allocA=[3,2] (n: 5=3+2),
	 * allocB=allocC=[8,8,8,8] (each 32=4·8) → 16×⟨3,8,8⟩=145 + 10×⟨2,8,8⟩=100.
	 *
	 * <p>The committed catalog held 3446 because the DEFAULT {@code rootPool}
	 * omits ⟨2,4,4⟩ as an outer base; the ⟨2,4,4⟩ base lives only in the
	 * derived-inclusive (extended) pool. This guards the mechanism: a regression
	 * that drops ⟨2,4,4⟩ from the extended pool, breaks 4-way ({@code [8,8,8,8]})
	 * allocations, or loses the ⟨3,8,8⟩/⟨2,8,8⟩ leaves would push this back to
	 * 3446 and fail. SOTA-or-better (≤), so a future improvement never breaks it.</p>
	 */
	@Test
	public void includeDerived_sweep_finds_5x32x32_3320_via_2x4x4() {
		List<BlockSplitSearch.NamedBase> pool = BlockSplitSearch.buildPool(PoolConfig.includeDerived());
		CitedBound sota = new CitedBound(lookup);
		// bound just above 3320 so the recombination B&B prunes hard and the test
		// stays fast, while still letting the 3320 route through.
		Optional<BlockSplitSearch.NonCubicStrategy> best = BlockSplitSearch.findBestStrategy(
				5, 32, 32, pool, sota, false,
				PoolConfig.UNBOUNDED_IMBALANCE, PoolConfig.UNBOUNDED_COMBINATIONS, 0, 3446L);
		assertThat(best).as("⟨5,32,32⟩ must resolve via the includeDerived pool").isPresent();
		assertThat(best.get().rank())
				.as("⟨5,32,32⟩ must reach FMM's 3320 or better (regression → 3446 = ⟨2,4,4⟩ base lost)")
				.isLessThanOrEqualTo(3320L);
		assertThat(best.get().recombination())
				.as("the 3320 route is a recombination, not concat/kronecker").isNotNull();
		assertThat(best.get().label())
				.as("the winning outer base must be ⟨2,4,4⟩").contains("2x4x4");
	}
}
