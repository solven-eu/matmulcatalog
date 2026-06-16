package eu.solven.matmul.search.flip;

import static org.assertj.core.api.Assertions.assertThat;

import org.junit.jupiter.api.BeforeAll;
import org.junit.jupiter.api.Test;

import eu.solven.matmul.NonCubicBilinearAlgorithm;
import eu.solven.matmul.Verifier;
import eu.solven.matmul.algebra.Field;
import eu.solven.matmul.catalog.FieldAwareLookup;

/**
 * Walk-level milestones. The classic validation (flip-graph literature): a
 * rank-descent walk seeded from the NAIVE ⟨2,2,2⟩ (r=8) must rediscover a
 * rank-7 scheme — over ternary integers, that is a Strassen-class vertex. The
 * structure-objective walks must never lose exactness, never exceed the seed
 * rank, never break the ternary cap, and never end below the seed's score.
 *
 * <p>All walks are deterministic (fixed RNG seeds); results are bounds, the
 * assertions are ≥/≤ so a genuinely better walk never breaks them.</p>
 */
public class TestFlipGraphWalk {

	private static FieldAwareLookup lookup;

	@BeforeAll
	static void setUp() {
		lookup = new FieldAwareLookup(Field.Z);
	}

	@Test
	public void naive_2x2x2_walk_rediscovers_rank_7() {
		FlipScheme seed = FlipScheme.of(NonCubicBilinearAlgorithm.naive(2, 2, 2));
		// A handful of deterministic seeds: the walk is stochastic, the milestone
		// is that SOME short ternary walk lands on r=7 (Strassen-class vertex).
		for (long rngSeed = 1; rngSeed <= 5; rngSeed++) {
			FlipGraphWalk.Config cfg = new FlipGraphWalk.Config(
					150_000, rngSeed, 1, 0.01, 1_000, 1, 0);
			FlipGraphWalk.Result r = FlipGraphWalk.walk(seed, FlipObjectives.minRank(), cfg);
			assertThat(Verifier.isExactNonCubic(r.best().toAlgorithm()))
					.as("every walk vertex must stay exact").isTrue();
			if (r.best().rank() <= 7) {
				assertThat(r.best().maxAbsCoefficient()).isLessThanOrEqualTo(1);
				return;
			}
		}
		org.assertj.core.api.Assertions.fail(
				"no ternary flip walk from naive ⟨2,2,2⟩ reached rank 7 — flip/reduce moves regressed");
	}

	@Test
	public void laderman_bud_walk_improves_or_holds_bud_score() {
		NonCubicBilinearAlgorithm laderman = lookup.find(3, 3, 3).orElseThrow();
		FlipScheme seed = FlipScheme.of(laderman);
		int seedBud = FlipObjectives.budScore(seed);
		FlipGraphWalk.Config cfg = new FlipGraphWalk.Config(
				20_000, 42, 1, 0.0, 2_000, 0, 0);
		FlipGraphWalk.Result r = FlipGraphWalk.walk(seed, FlipObjectives.maxBudScore(), cfg);
		assertThat(Verifier.isExactNonCubic(r.best().toAlgorithm())).isTrue();
		assertThat(r.best().rank())
				.as("rank-first objective must never trade rank away").isLessThanOrEqualTo(23);
		assertThat(FlipObjectives.budScore(r.best()))
				.as("best-found budScore must be ≥ seed's").isGreaterThanOrEqualTo(seedBud);
		assertThat(r.best().maxAbsCoefficient())
				.as("ternary cap must hold across the walk").isLessThanOrEqualTo(1);
	}

	/**
	 * The tradeoff regime the weighted objective exists for: rank is NOT
	 * lexicographic — at weights (rank·5, bud·1) the walk may pay rank+1
	 * whenever it buys ≥5 bud points. Splits are first-class cost-gated moves
	 * and merges stay cost-gated ({@code fullReduce=false}) so they cannot eat
	 * the doubly-proportional pairs the walk builds.
	 */
	@Test
	public void weighted_walk_can_trade_rank_for_structure() {
		NonCubicBilinearAlgorithm laderman = lookup.find(3, 3, 3).orElseThrow();
		FlipScheme seed = FlipScheme.of(laderman);
		FlipObjective obj = FlipObjectives.weighted(5, 1, 0);
		long seedCost = obj.cost(seed);
		FlipGraphWalk.Config cfg = new FlipGraphWalk.Config(
				30_000, 42, 1, 0.0, 2_000, 3, 0, 0.05, false);
		FlipGraphWalk.Result r = FlipGraphWalk.walk(seed, obj, cfg);
		assertThat(r.best().isExactIntTensor()).isTrue();
		assertThat(r.best().rank())
				.as("rank may rise, but only within the split budget")
				.isLessThanOrEqualTo(23 + 3);
		assertThat(r.bestCost())
				.as("best-found weighted cost must never end above the seed's")
				.isLessThanOrEqualTo(seedCost);
		assertThat(r.best().maxAbsCoefficient()).isLessThanOrEqualTo(1);
	}

	/**
	 * Phase-3a guard: the DIRECT serendipitous-cost objective (the real
	 * currency — predicted product rank under the catalog oracle) must never
	 * end above the seed's prediction. budScore maximization demonstrably
	 * WORSENS this metric (ProbeFlipBudHarvest 2026-06-11), which is exactly
	 * why this objective exists.
	 */
	@Test
	public void direct_serendipitous_walk_never_worsens_predicted_cost() {
		FieldAwareLookup qLookup = new FieldAwareLookup(Field.Q);
		NonCubicBilinearAlgorithm laderman = lookup.find(3, 3, 3).orElseThrow();
		FlipScheme seed = FlipScheme.of(laderman);
		FlipObjective obj = FlipObjectives.serendipitous(qLookup, 2, 2, 2);
		long seedCost = obj.cost(seed);
		FlipGraphWalk.Config cfg = new FlipGraphWalk.Config(
				4_000, 42, 1, 0.0, 1_000, 4, 0, 0.05, false);
		FlipGraphWalk.Result r = FlipGraphWalk.walk(seed, obj, cfg);
		assertThat(r.best().isExactIntTensor()).isTrue();
		assertThat(r.bestCost())
				.as("predicted ⟨6,6,6⟩-product cost must be ≤ the seed's prediction")
				.isLessThanOrEqualTo(seedCost);
	}

	@Test
	public void strassen_margin_walk_improves_or_holds_margin() {
		NonCubicBilinearAlgorithm strassen = lookup.find(2, 2, 2).orElseThrow();
		FlipScheme seed = FlipScheme.of(strassen);
		int seedMargin = FlipObjectives.projectionMargin(seed);
		FlipGraphWalk.Config cfg = new FlipGraphWalk.Config(
				5_000, 7, 1, 0.0, 500, 0, 0);
		FlipGraphWalk.Result r =
				FlipGraphWalk.walk(seed, FlipObjectives.maxProjectionMargin(), cfg);
		assertThat(Verifier.isExactNonCubic(r.best().toAlgorithm())).isTrue();
		assertThat(r.best().rank()).isLessThanOrEqualTo(7);
		assertThat(FlipObjectives.projectionMargin(r.best()))
				.isGreaterThanOrEqualTo(seedMargin);
	}
}
