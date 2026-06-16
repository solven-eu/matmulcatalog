package eu.solven.matmul.search.flip;

import static org.assertj.core.api.Assertions.assertThat;

import java.util.stream.Stream;

import org.junit.jupiter.api.BeforeAll;
import org.junit.jupiter.params.ParameterizedTest;
import org.junit.jupiter.params.provider.CsvSource;

import eu.solven.matmul.NonCubicBilinearAlgorithm;
import eu.solven.matmul.algebra.Field;
import eu.solven.matmul.catalog.FieldAwareLookup;
import eu.solven.matmul.catalog.ProjectionSearch;
import eu.solven.matmul.catalog.SerendipitousBudProduct;

/**
 * The walk's native integer metrics must agree with the catalog reference
 * implementations on real catalog schemes — a drift here would silently make
 * the flip search optimize a DIFFERENT score than the one the catalog ranks,
 * stamps and Pareto-selects by.
 */
public class TestFlipObjectives {

	private static FieldAwareLookup lookup;

	@BeforeAll
	static void setUp() {
		lookup = new FieldAwareLookup(Field.Z);
	}

	@ParameterizedTest(name = "⟨{0},{1},{2}⟩")
	@CsvSource({ "2, 2, 2", "3, 3, 3", "2, 3, 4", "4, 4, 4", "3, 3, 6" })
	public void native_scores_match_catalog_reference(int n, int m, int p) {
		NonCubicBilinearAlgorithm alg = lookup.find(n, m, p).orElseThrow();
		FlipScheme s = FlipScheme.of(alg);

		int expectedBud = Stream.of(SerendipitousBudProduct.independentClassSizes(alg))
				.flatMapToInt(java.util.Arrays::stream)
				.filter(size -> size >= 2)
				.sum();
		assertThat(FlipObjectives.budScore(s))
				.as("native budScore must match independentClassSizes Σ classes≥2")
				.isEqualTo(expectedBud);

		assertThat(FlipObjectives.projectionMargin(s))
				.as("native projection margin must match ProjectionSearch.projectionMargin")
				.isEqualTo(ProjectionSearch.projectionMargin(alg));
	}

	/**
	 * The clarified serendipity metric (Smith 2002 eq. 69 as in the paper):
	 * the self-objective must equal the explicit-inner objective at the
	 * scheme's own shape, and the per-axis savings triple must decompose the
	 * scalar saving exactly ({@code r·R(inner) − r_s = σ_U + σ_V + σ_W} for the
	 * same greedy decomposition).
	 */
	@ParameterizedTest(name = "⟨{0},{1},{2}⟩")
	@CsvSource({ "2, 2, 2", "3, 3, 3", "2, 3, 4" })
	public void self_serendipity_metric_is_consistent(int n, int m, int p) {
		FieldAwareLookup qLookup = new FieldAwareLookup(Field.Q);
		NonCubicBilinearAlgorithm alg = lookup.find(n, m, p).orElseThrow();
		FlipScheme s = FlipScheme.of(alg);
		assertThat(FlipObjectives.selfSerendipitous(qLookup).cost(s))
				.as("self objective == explicit-inner objective at own shape")
				.isEqualTo(FlipObjectives.serendipitous(qLookup, n, m, p).cost(s));
		long[] byAxis = FlipObjectives.serendipitySavingByAxis(s, qLookup, n, m, p);
		long inner = qLookup.findRank(n, m, p);
		long greedyCost = SerendipitousBudProduct.costOf(
				SerendipitousBudProduct.findBuds(alg), qLookup, n, m, p);
		assertThat(byAxis[0] + byAxis[1] + byAxis[2])
				.as("per-axis savings must decompose the greedy scalar saving")
				.isEqualTo((long) s.rank() * inner - greedyCost);
	}

	@ParameterizedTest(name = "⟨{0},{1},{2}⟩")
	@CsvSource({ "2, 2, 2", "3, 3, 3" })
	public void costs_are_rank_lexicographic(int n, int m, int p) {
		NonCubicBilinearAlgorithm alg = lookup.find(n, m, p).orElseThrow();
		FlipScheme s = FlipScheme.of(alg);
		// Any structure score is < SCALE, so the cost stays within the rank's
		// lexicographic band: a rank drop always dominates any structure gain.
		long band = (long) s.rank() * FlipObjectives.SCALE;
		assertThat(FlipObjectives.maxBudScore().cost(s))
				.isGreaterThan(band - FlipObjectives.SCALE).isLessThanOrEqualTo(band);
		assertThat(FlipObjectives.maxProjectionMargin().cost(s))
				.isGreaterThan(band - FlipObjectives.SCALE).isLessThanOrEqualTo(band);
		assertThat(FlipObjectives.minRank().cost(s)).isEqualTo(s.rank());
	}

	/** Deterministic projection currency: naive ⟨2,2,2⟩ (r=8) dropping one n-row
	 *  kills exactly the 4 products supported on it → projected rank 4 =
	 *  R⟨1,2,2⟩, i.e. margin 4 on the n axis. */
	@org.junit.jupiter.api.Test
	public void projected_cost_is_exact_on_naive() {
		NonCubicBilinearAlgorithm naive = NonCubicBilinearAlgorithm.naive(2, 2, 2);
		assertThat(ProjectionSearch.projectedRank(naive, 1, 2, 2, 2)).isEqualTo(4);
		assertThat(ProjectionSearch.axisMargins(naive)).containsExactly(4, 4, 4);
		FlipScheme s = FlipScheme.of(naive);
		assertThat(FlipObjectives.projectedTo(1, 2, 2, 2).cost(s)).isEqualTo(4);
		// Unreachable target → huge cost, not an exception.
		assertThat(FlipObjectives.projectedTo(3, 2, 2, 2).cost(s))
				.isGreaterThan(1_000_000L);
	}

	/** The scalar stamped margin must stay the max of the per-axis triple. */
	@ParameterizedTest(name = "⟨{0},{1},{2}⟩")
	@CsvSource({ "2, 2, 2", "3, 3, 3", "2, 4, 3" })
	public void axis_margins_max_is_scalar_margin(int n, int m, int p) {
		NonCubicBilinearAlgorithm alg = lookup.find(n, m, p).orElseThrow();
		int[] tri = ProjectionSearch.axisMargins(alg);
		assertThat(Math.max(tri[0], Math.max(tri[1], tri[2])))
				.isEqualTo(ProjectionSearch.projectionMargin(alg));
	}
}
