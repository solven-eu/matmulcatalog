package eu.solven.matmul.search.flip;

import static org.assertj.core.api.Assertions.assertThat;

import java.util.Map;

import org.junit.jupiter.api.BeforeAll;
import org.junit.jupiter.api.Test;
import org.junit.jupiter.params.ParameterizedTest;
import org.junit.jupiter.params.provider.EnumSource;

import eu.solven.matmul.NonCubicBilinearAlgorithm;
import eu.solven.matmul.algebra.Field;
import eu.solven.matmul.catalog.FieldAwareLookup;

/**
 * Meta-move soundness (Kauers–Wood meta flip graph): extend and project must
 * land on the EXACT matmul tensor of the neighbouring format — checked in
 * exact integer arithmetic — and the cross-format walk must only ever record
 * exact vertices.
 */
public class TestMetaFlipMoves {

	private static NonCubicBilinearAlgorithm strassen;

	@BeforeAll
	static void setUp() {
		strassen = new FieldAwareLookup(Field.Z).find(2, 2, 2).orElseThrow();
	}

	@ParameterizedTest(name = "extend {0}")
	@EnumSource(FlipScheme.Axis.class)
	public void extend_is_exact_on_grown_format(FlipScheme.Axis axis) {
		FlipScheme s = FlipScheme.of(strassen);
		FlipScheme grown = s.extendAxis(axis);
		int expectedAdded = switch (axis) {
			case N -> s.m * s.p;
			case M -> s.n * s.p;
			case P -> s.n * s.m;
		};
		assertThat(grown.n + grown.m + grown.p).isEqualTo(s.n + s.m + s.p + 1);
		assertThat(grown.rank()).isEqualTo(s.rank() + expectedAdded);
		assertThat(grown.isExactIntTensor())
				.as("extend(%s) must compute the grown format exactly", axis).isTrue();
		// The original is untouched (meta-moves return new schemes).
		assertThat(s.rank()).isEqualTo(strassen.r);
		assertThat(s.isExactIntTensor()).isTrue();
	}

	@ParameterizedTest(name = "project {0}")
	@EnumSource(FlipScheme.Axis.class)
	public void project_is_exact_on_shrunk_format(FlipScheme.Axis axis) {
		FlipScheme s = FlipScheme.of(NonCubicBilinearAlgorithm.naive(3, 3, 3));
		FlipScheme shrunk = s.projectAxis(axis, 1);
		assertThat(shrunk.n + shrunk.m + shrunk.p).isEqualTo(s.n + s.m + s.p - 1);
		// Naive products are single-index on every axis: dropping one index DCEs
		// exactly one 9-product slice.
		assertThat(shrunk.rank()).isEqualTo(18);
		assertThat(shrunk.isExactIntTensor())
				.as("project(%s) must compute the shrunk format exactly", axis).isTrue();
	}

	@Test
	public void extend_then_project_round_trips() {
		FlipScheme s = FlipScheme.of(strassen);
		// Grow P (adds the naive column products), then drop exactly that new
		// column: DCE removes exactly the added products — back to rank 7.
		FlipScheme back = s.extendAxis(FlipScheme.Axis.P).projectAxis(FlipScheme.Axis.P, 2);
		assertThat(back.n).isEqualTo(2);
		assertThat(back.m).isEqualTo(2);
		assertThat(back.p).isEqualTo(2);
		assertThat(back.rank()).isEqualTo(s.rank());
		assertThat(back.isExactIntTensor()).isTrue();
	}

	@Test
	public void meta_walk_visits_neighbour_shapes_and_stays_exact() {
		FlipScheme seed = FlipScheme.of(strassen);
		MetaFlipWalk.Config cfg = new MetaFlipWalk.Config(
				30_000, 42, 1, 0.01, 0.01, 0.01, 0.02, 2, 3, 80, true, 0);
		MetaFlipWalk.Result r = MetaFlipWalk.walk(seed, FlipObjectives.minRank(), cfg);
		assertThat(r.extendsTaken()).as("walk must exercise the extend meta-move").isPositive();
		assertThat(r.projectsTaken()).as("walk must exercise the project meta-move").isPositive();
		assertThat(r.bestByShape().size())
				.as("walk must visit neighbour formats, not just the seed's").isGreaterThan(1);
		for (Map.Entry<String, MetaFlipWalk.ShapeBest> e : r.bestByShape().entrySet()) {
			assertThat(e.getValue().scheme().isExactIntTensor())
					.as("per-shape best at %s must be exact", e.getKey()).isTrue();
		}
		// The seed's own format must never regress past its seed rank.
		MetaFlipWalk.ShapeBest home = r.bestByShape().get("2x2x2");
		assertThat(home).isNotNull();
		assertThat(home.scheme().rank()).isLessThanOrEqualTo(7);
	}
}
