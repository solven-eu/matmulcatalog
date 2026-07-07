package eu.solven.matmul.docs;

import static org.assertj.core.api.Assertions.assertThat;

import java.util.List;

import org.junit.jupiter.api.BeforeAll;
import org.junit.jupiter.api.Test;

import eu.solven.matmul.catalog.FieldAwareLookup;
import eu.solven.matmul.catalog.Lineage;
import eu.solven.matmul.recombination.BlockSplitSearch;

/**
 * Guards the 2026-07-07 silent bug: {@code SchemeSweep.userBasePool} ({@code --base})
 * created pool entries WITHOUT an {@code originLineage}, so recombination stubs
 * recorded their outer base as the display label ({@code base<3x4x7>=63}) — a
 * best-at-shape cited-bound ref that re-resolves against future catalogs. No crash;
 * the stubs verified fine but were {@code explicitable:false} and their recorded
 * ranks were hostage to whichever representative is catalog-best at replay time.
 * Every user-base pool entry must pin its base as {@code {shape}@{contentHash}}
 * (optionally wrapped in {@code OrientAs} for non-native orientations).
 */
public class TestSchemeSweepUserBasePool {

	private static FieldAwareLookup lookup;

	@BeforeAll
	static void setUp() {
		lookup = new FieldAwareLookup("Q");
	}

	@Test
	public void user_base_entries_carry_pinned_origin_lineage() {
		List<BlockSplitSearch.NamedBase> pool =
				SchemeSweep.userBasePool(List.of(new int[] { 3, 4, 7 }), lookup);
		assertThat(pool).isNotEmpty();
		for (BlockSplitSearch.NamedBase nb : pool) {
			assertThat(nb.originLineage())
					.as("pool entry %s must carry an origin lineage (unpinned label = cited bound)", nb.label())
					.isNotNull();
			Lineage.Atom atom = leafAtom(nb.originLineage());
			assertThat(atom.ref())
					.as("pool entry %s must pin its base as shape@contentHash", nb.label())
					.matches("\\d+x\\d+x\\d+@[0-9a-f]{7,}");
		}
	}

	@Test
	public void width1_axis_base_pins_as_naive() {
		List<BlockSplitSearch.NamedBase> pool =
				SchemeSweep.userBasePool(List.of(new int[] { 1, 2, 2 }), lookup);
		assertThat(pool).isNotEmpty();
		for (BlockSplitSearch.NamedBase nb : pool) {
			assertThat(nb.originLineage()).isNotNull();
			assertThat(leafAtom(nb.originLineage()).ref()).startsWith("naive-");
		}
	}

	private static Lineage.Atom leafAtom(Lineage.Node node) {
		if (node instanceof Lineage.Atom a) {
			return a;
		}
		if (node instanceof Lineage.OrientAs o) {
			return leafAtom(o.child());
		}
		throw new AssertionError("unexpected origin lineage node: " + node);
	}
}
