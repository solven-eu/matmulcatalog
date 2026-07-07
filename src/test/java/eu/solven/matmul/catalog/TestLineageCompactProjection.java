package eu.solven.matmul.catalog;

import static org.assertj.core.api.Assertions.assertThat;

import java.util.stream.IntStream;

import org.junit.jupiter.api.Test;

/**
 * Compact-lineage projection rendering (user 2026-07-07): a projection usually
 * drops a FEW coordinates and keeps many, so the compact form shows the dropped
 * complement ({@code ∖[dropped]}) when the child's source shape is derivable and
 * the drop-list is shorter — instead of the unreadable full keep-list
 * ({@code 29x29x29 ↓[0,1,…,27|0,…,28|0,…,28]}). Falls back to the keep-form when
 * the source shape is unknown or keeping is rarer than dropping.
 */
public class TestLineageCompactProjection {

	private static int[] range(int n) {
		return IntStream.range(0, n).toArray();
	}

	@Test
	public void drop_few_renders_set_minus_form() {
		Lineage.Node pr = new Lineage.Project(new Lineage.Atom("29x29x29"),
				range(28), range(29), range(29));
		assertThat(Lineage.prettyCompact(pr)).isEqualTo("29x29x29 ∖[28||]");
	}

	@Test
	public void dropped_middle_indices_are_listed() {
		int[] keepM = { 0, 1, 3, 4 };  // drops m-index 2
		Lineage.Node pr = new Lineage.Project(new Lineage.Atom("3x5x4@abc1234"),
				range(3), keepM, range(4));
		assertThat(Lineage.prettyCompact(pr)).isEqualTo("3x5x4 ∖[|2|]");
	}

	@Test
	public void keep_few_falls_back_to_keep_form() {
		Lineage.Node pr = new Lineage.Project(new Lineage.Atom("9x10x13"),
				new int[] { 0 }, new int[] { 1 }, new int[] { 2 });
		assertThat(Lineage.prettyCompact(pr)).isEqualTo("9x10x13 ↓[0|1|2]");
	}

	@Test
	public void unknown_child_shape_falls_back_to_keep_form() {
		Lineage.Node pr = new Lineage.Project(new Lineage.Atom("TA_lita(n=28)"),
				range(27), range(28), range(28));
		assertThat(Lineage.prettyCompact(pr)).startsWith("TA_lita(n=28) ↓[0,");
	}

	@Test
	public void source_shape_derivable_through_wrappers() {
		// OrientAs carries its shape explicitly; the projection drops one p-column.
		Lineage.Node child = new Lineage.OrientAs(new Lineage.Atom("weird-ref"), 4, 7, 5, null);
		Lineage.Node pr = new Lineage.Project(child, range(4), range(7), range(4));
		assertThat(Lineage.prettyCompact(pr)).endsWith(" ∖[||4]");
	}
}
