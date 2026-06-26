package eu.solven.matmul.catalog;

import static org.assertj.core.api.Assertions.assertThat;

import java.util.List;

import org.junit.jupiter.api.Test;

/**
 * Guards {@link BaseUsageStats} — the per-base "used as a building block N times by
 * recombination / projection / kronecker / …" accounting stamped into the catalog.
 */
public class TestBaseUsageStats {

	@Test
	public void baseKey_normalises_pinned_bare_and_skips_primitives() {
		assertThat(BaseUsageStats.baseKey("7x7x15@89c37dfaaaaaaaaaaaaaaaaaaaaaaaaa"))
				.isEqualTo("7x7x15@89c37df");
		assertThat(BaseUsageStats.baseKey("3x4x8-r73-perminov_Z-89c37df"))
				.isEqualTo("3x4x8@89c37df");
		assertThat(BaseUsageStats.baseKey("7x7x15")).isEqualTo("7x7x15");
		assertThat(BaseUsageStats.baseKey("naive-1x1x1")).isNull();
		assertThat(BaseUsageStats.baseKey("@ref?:L0")).isNull();
	}

	@Test
	public void recombination_attributes_base_and_leaves() {
		// RecombinationN(base=2x2x2@aaaaaaa, leaves=[7x7x15@bbbbbbb, naive-1x1x1])
		Lineage.Node lin = new Lineage.RecombinationN(
				new Lineage.Atom("2x2x2@aaaaaaa1111111"),
				new int[] { 1 }, new int[] { 1 }, new int[] { 1 },
				List.of(new Lineage.Atom("7x7x15@bbbbbbb2222222"),
						new Lineage.Atom("naive-1x1x1")));

		BaseUsageStats stats = new BaseUsageStats();
		stats.accumulate(lin);

		assertThat(stats.forKey("2x2x2@aaaaaaa")).containsEntry("recombination", 1);
		assertThat(stats.forKey("7x7x15@bbbbbbb")).containsEntry("recombination", 1);
		// the naïve primitive is not a catalog base → not counted
		assertThat(stats.all()).doesNotContainKey("naive-1x1x1");
	}

	@Test
	public void transform_wrappers_pass_the_using_op_through() {
		// Project(OrientAs(Atom 9x9x9)) → 9x9x9 is used BY projection (OrientAs is a
		// transparent wrapper, not a use of its own).
		Lineage.Node lin = new Lineage.Project(
				new Lineage.OrientAs(new Lineage.Atom("9x9x9@ccccccc3333333"),
						9, 9, 9, new int[] { 0, 1, 2 }),
				new int[] { 0, 1, 2 }, new int[] { 0, 1, 2 }, new int[] { 0, 1, 2 });

		BaseUsageStats stats = new BaseUsageStats();
		stats.accumulate(lin);

		assertThat(stats.forKey("9x9x9@ccccccc")).containsEntry("projection", 1);
	}

	@Test
	public void aggregates_across_schemes_by_op() {
		BaseUsageStats stats = new BaseUsageStats();
		// base 2x2x2@aaaaaaa used once by recombination, once by kronecker
		stats.accumulate(new Lineage.RecombinationN(
				new Lineage.Atom("2x2x2@aaaaaaa1"), new int[] { 1 }, new int[] { 1 },
				new int[] { 1 }, List.of()));
		stats.accumulate(new Lineage.KronProduct(
				new Lineage.Atom("2x2x2@aaaaaaa1"), new Lineage.Atom("3x3x3@ddddddd4")));

		assertThat(stats.forKey("2x2x2@aaaaaaa"))
				.containsEntry("recombination", 1)
				.containsEntry("kronecker", 1);
		assertThat(stats.forKey("3x3x3@ddddddd")).containsEntry("kronecker", 1);
	}
}
