package eu.solven.matmul.catalog;

import static org.assertj.core.api.Assertions.assertThat;

import org.junit.jupiter.api.Test;

/**
 * Key order is part of the canonical format: two writers that build the same lineage
 * node with keys inserted in different orders must produce the same bytes. Guards the
 * {@code OrientAs} {@code axisMap}/{@code child} swap that made every stub restamp in
 * PR #26 show a diff with no content change.
 */
public class TestMatrixJsonFormatterKeyOrder {

	private static final String CHILD_FIRST = """
			{"lineage": {"op": "OrientAs", "n": 12, "m": 8, "p": 10,
			  "child": {"op": "Atom", "ref": "8x10x12@c5b7abd"}, "axisMap": "pnm"}, "m": 624}""";
	private static final String AXISMAP_FIRST = """
			{"lineage": {"op": "OrientAs", "n": 12, "m": 8, "p": 10,
			  "axisMap": "pnm", "child": {"op": "Atom", "ref": "8x10x12@c5b7abd"}}, "m": 624}""";

	@Test
	public void orientAs_key_order_does_not_depend_on_the_writer() throws Exception {
		assertThat(MatrixJsonFormatter.format(CHILD_FIRST)).isEqualTo(MatrixJsonFormatter.format(AXISMAP_FIRST));
	}

	@Test
	public void lineage_node_puts_op_first_and_child_last() throws Exception {
		String out = MatrixJsonFormatter.format("""
				{"child": {"ref": "2x2x2", "op": "Atom"}, "perm": [1, 0, 2], "op": "Transpose"}""");
		assertThat(out).isEqualTo("""
				{
				  "op": "Transpose",
				  "perm": [1, 0, 2],
				  "child": {"op": "Atom", "ref": "2x2x2"}
				}
				""");
	}

	@Test
	public void non_lineage_objects_keep_insertion_order() throws Exception {
		// Top-level scheme keys (no "op") are not reordered: the policy is lineage-only.
		String out = MatrixJsonFormatter.format("{\"m\": 7, \"child\": [1, 2], \"n\": [2, 2, 2]}");
		assertThat(out).isEqualTo("""
				{
				  "m": 7,
				  "child": [1, 2],
				  "n": [2, 2, 2]
				}
				""");
	}

	@Test
	public void format_is_idempotent_after_reordering() throws Exception {
		String once = MatrixJsonFormatter.format(CHILD_FIRST);
		assertThat(MatrixJsonFormatter.format(once)).isEqualTo(once);
	}
}
