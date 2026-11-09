package eu.solven.matmul.docs.verify;

import static org.assertj.core.api.Assertions.assertThat;

import java.util.List;

import org.junit.jupiter.api.Test;

import tools.jackson.databind.node.ObjectNode;

/**
 * The FMM-Lille index parse ({@link SyncReferenceCatalogs#parseFmmCatalog}) against a fixture
 * cut from the live table (2026-11-09). Two silent drifts hid here: the recipe cell's
 * attribute went {@code content="description"} → {@code data-col="description"} and every
 * digest row's {@code description} became null (so a gap could not be told closure from
 * base), and the naive-rank regex matched the {@code #view} (fast) cell instead of
 * {@code #notview} (classical), so {@code naive_rank == rank} on every row. The rank index
 * kept parsing throughout, so nothing failed.
 */
public class TestSyncReferenceCatalogs {

	/** A derived row: the recipe cell is a projection of another format. */
	private static final String ROW_9x17x17 = """
			<tr><td data-col="format"><a href="9x17x17.html">&lang;9&times;17&times;17:1574&rang;</a>
			</td>
			<td data-col="fast_tensor_rank">
			<a href="#view9x9x9">1574</a>
			</td>
			<td data-col="classical_tensor_rank">
			<a href="#notview9x9x9">2601</a>
			</td>
			<td data-col="description">
			<a href="projalongmode.html">ProjOfAlongMode</a>(<a href="#9x17x18">&lang;9&times;17&times;18:1591&rang;</a>,2)</td>
			<td data-col="ratio">.6051518647</td>
			<td data-col="complexity">2.808380680</td>
			<td data-col="stabilizer">C<sub>1</sub>
			</td>
			</tr>""";

	/** An imported row: the recipe cell is a citation. */
	private static final String ROW_2x2x2 = """
			<tr><td data-col="format"><a href="2x2x2.html">&lang;2&times;2&times;2:7&rang;</a></td>
			<td data-col="fast_tensor_rank"><a href="#view2x2x2">7</a></td>
			<td data-col="classical_tensor_rank"><a href="#notview2x2x2">8</a></td>
			<td data-col="description"><a href="http://fmm.univ-lille.fr/biblio.html#strassen:1969">Strassen 1969</a> [<a href=http://dx.doi.org/10.1007/BF02165411>DOI</a>]
			</td>
			<td data-col="ratio">.8750000000</td>
			</tr>""";

	/** The pre-2026 markup (attribute {@code content=}) with a Kronecker recipe. */
	private static final String ROW_LEGACY = "<a href=\"2x3x6.html\">&lang;2&times;3&times;6:30&rang;</a>"
			+ "<td content=\"description\">&lang;1&times;1&times;2:2&rang; &otimes; &lang;2&times;3&times;3:15&rang;</td>";

	@Test
	public void index_rows_carry_rank_naive_rank_references_and_the_recipe_cell() {
		List<ObjectNode> rows = SyncReferenceCatalogs.parseFmmCatalog("<table>" + ROW_9x17x17 + ROW_2x2x2 + "</table>");
		assertThat(rows).hasSize(2);

		ObjectNode proj = rows.get(0);
		assertThat(proj.get("format").toString()).isEqualTo("[9,17,17]");
		assertThat(proj.get("rank").asInt()).isEqualTo(1574);
		assertThat(proj.get("naive_rank").asInt()).as("classical rank, not the fast one").isEqualTo(2601);
		assertThat(proj.get("description").asText()).isEqualTo("ProjOfAlongMode(⟨9×17×18:1591⟩,2)");
		assertThat(proj.get("references").size()).isZero();
		assertThat(proj.get("details_url").asText()).endsWith("/9x17x17.html");

		ObjectNode strassen = rows.get(1);
		assertThat(strassen.get("rank").asInt()).isEqualTo(7);
		assertThat(strassen.get("naive_rank").asInt()).isEqualTo(8);
		assertThat(strassen.get("description").asText()).isEqualTo("Strassen 1969 [DOI]");
		assertThat(strassen.get("references").get(0).asText()).isEqualTo("strassen:1969");
	}

	@Test
	public void legacy_content_attribute_and_kronecker_entity_still_parse() {
		List<ObjectNode> rows = SyncReferenceCatalogs.parseFmmCatalog(ROW_LEGACY + "</table>");
		assertThat(rows).hasSize(1);
		assertThat(rows.get(0).get("rank").asInt()).isEqualTo(30);
		assertThat(rows.get(0).get("description").asText()).isEqualTo("⟨1×1×2:2⟩ ⊗ ⟨2×3×3:15⟩");
	}

	@Test
	public void a_row_without_a_recipe_cell_reads_as_null_not_as_a_crash() {
		String bare = "<a href=\"4x4x4.html\">&lang;4&times;4&times;4:48&rang;</a><td>x</td></table>";
		List<ObjectNode> rows = SyncReferenceCatalogs.parseFmmCatalog(bare);
		assertThat(rows).hasSize(1);
		assertThat(rows.get(0).get("description").isNull()).isTrue();
		assertThat(rows.get(0).get("naive_rank").isNull()).isTrue();
	}
}
