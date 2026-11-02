package eu.solven.matmul.additions;

import static org.assertj.core.api.Assertions.assertThat;

import java.io.File;

import org.junit.jupiter.api.Test;

import eu.solven.matmul.NonCubicBilinearAlgorithm;
import eu.solven.matmul.additions.SlpJson.SchemeSlp;
import eu.solven.matmul.catalog.SchemeIO;

import tools.jackson.databind.JsonNode;
import tools.jackson.databind.json.JsonMapper;
import tools.jackson.databind.node.ObjectNode;

/**
 * A stored SLP is a claim about additions; these guards make it a replayed fact.
 * Fixture: Perminov's re-encoding of Strassen–Winograd as a 15-addition reduced
 * program (the legacy {@code _reduced} file on disk, 24 flat → 15 scheduled).
 */
public class TestSlpJson {
	private static final JsonMapper MAPPER = JsonMapper.builder().build();
	private static final File STRASSEN_REDUCED = new File(
			"src/main/resources/schemes/known/section2/2x2x2-r7-perminov_cr15_cn24_ZT_reduced-e498eb7.json");

	@Test
	public void fresh_var_layout_of_the_15_addition_strassen_replays_to_its_matrices() throws Exception {
		JsonNode root = SchemeIO.parseJson(STRASSEN_REDUCED);
		assertThat(SlpJson.isFreshVarLayout(root)).isTrue();
		SchemeSlp slp = SlpJson.fromFreshVars(root);
		NonCubicBilinearAlgorithm alg = SchemeIO.readBilinear(STRASSEN_REDUCED);

		assertThat(slp.reconstructs(alg)).as("the program replays to the expanded U/V/W").isTrue();
		assertThat(slp.additions()).as("the contributor's reduced count is what the program costs")
				.isEqualTo(root.get("complexity").get("reduced").asInt()).isEqualTo(15);
		// Perminov writes the sixth right factor as the lone negated entry −b5: a sign flip
		// (one scalar multiplication), never an addition — the count convention of both sides.
		assertThat(slp.scalarMults()).isEqualTo(1);
	}

	@Test
	public void object_layout_round_trips_and_replays() throws Exception {
		JsonNode root = SchemeIO.parseJson(STRASSEN_REDUCED);
		SchemeSlp slp = SlpJson.fromFreshVars(root);
		ObjectNode json = SlpJson.toJson(MAPPER, slp, SlpJson.OPTIMALITY_BOUND, "test", "Perminov 2023");

		SchemeSlp back = SlpJson.readScheme(json).orElseThrow();
		assertThat(back.additions()).isEqualTo(15);
		assertThat(back.a().ops()).isEqualTo(slp.a().ops());
		assertThat(back.out().formResult()).isEqualTo(slp.out().formResult());
		assertThat(back.reconstructs(SchemeIO.readBilinear(STRASSEN_REDUCED))).isTrue();
		assertThat(json.get("optimality").asString()).isEqualTo("bound");
		assertThat(json.get("source").asString()).isEqualTo("Perminov 2023");
		// Array-form (string) programs and missing sides are not this layout.
		assertThat(SlpJson.readScheme(MAPPER.createArrayNode())).isEmpty();
		assertThat(SlpJson.readScheme(null)).isEmpty();
	}

	@Test
	public void a_flipped_sign_in_the_program_is_caught_by_the_replay() throws Exception {
		ObjectNode root = (ObjectNode) SchemeIO.parseJson(STRASSEN_REDUCED);
		NonCubicBilinearAlgorithm alg = SchemeIO.readBilinear(STRASSEN_REDUCED);
		// Flip one coefficient of the first w_fresh definition.
		ObjectNode term = (ObjectNode) root.get("w_fresh").get(0).get(0);
		term.put("value", -term.get("value").asInt());
		assertThat(SlpJson.fromFreshVars(root).reconstructs(alg)).isFalse();
	}

	/**
	 * Term shapes the Strassen fixture does not exercise: a three-term fresh variable
	 * (two additions), a fresh variable made only of negated terms (one addition, one
	 * free negation), a row that is a single negated reference (one negation, no
	 * addition), and column-major {@code w} rows on a non-square output.
	 */
	@Test
	public void counting_convention_one_addition_per_extra_term_negation_free() throws Exception {
		String json = """
				{"n": [2, 3, 1], "m": 2,
				 "u_fresh": [[{"index": 0, "value": -1}, {"index": 4, "value": -1}]],
				 "u": [[{"index": 6, "value": 1}], [{"index": 1, "value": 1}, {"index": 5, "value": 1}]],
				 "v_fresh": [[{"index": 0, "value": 1}, {"index": 1, "value": -1}, {"index": 2, "value": 1}]],
				 "v": [[{"index": 3, "value": 1}], [{"index": 2, "value": -1}]],
				 "w_fresh": [],
				 "w": [[{"index": 0, "value": 1}, {"index": 1, "value": 1}], [{"index": 1, "value": -1}]]}
				""";
		JsonNode root = MAPPER.readTree(json);
		SchemeSlp slp = SlpJson.fromFreshVars(root);
		assertThat(slp.a().additions()).as("fresh −a0−a4 (1) + row a1+a5 (1)").isEqualTo(2);
		assertThat(slp.a().scalarMults()).as("the all-negative fresh variable is negated once").isEqualTo(1);
		assertThat(slp.b().additions()).as("three-term fresh variable").isEqualTo(2);
		assertThat(slp.b().scalarMults()).as("−b2 is a negation, not an addition").isEqualTo(1);
		assertThat(slp.out().additions()).isEqualTo(1);
		assertThat(slp.additions()).isEqualTo(5);
		// Linear forms agree with SchemeIO's own expansion of the same layout (w rows
		// are column-major C positions: ⟨2,3,1⟩ has a 2×1 output, so the orders coincide
		// only up to the rowMap — both readers must apply it).
		assertThat(slp.reconstructs(SchemeIO.readReduced(root))).isTrue();
	}

	@Test
	public void column_major_w_rows_are_mapped_on_a_wide_output() throws Exception {
		// ⟨1,1,2⟩ = 2 products, output 1×2: Perminov row p = j·n1 + i → with n1 = 1, p = j.
		// ⟨2,1,1⟩ … use ⟨2,1,2⟩ (4 outputs, n1 = 2, n3 = 2): row 1 is (i=1, j=0) → ours 1·2+0 = 2.
		String json = """
				{"n": [2, 1, 2], "m": 4,
				 "u": [[{"index": 0, "value": 1}], [{"index": 1, "value": 1}], [{"index": 0, "value": 1}], [{"index": 1, "value": 1}]],
				 "v": [[{"index": 0, "value": 1}], [{"index": 0, "value": 1}], [{"index": 1, "value": 1}], [{"index": 1, "value": 1}]],
				 "w": [[{"index": 0, "value": 1}], [{"index": 1, "value": 1}], [{"index": 2, "value": 1}], [{"index": 3, "value": 1}]]}
				""";
		JsonNode root = MAPPER.readTree(json);
		SchemeSlp slp = SlpJson.fromFreshVars(root);
		NonCubicBilinearAlgorithm alg = SchemeIO.readReduced(root);
		assertThat(slp.reconstructs(alg)).isTrue();
		assertThat(slp.out().formResult()).as("form 2 (ours c10) is Perminov's row 1 → product 1").containsExactly(0, 2, 1, 3);
	}
}
