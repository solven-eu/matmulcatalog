package eu.solven.matmul.catalog;

import static org.assertj.core.api.Assertions.assertThat;

import java.nio.file.Files;
import java.nio.file.Path;
import java.util.ArrayList;
import java.util.List;
import java.util.stream.Stream;

import org.junit.jupiter.api.Tag;
import org.junit.jupiter.api.Test;

import eu.solven.matmul.NonCubicBilinearAlgorithm;
import tools.jackson.databind.JsonNode;
import tools.jackson.databind.json.JsonMapper;

/**
 * Guards {@link SchemeIO#fieldsContradictedByCoefficients} — the cheap, EXACT
 * necessary-condition gate behind {@code VerifySchemes --check=coefficients} and
 * the {@code [skip ci]} import workflow.
 *
 * <p>Regression for the recurring field-discipline over-claim: a {@code 1/8}
 * coefficient scheme stamped {@code fields=[F2,F3,Z,Q,R,C]} (the materialiser's
 * "integer base ⇒ all fields" rule applied to a scheme that picked up rational
 * coefficients). The narrowest correct field is {@code Q}; {@code Z} and
 * {@code F2} are both impossible, while {@code F3} survives (8 is coprime to 3).</p>
 */
public class TestFieldCoefficientOverclaim {

	private static final JsonMapper MAPPER = JsonMapper.builder().build();

	private static JsonNode node(String fields, String coef) {
		String json = "{\"n\":[2,2,2],\"m\":1,"
				+ "\"u_sparse\":{\"0\":{\"i\":[0],\"c\":[" + coef + "]}},"
				+ "\"v_sparse\":{\"0\":{\"i\":[0],\"c\":[1]}},"
				+ "\"w_sparse\":{\"0\":{\"i\":[0],\"c\":[1]}},"
				+ "\"fields\":" + fields + "}";
		return MAPPER.readTree(json);
	}

	@Test
	public void rational_eighth_contradicts_Z_and_F2_but_not_F3() {
		// The exact shape of the recurring bug.
		List<String> bad = SchemeIO.fieldsContradictedByCoefficients(
				node("[\"F2\",\"F3\",\"Z\",\"Q\",\"R\",\"C\"]", "\"1/8\""));
		assertThat(bad).anySatisfy(s -> assertThat(s).contains("Z over-claim"));
		assertThat(bad).anySatisfy(s -> assertThat(s).contains("F2 over-claim"));
		// 8 is coprime to 3 → F3 is NOT contradicted by the denominator.
		assertThat(bad).noneSatisfy(s -> assertThat(s).contains("F3 over-claim"));
	}

	@Test
	public void all_integer_scheme_claiming_Z_is_clean() {
		assertThat(SchemeIO.fieldsContradictedByCoefficients(
				node("[\"F2\",\"F3\",\"Z\",\"Q\",\"R\",\"C\"]", "-1")))
				.isEmpty();
	}

	@Test
	public void third_contradicts_F3_but_not_F2() {
		// 1/3: denominator divisible by 3 → not representable mod 3; coprime to 2 → fine mod 2.
		List<String> badF3 = SchemeIO.fieldsContradictedByCoefficients(
				node("[\"F3\",\"Q\",\"R\",\"C\"]", "\"1/3\""));
		assertThat(badF3).anySatisfy(s -> assertThat(s).contains("F3 over-claim"));

		List<String> okF2 = SchemeIO.fieldsContradictedByCoefficients(
				node("[\"F2\",\"Q\",\"R\",\"C\"]", "\"1/3\""));
		assertThat(okF2).isEmpty();
	}

	/** A ⟨1,1,1⟩ rank-1 alg carrying coefficient {@code uCoef} (identity irrelevant —
	 *  the field-narrowing predicate inspects only coefficient denominators). */
	private static NonCubicBilinearAlgorithm algWithCoef(double uCoef) {
		return new NonCubicBilinearAlgorithm(1, 1, 1,
				new double[][] { { uCoef } }, new double[][] { { 1.0 } }, new double[][] { { 1.0 } });
	}

	@Test
	public void production_narrows_rational_alg_dropping_Z_and_F2_keeping_F3() {
		// The write-time guard: lineage inferred all six fields from integer leaves,
		// but the composed matrices hold 1/8 → drop Z and F2, keep F3/Q/R/C.
		List<String> inferred = List.of("F2", "F3", "Z", "Q", "R", "C");
		assertThat(SchemeIO.narrowFieldsToCoefficients(algWithCoef(0.125), inferred))
				.containsExactly("F3", "Q", "R", "C");
	}

	@Test
	public void production_leaves_integer_alg_untouched() {
		List<String> inferred = List.of("F2", "F3", "Z", "Q", "R", "C");
		assertThat(SchemeIO.narrowFieldsToCoefficients(algWithCoef(-2.0), inferred))
				.isEqualTo(inferred);
	}

	@Test
	public void rational_not_claiming_the_contradicted_field_is_clean() {
		// A 1/8 scheme that honestly declares only [Q,R,C] has nothing to flag.
		assertThat(SchemeIO.fieldsContradictedByCoefficients(
				node("[\"Q\",\"R\",\"C\"]", "\"1/8\"")))
				.isEmpty();
	}

	/**
	 * The published catalog on disk must carry ZERO coefficient/field over-claims.
	 * Cheap (raw-token scan, NO dense materialization) but iterates every scheme,
	 * so it is tagged {@code catalog-iterating}. This is the surefire mirror of the
	 * {@code VerifySchemes --check=coefficients} CI gate — it would have caught the
	 * 59 {@code 1/8}-tagged-Z/F2 derived schemes.
	 */
	@Test
	@Tag("catalog-iterating")
	public void disk_catalog_has_no_coefficient_field_overclaim() throws Exception {
		Path root = Path.of("src/main/resources/schemes");
		List<String> violations = new ArrayList<>();
		try (Stream<Path> walk = Files.walk(root)) {
			for (Path f : (Iterable<Path>) walk.filter(p -> p.toString().endsWith(".json"))::iterator) {
				JsonNode r = SchemeIO.parseJson(f.toFile());
				for (String why : SchemeIO.fieldsContradictedByCoefficients(r)) {
					violations.add(f.getFileName() + ": " + why);
				}
			}
		}
		assertThat(violations)
				.as("schemes whose coefficient denominators contradict their declared fields[]")
				.isEmpty();
	}
}
