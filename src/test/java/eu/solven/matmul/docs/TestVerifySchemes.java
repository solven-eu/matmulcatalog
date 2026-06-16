package eu.solven.matmul.docs;

import static org.assertj.core.api.Assertions.assertThat;

import java.nio.file.Path;
import java.util.EnumSet;
import java.util.List;
import java.util.Optional;

import org.junit.jupiter.api.BeforeAll;
import org.junit.jupiter.api.Test;

import eu.solven.matmul.catalog.FieldAwareLookup;
import eu.solven.matmul.catalog.SchemeIO;
import eu.solven.matmul.docs.verify.VerifyLineageFieldCompat;

import tools.jackson.databind.JsonNode;
import tools.jackson.databind.node.ArrayNode;
import tools.jackson.databind.node.ObjectNode;

/**
 * Regression guard for the derived ⟨3,3,8⟩=55 field over-claim: the scheme was
 * reported over <b>Z</b> even though it concatenates a Q-only {@code smirnov13}
 * ⟨3,3,6⟩ leaf (which carries {@code 1/8} coefficients), so it lives in Q, not Z.
 *
 * <p>The bug had two faces, both guarded here: (1) the on-disk {@code fields[]}
 * must not name a char-0 field narrower than the lineage supports, and (2) the
 * {@code VerifyLineageFieldCompat} check that backstops it must actually have
 * teeth — a fabricated Z claim over the same lineage must be flagged.</p>
 */
public class TestVerifySchemes {

	private static FieldAwareLookup lookup;

	private static final Path DERIVED_338 = Path.of(
			"src/main/resources/schemes/derived/section8/derived_recursive-3x3x8_m55_a920_b0.json");

	@BeforeAll
	static void setup() {
		lookup = new FieldAwareLookup("Q");
	}

	@Test
	void derived_3x3x8_does_not_overclaim_Z() throws Exception {
		JsonNode root = SchemeIO.parseJson(DERIVED_338.toFile());
		List<String> declared = SchemeIO.fieldTags(root);
		// It depends on the Q-only smirnov13 ⟨3,3,6⟩ leaf and holds 1/8 coefficients:
		// its true field set is [F3,Q,R,C] — emphatically NOT Z.
		assertThat(declared).doesNotContain("Z");
		assertThat(declared).contains("Q", "R", "C");
	}

	@Test
	void lineage_check_passes_for_the_fixed_scheme() throws Exception {
		JsonNode root = SchemeIO.parseJson(DERIVED_338.toFile());
		Optional<VerifyLineageFieldCompat.Violation> v =
				VerifyLineageFieldCompat.check(DERIVED_338, root, lookup);
		assertThat(v).as("fixed ⟨3,3,8⟩ must not trip the char-0 floor check").isEmpty();
	}

	@Test
	void full_check_pipeline_passes_for_the_fixed_scheme() throws Exception {
		JsonNode root = SchemeIO.parseJson(DERIVED_338.toFile());
		List<String> declared = SchemeIO.fieldTags(root);
		// Run the real three checks (lineage + coefficients + full identity).
		List<String> failures = VerifySchemes.checkOne(DERIVED_338, root, declared,
				EnumSet.allOf(VerifySchemes.Check.class), lookup);
		assertThat(failures).as("all VerifySchemes checks must pass on the fixed ⟨3,3,8⟩").isEmpty();
	}

	@Test
	void lineage_check_has_teeth_against_a_fabricated_Z_claim() throws Exception {
		// Take the real (Q-floor) lineage but FABRICATE a Z over-claim in fields[].
		ObjectNode root = (ObjectNode) SchemeIO.parseJson(DERIVED_338.toFile());
		ArrayNode fields = root.arrayNode();
		fields.add("Z");
		fields.add("Q");
		fields.add("R");
		fields.add("C");
		root.set("fields", fields);

		Optional<VerifyLineageFieldCompat.Violation> v =
				VerifyLineageFieldCompat.check(DERIVED_338, root, lookup);
		assertThat(v).as("a fabricated Z over a Q-floor lineage MUST be flagged").isPresent();
		assertThat(v.get().declaredNarrowest()).isEqualTo("Z");
		assertThat(v.get().floor()).isEqualTo("Q");
	}
}
