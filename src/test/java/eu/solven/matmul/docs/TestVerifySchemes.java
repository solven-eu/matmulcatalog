package eu.solven.matmul.docs;

import static org.assertj.core.api.Assertions.assertThat;

import java.nio.file.Files;
import java.nio.file.Path;
import java.util.EnumSet;
import java.util.List;
import java.util.Optional;

import org.junit.jupiter.api.BeforeAll;
import org.junit.jupiter.api.Test;
import org.junit.jupiter.api.io.TempDir;

import eu.solven.matmul.NonCubicBilinearAlgorithm;
import eu.solven.matmul.catalog.FieldAwareLookup;
import eu.solven.matmul.catalog.Lineage;
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

	// The fixed ⟨3,3,8⟩=55 scheme (fields [F3,Q,R,C], additions 920) — resolved by
	// shape+token so the 2026 filename migration ('derived_recursive-3x3x8_m55_a920_b0'
	// → '3x3x8-r55-derived-{hash7}') and future re-hashes cannot orphan this test.
	private static final Path DERIVED_338 = eu.solven.matmul.catalog.SchemeResolver
			.byHint("src/main/resources/schemes/derived/section8/derived-3x3x8_m55.json").toPath();

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

	/**
	 * 2026-10-14: two projection stubs recorded a rank 3–4 BELOW what their lineage
	 * replays today (the parent's leaves had moved since the stub was written) — false
	 * rank claims that the daily verification never looked at. The {@code rank} check
	 * must flag a stub whose recorded {@code m} is below its replay, accept an honest
	 * one, and accept a merely loose one (recorded above the replay).
	 */
	@Test
	void rank_check_flags_a_stub_recorded_below_its_replay(@TempDir Path dir) throws Exception {
		// A real Project lineage on the catalog's pinned ⟨2,2,2⟩: keep both A-rows and
		// inner indices, one B-column → ⟨2,2,1⟩, DCE'ing the products that only touch
		// the dropped column (5 of Strassen's 7).
		// A hash-stamped explicit rank-7 ⟨2,2,2⟩ (the lookup's first file may be a legacy
		// record without a `hash`, e.g. the Waksman commutative one).
		Path strassen = lookup.findFiles(2, 2, 2).stream().filter(p -> {
			try {
				JsonNode r = SchemeIO.parseJson(p.toFile());
				return r.has("hash") && r.path("m").asInt() == 7 && !SchemeIO.isStub(r) && !SchemeIO.isNonBilinear(r)
						&& !SchemeIO.isComplex(r);
			} catch (Exception e) {
				return false;
			}
		}).findFirst().orElseThrow();
		String hash = SchemeIO.parseJson(strassen.toFile()).get("hash").asString();
		Lineage.Node lineage = new Lineage.Project(new Lineage.Atom("2x2x2@" + hash), new int[] { 0, 1 },
				new int[] { 0, 1 }, new int[] { 0 });
		NonCubicBilinearAlgorithm honest = VerifySchemes.replayerFor(lookup).replay(lineage);
		assertThat(honest.n).isEqualTo(2);
		assertThat(honest.p).isEqualTo(1);

		Path stub = dir.resolve("2x2x1-stub.json");
		SchemeIO.writeStub(honest, stub.toFile(), lineage, List.of("Z", "Q", "R", "C"));
		assertThat(VerifySchemes.rankViolation(stub, SchemeIO.parseJson(stub.toFile()), lookup))
				.as("an honest stub passes").isNull();

		Path loose = dir.resolve("2x2x1-loose.json");
		Files.writeString(loose, Files.readString(stub).replace("\"m\": " + honest.r + ",", "\"m\": " + (honest.r + 1) + ","));
		assertThat(VerifySchemes.rankViolation(loose, SchemeIO.parseJson(loose.toFile()), lookup))
				.as("a recorded rank ABOVE the replay is loose, not false").isNull();

		Path over = dir.resolve("2x2x1-over.json");
		Files.writeString(over, Files.readString(stub).replace("\"m\": " + honest.r + ",", "\"m\": " + (honest.r - 1) + ","));
		assertThat(VerifySchemes.rankViolation(over, SchemeIO.parseJson(over.toFile()), lookup))
				.as("a recorded rank BELOW the replay is a false claim").contains("false rank claim");
		assertThat(VerifySchemes.checkOne(over, SchemeIO.parseJson(over.toFile()), List.of("Q"),
				EnumSet.of(VerifySchemes.Check.RANK), lookup)).hasSize(1);
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
