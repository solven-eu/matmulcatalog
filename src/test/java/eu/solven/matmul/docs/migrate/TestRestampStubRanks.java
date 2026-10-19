package eu.solven.matmul.docs.migrate;

import static org.assertj.core.api.Assertions.assertThat;

import java.nio.file.Files;
import java.nio.file.Path;
import java.util.List;
import java.util.Optional;

import org.junit.jupiter.api.BeforeAll;
import org.junit.jupiter.api.Test;
import org.junit.jupiter.api.io.TempDir;

import eu.solven.matmul.NonCubicBilinearAlgorithm;
import eu.solven.matmul.catalog.FieldAwareLookup;
import eu.solven.matmul.catalog.Lineage;
import eu.solven.matmul.catalog.SchemeIO;
import eu.solven.matmul.search.LineageReplayer;
import tools.jackson.databind.JsonNode;

/**
 * {@link RestampStubRanks} repairs a stub whose recorded rank drifted from its replay:
 * the rewritten stub records the replayed rank, keeps the lineage, fields and
 * provenance keys, is renamed to the new rank/hash, and the stale file is gone.
 */
public class TestRestampStubRanks {
	private static FieldAwareLookup lookup;
	private static LineageReplayer replayer;

	@BeforeAll
	static void setup() {
		lookup = new FieldAwareLookup("Q");
		replayer = LineageReplayer.withDefaultPool(lookup);
	}

	private static Lineage.Node projectionOfStrassen() throws Exception {
		// A hash-stamped explicit rank-7 ⟨2,2,2⟩ (the lookup's first file may be a legacy
		// record without a `hash`).
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
		return new Lineage.Project(new Lineage.Atom("2x2x2@" + hash), new int[] { 0, 1 }, new int[] { 0, 1 },
				new int[] { 0 });
	}

	@Test
	public void an_over_claiming_stub_is_rewritten_at_its_replayed_rank(@TempDir Path dir) throws Exception {
		Lineage.Node lineage = projectionOfStrassen();
		NonCubicBilinearAlgorithm honest = replayer.replay(lineage);
		Path stub = dir.resolve("2x2x1-r" + (honest.r - 1) + "-derived-0000000.json");
		SchemeIO.writeStub(honest, stub.toFile(), lineage, List.of("Z", "Q", "R", "C"));
		String overClaimed = Files.readString(stub).replace("\"m\": " + honest.r + ",", "\"m\": " + (honest.r - 1) + ",");
		Files.writeString(stub, overClaimed);
		SchemeIO.addFields(stub.toFile(), java.util.Map.of("source", "test-provenance", "year", 2026), true);

		Optional<RestampStubRanks.Drift> drift = RestampStubRanks.driftOf(stub, replayer);
		assertThat(drift).isPresent();
		assertThat(drift.get().overClaim()).isTrue();
		assertThat(drift.get().replayed()).isEqualTo(honest.r);

		Path out = RestampStubRanks.restamp(stub, replayer);

		assertThat(Files.exists(stub)).as("the stale file is deleted").isFalse();
		assertThat(out.getFileName().toString()).startsWith("2x2x1-r" + honest.r + "-derived-");
		JsonNode fixed = SchemeIO.parseJson(out.toFile());
		assertThat(fixed.get("m").asInt()).isEqualTo(honest.r);
		assertThat(SchemeIO.isStub(fixed)).isTrue();
		assertThat(SchemeIO.fieldTags(fixed)).containsExactly("Z", "Q", "R", "C");
		assertThat(fixed.get("source").asString()).as("provenance keys carried over").isEqualTo("test-provenance");
		assertThat(fixed.get("year").asInt()).isEqualTo(2026);
		assertThat(RestampStubRanks.driftOf(out, replayer)).as("honest after the rewrite").isEmpty();
	}

	@Test
	public void note_token_is_taken_from_the_old_name() {
		assertThat(RestampStubRanks.noteOf("13x13x17-r1876-derived-81d3a84.json")).isEqualTo("derived");
		assertThat(RestampStubRanks.noteOf("16x20x28-r4960-derived_serendip-818b551.json")).isEqualTo("derived_serendip");
		assertThat(RestampStubRanks.noteOf("legacy-name.json")).isEqualTo("derived");
	}
}
