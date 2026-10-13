package eu.solven.matmul.docs.migrate;

import static org.assertj.core.api.Assertions.assertThat;

import java.io.File;
import java.nio.charset.StandardCharsets;
import java.nio.file.Files;
import java.nio.file.Path;
import java.util.List;
import java.util.Map;

import org.junit.jupiter.api.Test;
import org.junit.jupiter.api.io.TempDir;

import eu.solven.matmul.NonCubicBilinearAlgorithm;
import eu.solven.matmul.catalog.SchemeIO;

import tools.jackson.databind.JsonNode;
import tools.jackson.databind.json.JsonMapper;

/**
 * The audit step of a contributed import: {@code discovery: "TBD"} → a verdict, taken
 * from the contribution's meta file and written through the canonical writer. The
 * trap is {@code SchemeIO.updateFields}, which never clobbers an existing key — a
 * re-stamp that only "adds" leaves every file at TBD while reporting success
 * (2026-10-12, first run on issue #12).
 */
public class TestRestampContributedProvenance {

	@Test
	public void provenance_is_read_from_the_meta_with_typed_discovery() throws Exception {
		JsonNode meta = JsonMapper.builder().build().readTree(
				"{\"source\":\"X 2026\",\"year\":2026,\"source_url\":\"https://x\",\"discovery\":\"true\","
						+ "\"attribution_for_rank\":null,\"discovery_note\":\"n\",\"notes\":\"m\"}");
		Map<String, Object> p = RestampContributedProvenance.provenance(meta);
		assertThat(p).containsEntry("discovery", Boolean.TRUE).containsEntry("year", 2026)
				.containsEntry("source_url", "https://x").containsEntry("discovery_note", "n").containsEntry("notes", "m")
				.doesNotContainKey("attribution_for_rank");
		assertThat(RestampContributedProvenance.provenance(JsonMapper.builder().build().readTree("{\"discovery\":\"TBD\"}")))
				.containsEntry("discovery", "TBD");
	}

	@Test
	public void an_existing_verdict_is_overwritten_not_kept(@TempDir Path tmp) throws Exception {
		File f = tmp.resolve("1x1x1-r1-x_2026-0000000.json").toFile();
		SchemeIO.write(NonCubicBilinearAlgorithm.naive(1, 1, 1), f);
		SchemeIO.addFields(f, Map.of("source", "X 2026", "discovery", "TBD", "discovery_note", "old"), true);
		assertThat(Files.readString(f.toPath(), StandardCharsets.UTF_8)).contains("\"discovery\": \"TBD\"");

		Map<String, Object> stamp = Map.of("discovery", Boolean.TRUE, "discovery_note", "audited");
		// The exact call the driver makes: remove-then-add, because updateFields never clobbers.
		SchemeIO.updateFields(f, stamp, stamp.keySet(), true);

		JsonNode root = SchemeIO.parseJson(f);
		assertThat(root.get("discovery").isBoolean() && root.get("discovery").booleanValue()).isTrue();
		assertThat(root.get("discovery_note").asString()).isEqualTo("audited");
		assertThat(root.get("source").asString()).as("untouched keys stay").isEqualTo("X 2026");
		assertThat(SchemeIO.readBilinear(f).r).isEqualTo(1);
	}

	@Test
	public void keys_are_the_ones_the_importer_stamps() {
		assertThat(RestampContributedProvenance.KEYS)
				.containsExactly(List.of("year", "source_url", "discovery", "attribution_for_rank", "discovery_note", "notes")
						.toArray(new String[0]));
	}
}
