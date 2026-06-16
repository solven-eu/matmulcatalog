package eu.solven.matmul.catalog;

import static org.assertj.core.api.Assertions.assertThat;

import java.io.File;
import java.nio.file.Path;
import java.util.List;
import java.util.Map;

import org.junit.jupiter.api.Test;
import org.junit.jupiter.api.io.TempDir;

import eu.solven.matmul.NonCubicBilinearAlgorithm;
import tools.jackson.databind.JsonNode;

/**
 * Regression guard for a SILENT metadata bug (2026-06-11): {@code addFields}
 * stringified List values, so a stamped {@code fields: ["F2",…]} landed as the
 * single string {@code "[F2, F3, …]"} — and the content-driven index treated
 * the scheme as field-less (invisible to every field-chain lookup). Lists must
 * round-trip as JSON ARRAYS.
 */
public class TestSchemeIOAddFieldsArray {

	@TempDir
	Path tmp;

	@Test
	public void list_values_become_json_arrays() throws Exception {
		File f = tmp.resolve("naive-1x1x1.json").toFile();
		SchemeIO.write(NonCubicBilinearAlgorithm.naive(1, 1, 1), f);
		SchemeIO.addFields(f, Map.of(
				"fields", List.of("F2", "Z"),
				"self_serendipity_savings", List.of(18L, 12L, 12L),
				"source", "metaflip"), true);
		JsonNode root = SchemeIO.parseJson(f);
		assertThat(root.get("fields").isArray())
				.as("fields must be a JSON array, not a stringified list").isTrue();
		assertThat(root.get("fields").get(0).asString()).isEqualTo("F2");
		assertThat(root.get("self_serendipity_savings").isArray()).isTrue();
		assertThat(root.get("self_serendipity_savings").get(0).asLong()).isEqualTo(18);
		assertThat(root.get("source").asString()).isEqualTo("metaflip");
		assertThat(SchemeIO.fieldTags(root)).contains("F2", "Z");
	}
}
