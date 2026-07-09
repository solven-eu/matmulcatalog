package eu.solven.matmul.docs.verify;

import static org.assertj.core.api.Assertions.assertThat;

import java.io.File;

import org.junit.jupiter.api.Test;

import tools.jackson.databind.JsonNode;
import tools.jackson.databind.ObjectMapper;

/**
 * Curation guard for {@code references/fmm-artifact-audit.json} — the manual
 * per-shape audit {@link FmmCrossCheck} consumes to split WORSE into
 * artifact-backed vs upstream-unverified (2026-07-09). A malformed entry
 * would silently drop the split (rows misclassified back into WORSE) or
 * crash the generator.
 */
public class TestFmmArtifactAudit {

	@Test
	public void audit_json_is_well_formed() throws Exception {
		File f = new File("references/fmm-artifact-audit.json");
		assertThat(f).exists();
		JsonNode root = new ObjectMapper().readTree(f);
		assertThat(root.get("audited")).isNotNull();
		JsonNode shapes = root.get("shapes");
		assertThat(shapes).isNotNull();
		assertThat(shapes.size()).isGreaterThan(0);
		shapes.properties().forEach(e -> {
			String[] d = e.getKey().split("x");
			assertThat(d).as("shape key " + e.getKey()).hasSize(3);
			int n = Integer.parseInt(d[0]), m = Integer.parseInt(d[1]), p = Integer.parseInt(d[2]);
			// FmmCrossCheck matches rows by their canonical (sorted) shape — an
			// unsorted key would never match and the row would silently stay WORSE.
			assertThat(n <= m && m <= p).as("key must be sorted n≤m≤p: " + e.getKey()).isTrue();
			String cls = e.getValue().get("class").asText();
			assertThat(cls).isIn("index_only", "placeholder");
		});
	}
}
