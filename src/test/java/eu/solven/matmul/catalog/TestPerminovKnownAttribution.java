package eu.solven.matmul.catalog;

import static org.assertj.core.api.Assertions.assertThat;

import java.nio.file.Files;
import java.nio.file.Path;
import java.util.List;
import java.util.stream.Collectors;

import org.junit.jupiter.api.Test;

import eu.solven.matmul.catalog.PerminovKnownAttribution.Attribution;
import eu.solven.matmul.catalog.PerminovKnownAttribution.Disposition;

import tools.jackson.databind.JsonNode;

/**
 * Guards the Perminov re-attribution rule: a scheme that Perminov files under
 * {@code schemes/known/<sub>} originates ELSEWHERE and must NOT be credited to
 * Perminov. The canonical regression is {@code ⟨2,7,7⟩=76} (Kauers & Wood 2025,
 * not Perminov). Only {@code schemes/results/*} stays Perminov.
 */
class TestPerminovKnownAttribution {

	@Test
	void metaFlipGraph_is_kauers_wood_2025() {
		// The ⟨2,7,7⟩=76 case the user flagged: known/meta_flip_graph → Kauers & Wood 2025.
		Attribution a = PerminovKnownAttribution
				.forPath("schemes/known/meta_flip_graph/277/k36da0c1bbdd8fb9.m").orElseThrow();
		assertThat(a.source()).isEqualTo("Kauers & Wood 2025");
		assertThat(a.disposition()).isEqualTo(Disposition.EXTERNAL);
		assertThat(a.isPerminovOwn()).isFalse();
	}

	@Test
	void known_subsources_map_to_their_true_origin() {
		assertThat(PerminovKnownAttribution.forPath("schemes/known/alpha_tensor/x.m").orElseThrow().source())
				.isEqualTo("AlphaTensor 2022");
		assertThat(PerminovKnownAttribution.forPath("schemes/known/alpha_evolve/x.m").orElseThrow().source())
				.isEqualTo("AlphaEvolve 2025");
		assertThat(PerminovKnownAttribution.forPath("schemes/known/a_60_addition/x.m").orElseThrow().source())
				.isEqualTo("Stapleton 2025 (a=60)");
		assertThat(PerminovKnownAttribution.forPath("schemes/known/jakobmoosbauer_flips/x.m").orElseThrow().source())
				.isEqualTo("Kauers-Moosbauer 2023");
		assertThat(PerminovKnownAttribution.forPath("schemes/known/jakobmoosbauer_symmetric_flips/x.m")
				.orElseThrow().source()).isEqualTo("Moosbauer-Poole 2025");
	}

	@Test
	void classic_attributes_per_filename_author() {
		assertThat(PerminovKnownAttribution.forPath("schemes/known/classic/Strassen-222-7-18.m")
				.orElseThrow().source()).isEqualTo("Strassen 1969");
		assertThat(PerminovKnownAttribution.forPath("schemes/known/classic/Laderman-333-23.m")
				.orElseThrow().source()).isEqualTo("Laderman 1976");
		assertThat(PerminovKnownAttribution.forPath("schemes/known/classic/Smirnov-333-23-139.m")
				.orElseThrow().source()).isEqualTo("Smirnov 2013");
	}

	@Test
	void tensor_and_matmulcatalog_are_skip_fresh_import_mirrors() {
		Attribution tensor = PerminovKnownAttribution
				.forPath("schemes/known/tensor/2x4x9_tensor.mpl").orElseThrow();
		assertThat(tensor.source()).isEqualTo("FMM-Lille");
		assertThat(tensor.disposition()).isEqualTo(Disposition.SKIP_FRESH_IMPORT);

		Attribution self = PerminovKnownAttribution
				.forPath("schemes/known/matmulcatalog/whatever.json").orElseThrow();
		assertThat(self.disposition()).isEqualTo(Disposition.SKIP_FRESH_IMPORT);
	}

	@Test
	void results_subtree_stays_perminov_own() {
		Attribution a = PerminovKnownAttribution
				.forPath("schemes/results/Z/2x7x7_m76_cr320.json").orElseThrow();
		assertThat(a.source()).isEqualTo("Perminov 2023");
		assertThat(a.isPerminovOwn()).isTrue();
	}

	/**
	 * The catalog-wide invariant after {@link eu.solven.matmul.docs.migrate.AttributePerminovKnown}:
	 * NO on-disk scheme whose {@code original_source_path} sits under {@code known/}
	 * may still be credited to Perminov. Fast — reads JSON metadata only, never
	 * expands matrices.
	 */
	@Test
	void no_known_scheme_remains_credited_to_perminov() throws Exception {
		Path root = Path.of("src/main/resources/schemes");
		if (!Files.isDirectory(root)) {
			return; // resources not present in this run — nothing to assert
		}
		List<Path> offenders;
		try (var s = Files.walk(root)) {
			offenders = s.filter(p -> p.toString().endsWith(".json")).filter(p -> {
				try {
					JsonNode n = SchemeIO.parseJson(p.toFile());
					String src = n.has("source") ? n.get("source").asString().toLowerCase() : "";
					String osp = n.has("original_source_path") ? n.get("original_source_path").asString() : "";
					return src.contains("perminov") && osp.contains("known/");
				} catch (Exception e) {
					return false;
				}
			}).collect(Collectors.toList());
		}
		assertThat(offenders)
				.as("schemes under Perminov's known/ subtree must be attributed to their true origin, not Perminov")
				.isEmpty();
	}
}
