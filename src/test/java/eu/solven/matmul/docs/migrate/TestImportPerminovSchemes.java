package eu.solven.matmul.docs.migrate;

import static org.assertj.core.api.Assertions.assertThat;

import java.util.List;

import org.junit.jupiter.api.Test;

/**
 * The Perminov importer only pulls {@code schemes/results/*}. A third-party source he
 * starts mirroring under {@code schemes/known/<new-sub>/} is therefore invisible to it
 * — which is how {@code MerlijnW70_fmm_schemes} (a new best over every ring at
 * ⟨11,13,15⟩) and {@code lita} sat upstream, un-imported, behind a green sync job
 * (2026-09-29/30). The importer must name every sub-folder it has no mapping for.
 */
public class TestImportPerminovSchemes {

	@Test
	public void an_unmapped_known_subtree_is_reported_with_its_file_count() {
		List<String> upstream = List.of(
				"README.md",
				"schemes/status.json",
				"schemes/results/ZT/2x2x2_m7_ZT.json",
				"schemes/known/tensor/2x4x9_tensor.mpl",
				"schemes/known/meta_flip_graph/277/k36da0c1bbdd8fb9.m",
				"schemes/known/MerlijnW70_fmm_schemes/11x13x15_m1364_ZT.json",
				"schemes/known/lita/13x13x13_m1379_Q.json",
				"schemes/known/someone_new_2027/3x3x3_m23_ZT.json",
				"schemes/known/someone_new_2027/4x4x4_m48_Q.json");

		assertThat(ImportPerminovSchemes.unmappedKnownSubtrees(upstream))
				.containsExactly(java.util.Map.entry("someone_new_2027", 2));
	}

	@Test
	public void the_current_upstream_layout_has_no_unmapped_subtree() {
		// Every sub-folder present upstream on 2026-09-30.
		List<String> upstream = List.of("MerlijnW70_fmm_schemes", "a_60_addition", "alpha_evolve", "alpha_tensor", "classic",
				"fmm_add_reduction", "jakobmoosbauer_flips", "jakobmoosbauer_symmetric_flips", "lita", "matmulcatalog",
				"meta_flip_graph", "tensor").stream().map(sub -> "schemes/known/" + sub + "/x.json").toList();

		assertThat(ImportPerminovSchemes.unmappedKnownSubtrees(upstream)).isEmpty();
	}
}
