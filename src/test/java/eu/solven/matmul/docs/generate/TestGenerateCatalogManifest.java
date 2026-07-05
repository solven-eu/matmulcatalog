package eu.solven.matmul.docs.generate;

import static org.assertj.core.api.Assertions.assertThat;

import org.junit.jupiter.api.Test;

/**
 * Guards the lineage→source attribution of source-less derived stubs. The
 * silent failure mode: a formula-derived stub (e.g. the projected
 * {@code DIS09Lemma4(n=30) ↓[…]} behind ⟨29,30,30⟩) carries no on-disk
 * {@code source}, so the manifest emitted {@code "unknown"} and the SPA showed
 * a literal "unknown" source — misplacing the DIS09 formula credit for 150+
 * rows without any crash.
 */
public class TestGenerateCatalogManifest {

	@Test
	public void dis09_lemma4_formula_is_attributed() {
		assertThat(GenerateCatalogManifest.attributeSourceFromLineage("DIS09Lemma4(n=30)"))
				.isEqualTo("DIS09");
	}

	@Test
	public void projected_dis09_lemma4_keeps_the_formula_credit() {
		// The ⟨29,30,30⟩ case: formula call + projection suffix. startsWith must
		// still attribute the constructive head.
		assertThat(GenerateCatalogManifest.attributeSourceFromLineage(
				"DIS09Lemma4(n=30) ↓[0,1,2|0,1,2,3|0,1,2,3]"))
				.isEqualTo("DIS09");
	}

	@Test
	public void ta_lita_formula_is_attributed() {
		assertThat(GenerateCatalogManifest.attributeSourceFromLineage("TA_lita(n=19)"))
				.isEqualTo("Khoruzhii, Gelß & Pokutta 2026 (LITA)");
	}

	@Test
	public void non_formula_lineages_stay_unknown() {
		// A bare scheme-ref projection is OUR derivation of a catalogued base —
		// no paper to credit; a composition likewise.
		assertThat(GenerateCatalogManifest.attributeSourceFromLineage("4x4x11 ↓[0,1,2,3|0,1,2,3|0,1,2]"))
				.isEqualTo("unknown");
		assertThat(GenerateCatalogManifest.attributeSourceFromLineage("2x2x2 ⊗ 3x3x3"))
				.isEqualTo("unknown");
		assertThat(GenerateCatalogManifest.attributeSourceFromLineage(null))
				.isEqualTo("unknown");
	}
}
