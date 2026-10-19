package eu.solven.matmul.catalog;

import static org.assertj.core.api.Assertions.assertThat;

import java.nio.file.Path;
import java.util.regex.Pattern;

import org.junit.jupiter.api.Test;

import eu.solven.matmul.NonCubicBilinearAlgorithm;
import eu.solven.matmul.search.LineageReplayer;

import tools.jackson.databind.JsonNode;

/**
 * {@link Lineage#atomFromFilename} is what the materialiser emits as the leaf for a
 * disk hit without its own lineage. For a canonical {@code {shape}-r{rank}-{note}-{hash7}}
 * name it must be the PINNED ref {@code shape@hash7} — exact content, replayed as such —
 * not the whole stem, which the replayer can only resolve by shape (the catalog-best of
 * the day). The unpinned stems were the mechanism behind stub rank drift: 419 distinct
 * ones sat in the catalog's lineages on 2026-10-18, and 510 stubs had moved away from
 * their recorded rank.
 */
public class TestLineageAtomFromFilename {

	private static final Pattern CANONICAL = Pattern.compile("^\\d+x\\d+x\\d+-r\\d+-.+-([0-9a-f]{7})\\.json$");

	@Test
	void canonical_filename_becomes_a_pinned_ref() {
		assertThat(Lineage.atomFromFilename("9x11x15-r960-derived-55c1918.json").ref()).isEqualTo("9x11x15@55c1918");
		assertThat(Lineage.atomFromFilename("known/section4/4x4x4-r47-alphatensor_F2-258e5b7.json").ref())
				.isEqualTo("4x4x4@258e5b7");
		assertThat(Lineage.atomFromFilename("16x20x28-r4960-derived_serendip-818b551").ref())
				.isEqualTo("16x20x28@818b551");
	}

	@Test
	void legacy_names_keep_their_canonical_key() {
		// Pre-existing behaviour, unchanged: an underscore-joined legacy name reduces to its
		// {shape}_m{rank}_a{adds} key; a hyphen-joined one (and a name without the a-count)
		// stays a whole stem that the replayer resolves by shape.
		assertThat(Lineage.atomFromFilename("alphatensor_Z_2x3x3_m15_a58.json").ref()).isEqualTo("2x3x3_m15_a58");
		assertThat(Lineage.atomFromFilename("fmm_lille-3x3x6_m40_a862.json").ref()).isEqualTo("fmm_lille-3x3x6_m40_a862");
		assertThat(Lineage.atomFromFilename("rosowski_2019_alg1-3x3x3_m21.json").ref())
				.isEqualTo("rosowski_2019_alg1-3x3x3_m21");
	}

	@Test
	void the_pinned_ref_of_a_catalog_file_replays_that_exact_file() throws Exception {
		FieldAwareLookup q = new FieldAwareLookup("Q");
		LineageReplayer replayer = LineageReplayer.withDefaultPool(q);
		Path file = q.findFiles(2, 2, 2).stream().filter(f -> {
			if (!CANONICAL.matcher(f.getFileName().toString()).matches()) return false;
			try {
				JsonNode r = SchemeIO.parseJson(f.toFile());
				return !SchemeIO.isStub(r) && !SchemeIO.isComplex(r) && !SchemeIO.isNonBilinear(r);
			} catch (Exception e) {
				return false;
			}
		}).findFirst().orElseThrow();

		Lineage.Atom atom = Lineage.atomFromFilename(file.getFileName().toString());
		assertThat(atom.ref()).contains("@");
		assertThat(replayer.classifyLeaf(atom.ref()).cls()).isEqualTo(LineageReplayer.LeafClass.PINNED);
		NonCubicBilinearAlgorithm a = replayer.replay(atom);
		assertThat(SchemeIO.contentHash(a)).startsWith(atom.ref().substring(atom.ref().indexOf('@') + 1));
		assertThat(SchemeIO.contentHash(a)).isEqualTo(SchemeIO.contentHash(SchemeIO.readBilinear(file.toFile())));
	}
}
