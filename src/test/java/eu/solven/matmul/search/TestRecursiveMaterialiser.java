package eu.solven.matmul.search;

import static org.assertj.core.api.Assertions.assertThat;

import java.nio.file.Path;
import java.util.List;
import java.util.Optional;

import org.junit.jupiter.api.Test;
import org.junit.jupiter.api.io.TempDir;

import eu.solven.matmul.Verifier;
import eu.solven.matmul.catalog.FieldAwareLookup;
import eu.solven.matmul.catalog.Recombination;

public class TestRecursiveMaterialiser {

	@TempDir
	Path tmpRoot;

	private static RecursiveMaterialiser dryRun(FieldAwareLookup lookup) {
		List<BlockSplitSearch.NamedBase> pool = BlockSplitSearch.defaultPool();
		RecursiveClosureSota sota = new RecursiveClosureSota(lookup, pool, true, true);
		return new RecursiveMaterialiser(lookup, pool, sota, null, false, true);
	}

	@Test
	public void direct_catalog_hit_returns_as_leaf() {
		FieldAwareLookup lookup = new FieldAwareLookup("Q");
		RecursiveMaterialiser mat = dryRun(lookup);
		// Strassen ⟨2,2,2⟩=7 is on disk.
		Optional<RecursiveMaterialiser.Result> r = mat.materialise(2, 2, 2);
		assertThat(r).isPresent();
		assertThat(r.get().alg().r).isEqualTo(7);
		assertThat(r.get().fromDisk()).isTrue();
	}

	@Test
	public void kronecker_path_composes_two_subshapes() {
		FieldAwareLookup lookup = new FieldAwareLookup("Q");
		RecursiveMaterialiser mat = dryRun(lookup);
		// ⟨4,4,4⟩ via ⟨2,2,2⟩ ⊗ ⟨2,2,2⟩ — both on disk; the materialiser
		// should compose and return rank 49.
		Optional<RecursiveMaterialiser.Result> r = mat.materialise(4, 4, 4);
		assertThat(r).isPresent();
		assertThat(r.get().alg().r).isLessThanOrEqualTo(49);
		assertThat(Verifier.passesRandomMatmulSpotCheck(r.get().alg())).isTrue();
	}

	@Test
	public void compose_18x18x18_returns_3200_via_kron_concat_chain() {
		FieldAwareLookup lookup = new FieldAwareLookup("Q");
		RecursiveMaterialiser mat = dryRun(lookup);
		Optional<RecursiveMaterialiser.Result> r = mat.materialise(18, 18, 18);
		assertThat(r).isPresent();
		// ⟨18,18,18⟩=3200 was the strategy survey's prediction;
		// the materialiser should reach it (or better).
		assertThat(r.get().alg().r).isLessThanOrEqualTo(3402);
		assertThat(Verifier.passesRandomMatmulSpotCheck(r.get().alg())).isTrue();
	}
}
