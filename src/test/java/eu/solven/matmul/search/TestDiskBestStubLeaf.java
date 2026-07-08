package eu.solven.matmul.search;

import static org.assertj.core.api.Assertions.assertThat;

import java.nio.file.Path;
import java.util.Optional;

import org.junit.jupiter.api.Test;

import eu.solven.matmul.catalog.FieldAwareLookup;
import eu.solven.matmul.catalog.Lineage;
import eu.solven.matmul.recombination.BlockSplitSearch;

/**
 * Regression guard (fmm-gap 2026-07-08, ⟨7,14,32⟩): the composed-build leaf
 * fetch ({@code RecursiveLookup.find} → {@code diskBest}) used stub-blind
 * {@code findWithSource}, so a recombination/concat whose scored leaf existed
 * only as a lineage stub (⟨7,14,30⟩=1865) threw "construct: missing
 * sub-algorithm" at build time and the strategy was silently discarded —
 * FMM's plain concat ⟨7,14,2⟩:152 + ⟨7,14,30⟩:1865 = 2017 looked unreachable.
 * {@code diskBest} must resolve a stub-only-best shape by REPLAY (via
 * resolveParentHit) and pin the leaf durably.
 */
public class TestDiskBestStubLeaf {

	@Test
	public void disk_best_resolves_stub_only_leaf_by_replay() {
		FieldAwareLookup lookup = new FieldAwareLookup("Q");
		RecursiveMaterialiser mat = new RecursiveMaterialiser(
				lookup, BlockSplitSearch.defaultPool(), new RecursiveClosureSota(
						lookup, BlockSplitSearch.defaultPool(), true, true),
				Path.of("src/main/resources/schemes"), false, false, false, false);
		// ⟨7,14,30⟩'s best is the r=1865 lineage-only stub (2026-07-08 harvest);
		// findWithSource alone cannot see it.
		int claimed = lookup.findRank(7, 14, 30);
		assertThat(claimed).as("precondition: stub rank visible to findRank").isLessThanOrEqualTo(1865);
		Optional<RecursiveMaterialiser.Result> r = mat.diskBest(7, 14, 30);
		assertThat(r).as("diskBest must resolve the stub-only ⟨7,14,30⟩").isPresent();
		assertThat(r.get().alg().r).as("resolved at the stub's claimed rank, not a worse dense sibling")
				.isLessThanOrEqualTo(claimed);
		// The leaf must be durably pinned (hash ref or explicit naive) — never a
		// bare best-at-shape label that replays to arbitrary future catalog state.
		Lineage.Node leaf = r.get().lineage();
		assertThat(leaf).isNotNull();
		String s = Lineage.prettyCompact(leaf);
		assertThat(s).as("durable pin, got: " + s).contains("7x14x30");
	}
}
