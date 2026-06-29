package eu.solven.matmul.search;

import eu.solven.matmul.recombination.BlockSplitSearch;

import static org.assertj.core.api.Assertions.assertThat;

import java.nio.file.Path;
import java.util.List;
import java.util.Optional;

import org.junit.jupiter.api.Test;
import org.junit.jupiter.api.io.TempDir;

import eu.solven.matmul.Verifier;
import eu.solven.matmul.catalog.FieldAwareLookup;
import eu.solven.matmul.catalog.Lineage;
import eu.solven.matmul.recombination.Recombination;

public class TestRecursiveMaterialiser {

	@TempDir
	Path tmpRoot;

	private static RecursiveMaterialiser dryRun(FieldAwareLookup lookup) {
		List<BlockSplitSearch.NamedBase> pool = BlockSplitSearch.defaultPool();
		RecursiveClosureSota sota = new RecursiveClosureSota(lookup, pool, true, true);
		return new RecursiveMaterialiser(lookup, pool, sota, null, false, true);
	}

	@Test
	public void write_guard_rejects_same_shape_atom_lineage() {
		// The strengthened write-time guard must REFUSE a lineage that builds a shape from an
		// Atom of its OWN ordered shape (degenerate self-derivation: "⟨2,2,3⟩ from ⟨2,2,3⟩").
		// This terminates (the atom is acyclic), so the @hash-cycle DFS alone would not catch it.
		RecursiveMaterialiser mat = dryRun(new FieldAwareLookup("Q"));

		Lineage.Node selfShape = new Lineage.ConcatCols(
				new Lineage.Atom("naive-2x2x1"), new Lineage.Atom("2x2x3@deadbee"));
		assertThat(mat.lineageCorruption(selfShape, 2, 2, 3, "abc1234"))
				.as("same-ordered-shape atom must be flagged SELF-SHAPE")
				.contains("SELF-SHAPE");

		// An orientation/transpose atom (same multiset, DIFFERENT ordered shape) is legitimate
		// — it must NOT be flagged SELF-SHAPE. Bare refs resolve to catalog-best (no @hash edge
		// followed), so the guard returns null for this clean orientation lineage.
		Lineage.Node orientation = new Lineage.ConcatCols(
				new Lineage.Atom("naive-2x2x1"), new Lineage.Atom("2x3x2"));
		String r = mat.lineageCorruption(orientation, 2, 2, 3, "abc1234");
		assertThat(r == null || !r.contains("SELF-SHAPE"))
				.as("orientation/transpose atom must NOT be SELF-SHAPE (got: %s)", r)
				.isTrue();
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
	public void sota_18x18x18_is_3200_in_catalog_and_kron_path_is_exact() {
		FieldAwareLookup lookup = new FieldAwareLookup("Q");
		// Catalog SOTA ⟨18,18,18⟩=3200 — a deep ConcatCols(KronProduct(Project(⟨3,3,7⟩…)))
		// chain (maxDim>16 ⇒ a lineage stub). findRank reads it straight from the index
		// (instant) — the right way to assert SOTA without a multi-minute compose sweep.
		assertThat(lookup.findRank(18, 18, 18)).isLessThanOrEqualTo(3200);

		// The materialiser's FAST always-on path (kron + concat only; the expensive
		// recombination/serendipitous/projection strategies disabled) must still compose an
		// EXACT scheme. It reaches 3306 here (the kron/concat optimum); 3200 needs the deep
		// chain above, which is the slow path we deliberately exclude. ≤ ⇒ SOTA-or-better.
		RecursiveMaterialiser mat = dryRun(lookup);
		mat.setStrategies(java.util.Set.of());
		Optional<RecursiveMaterialiser.Result> r = mat.materialise(18, 18, 18);
		assertThat(r).isPresent();
		assertThat(r.get().alg().r).isLessThanOrEqualTo(3306);
		assertThat(Verifier.passesRandomMatmulSpotCheck(r.get().alg())).isTrue();
	}
}
