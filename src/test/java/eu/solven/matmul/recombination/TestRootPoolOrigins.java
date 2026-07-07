package eu.solven.matmul.recombination;

import static org.assertj.core.api.Assertions.assertThat;

import java.util.List;

import org.junit.jupiter.api.Test;

import eu.solven.matmul.AxisSplitBases;
import eu.solven.matmul.NonCubicBilinearAlgorithm;
import eu.solven.matmul.catalog.SchemeIO;
import eu.solven.matmul.search.RecombinationPoolConfig;

/**
 * Going-forward pinning guard (user 2026-07-07): every stock recombination pool
 * entry must carry an {@code originLineage}, so recombination stubs record their
 * outer base as a pinned, replay-exact ref ({@code shape@contentHash} /
 * {@code naive-NxMxP} / a Transpose-or-AxisPermute wrapper thereof) — never the
 * display label. A null origin makes {@code RecursiveMaterialiser} fall back to
 * {@code Atom(label)}: a best-at-shape cited-bound ref that re-resolves against
 * future catalogs and marks the stub {@code explicitable:false}. The historical
 * ~2.1k unpinned derived stubs all trace back to null-origin pool entries.
 */
public class TestRootPoolOrigins {

	@Test
	public void every_root_pool_entry_carries_an_origin() {
		for (BlockSplitSearch.NamedBase nb : BlockSplitSearch.rootPool()) {
			assertThat(nb.originLineage())
					.as("rootPool entry '%s' must carry an origin lineage", nb.label())
					.isNotNull();
		}
	}

	@Test
	public void every_default_pool_entry_carries_an_origin() {
		for (BlockSplitSearch.NamedBase nb : BlockSplitSearch.defaultPool()) {
			assertThat(nb.originLineage())
					.as("defaultPool entry '%s' must carry an origin lineage", nb.label())
					.isNotNull();
		}
	}

	@Test
	public void every_thorough_pool_entry_carries_an_origin() {
		List<BlockSplitSearch.NamedBase> pool =
				BlockSplitSearch.buildPool(RecombinationPoolConfig.thorough(), "Q");
		assertThat(pool).isNotEmpty();
		for (BlockSplitSearch.NamedBase nb : pool) {
			assertThat(nb.originLineage())
					.as("thorough-pool entry '%s' must carry an origin lineage", nb.label())
					.isNotNull();
		}
	}

	/**
	 * The AxisSplit bases are pinned as {@code naive-NxMxP}: that ref must replay to
	 * the EXACT same scheme (same product order — the recombination allocation maps
	 * products by index, so a reordered replay would silently change the composition).
	 */
	@Test
	public void axis_split_bases_are_content_identical_to_naive() {
		assertThat(SchemeIO.contentHash(AxisSplitBases.mul211()))
				.isEqualTo(SchemeIO.contentHash(NonCubicBilinearAlgorithm.naive(2, 1, 1)));
		assertThat(SchemeIO.contentHash(AxisSplitBases.mul121()))
				.isEqualTo(SchemeIO.contentHash(NonCubicBilinearAlgorithm.naive(1, 2, 1)));
		assertThat(SchemeIO.contentHash(AxisSplitBases.mul112()))
				.isEqualTo(SchemeIO.contentHash(NonCubicBilinearAlgorithm.naive(1, 1, 2)));
	}
}
