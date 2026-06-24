package eu.solven.matmul.search;

import eu.solven.matmul.recombination.BlockSplitSearch;

import static org.assertj.core.api.Assertions.assertThat;

import java.util.List;

import org.junit.jupiter.api.Test;

/**
 * Regression guard: the char-0 ("Z-arithmetic NC") outer-base pool must NOT contain
 * field-restricted bases (F₂-only / C-only). The field filter
 * ({@code BlockSplitSearch.isLeafZArithmeticNC}) historically checked the obsolete
 * singular {@code "field"} key; after the catalog migrated to a {@code fields[]} ARRAY
 * the check silently matched NOTHING, so every F₂-only base leaked in — e.g. AlphaTensor
 * ⟨4,4,4⟩=47/F₂ used as a Q-sweep outer base produced a phantom ⟨22,28,28⟩=9316 that the
 * materialise spot-check then rejected (a "best scheme whatever the field" footgun).
 */
public class TestPoolFieldFilter {

	private static boolean hasFourCube47(List<BlockSplitSearch.NamedBase> pool) {
		// AlphaTensor ⟨4,4,4⟩=47 is F₂-only (the Z/char-0 optimum is 49), so a ⟨4,4,4⟩=47
		// entry is necessarily the F₂ scheme.
		return pool.stream().anyMatch(nb -> {
			var b = nb.base();
			return b.n == 4 && b.m == 4 && b.p == 4 && b.r == 47;
		});
	}

	@Test
	public void char0_pool_excludes_f2_only_bases() {
		List<BlockSplitSearch.NamedBase> pool = BlockSplitSearch.extendedPool(5);  // null = char-0
		assertThat(pool).as("the extended pool must not be empty").isNotEmpty();
		assertThat(hasFourCube47(pool))
				.as("F₂-only ⟨4,4,4⟩=47 must NOT be in the char-0 pool (the fields[] filter bug)")
				.isFalse();
		assertThat(BlockSplitSearch.extendedPool(5, "Q").stream().anyMatch(nb -> {
			var b = nb.base();
			return b.n == 2 && b.m == 2 && b.p == 2 && b.r == 7;
		})).as("a char-0 base (Strassen ⟨2,2,2⟩=7) must remain in the Q pool").isTrue();
	}

	@Test
	public void q_pool_excludes_but_f2_pool_includes_the_f2_native_base() {
		// Field-coherent pools: a Q sweep must NOT see the F₂-only ⟨4,4,4⟩=47, but an F₂
		// sweep MUST (it's a legitimate F₂-native base — the point of an explicit field).
		assertThat(hasFourCube47(BlockSplitSearch.extendedPool(5, "Q")))
				.as("Q pool excludes F₂-only ⟨4,4,4⟩=47").isFalse();
		assertThat(hasFourCube47(BlockSplitSearch.extendedPool(5, "F2")))
				.as("F₂ pool INCLUDES the F₂-native ⟨4,4,4⟩=47").isTrue();
	}
}
