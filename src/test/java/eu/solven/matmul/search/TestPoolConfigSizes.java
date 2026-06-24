package eu.solven.matmul.search;

import eu.solven.matmul.recombination.BlockSplitSearch;

import org.junit.jupiter.api.Tag;


import static org.assertj.core.api.Assertions.assertThat;

import org.junit.jupiter.api.Test;

import eu.solven.matmul.SymmetryTransforms.InternalOrbitMode;

/**
 * Print the pool size for each of the 8 meaningful PoolConfig points
 * (2³ binary axes × canonical orbit; permutation only sampled on the
 * tractable maxBaseDim=3 corner). Surfaces the "fast → slow" cost
 * gradient up-front so the runner can size compute budgets.
 */
@Tag("slow")
class TestPoolConfigSizes {

	@Test
	void enumerateAllConfigSizes() {
		System.out.println();
		System.out.println("=== Pool sizes across PoolConfig ===");
		System.out.println();
		System.out.printf("  %-30s  %-50s  %-6s%n", "preset", "config", "size");
		System.out.println("  " + "-".repeat(95));
		print("simple()", PoolConfig.simple());
		print("auditAxisFlip()", PoolConfig.auditAxisFlip());
		print("auditPermutation()", PoolConfig.auditPermutation());
		print("rectangular()", PoolConfig.rectangular());
		print("includeDerived()", PoolConfig.includeDerived());
		print("thorough()", PoolConfig.thorough());
		// Custom corners not covered by presets:
		print("(cub, root, perm, maxDim=3)",
				new PoolConfig(true, true, InternalOrbitMode.PERMUTATION_BOUNDED, 3, 216));
		print("(any, all, axis-flip, maxDim=5)",
				new PoolConfig(false, false, InternalOrbitMode.AXIS_FLIP, 5, 8));
		System.out.println();
		// Sanity: simple is the smallest.
		assertThat(BlockSplitSearch.buildPool(PoolConfig.simple()).size())
				.isLessThan(BlockSplitSearch.buildPool(PoolConfig.thorough()).size());
	}

	private static void print(String preset, PoolConfig cfg) {
		int size = BlockSplitSearch.buildPool(cfg).size();
		System.out.printf("  %-30s  %-50s  %6d%n", preset, cfg.shortLabel(), size);
	}
}
