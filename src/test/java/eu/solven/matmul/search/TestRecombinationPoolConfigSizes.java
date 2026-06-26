package eu.solven.matmul.search;

import eu.solven.matmul.recombination.BlockSplitSearch;

import org.junit.jupiter.api.Tag;


import static org.assertj.core.api.Assertions.assertThat;

import org.junit.jupiter.api.Test;

import eu.solven.matmul.SymmetryTransforms.InternalOrbitMode;

/**
 * Print the pool size for each of the 8 meaningful RecombinationPoolConfig points
 * (2³ binary axes × canonical orbit; permutation only sampled on the
 * tractable maxBaseDim=3 corner). Surfaces the "fast → slow" cost
 * gradient up-front so the runner can size compute budgets.
 */
@Tag("slow")
class TestRecombinationPoolConfigSizes {

	@Test
	void enumerateAllConfigSizes() {
		System.out.println();
		System.out.println("=== Pool sizes across RecombinationPoolConfig ===");
		System.out.println();
		System.out.printf("  %-30s  %-50s  %-6s%n", "preset", "config", "size");
		System.out.println("  " + "-".repeat(95));
		print("simple()", RecombinationPoolConfig.simple());
		print("auditAxisFlip()", RecombinationPoolConfig.auditAxisFlip());
		print("auditPermutation()", RecombinationPoolConfig.auditPermutation());
		print("rectangular()", RecombinationPoolConfig.rectangular());
		print("includeDerived()", RecombinationPoolConfig.includeDerived());
		print("thorough()", RecombinationPoolConfig.thorough());
		// Custom corners not covered by presets:
		print("(cub, root, perm, maxDim=3)",
				new RecombinationPoolConfig(true, true, InternalOrbitMode.PERMUTATION_BOUNDED, 3, 216));
		print("(any, all, axis-flip, maxDim=5)",
				new RecombinationPoolConfig(false, false, InternalOrbitMode.AXIS_FLIP, 5, 8));
		System.out.println();
		// Sanity: simple is the smallest.
		assertThat(BlockSplitSearch.buildPool(RecombinationPoolConfig.simple()).size())
				.isLessThan(BlockSplitSearch.buildPool(RecombinationPoolConfig.thorough()).size());
	}

	private static void print(String preset, RecombinationPoolConfig cfg) {
		int size = BlockSplitSearch.buildPool(cfg).size();
		System.out.printf("  %-30s  %-50s  %6d%n", preset, cfg.shortLabel(), size);
	}
}
