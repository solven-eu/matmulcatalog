package eu.solven.matmul.search;

import eu.solven.matmul.recombination.BlockSplitSearch;

import static org.assertj.core.api.Assertions.assertThat;

import org.junit.jupiter.api.Test;

/**
 * Regression guard for the field-aware base-pool inclusion (Z⊂Q⊂R⊂C).
 *
 * <p>Silent bug: {@code isFieldValidLeafNC} filtered the extended pool by EXACT
 * {@code fieldTag ∈ fields[]}, so the ~98 base files stamped {@code fields:["Z"]}
 * (integer schemes) were rejected from an R/Q/C sweep — collapsing the
 * {@code includeDerived} R pool from 103 to 5 (Strassen-only). The search then
 * silently couldn't reach ⟨2,4,4⟩/⟨2,5,5⟩-based recombinations on R sweeps, and
 * nobody noticed because the curated {@code --base} pool masked it. The fix makes
 * the filter inclusion-aware via {@code Field.fallbackChain()}.
 *
 * <p>Invariants asserted (pool monotonicity under field inclusion + anti-collapse):
 * a scheme valid over a sub-field is valid over the larger field, so
 * {@code |Z| ≤ |Q| ≤ |R| ≤ |C|}, and none of the char-0 pools may collapse to the
 * bare root templates.</p>
 */
public class TestPoolFieldInclusion {

	private static int poolSize(String fieldTag) {
		return BlockSplitSearch.buildPool(PoolConfig.includeDerived(), fieldTag).size();
	}

	@Test
	public void char0PoolsAreInclusionMonotoneAndNotStarved() {
		int z = poolSize("Z");
		int q = poolSize("Q");
		int r = poolSize("R");
		int c = poolSize("C");

		// Z⊂Q⊂R⊂C ⇒ each larger field admits at least as many bases.
		assertThat(q).as("|Q| ≥ |Z|").isGreaterThanOrEqualTo(z);
		assertThat(r).as("|R| ≥ |Q|").isGreaterThanOrEqualTo(q);
		assertThat(c).as("|C| ≥ |R|").isGreaterThanOrEqualTo(r);

		// Anti-collapse: the pre-fix bug left R at the ~5 bare root templates. Any
		// healthy char-0 pool over the real catalog is far richer; 40 is a safe floor
		// well below the ~100 observed and well above the collapsed 5.
		assertThat(r).as("R pool must not collapse to the root templates").isGreaterThan(40);
		assertThat(z).as("Z pool not starved").isGreaterThan(40);
	}
}
