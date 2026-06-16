package eu.solven.matmul.search;

import static org.assertj.core.api.Assertions.assertThat;
import static org.assertj.core.api.Assertions.assertThatCode;
import static org.assertj.core.api.Assertions.assertThatThrownBy;

import org.junit.jupiter.api.Test;

/**
 * Regression guard for the silent predict/build-divergence bug: a rebuilt
 * lineage whose replayed rank was WORSE than the evaluated rank used to be
 * logged as a {@code [replay-diag]} WARN and silently discarded — letting
 * over-claiming stubs survive in the catalog at their inflated rank (e.g.
 * ⟨4,19,20⟩ evaluated 1000 but replaying to 1016). The fix makes
 * {@link RecursiveMaterialiser#assertRebuildNotWorse} FAIL LOUD on a worse
 * rebuild, while still accepting a BETTER one (projection DCE legitimately
 * drops products). See {@code feedback_fail_loud_dont_swallow}.
 */
public class TestRebuildNotWorseGuard {

	@Test
	public void worse_rebuild_throws() {
		// evaluated 1000, replayed 1016 — the canonical ⟨4,19,20⟩ phantom.
		assertThatThrownBy(() -> RecursiveMaterialiser.assertRebuildNotWorse(4, 19, 20, 1000, 1016))
				.isInstanceOf(IllegalStateException.class)
				.hasMessageContaining("rebuild WORSE than evaluated")
				.hasMessageContaining("evaluated r=1000")
				.hasMessageContaining("replayed to r=1016");
	}

	@Test
	public void equal_rebuild_is_accepted() {
		assertThatCode(() -> RecursiveMaterialiser.assertRebuildNotWorse(4, 19, 20, 1000, 1000))
				.doesNotThrowAnyException();
	}

	@Test
	public void better_rebuild_is_accepted() {
		// DCE in projection can drop products → a strictly better replay is fine.
		assertThatCode(() -> RecursiveMaterialiser.assertRebuildNotWorse(4, 19, 20, 1000, 996))
				.doesNotThrowAnyException();
	}

	@Test
	public void off_by_one_worse_still_throws() {
		// The asymmetry is exact: even +1 worse is a bug, not rounding tolerance.
		assertThatThrownBy(() -> RecursiveMaterialiser.assertRebuildNotWorse(2, 2, 2, 7, 8))
				.isInstanceOf(IllegalStateException.class);
		assertThatCode(() -> RecursiveMaterialiser.assertRebuildNotWorse(2, 2, 2, 8, 7))
				.doesNotThrowAnyException();
	}

	@Test
	public void message_names_the_shape() {
		assertThatThrownBy(() -> RecursiveMaterialiser.assertRebuildNotWorse(13, 19, 29, 4248, 4250))
				.hasMessageContaining("⟨13,19,29⟩");
		// sanity: a correct rebuild path is silent
		assertThat(true).isTrue();
	}
}
