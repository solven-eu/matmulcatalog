package eu.solven.matmul.papers.khoruzhii2026;

import static org.assertj.core.api.Assertions.assertThat;
import static org.assertj.core.api.Assertions.assertThatThrownBy;

import org.junit.jupiter.api.Test;

/**
 * LITA cubic closed form (Khoruzhii–Gelß–Pokutta 2026) must reproduce the
 * values FMM-Lille's index reports for the large cubic formats. These doubled
 * as the "WORSE" cubic band in the FMM cross-check before the formula was wired
 * in, so this guard pins each one.
 */
class TestLitaTrilinearAggregation {

	@Test
	void odd_N_matches_FMM_index_values() {
		// R⟨N,N,N⟩ per FMM-Lille index (= LITA odd closed form).
		assertThat(LitaTrilinearAggregation.cubicRank(21)).isEqualTo(5198L);
		assertThat(LitaTrilinearAggregation.cubicRank(23)).isEqualTo(6586L);
		assertThat(LitaTrilinearAggregation.cubicRank(25)).isEqualTo(8196L);
		assertThat(LitaTrilinearAggregation.cubicRank(27)).isEqualTo(10045L);
		assertThat(LitaTrilinearAggregation.cubicRank(29)).isEqualTo(12147L);
		assertThat(LitaTrilinearAggregation.cubicRank(31)).isEqualTo(14519L);
	}

	@Test
	void even_N_matches_FMM_index_values_past_crossover() {
		// Even crossover is N≥26; FMM index matches the LITA even closed form there.
		assertThat(LitaTrilinearAggregation.cubicRank(26)).isEqualTo(8652L);
		assertThat(LitaTrilinearAggregation.cubicRank(28)).isEqualTo(10535L);
		assertThat(LitaTrilinearAggregation.cubicRank(30)).isEqualTo(12672L);
	}

	@Test
	void formula_is_integer_valued_across_the_valid_domain() {
		// Both branches must divide exactly by 12 — cubicRank throws otherwise.
		// Domain starts at MIN_N (=19); the construction is undefined below.
		for (int n = LitaTrilinearAggregation.MIN_N; n <= 64; n++) {
			assertThat(LitaTrilinearAggregation.cubicRank(n)).isPositive();
		}
	}

	@Test
	void rejects_N_below_19() {
		// TA_lita requires N > 18 — no construction exists below, so the closed
		// form must not be evaluated (let alone claimed as a bound) there.
		for (int n : new int[] { 0, 1, 4, 17, 18 }) {
			assertThatThrownBy(() -> LitaTrilinearAggregation.cubicRank(n))
					.isInstanceOf(IllegalArgumentException.class);
		}
		assertThat(LitaTrilinearAggregation.cubicRank(19)).isEqualTo(4016L);
	}
}
