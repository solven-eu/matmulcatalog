package eu.solven.matmul.papers.ta;

import static org.assertj.core.api.Assertions.assertThat;

import java.util.OptionalLong;

import org.junit.jupiter.api.Test;

/**
 * The TA-family registry: per-member domains, the Pan-family shim parity with the
 * historical {@code bestPanTaBound}, and the broader {@code bestRank} that admits
 * LITA (which wins large-even + all odd ≥19, but NOT every N — DIS still wins at 22).
 */
class TestTrilinearAggregations {

	@Test
	void pan_family_min_matches_historical_bestPanTaBound() {
		// Same table asserted by TestPanTrilinearAggregationBounds — the shim must not drift.
		assertThat(TrilinearAggregations.bestPanFamilyRank(4)).hasValue(62L);
		assertThat(TrilinearAggregations.bestPanFamilyRank(9)).hasValue(688L);  // odd → DIS
		assertThat(TrilinearAggregations.bestPanFamilyRank(18)).hasValue(3306L);
		assertThat(TrilinearAggregations.bestPanFamilyRank(22)).hasValue(5566L); // DIS even
		assertThat(TrilinearAggregations.bestPanFamilyRank(28)).hasValue(10550L); // SZ
		assertThat(TrilinearAggregations.bestPanFamilyRank(32)).hasValue(15096L); // SZ
		assertThat(TrilinearAggregations.bestPanFamilyRank(44)).hasValue(36110L); // SZ
	}

	@Test
	void bestRank_admits_lita_where_it_wins() {
		// Odd ≥19 and large even: LITA wins (these match FMM's index + the repo's npz).
		assertThat(TrilinearAggregations.best(19).orElseThrow().method().tag()).isEqualTo("TA_lita");
		assertThat(TrilinearAggregations.bestRank(19)).hasValue(4016L);
		assertThat(TrilinearAggregations.bestRank(21)).hasValue(5198L);
		assertThat(TrilinearAggregations.bestRank(23)).hasValue(6586L);
		assertThat(TrilinearAggregations.best(28).orElseThrow().method().tag()).isEqualTo("TA_lita");
		assertThat(TrilinearAggregations.bestRank(28)).hasValue(10535L);
		assertThat(TrilinearAggregations.bestRank(32)).hasValue(15079L);
		assertThat(TrilinearAggregations.bestRank(44)).hasValue(36087L);
	}

	@Test
	void lita_does_not_win_everywhere() {
		// At n=22, DIS (Islam even, 5566) beats LITA (5584) — the crossover isn't monotone.
		assertThat(TrilinearAggregations.best(22).orElseThrow().method().tag()).isEqualTo("TA_dis");
		assertThat(TrilinearAggregations.bestRank(22)).hasValue(5566L);
	}

	@Test
	void per_member_domains() {
		// HS / SZ: even, excluded at the n=16 pole.
		assertThat(TrilinearAggregations.HS.cubicRank(16)).isEmpty();
		assertThat(TrilinearAggregations.SZ.cubicRank(16)).isEmpty();
		assertThat(TrilinearAggregations.HS.cubicRank(18)).isNotEmpty();
		// PAN: even only.
		assertThat(TrilinearAggregations.PAN.cubicRank(9)).isEmpty();
		// LITA: n ≥ 19.
		assertThat(TrilinearAggregations.LITA.cubicRank(18)).isEmpty();
		assertThat(TrilinearAggregations.LITA.cubicRank(19)).hasValue(4016L);
		// DIS: both parities, no exclusions.
		assertThat(TrilinearAggregations.DIS.cubicRank(16)).isNotEmpty();
		assertThat(TrilinearAggregations.DIS.cubicRank(17)).isNotEmpty();
	}

	@Test
	void buildable_member_actually_builds_and_matches_rank() {
		// DIS is buildable for all n≥2; the built rank must equal the formula.
		var built = TrilinearAggregations.DIS.build(8).orElseThrow();
		assertThat((long) built.r).isEqualTo(TrilinearAggregations.DIS.cubicRank(8).getAsLong());
		// LITA now builds (LitaTaConstruction) for n≥19 — capability is cheap to
		// check; the actual build+verify is covered by TestLitaTaConstruction.
		assertThat(TrilinearAggregations.LITA.canBuild(18)).isFalse();
		assertThat(TrilinearAggregations.LITA.canBuild(21)).isTrue();
		// SZ now builds (TaNew25Construction) for even n≠16 — the built rank matches
		// the bound; exact tensor-verify is covered by TestTaNew25Construction.
		assertThat(TrilinearAggregations.SZ.canBuild(9)).isFalse();   // odd
		assertThat(TrilinearAggregations.SZ.canBuild(16)).isFalse();  // n=16 pole
		assertThat(TrilinearAggregations.SZ.canBuild(8)).isTrue();
		var sz = TrilinearAggregations.SZ.build(8).orElseThrow();
		assertThat((long) sz.r).isEqualTo(TrilinearAggregations.SZ.cubicRank(8).getAsLong());
		// The best BUILDABLE TA at a LITA-winning even N is now LITA (10535 < DIS 10556).
		assertThat(TrilinearAggregations.bestBuildable(28).orElseThrow().method().tag()).isEqualTo("TA_lita");
	}

	@Test
	void best_is_empty_below_domain() {
		// Below n=2 nothing applies.
		assertThat(TrilinearAggregations.bestRank(1)).isEqualTo(OptionalLong.empty());
	}
}
