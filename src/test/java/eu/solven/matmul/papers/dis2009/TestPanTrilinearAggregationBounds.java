package eu.solven.matmul.papers.dis2009;

import static org.assertj.core.api.Assertions.assertThat;

import org.junit.jupiter.api.Test;

/**
 * Unit tests for the closed-form Pan-TA bound formulas in
 * {@link PanTrilinearAggregation}. Anchors are pulled from:
 * <ul>
 *   <li>DIS09 Table 3 (Islam 2009 / Drevet–Islam–Schost 2009),</li>
 *   <li>Schwartz–Zwecher 2025 (arXiv:2508.01748) Table 1 / Thm 3.4,</li>
 *   <li>Hand-verification of Pan 1980 (SIAM J. Comput. 9(2)) and the
 *   Hadas–Schwartz 1982 formalisation at n=18.</li>
 * </ul>
 */
class TestPanTrilinearAggregationBounds {

	// ---------------------------------------------------------------
	// Anchor: Pan 1980 at n=18 → 3618.
	// (n³ + 9n²/2 − 3n)/2 = (5832 + 1458 − 54)/2 = 7236/2 = 3618.
	// ---------------------------------------------------------------
	@Test
	void panSiam1980_n18() {
		assertThat(PanTrilinearAggregation.panSiam1980Bound(18)).isEqualTo(3618L);
	}

	@Test
	void panSiam1980_oddReturnsMinusOne() {
		assertThat(PanTrilinearAggregation.panSiam1980Bound(17)).isEqualTo(-1L);
		assertThat(PanTrilinearAggregation.panSiam1980Bound(23)).isEqualTo(-1L);
	}

	// ---------------------------------------------------------------
	// Anchor: Hadas–Schwartz 1982 at n=18 → 3360.
	// n³/3 + 15n²/4 + 32n/3 + 9 = 1944 + 1215 + 192 + 9 = 3360.
	// ---------------------------------------------------------------
	@Test
	void panHadasSchwartz1982_n18() {
		assertThat(PanTrilinearAggregation.panHadasSchwartz1982Bound(18)).isEqualTo(3360L);
	}

	@Test
	void panHadasSchwartz1982_n16NotApplicable() {
		// n=16 is excluded — construction divides by (1 − 9/(n/2+1)) = 0.
		assertThat(PanTrilinearAggregation.panHadasSchwartz1982Bound(16)).isEqualTo(-1L);
	}

	@Test
	void panHadasSchwartz1982_oddReturnsMinusOne() {
		assertThat(PanTrilinearAggregation.panHadasSchwartz1982Bound(17)).isEqualTo(-1L);
	}

	// ---------------------------------------------------------------
	// Anchor: Schwartz–Zwecher 2025 at n=18 → 3350.
	// n³/3 + 15n²/4 + 61n/6 + 8 = 1944 + 1215 + 183 + 8 = 3350.
	// ---------------------------------------------------------------
	@Test
	void schwartzZwecher2025_n18() {
		assertThat(PanTrilinearAggregation.schwartzZwecher2025Bound(18)).isEqualTo(3350L);
	}

	// SZ 2025 Table 1 anchor.
	@Test
	void schwartzZwecher2025_n44() {
		assertThat(PanTrilinearAggregation.schwartzZwecher2025Bound(44)).isEqualTo(36110L);
	}

	@Test
	void schwartzZwecher2025_n16NotApplicable() {
		assertThat(PanTrilinearAggregation.schwartzZwecher2025Bound(16)).isEqualTo(-1L);
	}

	@Test
	void schwartzZwecher2025_oddReturnsMinusOne() {
		assertThat(PanTrilinearAggregation.schwartzZwecher2025Bound(17)).isEqualTo(-1L);
	}

	// ---------------------------------------------------------------
	// Existing Islam 2009 anchors (sanity check that the new methods
	// don't interfere with the existing cubicBound entry point).
	// ---------------------------------------------------------------
	@Test
	void islam2009Even_n18_isCurrentBestAtThisN() {
		// (n³ + 12n² + 11n)/3 = (5832 + 3888 + 198)/3 = 9918/3 = 3306.
		assertThat(PanTrilinearAggregation.cubicBound(18)).isEqualTo(3306L);
	}

	@Test
	void islam2009Odd_n17() {
		// (n³ + 15n² + 14n − 6)/3 = (4913 + 4335 + 238 − 6)/3 = 9480/3 = 3160.
		assertThat(PanTrilinearAggregation.cubicBound(17)).isEqualTo(3160L);
	}

	// ---------------------------------------------------------------
	// SZ 2025 STRICTLY dominates HS 1982 for all valid n.
	// ---------------------------------------------------------------
	@Test
	void sz2025StrictlyDominatesHadasSchwartz1982() {
		for (int n = 2; n <= 60; n += 2) {
			if (n == 16) continue;
			long hs = PanTrilinearAggregation.panHadasSchwartz1982Bound(n);
			long sz = PanTrilinearAggregation.schwartzZwecher2025Bound(n);
			assertThat(sz).as("SZ2025 < HS1982 at n=" + n).isLessThan(hs);
		}
	}

	// ---------------------------------------------------------------
	// bestPanTaBound — assert MIN matches the hand-computed table.
	// Table covers each n ∈ {4, 6, 8, 9, 12, 17, 18, 22, 28, 32, 44}.
	// At n=18 the winner is Islam even (3306) — sharper than the
	// "newer" Pan-family formulas. At n=28/32 SZ 2025 starts winning
	// over Islam even.
	// ---------------------------------------------------------------
	@Test
	void bestPanTaBound_smallEven_pan1980Wins() {
		// n=4: Pan 1980 = 62 vs Islam = 100. MIN = 62.
		assertThat(PanTrilinearAggregation.bestPanTaBound(4)).isEqualTo(62L);
		// n=6: Pan 1980 = 180 vs Islam = 238. MIN = 180.
		assertThat(PanTrilinearAggregation.bestPanTaBound(6)).isEqualTo(180L);
		// n=8: Pan 1980 = 388 vs Islam = 456. MIN = 388.
		assertThat(PanTrilinearAggregation.bestPanTaBound(8)).isEqualTo(388L);
		// n=12: Pan 1980 = 1170 vs Islam = 1196. MIN = 1170.
		assertThat(PanTrilinearAggregation.bestPanTaBound(12)).isEqualTo(1170L);
	}

	@Test
	void bestPanTaBound_odd_islamWins() {
		// Only Islam's odd formula applies — bestPanTaBound = cubicBound.
		assertThat(PanTrilinearAggregation.bestPanTaBound(9)).isEqualTo(688L);
		assertThat(PanTrilinearAggregation.bestPanTaBound(17)).isEqualTo(3160L);
	}

	@Test
	void bestPanTaBound_n18_islamEvenWins() {
		// Per the literature finding — Islam even at n=18 beats Pan 1980,
		// HS 1982, AND SZ 2025.
		assertThat(PanTrilinearAggregation.bestPanTaBound(18)).isEqualTo(3306L);
	}

	@Test
	void bestPanTaBound_n22_islamEvenWins() {
		// Islam=5566, Pan80=6380, HS=5608, SZ=5596 → MIN=5566.
		assertThat(PanTrilinearAggregation.bestPanTaBound(22)).isEqualTo(5566L);
	}

	@Test
	void bestPanTaBound_n28_sz2025Wins() {
		// Crossover region: Islam=10556, SZ=10550. MIN=10550.
		assertThat(PanTrilinearAggregation.bestPanTaBound(28)).isEqualTo(10550L);
	}

	@Test
	void bestPanTaBound_n32_sz2025Wins() {
		// Islam=15136, SZ=15096 → MIN=15096.
		assertThat(PanTrilinearAggregation.bestPanTaBound(32)).isEqualTo(15096L);
	}

	@Test
	void bestPanTaBound_n44_sz2025Wins() {
		// Islam=36300, SZ=36110 → MIN=36110.
		assertThat(PanTrilinearAggregation.bestPanTaBound(44)).isEqualTo(36110L);
	}
}
