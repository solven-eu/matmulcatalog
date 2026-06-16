package eu.solven.matmul.commutative;

import static org.assertj.core.api.Assertions.assertThat;

import org.junit.jupiter.api.Test;

import eu.solven.matmul.papers.waksman1970.WaksmanBound;

/**
 * Verifies {@link WaksmanBound#bestCubic} reproduces the
 * "Waksman" cells in DIS09 Table 4.
 */
public class TestWaksmanBound {

	@Test
	public void waksman_cubic_matches_dis09_table4() {
		// n=4 even: 4·(16+4+4-1)/2 = 4·23/2 = 46
		assertThat(WaksmanBound.bestCubic(4)).isEqualTo(46L);
		// n=5 odd: (5-1)·(25+5+5-1)/2 + 25 = 4·34/2 + 25 = 68 + 25 = 93
		assertThat(WaksmanBound.bestCubic(5)).isEqualTo(93L);
		// n=6 even: 6·(36+6+6-1)/2 = 6·47/2 = 141
		assertThat(WaksmanBound.bestCubic(6)).isEqualTo(141L);
		// n=7 odd: 6·(49+7+7-1)/2 + 49 = 6·62/2 + 49 = 186 + 49 = 235
		assertThat(WaksmanBound.bestCubic(7)).isEqualTo(235L);
		// n=8 even: 8·(64+8+8-1)/2 = 8·79/2 = 316
		assertThat(WaksmanBound.bestCubic(8)).isEqualTo(316L);
		// n=10 even: 10·(100+10+10-1)/2 = 10·119/2 = 595
		assertThat(WaksmanBound.bestCubic(10)).isEqualTo(595L);
	}

	@Test
	public void waksman_handles_non_cubic_shapes() {
		// ⟨3,3,3⟩: (3-1)·(9+3+3-1)/2 + 9 = 14 + 9 = 23
		assertThat(WaksmanBound.forShape(3, 3, 3)).isEqualTo(23L);
		// ⟨2,2,2⟩: 2·(4+2+2-1)/2 = 7
		assertThat(WaksmanBound.forShape(2, 2, 2)).isEqualTo(7L);
	}
}
