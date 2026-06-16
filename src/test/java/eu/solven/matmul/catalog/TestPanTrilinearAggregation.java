package eu.solven.matmul.catalog;

import static org.assertj.core.api.Assertions.assertThat;

import java.util.Map;

import org.junit.jupiter.api.Test;

import eu.solven.matmul.papers.dis2009.PanTrilinearAggregation;

/**
 * Verifies {@link PanTrilinearAggregation#cubicBound} reproduces the
 * "TA" column of DIS09 Table 3 exactly.
 */
public class TestPanTrilinearAggregation {

	/** DIS09 Table 3 rows attributed to "TA" — n → reported bound. */
	private static final Map<Integer, Long> DIS09_TA = Map.ofEntries(
			Map.entry(18, 3306L),
			Map.entry(20, 4340L),
			Map.entry(22, 5566L),
			Map.entry(23, 6806L),
			Map.entry(24, 7000L),
			Map.entry(25, 8448L),
			Map.entry(26, 8658L),
			Map.entry(27, 10330L),
			Map.entry(28, 10556L),
			Map.entry(29, 12468L),
			Map.entry(30, 12710L));

	@Test
	public void reproduces_dis09_ta_column_exactly() {
		for (Map.Entry<Integer, Long> e : DIS09_TA.entrySet()) {
			int n = e.getKey();
			long expected = e.getValue();
			long got = PanTrilinearAggregation.cubicBound(n);
			assertThat(got).as("Pan TA at n=" + n).isEqualTo(expected);
		}
	}

	@Test
	public void parity_branches_documented() {
		assertThat(PanTrilinearAggregation.branchLabel(18)).contains("even");
		assertThat(PanTrilinearAggregation.branchLabel(23)).contains("odd");
	}
}
