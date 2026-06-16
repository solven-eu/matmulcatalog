package eu.solven.matmul.papers.hopcroftkerr1971;

import eu.solven.matmul.papers.hopcroftkerr1971.LemmaOneAugmentation;

import static org.assertj.core.api.Assertions.assertThat;

import org.junit.jupiter.api.Test;

public class TestLemmaOneAugmentation {

	@Test
	public void identity_when_n_equals_p() {
		for (int p = 1; p <= 8; p++) {
			int[][] M = LemmaOneAugmentation.build(p, p);
			assertThat(M.length).isEqualTo(p);
			for (int i = 0; i < p; i++) {
				for (int j = 0; j < p; j++) {
					assertThat(M[i][j]).isEqualTo(i == j ? 1 : 0);
				}
			}
		}
	}

	@Test
	public void identity_prefix_preserved_for_all_n_p() {
		for (int p : new int[] { 2, 3, 5, 10, 12 }) {
			for (int extra = 1; extra <= 6; extra++) {
				int n = p + extra;
				int[][] M = LemmaOneAugmentation.build(p, n);
				assertThat(M.length).as("rows for p=" + p + ", n=" + n).isEqualTo(n);
				for (int i = 0; i < p; i++) {
					for (int j = 0; j < p; j++) {
						assertThat(M[i][j])
								.as("identity row " + i + " col " + j)
								.isEqualTo(i == j ? 1 : 0);
					}
				}
			}
		}
	}

	@Test
	public void every_cyclic_window_invertible_p2_n5() {
		int p = 2, n = 5;
		int[][] M = LemmaOneAugmentation.build(p, n);
		for (int start = 0; start < n; start++) {
			int[][] window = new int[p][p];
			for (int r = 0; r < p; r++) {
				int src = (start + r) % n;
				System.arraycopy(M[src], 0, window[r], 0, p);
			}
			long det = LemmaOneAugmentation.detMatrix(window);
			assertThat(det).as("cyclic window starting at row " + start).isNotEqualTo(0);
		}
	}

	@Test
	public void every_cyclic_window_invertible_for_our_target_shapes() {
		// Asymmetric HK target shapes for the 3 narrow FMM gaps:
		// ⟨10, 2, 15⟩, ⟨10, 2, 16⟩, ⟨12, 2, 16⟩
		int[][] cases = { {10, 15}, {10, 16}, {12, 16} };
		for (int[] c : cases) {
			int p = c[0], n = c[1];
			int[][] M = LemmaOneAugmentation.build(p, n);
			for (int start = 0; start < n; start++) {
				int[][] window = new int[p][p];
				for (int r = 0; r < p; r++) {
					int src = (start + r) % n;
					System.arraycopy(M[src], 0, window[r], 0, p);
				}
				long det = LemmaOneAugmentation.detMatrix(window);
				assertThat(det)
						.as("p=" + p + ", n=" + n + ", window-start=" + start)
						.isNotEqualTo(0);
			}
		}
	}

	@Test
	public void rejects_bad_dimensions() {
		org.junit.jupiter.api.Assertions.assertThrows(IllegalArgumentException.class,
				() -> LemmaOneAugmentation.build(0, 1));
		org.junit.jupiter.api.Assertions.assertThrows(IllegalArgumentException.class,
				() -> LemmaOneAugmentation.build(5, 4));   // n < p
	}
}
