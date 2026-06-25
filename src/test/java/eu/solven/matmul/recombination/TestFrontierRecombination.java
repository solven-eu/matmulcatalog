package eu.solven.matmul.recombination;

import static org.assertj.core.api.Assertions.assertThat;

import java.io.File;
import java.util.List;

import org.junit.jupiter.api.Test;

import eu.solven.matmul.NonCubicBilinearAlgorithm;
import eu.solven.matmul.Verifier;
import eu.solven.matmul.algebra.Algebra;
import eu.solven.matmul.algebra.Field;
import eu.solven.matmul.catalog.SchemeIO;
import eu.solven.matmul.recombination.Recombination.SotaResolver;

/**
 * b-2a guard: the frontier-recombination engine reproduces the orbit×orientation wins (⟨2,9,10⟩→160,
 * ⟨3,10,10⟩→275 via the "3"-axis re-orientation), and its winner materialises to a valid matmul.
 * SOTA-or-better assertions (≤) so a real improvement never breaks the test.
 */
public class TestFrontierRecombination {

	private static final String AT_223 =
			"src/main/resources/schemes/known/section3/2x2x3-r11-alphatensor_Z-682e003.json";

	@Test
	public void orientation_win_2x9x10() {
		assertWin(2, 9, 10, 160);
	}

	@Test
	public void orientation_win_3x10x10() {
		assertWin(3, 10, 10, 275);
	}

	private void assertWin(int N, int M, int P, long expected) {
		NonCubicBilinearAlgorithm base = read(AT_223);
		SotaResolver sota = Recombination.catalogResolver(Algebra.nonCommutative(Field.R));
		FrontierRecombination.Best best = FrontierRecombination.bestFor(N, M, P, List.of(base), sota);
		assertThat(best).isNotNull();
		// SOTA-or-better: the engine must reach at least the known orbit win.
		assertThat(best.rank()).as("⟨%d,%d,%d⟩ frontier rank", N, M, P).isLessThanOrEqualTo(expected);

		// the winning scheme rebuilds and still computes the base's matmul
		NonCubicBilinearAlgorithm winner = best.materialiseWinner();
		assertThat(winner.r).isEqualTo(base.r);
		assertThat(Verifier.isExactNonCubic(winner)).as("winner computes matmul").isTrue();
	}

	private static NonCubicBilinearAlgorithm read(String p) {
		try {
			File f = new File(p);
			assertThat(f).exists();
			return SchemeIO.read(f);
		} catch (Exception e) {
			throw new AssertionError(e);
		}
	}
}
