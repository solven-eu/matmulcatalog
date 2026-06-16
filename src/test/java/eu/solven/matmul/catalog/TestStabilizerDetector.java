package eu.solven.matmul.catalog;

import static org.assertj.core.api.Assertions.assertThat;

import org.junit.jupiter.api.Test;

import eu.solven.matmul.papers.laderman1976.Laderman23;
import eu.solven.matmul.NonCubicBilinearAlgorithm;
import eu.solven.matmul.papers.strassen1969.Strassen7;
import eu.solven.matmul.catalog.StabilizerDetector;

/**
 * Validates the {@link StabilizerDetector} on the two algorithms whose
 * stabilizers are independently known:
 *
 * <ul>
 *   <li>Strassen `⟨2,2,2⟩=7` — has Z/3 cyclic stabilizer (de Groote 1978).</li>
 *   <li>Laderman `⟨3,3,3⟩=23` — no Z/3 (the multiplications break the
 *       trilinear cyclic symmetry; widely-stated fact in the matmul-rank
 *       literature).</li>
 * </ul>
 */
public class TestStabilizerDetector {

	@Test
	public void strassen_has_z3_stabilizer() {
		NonCubicBilinearAlgorithm alg = NonCubicBilinearAlgorithm.fromCubic(Strassen7.get());
		assertThat(StabilizerDetector.detectZ3(alg)).isEqualTo("Z/3");
	}

	@Test
	public void laderman_no_z3() {
		NonCubicBilinearAlgorithm alg = NonCubicBilinearAlgorithm.fromCubic(Laderman23.get());
		assertThat(StabilizerDetector.detectZ3(alg)).isEqualTo("trivial");
	}

	@Test
	public void non_cubic_returns_na() {
		// Build a fake non-cubic ⟨2,3,3⟩ (use a trivial dummy).
		double[][] U = new double[6][1]; U[0][0] = 1.0;
		double[][] V = new double[9][1]; V[0][0] = 1.0;
		double[][] W = new double[6][1]; W[0][0] = 1.0;
		NonCubicBilinearAlgorithm alg = new NonCubicBilinearAlgorithm(2, 3, 3, U, V, W);
		assertThat(StabilizerDetector.detectZ3(alg)).isEqualTo("n/a");
	}
}
