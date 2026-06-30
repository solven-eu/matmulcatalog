package eu.solven.matmul;

import static org.assertj.core.api.Assertions.assertThat;
import static org.assertj.core.api.Assertions.within;

import org.junit.jupiter.api.Test;

import eu.solven.matmul.BilinearAlgorithm;
import eu.solven.matmul.papers.laderman1976.Laderman23;
import eu.solven.matmul.verifiers.NaiveMatMul;
import eu.solven.matmul.papers.strassen1969.Strassen7;
import eu.solven.matmul.verifiers.Verifier;

public class TestVerifier {

	@Test
	public void naiveMatMul2x2() {
		BilinearAlgorithm naive = NaiveMatMul.of(2);
		assertThat(naive.r).isEqualTo(8);
		assertThat(Verifier.residual(naive)).isCloseTo(0.0, within(1e-12));
	}

	@Test
	public void naiveMatMul3x3() {
		BilinearAlgorithm naive = NaiveMatMul.of(3);
		assertThat(naive.r).isEqualTo(27);
		assertThat(Verifier.residual(naive)).isCloseTo(0.0, within(1e-12));
	}

	@Test
	public void strassen7() {
		BilinearAlgorithm strassen = Strassen7.get();
		assertThat(strassen.r).isEqualTo(7);
		assertThat(Verifier.residual(strassen)).isCloseTo(0.0, within(1e-12));
	}

	@Test
	public void laderman23() {
		BilinearAlgorithm laderman = Laderman23.get();
		assertThat(laderman.n).isEqualTo(3);
		assertThat(laderman.r).isEqualTo(23);
		assertThat(Verifier.residual(laderman)).isCloseTo(0.0, within(1e-12));
	}

	@Test
	public void strassenTransposesToTrilin() {
		BilinearAlgorithm s = Strassen7.get();
		assertThat(Verifier.residual(s)).isCloseTo(0.0, within(1e-12));
		BilinearAlgorithm sTrilin = Verifier.transposeW(s);
		double[][][] trilin = Verifier.trilinTensor(2);
		assertThat(Verifier.residualAgainst(sTrilin, trilin)).isCloseTo(0.0, within(1e-12));
		// And the conversion is self-inverse: transpose twice = original.
		BilinearAlgorithm sBack = Verifier.transposeW(sTrilin);
		assertThat(Verifier.residual(sBack)).isCloseTo(0.0, within(1e-12));
	}

	@Test
	public void ladermanTransposesToTrilin() {
		BilinearAlgorithm l = Laderman23.get();
		assertThat(Verifier.residual(l)).isCloseTo(0.0, within(1e-12));
		BilinearAlgorithm lTrilin = Verifier.transposeW(l);
		double[][][] trilin = Verifier.trilinTensor(3);
		assertThat(Verifier.residualAgainst(lTrilin, trilin)).isCloseTo(0.0, within(1e-12));
	}

	@Test
	public void perturbingStrassenBreaksIt() {
		BilinearAlgorithm strassen = Strassen7.get();
		strassen.U[0][0] += 0.1; // corrupt one entry
		assertThat(Verifier.residual(strassen)).isGreaterThan(1e-3);
	}
}
