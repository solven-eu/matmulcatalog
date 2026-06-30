package eu.solven.matmul;

import static org.assertj.core.api.Assertions.assertThat;

import org.junit.jupiter.api.Test;

import eu.solven.matmul.AxisSplitBases;
import eu.solven.matmul.NonCubicBilinearAlgorithm;
import eu.solven.matmul.verifiers.Verifier;

public class TestAxisSplitBases {

	@Test
	public void mul211_verifies() {
		NonCubicBilinearAlgorithm alg = AxisSplitBases.mul211();
		assertThat(alg.n).isEqualTo(2);
		assertThat(alg.m).isEqualTo(1);
		assertThat(alg.p).isEqualTo(1);
		assertThat(alg.r).isEqualTo(2);
		assertThat(Verifier.isExactNonCubic(alg)).isTrue();
	}

	@Test
	public void mul121_verifies() {
		NonCubicBilinearAlgorithm alg = AxisSplitBases.mul121();
		assertThat(alg.n).isEqualTo(1);
		assertThat(alg.m).isEqualTo(2);
		assertThat(alg.p).isEqualTo(1);
		assertThat(alg.r).isEqualTo(2);
		assertThat(Verifier.isExactNonCubic(alg)).isTrue();
	}

	@Test
	public void mul112_verifies() {
		NonCubicBilinearAlgorithm alg = AxisSplitBases.mul112();
		assertThat(alg.n).isEqualTo(1);
		assertThat(alg.m).isEqualTo(1);
		assertThat(alg.p).isEqualTo(2);
		assertThat(alg.r).isEqualTo(2);
		assertThat(Verifier.isExactNonCubic(alg)).isTrue();
	}
}
