package eu.solven.matmul.papers.dis2009;

import static org.assertj.core.api.Assertions.assertThat;

import org.junit.jupiter.api.Test;

import eu.solven.matmul.NonCubicBilinearAlgorithm;
import eu.solven.matmul.verifiers.Verifier;

class TestPanTrilinearAggregationBuild {

	@Test
	void evenN4() {
		NonCubicBilinearAlgorithm alg = PanTrilinearAggregation.build(4);
		assertThat(alg.n).isEqualTo(4);
		assertThat(alg.m).isEqualTo(4);
		assertThat(alg.p).isEqualTo(4);
		assertThat(alg.r).isEqualTo((int) PanTrilinearAggregation.cubicBound(4));
		assertThat(alg.r).isEqualTo(100);
		assertThat(Verifier.isExactNonCubic(alg)).isTrue();
	}

	@Test
	void evenN6() {
		NonCubicBilinearAlgorithm alg = PanTrilinearAggregation.build(6);
		assertThat(alg.r).isEqualTo((int) PanTrilinearAggregation.cubicBound(6));
		assertThat(alg.r).isEqualTo(238);
		assertThat(Verifier.isExactNonCubic(alg)).isTrue();
	}

	@Test
	void oddN3() {
		NonCubicBilinearAlgorithm alg = PanTrilinearAggregation.build(3);
		assertThat(alg.r).isEqualTo((int) PanTrilinearAggregation.cubicBound(3));
		assertThat(alg.r).isEqualTo(66);
		assertThat(Verifier.isExactNonCubic(alg)).isTrue();
	}

	@Test
	void oddN5() {
		NonCubicBilinearAlgorithm alg = PanTrilinearAggregation.build(5);
		assertThat(alg.r).isEqualTo((int) PanTrilinearAggregation.cubicBound(5));
		assertThat(alg.r).isEqualTo(188);
		assertThat(Verifier.isExactNonCubic(alg)).isTrue();
	}
}
