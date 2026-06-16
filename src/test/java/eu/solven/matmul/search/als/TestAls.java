package eu.solven.matmul.search.als;

import static org.assertj.core.api.Assertions.assertThat;
import static org.assertj.core.api.Assertions.within;

import java.util.SplittableRandom;

import org.junit.jupiter.api.Test;

import eu.solven.matmul.search.als.Als;
import eu.solven.matmul.BilinearAlgorithm;
import eu.solven.matmul.papers.laderman1976.Laderman23;
import eu.solven.matmul.LinAlg;
import eu.solven.matmul.papers.strassen1969.Strassen7;
import eu.solven.matmul.Verifier;

public class TestAls {

	@Test
	public void linAlgInvertIdentity() {
		double[][] I = { { 1, 0, 0 }, { 0, 1, 0 }, { 0, 0, 1 } };
		double[][] inv = LinAlg.invert(I);
		for (int i = 0; i < 3; i++) {
			for (int j = 0; j < 3; j++) {
				assertThat(inv[i][j]).isEqualTo(I[i][j]);
			}
		}
	}

	@Test
	public void linAlgInvertRoundTrip() {
		double[][] A = { { 4, 3, 2 }, { 1, 5, 7 }, { 2, 6, 9 } };
		double[][] inv = LinAlg.invert(A);
		double[][] prod = LinAlg.mul(A, inv);
		for (int i = 0; i < 3; i++) {
			for (int j = 0; j < 3; j++) {
				double expected = (i == j) ? 1.0 : 0.0;
				assertThat(prod[i][j]).isCloseTo(expected, within(1e-10));
			}
		}
	}

	/**
	 * Implementation correctness: starting from Strassen ± tiny noise, ALS at r=7
	 * snaps right back to a residual-0 decomposition. This isolates the ALS update
	 * rule from the (genuinely hard) global optimization problem.
	 */
	@Test
	public void alsConvergesNearStrassen() {
		BilinearAlgorithm s = Strassen7.get();
		SplittableRandom rng = new SplittableRandom(123L);
		double[][] U = perturb(s.U, rng, 0.01);
		double[][] V = perturb(s.V, rng, 0.01);
		double[][] W = perturb(s.W, rng, 0.01);
		Als als = new Als(2, 7);
		Als.Result result = als.fitFrom(U, V, W, 1000, 1e-10);
		assertThat(result.converged).isTrue();
		assertThat(result.residual).isLessThan(1e-10);
		assertThat(Verifier.isExact(result.algorithm)).isTrue();
	}

	/** Same correctness check for Laderman at r=23. */
	@Test
	public void alsConvergesNearLaderman() {
		BilinearAlgorithm l = Laderman23.get();
		SplittableRandom rng = new SplittableRandom(456L);
		double[][] U = perturb(l.U, rng, 0.01);
		double[][] V = perturb(l.V, rng, 0.01);
		double[][] W = perturb(l.W, rng, 0.01);
		Als als = new Als(3, 23);
		Als.Result result = als.fitFrom(U, V, W, 2000, 1e-10);
		assertThat(result.converged).isTrue();
		assertThat(result.residual).isLessThan(1e-8);
	}

	/**
	 * Random-init ALS at over-rank for 2×2 — known to suffer the ALS swamp.
	 * With a handful of restarts we usually find at least one good basin reaching
	 * residual ≪ 1, far below the initial ‖T‖ = √8 ≈ 2.83.
	 */
	@Test
	public void alsRandomInitMakesProgress2x2() {
		Als als = new Als(2, 8);
		Als.Result result = als.fitWithRestarts(20, 500, 1e-8, 0L);
		assertThat(result.residual).isLessThan(0.1);
	}

	private static double[][] perturb(double[][] M, SplittableRandom rng, double sigma) {
		double[][] out = new double[M.length][M[0].length];
		for (int i = 0; i < M.length; i++) {
			for (int j = 0; j < M[0].length; j++) {
				out[i][j] = M[i][j] + sigma * rng.nextGaussian();
			}
		}
		return out;
	}
}
