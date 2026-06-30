package eu.solven.matmul.catalog;

import static org.assertj.core.api.Assertions.assertThat;

import java.io.File;
import java.nio.file.Files;
import java.nio.file.Path;

import org.junit.jupiter.api.Test;

import eu.solven.matmul.NonCubicBilinearAlgorithm;
import eu.solven.matmul.verifiers.Verifier;
import eu.solven.matmul.papers.dis2009.PanTrilinearAggregation;

/**
 * Coefficients must be stored <strong>exactly</strong>: a rational like {@code 1/3}
 * is written as the string {@code "1/3"}, never the rounded decimal
 * {@code 0.3333333333333333}. {@link PanTrilinearAggregation#build(int)} for n=4
 * carries {@code 1/(n/2+1) = 1/3}, so it is a natural probe.
 */
public class TestRationalCoef {

	@Test
	public void rational_coeffs_round_trip_exactly() throws Exception {
		NonCubicBilinearAlgorithm alg = PanTrilinearAggregation.build(4); // ⟨4,4,4⟩=100, has 1/3
		Path tmp = Files.createTempFile("rationalcoef", ".json");
		try {
			SchemeIO.write(alg, tmp.toFile(), null);
			String json = Files.readString(tmp);

			// Exact fraction strings present; NO rounded 1/3 decimal anywhere.
			assertThat(json).contains("1/3");
			assertThat(json).doesNotContain("0.3333");
			assertThat(json).doesNotContain("0.6666");
			assertThat(json).doesNotContain("0.33333333");

			// Re-read reproduces a verifying ⟨4,4,4⟩=100 scheme.
			NonCubicBilinearAlgorithm back = SchemeIO.readBilinear(tmp.toFile());
			assertThat(back.n).isEqualTo(4);
			assertThat(back.r).isEqualTo(100);
			assertThat(Verifier.passesRandomMatmulSpotCheck(back)).isTrue();

			// Coefficients equal the TRUE rationals (not bit-exact to the build's
			// double: the constructor's 1−1/3 yields 0.666…667, one ULP off 2/3;
			// rationalisation correctly stores the exact 2/3 → 0.666…666 on re-read).
			double[][] u0 = alg.denseU(), u1 = back.denseU();
			for (int i = 0; i < u0.length; i++) {
				for (int k = 0; k < u0[i].length; k++) {
					assertThat(u1[i][k]).isCloseTo(u0[i][k], org.assertj.core.api.Assertions.within(1e-9));
				}
			}

			// Idempotence: once stored as exact rationals, re-emitting is byte-stable
			// (no drift) — the load→store cycle has reached its fixed point.
			Path tmp2 = Files.createTempFile("rationalcoef2", ".json");
			try {
				SchemeIO.write(back, tmp2.toFile(), null);
				assertThat(Files.readString(tmp2)).isEqualTo(json);
			} finally {
				Files.deleteIfExists(tmp2);
			}
		} finally {
			Files.deleteIfExists(tmp);
		}
	}
}
