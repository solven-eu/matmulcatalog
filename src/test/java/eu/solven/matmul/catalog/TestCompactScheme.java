package eu.solven.matmul.catalog;

import static org.assertj.core.api.Assertions.assertThat;

import java.io.File;

import org.junit.jupiter.api.Test;

import eu.solven.matmul.NonCubicBilinearAlgorithm;
import eu.solven.matmul.verifiers.Verifier;

public class TestCompactScheme {

	@Test
	public void strassen_round_trip_identical() throws Exception {
		NonCubicBilinearAlgorithm orig = SchemeIO.read(
				eu.solven.matmul.catalog.SchemeResolver.byHint("src/main/resources/schemes/known/section2/strassen-2x2x2_m7_a18.json"));
		CompactScheme cs = CompactScheme.of(orig);
		NonCubicBilinearAlgorithm round = cs.expand();
		assertThat(round.n).isEqualTo(orig.n);
		assertThat(round.r).isEqualTo(orig.r);
		double[][] origU = orig.denseU();
		double[][] origV = orig.denseV();
		double[][] origW = orig.denseW();
		double[][] roundU = round.denseU();
		double[][] roundV = round.denseV();
		double[][] roundW = round.denseW();
		// All entries equal (Strassen is ternary).
		for (int i = 0; i < origU.length; i++) {
			assertThat(roundU[i]).containsExactly(origU[i]);
		}
		for (int i = 0; i < origV.length; i++) {
			assertThat(roundV[i]).containsExactly(origV[i]);
		}
		for (int i = 0; i < origW.length; i++) {
			assertThat(roundW[i]).containsExactly(origW[i]);
		}
		assertThat(Verifier.passesRandomMatmulSpotCheck(round)).isTrue();
	}

	@Test
	public void byte_size_reports_reasonable_compression() throws Exception {
		NonCubicBilinearAlgorithm orig = SchemeIO.read(
				eu.solven.matmul.catalog.SchemeResolver.byHint("src/main/resources/schemes/known/section2/strassen-2x2x2_m7_a18.json"));
		CompactScheme cs = CompactScheme.of(orig);
		long denseBytes = (long) (orig.n * orig.m + orig.m * orig.p + orig.n * orig.p) * orig.r * 8L;
		// Compact representation should be much smaller.
		assertThat(cs.byteSize()).isLessThan(denseBytes);
	}

	@Test
	public void perminov_q_rational_round_trips_via_dict_encoding() throws Exception {
		// Find any Perminov Q-rational scheme (has values like ±1/2, ±1/4).
		java.io.File[] qFiles;
		try (var stream = java.nio.file.Files.walk(java.nio.file.Path.of("src/main/resources/schemes"))) {
			qFiles = stream.filter(p -> p.getFileName().toString().startsWith("perminov-Q_"))
					.limit(1).map(java.nio.file.Path::toFile).toArray(java.io.File[]::new);
		}
		if (qFiles.length == 0) return;  // skip if no Q schemes
		java.io.File qFile = qFiles[0];
		NonCubicBilinearAlgorithm orig = SchemeIO.read(qFile);
		CompactScheme cs = CompactScheme.of(orig);
		NonCubicBilinearAlgorithm round = cs.expand();
		double[][] origU = orig.denseU();
		double[][] origV = orig.denseV();
		double[][] origW = orig.denseW();
		double[][] roundU = round.denseU();
		double[][] roundV = round.denseV();
		double[][] roundW = round.denseW();
		// Round-trip must preserve every U/V/W entry exactly.
		for (int i = 0; i < origU.length; i++) assertThat(roundU[i]).containsExactly(origU[i]);
		for (int i = 0; i < origV.length; i++) assertThat(roundV[i]).containsExactly(origV[i]);
		for (int i = 0; i < origW.length; i++) assertThat(roundW[i]).containsExactly(origW[i]);
		assertThat(Verifier.passesRandomMatmulSpotCheck(round)).isTrue();
		// Compact form should be much smaller than dense doubles.
		long denseBytes = (long) (orig.n * orig.m + orig.m * orig.p + orig.n * orig.p) * orig.r * 8L;
		System.out.printf("Q-rational ⟨%d,%d,%d⟩=r%d  dense=%d B  compact=%d B  (%.1f%% of dense)%n",
				orig.n, orig.m, orig.p, orig.r, denseBytes, cs.byteSize(),
				100.0 * cs.byteSize() / denseBytes);
		assertThat(cs.byteSize()).isLessThan(denseBytes);
	}

	@Test
	public void laderman_round_trip_identical() throws Exception {
		NonCubicBilinearAlgorithm orig = SchemeIO.read(
				eu.solven.matmul.catalog.SchemeResolver.byHint("src/main/resources/schemes/known/section3/laderman_1976-3x3x3_m23_a98.json"));
		CompactScheme cs = CompactScheme.of(orig);
		NonCubicBilinearAlgorithm round = cs.expand();
		assertThat(round.r).isEqualTo(orig.r);
		double[][] origU = orig.denseU();
		double[][] origV = orig.denseV();
		double[][] origW = orig.denseW();
		double[][] roundU = round.denseU();
		double[][] roundV = round.denseV();
		double[][] roundW = round.denseW();
		for (int i = 0; i < origU.length; i++) assertThat(roundU[i]).containsExactly(origU[i]);
		for (int i = 0; i < origV.length; i++) assertThat(roundV[i]).containsExactly(origV[i]);
		for (int i = 0; i < origW.length; i++) assertThat(roundW[i]).containsExactly(origW[i]);
		assertThat(Verifier.passesRandomMatmulSpotCheck(round)).isTrue();
	}
}
