package eu.solven.matmul.catalog;

import static org.assertj.core.api.Assertions.assertThat;

import java.io.IOException;
import java.nio.file.Files;
import java.nio.file.Path;

import org.junit.jupiter.api.Test;

import eu.solven.matmul.NonCubicBilinearAlgorithm;

/**
 * Sanity test for the row-oriented {@code u_sparse}/{@code v_sparse}/{@code w_sparse}
 * format introduced 2026-06-03 (task #174). Confirms:
 *
 * <ul>
 *   <li>The legacy array-of-per-key-map format still reads correctly.</li>
 *   <li>The new writer emits the row-oriented map shape
 *       {@code {"k":{"i":[...],"c":[...]}}}.</li>
 *   <li>Round-trip (read legacy → write new → re-read) preserves U/V/W exactly.</li>
 * </ul>
 */
public class TestSparseFactorRoundTrip {

	private static final Path LEGACY =
			Path.of("src/main/resources/schemes/known/section4/4x4x4-r48-dumas_pernet_sedoglavic_2025-929db4e.json");

	@Test
	public void roundTripLegacyToNew() throws IOException {
		String legacyJson = Files.readString(LEGACY);
		NonCubicBilinearAlgorithm original = SchemeIO.read(legacyJson);

		String newJson = SchemeIO.toJsonSparse(original);
		// New format: outer container of u_sparse is an OBJECT, not an array.
		assertThat(newJson).contains("\"u_sparse\": {");
		assertThat(newJson).contains("\"v_sparse\": {");
		assertThat(newJson).contains("\"w_sparse\": {");
		// Per-product entries carry parallel i/c arrays.
		assertThat(newJson).contains("\"i\":");
		assertThat(newJson).contains("\"c\":");

		NonCubicBilinearAlgorithm reread = SchemeIO.read(newJson);
		assertThat(reread.n).isEqualTo(original.n);
		assertThat(reread.m).isEqualTo(original.m);
		assertThat(reread.p).isEqualTo(original.p);
		assertThat(reread.r).isEqualTo(original.r);
		double[][] origU = original.denseU();
		double[][] origV = original.denseV();
		double[][] origW = original.denseW();
		double[][] rereadU = reread.denseU();
		double[][] rereadV = reread.denseV();
		double[][] rereadW = reread.denseW();
		assertThat(rereadU).isDeepEqualTo(origU);
		assertThat(rereadV).isDeepEqualTo(origV);
		assertThat(rereadW).isDeepEqualTo(origW);
	}
}
