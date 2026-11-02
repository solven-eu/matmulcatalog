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

	/**
	 * The legacy sparse layout has a second spelling upstream: per-row pair lists
	 * {@code [[pos, coef], …]} (mizopl / we*star) next to Perminov's per-row objects
	 * {@code {"pos": coef}}. Until 2026-11-02 the pair list parsed as an empty object — an
	 * all-zero scheme that failed the import gate as if the contributor's file were wrong
	 * (seven valid ZT records). Both spellings must read as the same scheme, and a row of
	 * any other shape must be refused, never zero-filled.
	 */
	@Test
	public void pair_list_sparse_rows_read_as_the_same_scheme_and_unknown_rows_are_refused() throws Exception {
		NonCubicBilinearAlgorithm strassen = eu.solven.matmul.catalog.SchemeIO.readBilinear(new java.io.File(
				"src/main/resources/schemes/known/section2/2x2x2-r7-perminov_cr15_cn24_ZT_reduced-e498eb7.json"));
		tools.jackson.databind.json.JsonMapper mapper = tools.jackson.databind.json.JsonMapper.builder().build();
		tools.jackson.databind.node.ObjectNode root = mapper.createObjectNode();
		root.putArray("n").add(strassen.n).add(strassen.m).add(strassen.p);
		root.put("m", strassen.r);
		double[][] U = strassen.denseU(), V = strassen.denseV(), W = strassen.denseW();
		tools.jackson.databind.node.ArrayNode us = root.putArray("u_sparse"), vs = root.putArray("v_sparse"),
				ws = root.putArray("w_sparse");
		for (int k = 0; k < strassen.r; k++) {
			tools.jackson.databind.node.ArrayNode ru = us.addArray(), rv = vs.addArray(), rw = ws.addArray();
			for (int v = 0; v < U.length; v++) if (U[v][k] != 0) ru.addArray().add(v).add((int) U[v][k]);
			for (int v = 0; v < V.length; v++) if (V[v][k] != 0) rv.addArray().add(v).add((int) V[v][k]);
			// W positions are Perminov's column-major C order (j·n1 + i), as in the object spelling.
			for (int i = 0; i < strassen.n; i++) for (int j = 0; j < strassen.p; j++) {
				double c = W[i * strassen.p + j][k];
				if (c != 0) rw.addArray().add(j * strassen.n + i).add((int) c);
			}
		}

		NonCubicBilinearAlgorithm back = eu.solven.matmul.catalog.SchemeIO.read(root);
		assertThat(eu.solven.matmul.verifiers.Verifier.isExactNonCubic(back)).isTrue();
		assertThat(eu.solven.matmul.catalog.SchemeIO.contentHash(back))
				.isEqualTo(eu.solven.matmul.catalog.SchemeIO.contentHash(strassen));

		// A row that is neither an object nor a pair list is an error, not a row of zeros.
		((tools.jackson.databind.node.ArrayNode) root.get("u_sparse")).set(0, mapper.createArrayNode().add("0:1"));
		org.assertj.core.api.Assertions.assertThatThrownBy(() -> eu.solven.matmul.catalog.SchemeIO.read(root))
				.isInstanceOf(IOException.class).hasMessageContaining("pairs");
	}
}
