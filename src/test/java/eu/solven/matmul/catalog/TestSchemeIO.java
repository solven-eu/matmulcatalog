package eu.solven.matmul.catalog;

import static org.assertj.core.api.Assertions.assertThat;

import java.io.File;
import java.io.IOException;
import java.nio.file.Files;
import java.nio.file.Path;

import org.junit.jupiter.api.Test;
import org.junit.jupiter.api.io.TempDir;
import java.nio.file.Path;

import eu.solven.matmul.BilinearAlgorithm;
import eu.solven.matmul.papers.laderman1976.Laderman23;
import eu.solven.matmul.NonCubicBilinearAlgorithm;
import eu.solven.matmul.papers.strassen1969.Strassen7;
import eu.solven.matmul.verifiers.Verifier;
import eu.solven.matmul.catalog.SchemeIO;

/**
 * Round-trip tests for the dronperminov scheme JSON format.
 */
public class TestSchemeIO {

	@TempDir
	public Path tmp;

	@Test
	public void strassen_round_trip_preserves_factors_and_verifies() throws IOException {
		BilinearAlgorithm original = Strassen7.get();
		File f = tmp.resolve("strassen.json").toFile();
		SchemeIO.write(original, f);
		NonCubicBilinearAlgorithm loaded = SchemeIO.read(f);

		assertThat(loaded.n).isEqualTo(2);
		assertThat(loaded.m).isEqualTo(2);
		assertThat(loaded.p).isEqualTo(2);
		assertThat(loaded.r).isEqualTo(7);
		assertThat(loaded.isCubic()).isTrue();
		assertThat(Verifier.isExactNonCubic(loaded)).isTrue();
		assertThat(Verifier.isExact(loaded.asCubic())).isTrue();

		// Bytewise structural check: U/V/W match within float tolerance.
		double[][] loadedU = loaded.denseU();
		double[][] loadedV = loaded.denseV();
		double[][] loadedW = loaded.denseW();
		assertThat(loadedU).isDeepEqualTo(original.U);
		assertThat(loadedV).isDeepEqualTo(original.V);
		assertThat(loadedW).isDeepEqualTo(original.W);
	}

	@Test
	public void laderman_round_trip_preserves_factors_and_verifies() throws IOException {
		BilinearAlgorithm original = Laderman23.get();
		File f = tmp.resolve("laderman.json").toFile();
		SchemeIO.write(original, f);
		NonCubicBilinearAlgorithm loaded = SchemeIO.read(f);

		assertThat(loaded.n).isEqualTo(3);
		assertThat(loaded.r).isEqualTo(23);
		assertThat(Verifier.isExactNonCubic(loaded)).isTrue();
	}

	@Test
	public void canonical_strassen_json_from_dronperminov_loads_and_verifies() throws IOException {
		// Verbatim from https://github.com/dronperminov/FastMatrixMultiplication
		String json = "{\n"
				+ "  \"n\": [2, 2, 2],\n"
				+ "  \"m\": 7,\n"
				+ "  \"z2\": false,\n"
				+ "  \"u\": [[1,0,0,1], [0,1,0,-1], [-1,0,1,0], [1,1,0,0], [1,0,0,0], [0,0,0,1], [0,0,1,1]],\n"
				+ "  \"v\": [[1,0,0,1], [0,0,1,1], [1,1,0,0], [0,0,0,1], [0,1,0,-1], [-1,0,1,0], [1,0,0,0]],\n"
				+ "  \"w\": [[1,0,0,1], [1,0,0,0], [0,0,0,1], [-1,0,1,0], [0,0,1,1], [1,1,0,0], [0,1,0,-1]]\n"
				+ "}\n";
		NonCubicBilinearAlgorithm loaded = SchemeIO.read(json);
		assertThat(loaded.r).isEqualTo(7);
		assertThat(loaded.isCubic()).isTrue();
		// This is a different Strassen variant from our Strassen7 but must be
		// a valid ⟨2,2,2⟩ rank-7 decomposition.
		assertThat(Verifier.isExactNonCubic(loaded)).isTrue();
	}

	@Test
	public void write_emits_n_array_and_correct_rank() throws IOException {
		String json = SchemeIO.toJson(
				NonCubicBilinearAlgorithm.fromCubic(Strassen7.get()),
				null, null);
		assertThat(json).contains("\"n\": [2, 2, 2]");
		assertThat(json).contains("\"m\": 7");
		assertThat(json).doesNotContain("\"z2\"");
	}

	/**
	 * The retired {@code "z2"} boolean must never reappear in emitted JSON —
	 * field membership lives solely in {@code fields[]}/{@code fields_not[]}.
	 * Guards both the dense and sparse writers.
	 */
	@Test
	public void writers_never_emit_the_retired_z2_key() {
		NonCubicBilinearAlgorithm alg = NonCubicBilinearAlgorithm.fromCubic(Strassen7.get());
		assertThat(SchemeIO.toJson(alg, null, null)).doesNotContain("z2");
		assertThat(SchemeIO.toJsonSparse(alg)).doesNotContain("z2");
	}

	/**
	 * {@code isZ2} (F₂-only) is decided purely from {@code fields[]}: present-but-
	 * char-0 fields ⇒ not F₂-only; a lone {@code F2} ⇒ F₂-only; no {@code fields[]}
	 * (the old {@code z2} fallback is gone) ⇒ not F₂-only.
	 */
	@Test
	public void isZ2_reads_fields_not_the_legacy_flag() throws IOException {
		tools.jackson.databind.json.JsonMapper mapper = tools.jackson.databind.json.JsonMapper.builder().build();
		// Integer scheme valid over many fields incl. F2 — NOT "z2".
		assertThat(SchemeIO.isZ2(mapper.readTree("{\"fields\":[\"F2\",\"F3\",\"Z\",\"Q\",\"R\",\"C\"]}"))).isFalse();
		// Genuine GF(2)-only scheme.
		assertThat(SchemeIO.isZ2(mapper.readTree("{\"fields\":[\"F2\"]}"))).isTrue();
		// No fields[] and a stale z2:true must NOT be honoured (legacy flag retired).
		assertThat(SchemeIO.isZ2(mapper.readTree("{\"z2\":true}"))).isFalse();
	}

	@Test
	public void file_written_then_re_read_byte_for_byte_factors_match() throws IOException {
		Path file = tmp.resolve("rt.json");
		NonCubicBilinearAlgorithm original = NonCubicBilinearAlgorithm.fromCubic(Strassen7.get());
		SchemeIO.write(original, file.toFile());
		String written = Files.readString(file);
		NonCubicBilinearAlgorithm reloaded = SchemeIO.read(written);
		double[][] origU = original.denseU();
		double[][] origV = original.denseV();
		double[][] origW = original.denseW();
		double[][] reloadedU = reloaded.denseU();
		double[][] reloadedV = reloaded.denseV();
		double[][] reloadedW = reloaded.denseW();
		assertThat(reloadedU).isDeepEqualTo(origU);
		assertThat(reloadedV).isDeepEqualTo(origV);
		assertThat(reloadedW).isDeepEqualTo(origW);
	}

	@Test
	public void nonbilinear_waksman_round_trip_via_materialiser_format() throws IOException {
		// MaterializeWaksman1970 writes the canonical non-bilinear JSON:
		// scheme_type "non_bilinear", sparse u_a/u_b/v_a/v_b/w as [[row, val], ...].
		// We round-trip the smallest case via that format and check
		// SchemeIO.readNonBilinear reconstructs the original factors.
		eu.solven.matmul.NonBilinearAlgorithm original =
				eu.solven.matmul.papers.waksman1970.Waksman1970.build(2);
		// Hand-write the matching JSON (same format as MaterializeWaksman1970.writeSparseJson).
		StringBuilder sb = new StringBuilder();
		sb.append("{\n");
		sb.append("  \"n\": [2, 2, 2],\n");
		sb.append("  \"scheme_type\": \"non_bilinear\",\n");
		sb.append("  \"rank\": ").append(original.r).append(",\n");
		appendSparse(sb, "u_a", original.Ua, original.r);
		appendSparse(sb, "u_b", original.Ub, original.r);
		appendSparse(sb, "v_a", original.Va, original.r);
		appendSparse(sb, "v_b", original.Vb, original.r);
		appendSparseLast(sb, "w", original.W, original.r);
		sb.append("}\n");
		Path file = tmp.resolve("waksman2.json");
		Files.writeString(file, sb.toString());

		tools.jackson.databind.JsonNode root = SchemeIO.parseJson(file.toFile());
		assertThat(SchemeIO.isNonBilinear(root)).isTrue();
		eu.solven.matmul.NonBilinearAlgorithm reloaded = SchemeIO.readNonBilinear(root);
		assertThat(reloaded.n).isEqualTo(original.n);
		assertThat(reloaded.m).isEqualTo(original.m);
		assertThat(reloaded.p).isEqualTo(original.p);
		assertThat(reloaded.r).isEqualTo(original.r);
		assertThat(reloaded.Ua).isDeepEqualTo(original.Ua);
		assertThat(reloaded.Ub).isDeepEqualTo(original.Ub);
		assertThat(reloaded.Va).isDeepEqualTo(original.Va);
		assertThat(reloaded.Vb).isDeepEqualTo(original.Vb);
		assertThat(reloaded.W).isDeepEqualTo(original.W);
		assertThat(Verifier.isExactNonBilinear(reloaded)).isTrue();
	}

	private static void appendSparse(StringBuilder sb, String key, double[][] M, int r) {
		appendSparseImpl(sb, key, M, r, ",");
	}

	private static void appendSparseLast(StringBuilder sb, String key, double[][] M, int r) {
		appendSparseImpl(sb, key, M, r, "");
	}

	private static void appendSparseImpl(StringBuilder sb, String key, double[][] M, int r, String trailing) {
		sb.append("  \"").append(key).append("\": [");
		for (int k = 0; k < r; k++) {
			if (k > 0) sb.append(", ");
			sb.append("[");
			boolean first = true;
			for (int row = 0; row < M.length; row++) {
				double v = M[row][k];
				if (v == 0.0) continue;
				if (!first) sb.append(", ");
				first = false;
				sb.append("[").append(row).append(", ");
				if (v == Math.rint(v)) sb.append((long) v);
				else sb.append(v);
				sb.append("]");
			}
			sb.append("]");
		}
		sb.append("]").append(trailing).append("\n");
	}
}
