package eu.solven.matmul.io;

import static org.assertj.core.api.Assertions.assertThat;
import static org.assertj.core.api.Assertions.assertThatThrownBy;

import java.io.IOException;

import org.junit.jupiter.api.Test;

import eu.solven.matmul.NonCubicBilinearAlgorithm;
import eu.solven.matmul.verifiers.Verifier;

import tools.jackson.databind.JsonNode;
import tools.jackson.databind.json.JsonMapper;
import tools.jackson.databind.node.ArrayNode;
import tools.jackson.databind.node.ObjectNode;

/**
 * The rank-row sparse layout of {@code gsidebottom/mat_mul} and
 * {@code heht24-cyber/fmm-oblique-records} (2026-11): {@code w} has one row per product,
 * unlike Perminov's {@code _reduced} layout, and neither repository states the output
 * index convention. Handing such a file to the Perminov reader would index {@code w}
 * the wrong way round — the guard here is that both spellings and both conventions read
 * back to an exact scheme, and that a file which is a scheme under neither is refused
 * rather than returned.
 */
public class TestRankRowSchemeJson {
	private static final JsonMapper MAPPER = JsonMapper.builder().build();

	/** The naive ⟨2,3,4⟩ algorithm (rank 24), with a rational scaling pushed through two factors. */
	private static NonCubicBilinearAlgorithm naive234() {
		int n = 2, m = 3, p = 4, r = n * m * p;
		double[][] u = new double[n * m][r], v = new double[m * p][r], w = new double[n * p][r];
		int k = 0;
		for (int i = 0; i < n; i++) {
			for (int j = 0; j < m; j++) {
				for (int l = 0; l < p; l++) {
					u[i * m + j][k] = 3;
					v[j * p + l][k] = 0.5;
					w[i * p + l][k] = 2.0 / 3.0;
					k++;
				}
			}
		}
		return new NonCubicBilinearAlgorithm(n, m, p, u, v, w);
	}

	/** Encode as rank rows; {@code transposeW}: positions {@code j·n + i} instead of {@code i·p + j}. */
	private static ObjectNode encode(NonCubicBilinearAlgorithm alg, boolean pairs, boolean transposeW) {
		ObjectNode root = MAPPER.createObjectNode();
		root.putArray("n").add(alg.n).add(alg.m).add(alg.p);
		root.put("m", alg.r);
		String suffix = pairs ? "_sparse" : "";
		rows(root.putArray("u" + suffix), alg.denseU(), alg.r, pairs, -1, -1);
		rows(root.putArray("v" + suffix), alg.denseV(), alg.r, pairs, -1, -1);
		rows(root.putArray("w" + suffix), alg.denseW(), alg.r, pairs, transposeW ? alg.n : -1, alg.p);
		return root;
	}

	private static void rows(ArrayNode out, double[][] f, int rank, boolean pairs, int n, int p) {
		for (int k = 0; k < rank; k++) {
			ArrayNode row = out.addArray();
			for (int pos = 0; pos < f.length; pos++) {
				double c = f[pos][k];
				if (c == 0) continue;
				int index = n < 0 ? pos : (pos % p) * n + pos / p;
				// integers as JSON numbers, everything else as a string fraction - as upstream does
				boolean integral = c == Math.rint(c);
				String frac = c == 0.5 ? "1/2" : "2/3";
				if (pairs) {
					ArrayNode e = row.addArray().add(index);
					if (integral) e.add((int) c);
					else e.add(frac);
				} else {
					ObjectNode e = row.addObject().put("index", index);
					if (integral) e.put("value", (int) c);
					else e.put("value", frac);
				}
			}
		}
	}

	@Test
	public void both_spellings_and_both_output_conventions_read_back_exact() throws IOException {
		NonCubicBilinearAlgorithm ref = naive234();
		assertThat(Verifier.isExactNonCubic(ref)).isTrue();
		for (boolean pairs : new boolean[] { false, true }) {
			for (boolean transposeW : new boolean[] { false, true }) {
				JsonNode root = encode(ref, pairs, transposeW);
				assertThat(RankRowSchemeJson.matches(root)).isTrue();
				NonCubicBilinearAlgorithm back = RankRowSchemeJson.read(root);
				assertThat(back.r).isEqualTo(24);
				assertThat(Verifier.isExactNonCubic(back))
						.as("pairs=%s transposeW=%s must read back as an exact <2,3,4> scheme", pairs, transposeW).isTrue();
				assertThat(back.denseW()).isDeepEqualTo(ref.denseW());
			}
		}
	}

	@Test
	public void a_file_that_is_a_scheme_under_neither_convention_is_refused() {
		ObjectNode root = encode(naive234(), false, false);
		// drop one product's output: no longer the matrix product, whichever way w is indexed
		((ArrayNode) root.get("w")).set(5, MAPPER.createArrayNode());
		assertThatThrownBy(() -> RankRowSchemeJson.read(root)).isInstanceOf(IOException.class)
				.hasMessageContaining("neither output index convention");
	}

	@Test
	public void an_unreadable_coefficient_is_refused_not_read_as_zero() {
		ObjectNode root = encode(naive234(), false, false);
		((ObjectNode) root.get("u").get(0).get(0)).put("value", "one half");
		assertThatThrownBy(() -> RankRowSchemeJson.read(root)).isInstanceOf(IOException.class)
				.hasMessageContaining("unreadable coefficient");
	}

	/** Shukla's spelling: {@code dimensions} / {@code rank} / {@code scale}, pair entries under u, v, w. */
	@Test
	public void shukla_keys_are_read_and_a_scale_other_than_one_is_refused() throws IOException {
		ObjectNode he = encode(naive234(), true, false);
		ObjectNode root = MAPPER.createObjectNode();
		root.set("dimensions", he.get("n"));
		root.put("rank", 24).put("scale", 1).put("layout", "row-major");
		root.set("u", he.get("u_sparse"));
		root.set("v", he.get("v_sparse"));
		root.set("w", he.get("w_sparse"));
		assertThat(RankRowSchemeJson.matches(root)).isTrue();
		assertThat(Verifier.isExactNonCubic(RankRowSchemeJson.read(root))).isTrue();

		root.put("scale", 2);
		assertThatThrownBy(() -> RankRowSchemeJson.read(root)).isInstanceOf(IOException.class)
				.hasMessageContaining("scale");
	}

	@Test
	public void perminov_reduced_layout_is_not_mistaken_for_rank_rows() {
		// Perminov's _reduced w: one row per OUTPUT POSITION (n*p = 8 rows), not per product (24).
		ObjectNode root = encode(naive234(), false, false);
		ArrayNode w = root.putArray("w");
		for (int pos = 0; pos < 8; pos++) w.addArray().addObject().put("index", pos).put("value", 1);
		assertThat(RankRowSchemeJson.matches(root)).isFalse();

		// Perminov's DENSE layout (rank rows of dim numbers) - what Sidebottom's schemes/ folder
		// uses beside the rank-row files of large/ and products/: the importer dispatches on this.
		ArrayNode dense = root.putArray("w");
		for (int k = 0; k < 24; k++) {
			ArrayNode row = dense.addArray();
			for (int pos = 0; pos < 8; pos++) row.add(pos == k % 8 ? 1 : 0);
		}
		assertThat(RankRowSchemeJson.matches(root)).isFalse();
	}
}
