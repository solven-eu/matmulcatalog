package eu.solven.matmul.io;

import java.io.IOException;

import eu.solven.matmul.NonCubicBilinearAlgorithm;
import eu.solven.matmul.verifiers.Verifier;

import tools.jackson.databind.JsonNode;

/**
 * Reader for the <b>rank-row sparse</b> scheme JSON two 2026 repositories publish
 * (gzipped): all three factors are lists of {@code rank} rows, each row the non-zero
 * entries of one product.
 *
 * <ul>
 *   <li><b>Sidebottom</b> ({@code gsidebottom/mat_mul}): {@code u} / {@code v} / {@code w},
 *       entries {@code {"index": i, "value": c}}.</li>
 *   <li><b>He</b> ({@code heht24-cyber/fmm-oblique-records}): {@code u_sparse} /
 *       {@code v_sparse} / {@code w_sparse}, entries {@code [i, c]}.</li>
 *   <li><b>Shukla</b> ({@code chandanshukla1989/fmm_decompositions}): {@code u} / {@code v} /
 *       {@code w}, entries {@code [i, c]}; shape under {@code dimensions}, rank under
 *       {@code rank}, and a {@code scale} divisor that is 1 in every published file.</li>
 * </ul>
 *
 * <p>A coefficient is a JSON integer or a string fraction ({@code "-5/6"}). Shape is
 * {@code n}, rank is {@code m} — the same keys as Perminov's layout, but this is
 * <em>not</em> Perminov's {@code _reduced} layout (there {@code w} has one row per
 * output position, here one per product), so {@code SchemeIO.readBilinear} must not
 * be handed these files: it would index {@code w} the wrong way round.</p>
 *
 * <p><b>The output index convention is decided by the identity, not assumed.</b> Neither
 * repository states whether a {@code w} position is {@code i·p + j} (row-major over the
 * {@code n × p} product) or {@code j·n + i} (Perminov's transposed flatten). Both are
 * built and the one that computes the matrix product on random inputs is returned; a
 * file where neither does is refused. The caller's gate still runs the exact proof.</p>
 */
public final class RankRowSchemeJson {
	private RankRowSchemeJson() {}

	/** Whether {@code root} carries the rank-row layout under either spelling. */
	public static boolean matches(JsonNode root) {
		JsonNode w = factor(root, "w");
		if (w == null || !w.isArray() || rank(root) < 0 || w.size() != rank(root)) return false;
		for (JsonNode row : w) {
			if (!row.isArray()) return false;
			if (row.isEmpty()) continue;
			JsonNode e = row.get(0);
			return e.isObject() ? e.has("index") && e.has("value") : e.isArray() && e.size() == 2;
		}
		return false;
	}

	public static NonCubicBilinearAlgorithm read(JsonNode root) throws IOException {
		JsonNode nArr = shape(root);
		int rank = rank(root);
		if (nArr == null || !nArr.isArray() || nArr.size() != 3 || rank < 0) {
			throw new IOException("rank-row scheme: expected n (or dimensions) = [n, m, p] and an integer m (or rank)");
		}
		JsonNode scale = root.get("scale");
		if (scale != null && !(scale.isInt() && scale.asInt() == 1)) {
			// Shukla's layout divides the output by `scale`; every published file has 1.
			throw new IOException("rank-row scheme: scale " + scale + " is not supported (only 1)");
		}
		int n = nArr.get(0).asInt(), m = nArr.get(1).asInt(), p = nArr.get(2).asInt();
		double[][] u = rows(root, "u", rank, n * m);
		double[][] v = rows(root, "v", rank, m * p);
		double[][] wRaw = rows(root, "w", rank, n * p);
		NonCubicBilinearAlgorithm rowMajor = new NonCubicBilinearAlgorithm(n, m, p, u, v, wRaw);
		if (Verifier.passesRandomMatmulSpotCheck(rowMajor)) return rowMajor;
		// Perminov's flatten: position j·n + i holds the coefficient of C[i][j].
		double[][] w = new double[n * p][];
		for (int i = 0; i < n; i++) {
			for (int j = 0; j < p; j++) {
				w[i * p + j] = wRaw[j * n + i];
			}
		}
		NonCubicBilinearAlgorithm transposed = new NonCubicBilinearAlgorithm(n, m, p, u, v, w);
		if (Verifier.passesRandomMatmulSpotCheck(transposed)) return transposed;
		throw new IOException("rank-row scheme <" + n + "," + m + "," + p + ":" + rank + ">: computes the matrix product"
				+ " under neither output index convention (i*p+j, j*n+i) - not a scheme we can read");
	}

	/** {@code n} (Perminov's key, used by Sidebottom and He) or {@code dimensions} (Shukla). */
	private static JsonNode shape(JsonNode root) {
		JsonNode n = root.get("n");
		return n != null ? n : root.get("dimensions");
	}

	/** {@code m} or, without it, {@code rank}; -1 when neither is an integer. */
	private static int rank(JsonNode root) {
		JsonNode m = root.get("m");
		if (m == null) m = root.get("rank");
		return m != null && m.isInt() ? m.asInt() : -1;
	}

	private static JsonNode factor(JsonNode root, String key) {
		JsonNode f = root.get(key);
		return f != null ? f : root.get(key + "_sparse");
	}

	/** {@code dim × rank} dense factor from {@code rank} sparse rows. */
	private static double[][] rows(JsonNode root, String key, int rank, int dim) throws IOException {
		JsonNode f = factor(root, key);
		if (f == null || !f.isArray() || f.size() != rank) {
			throw new IOException("rank-row scheme: factor '" + key + "' must be a list of " + rank + " rows, got "
					+ (f == null ? "nothing" : f.isArray() ? f.size() + " rows" : f.getNodeType().toString()));
		}
		double[][] out = new double[dim][rank];
		for (int k = 0; k < rank; k++) {
			for (JsonNode e : f.get(k)) {
				JsonNode idx, val;
				if (e.isObject()) {
					idx = e.get("index");
					val = e.get("value");
				} else if (e.isArray() && e.size() == 2) {
					idx = e.get(0);
					val = e.get(1);
				} else {
					throw new IOException("rank-row scheme: " + key + " row " + k + ": expected {index, value} or [index,"
							+ " value], got " + e);
				}
				if (idx == null || !idx.isIntegralNumber() || val == null) {
					throw new IOException("rank-row scheme: " + key + " row " + k + ": malformed entry " + e);
				}
				int i = idx.asInt();
				if (i < 0 || i >= dim) {
					throw new IOException("rank-row scheme: " + key + " row " + k + ": index " + i + " outside [0, " + dim + ")");
				}
				out[i][k] += coefficient(val);
			}
		}
		return out;
	}

	/** A JSON number, or a string integer / fraction ({@code "-5/6"}). Anything else is refused. */
	static double coefficient(JsonNode v) throws IOException {
		if (v.isNumber()) return v.asDouble();
		if (v.isString()) {
			String t = v.asString().trim();
			try {
				int slash = t.indexOf('/');
				return slash < 0 ? Double.parseDouble(t)
						: Double.parseDouble(t.substring(0, slash)) / Double.parseDouble(t.substring(slash + 1));
			} catch (NumberFormatException e) {
				throw new IOException("rank-row scheme: unreadable coefficient " + v);
			}
		}
		throw new IOException("rank-row scheme: unreadable coefficient " + v);
	}
}
