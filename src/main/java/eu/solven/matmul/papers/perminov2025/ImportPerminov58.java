package eu.solven.matmul.papers.perminov2025;

import java.io.File;

import eu.solven.matmul.NonCubicBilinearAlgorithm;
import eu.solven.matmul.Verifier;
import eu.solven.matmul.catalog.SchemeIO;
import lombok.extern.slf4j.Slf4j;
import tools.jackson.databind.JsonNode;
import tools.jackson.databind.json.JsonMapper;
import tools.jackson.databind.node.ArrayNode;
import tools.jackson.databind.node.ObjectNode;

/**
 * Import Perminov 2025 — "A 58-Addition, Rank-23 Scheme for General 3×3 Matrix
 * Multiplication" (arXiv:2512.21980). A new rank-23 ⟨3,3,3⟩ scheme with only
 * <strong>58 scalar additions</strong> (34 add + 24 sub), ternary {−1,0,1}, no
 * change of basis — the additive-complexity record (prev. 60, Stapleton 2025).
 *
 * <p>The 58 is the <em>scheduled</em> count: it reuses 4 A-combos (u₁..u₄), 8
 * B-combos (v₁..v₈) and 8 product-combos (w₁..w₈) as common subexpressions. We
 * reconstruct the bilinear factor matrices U/V/W by expanding the paper's SLP,
 * verify they compute matmul, then stamp {@code scheduled_additions:58} + the
 * SLP so the SPA surfaces the +58.</p>
 */
@Slf4j
public final class ImportPerminov58 {
	private ImportPerminov58() {}

	private static final int SZ = 9;     // 3·3 entries per factor
	private static final int R = 23;

	// a_{ij} → i*3+j (0-based); same for b_{jl}, c_{il}.
	private static double[] e(int i) {
		double[] v = new double[SZ];
		v[i] = 1;
		return v;
	}

	private static double[] add(double[] x, double[] y) {
		double[] o = x.clone();
		for (int i = 0; i < o.length; i++) o[i] += y[i];
		return o;
	}

	private static double[] sub(double[] x, double[] y) {
		double[] o = x.clone();
		for (int i = 0; i < o.length; i++) o[i] -= y[i];
		return o;
	}

	/** The paper's SLP, verbatim — the source of truth for both matrices and the slp[] field. */
	private static final String[] SLP = {
			"u1 = a31 + a33", "u2 = a21 + a22", "u3 = a13 + u1", "u4 = a32 - u2",
			"v1 = b22 + b32", "v2 = b31 - v1", "v3 = b12 + v2", "v4 = b11 - v3",
			"v5 = b33 + v1", "v6 = b21 - b23", "v7 = b12 - v5", "v8 = v4 - v6",
			"m1 = u1 * v5", "m2 = u2 * (v2 + v6)", "m3 = a32 * b23", "m4 = a31 * (b13 + v7)",
			"m5 = u3 * v7", "m6 = (a32 - a33) * b22", "m7 = a23 * b33", "m8 = (u1 - u2) * v2",
			"m9 = (a12 - a13) * b22", "m10 = (u3 - a21) * v3", "m11 = (a13 + a33) * (v1 - b12)",
			"m12 = a13 * b33", "m13 = (a31 + u4) * v4", "m14 = a11 * b13", "m15 = (a11 + u3) * b12",
			"m16 = (a13 + a23) * b31", "m17 = (a22 - a32) * (v4 - b21)", "m18 = (a11 + a21) * b11",
			"m19 = a12 * b23", "m20 = (a12 + a22) * b21", "m21 = a21 * (b13 - v8)",
			"m22 = (a22 - a23) * (b31 - b32)", "m23 = u4 * v8",
			"w1 = m5 + m12", "w2 = m1 + w1", "w3 = m8 + w2", "w4 = m3 - m23",
			"w5 = m2 + w4", "w6 = w3 + w5", "w7 = m10 - w6", "w8 = m17 + w7",
			"c11 = m18 + m20 + w8", "c12 = m9 + m15 - w2", "c13 = m12 + m14 + m19",
			"c21 = m16 - w8", "c22 = m16 + m22 + w3 - m10", "c23 = m7 + m21 + w4 - m17",
			"c31 = m11 + m13 + w6", "c32 = m6 + m11 + w2", "c33 = m3 + m4 - m11 - w1",
	};

	public static void main(String[] args) throws Exception {
		// A index: a11=0 a12=1 a13=2 a21=3 a22=4 a23=5 a31=6 a32=7 a33=8
		double[] u1 = add(e(6), e(8)), u2 = add(e(3), e(4)), u3 = add(e(2), u1), u4 = sub(e(7), u2);
		// B index: b11=0 b12=1 b13=2 b21=3 b22=4 b23=5 b31=6 b32=7 b33=8
		double[] v1 = add(e(4), e(7)), v2 = sub(e(6), v1), v3 = add(e(1), v2), v4 = sub(e(0), v3);
		double[] v5 = add(e(8), v1), v6 = sub(e(3), e(5)), v7 = sub(e(1), v5), v8 = sub(v4, v6);

		double[][] Acol = new double[R][], Bcol = new double[R][];
		Acol[0] = u1;                 Bcol[0] = v5;
		Acol[1] = u2;                 Bcol[1] = add(v2, v6);
		Acol[2] = e(7);               Bcol[2] = e(5);
		Acol[3] = e(6);               Bcol[3] = add(e(2), v7);
		Acol[4] = u3;                 Bcol[4] = v7;
		Acol[5] = sub(e(7), e(8));    Bcol[5] = e(4);
		Acol[6] = e(5);               Bcol[6] = e(8);
		Acol[7] = sub(u1, u2);        Bcol[7] = v2;
		Acol[8] = sub(e(1), e(2));    Bcol[8] = e(4);
		Acol[9] = sub(u3, e(3));      Bcol[9] = v3;
		Acol[10] = add(e(2), e(8));   Bcol[10] = sub(v1, e(1));
		Acol[11] = e(2);              Bcol[11] = e(8);
		Acol[12] = add(e(6), u4);     Bcol[12] = v4;
		Acol[13] = e(0);              Bcol[13] = e(2);
		Acol[14] = add(e(0), u3);     Bcol[14] = e(1);
		Acol[15] = add(e(2), e(5));   Bcol[15] = e(6);
		Acol[16] = sub(e(4), e(7));   Bcol[16] = sub(v4, e(3));
		Acol[17] = add(e(0), e(3));   Bcol[17] = e(0);
		Acol[18] = e(1);              Bcol[18] = e(5);
		Acol[19] = add(e(1), e(4));   Bcol[19] = e(3);
		Acol[20] = e(3);              Bcol[20] = sub(e(2), v8);
		Acol[21] = sub(e(4), e(5));   Bcol[21] = sub(e(6), e(7));
		Acol[22] = u4;                Bcol[22] = v8;

		// product-combos w₁..w₈ and C reconstruction, as vectors over m₁..m₂₃.
		double[] mw1 = madd(mc(5), mc(12));
		double[] mw2 = madd(mc(1), mw1);
		double[] mw3 = madd(mc(8), mw2);
		double[] mw4 = msub(mc(3), mc(23));
		double[] mw5 = madd(mc(2), mw4);
		double[] mw6 = madd(mw3, mw5);
		double[] mw7 = msub(mc(10), mw6);
		double[] mw8 = madd(mc(17), mw7);

		double[][] C = new double[SZ][];
		C[0] = madd(madd(mc(18), mc(20)), mw8);                 // c11
		C[1] = msub(madd(mc(9), mc(15)), mw2);                  // c12
		C[2] = madd(madd(mc(12), mc(14)), mc(19));              // c13
		C[3] = msub(mc(16), mw8);                               // c21
		C[4] = msub(madd(madd(mc(16), mc(22)), mw3), mc(10));   // c22
		C[5] = msub(madd(madd(mc(7), mc(21)), mw4), mc(17));    // c23
		C[6] = madd(madd(mc(11), mc(13)), mw6);                 // c31
		C[7] = madd(madd(mc(6), mc(11)), mw2);                  // c32
		C[8] = msub(msub(madd(mc(3), mc(4)), mc(11)), mw1);     // c33

		double[][] U = new double[SZ][R], V = new double[SZ][R], W = new double[SZ][R];
		for (int k = 0; k < R; k++) {
			for (int i = 0; i < SZ; i++) {
				U[i][k] = Acol[k][i];
				V[i][k] = Bcol[k][i];
			}
		}
		for (int row = 0; row < SZ; row++) {
			System.arraycopy(C[row], 0, W[row], 0, R);
		}

		NonCubicBilinearAlgorithm alg = new NonCubicBilinearAlgorithm(3, 3, 3, U, V, W);
		boolean exact = Verifier.isExactNonCubic(alg);
		int flat = Verifier.additionCount(alg);
		log.info("Perminov58 ⟨3,3,3⟩ r={} exact={} flatAdds={} scheduled=58", alg.r, exact, flat);
		if (!exact) {
			throw new IllegalStateException("reconstructed scheme does not compute matmul — check SLP transcription");
		}

		File dir = new File("src/main/resources/schemes/known/section3");
		File f = new File(dir, SchemeIO.canonicalName(alg, "perminov_2025"));
		SchemeIO.write(alg, f);

		// Stamp the SLP metadata: scheduled (CSE) additions, year, reference, slp[].
		JsonMapper mapper = JsonMapper.builder().build();
		JsonNode parsed = SchemeIO.parseJson(f);
		ObjectNode root = (ObjectNode) parsed;
		root.put("scheduled_additions", 58);
		root.put("year", 2025);
		root.put("reference", "Perminov 2025, arXiv:2512.21980 — A 58-Addition, Rank-23 Scheme for 3x3 MatMul");
		ArrayNode slp = mapper.createArrayNode();
		for (String line : SLP) slp.add(line);
		root.set("slp", slp);
		try (java.io.PrintWriter pw = new java.io.PrintWriter(new java.io.BufferedWriter(new java.io.FileWriter(f)))) {
			pw.print(eu.solven.matmul.catalog.MatrixJsonFormatter.format(mapper.writeValueAsString(root)));
		}
		log.info("wrote {} (scheduled_additions=58, slp {} lines)", f.getName(), SLP.length);
	}

	// ── m-vectors (over m₁..m₂₃, 1-based label) ──
	private static double[] mc(int kOneBased) {
		double[] v = new double[R];
		v[kOneBased - 1] = 1;
		return v;
	}

	private static double[] madd(double[] x, double[] y) {
		double[] o = x.clone();
		for (int i = 0; i < o.length; i++) o[i] += y[i];
		return o;
	}

	private static double[] msub(double[] x, double[] y) {
		double[] o = x.clone();
		for (int i = 0; i < o.length; i++) o[i] -= y[i];
		return o;
	}
}
