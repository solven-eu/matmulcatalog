package eu.solven.matmul.papers.pan1978;

import static org.assertj.core.api.Assertions.assertThat;

import org.junit.jupiter.api.Test;

/**
 * Rectangular Pan-TA (see references/RECTANGULAR_TA.md): cyclic-rotation pair
 * detection, the exact aggregation identity (3) and disjoint assembly, the
 * derived fused-rank formula {@code nrp+np+nr+rp} (= FMM's 2860 for ⟨26,3,26⟩),
 * and the exact construction {@code build(n,r,p)} verified bit-exact against the
 * disjoint-sum tensor. {@code build} is deterministic in {@code (n,r,p)} → the
 * lineage replays by re-running it.
 */
public class TestRectangularTrilinearAggregation {

	@Test
	public void detects_the_26x29x29_cross_rotation_pair() {
		// ⟨1,2,2⟩ symmetric peel of ⟨26,29,29⟩ → 4 leaves; the two cross-blocks
		// ⟨26,3,26⟩ and ⟨26,26,3⟩ are cyclic rotations of {26,3,26}.
		int[][] leaves = { { 26, 26, 26 }, { 26, 3, 26 }, { 26, 26, 3 }, { 26, 3, 3 } };
		var pairs = RectangularTrilinearAggregation.cyclicRotationPairs(leaves);
		assertThat(pairs).hasSize(1);
		assertThat(pairs.get(0)).containsExactly(1, 2);  // ⟨26,3,26⟩ × ⟨26,26,3⟩
	}

	@Test
	public void cyclic_rotation_semantics() {
		// (26,3,26) → rot1 (3,26,26) → rot2 (26,26,3)
		assertThat(RectangularTrilinearAggregation.isCyclicRotation(
				new int[] { 26, 3, 26 }, new int[] { 26, 26, 3 })).isTrue();
		assertThat(RectangularTrilinearAggregation.isCyclicRotation(
				new int[] { 26, 3, 26 }, new int[] { 3, 26, 26 })).isTrue();
		// identity is not a (distinct) rotation pair
		assertThat(RectangularTrilinearAggregation.isCyclicRotation(
				new int[] { 26, 3, 26 }, new int[] { 26, 3, 26 })).isFalse();
		// ⟨n,n,n⟩ rotations coincide → no distinct pair
		assertThat(RectangularTrilinearAggregation.isCyclicRotation(
				new int[] { 7, 7, 7 }, new int[] { 7, 7, 7 })).isFalse();
		// unrelated shapes
		assertThat(RectangularTrilinearAggregation.isCyclicRotation(
				new int[] { 2, 3, 4 }, new int[] { 2, 4, 3 })).isFalse();
	}

	@Test
	public void aggregation_identity_3_is_exact() {
		// Pan (3): a⊗b⊗c + a'⊗b'⊗c' == the 4-term aggregation. Verify the tensors
		// are EXACTLY equal over arbitrary small integer factor vectors.
		int[] a = { 1, -2, 3 }, b = { 2, 1 }, c = { -1, 0, 2, 1 };
		int[] a2 = { 0, 4, -1 }, b2 = { -3, 5 }, c2 = { 2, -2, 1, 3 };
		var agg = RectangularTrilinearAggregation.aggregatePair(a, b, c, a2, b2, c2);

		int du = a.length, dv = b.length, dw = c.length;
		long[][][] lhs = new long[du][dv][dw];  // a⊗b⊗c + a'⊗b'⊗c'
		long[][][] rhs = new long[du][dv][dw];  // the 4 aggregation terms
		addOuter(lhs, a, b, c);
		addOuter(lhs, a2, b2, c2);
		for (var t : agg) addOuter(rhs, t.u(), t.v(), t.w());
		for (int i = 0; i < du; i++)
			for (int j = 0; j < dv; j++)
				for (int k = 0; k < dw; k++)
					assertThat(rhs[i][j][k]).as("entry (%d,%d,%d)", i, j, k).isEqualTo(lhs[i][j][k]);
	}

	@Test
	public void disjoint_aggregation_computes_the_disjoint_sum_exactly() {
		// Two small rank-one decompositions over disjoint spaces (the cross-block
		// analogue). Aggregating them must reproduce p1 ⊕ p2 (block-diagonal) exactly.
		int du1 = 2, dv1 = 2, dw1 = 2, du2 = 3, dv2 = 2, dw2 = 2;
		var p1 = java.util.List.of(
				new RectangularTrilinearAggregation.Term(new int[] { 1, 0 }, new int[] { 2, -1 }, new int[] { 1, 1 }),
				new RectangularTrilinearAggregation.Term(new int[] { -1, 2 }, new int[] { 0, 3 }, new int[] { 2, -1 }));
		var p2 = java.util.List.of(
				new RectangularTrilinearAggregation.Term(new int[] { 1, 1, 0 }, new int[] { 1, 2 }, new int[] { -1, 1 }),
				new RectangularTrilinearAggregation.Term(new int[] { 0, -2, 1 }, new int[] { 3, 0 }, new int[] { 1, 2 }));

		var agg = RectangularTrilinearAggregation.aggregateDisjoint(p1, du1, dv1, dw1, p2, du2, dv2, dw2);

		int du = du1 + du2, dv = dv1 + dv2, dw = dw1 + dw2;
		long[][][] expected = new long[du][dv][dw];  // p1 in block 0, p2 offset
		for (var t : p1) addOuter(expected, embed(t.u(), du, 0), embed(t.v(), dv, 0), embed(t.w(), dw, 0));
		for (var t : p2) addOuter(expected, embed(t.u(), du, du1), embed(t.v(), dv, dv1), embed(t.w(), dw, dw1));
		long[][][] got = new long[du][dv][dw];
		for (var t : agg) addOuter(got, t.u(), t.v(), t.w());

		for (int i = 0; i < du; i++)
			for (int j = 0; j < dv; j++)
				for (int k = 0; k < dw; k++)
					assertThat(got[i][j][k]).as("(%d,%d,%d)", i, j, k).isEqualTo(expected[i][j][k]);
	}

	private static int[] embed(int[] v, int dim, int offset) {
		int[] r = new int[dim];
		System.arraycopy(v, 0, r, offset, v.length);
		return r;
	}

	@Test
	public void buildPeeledViaTa_assembles_a_correct_matmul() {
		// ⟨2,3,3⟩ = ⟨N=2, N+s=3, N+s=3⟩ with s=1: diag ⟨2,2,2⟩ + TA(2,1,2) + corner ⟨2,1,1⟩.
		// Rank = 8 + fusedRank(2,1,2)=12 + 2 = 22 (not optimal — verifying CORRECTNESS).
		var scheme = RectangularTrilinearAggregation.buildPeeledViaTa(
				2, 1, naive(2, 2, 2), naive(2, 1, 1));
		assertThat(scheme.n).isEqualTo(2);
		assertThat(scheme.m).isEqualTo(3);
		assertThat(scheme.p).isEqualTo(3);
		assertThat(scheme.r).isEqualTo(8 + (int) RectangularTrilinearAggregation.fusedRank(2, 1, 2) + 2);
		assertThat(eu.solven.matmul.verifiers.Verifier.isExactNonCubic(scheme))
				.as("buildPeeledViaTa must compute ⟨2,3,3⟩ exactly").isTrue();
	}

	@org.junit.jupiter.api.Test
	public void peeledViaTa_lineage_round_trips_for_replay(@org.junit.jupiter.api.io.TempDir java.io.File dir)
			throws Exception {
		// The lineage must carry ALL info needed for bit-exact replay: (n, s, cube, corner).
		var alg = RectangularTrilinearAggregation.buildPeeledViaTa(2, 1, naive(2, 2, 2), naive(2, 1, 1));
		eu.solven.matmul.catalog.Lineage.Node lin = new eu.solven.matmul.catalog.Lineage.PeeledViaTa(
				2, 1, new eu.solven.matmul.catalog.Lineage.Atom("2x2x2-cube"),
				new eu.solven.matmul.catalog.Lineage.Atom("2x1x1-corner"));
		java.io.File f = new java.io.File(dir, "ta.json");
		eu.solven.matmul.catalog.SchemeIO.write(alg, f, lin);
		var read = eu.solven.matmul.catalog.SchemeIO.readLineage(f);
		assertThat(read).isPresent();
		assertThat(read.get()).isInstanceOf(eu.solven.matmul.catalog.Lineage.PeeledViaTa.class);
		var pt = (eu.solven.matmul.catalog.Lineage.PeeledViaTa) read.get();
		assertThat(pt.n()).isEqualTo(2);
		assertThat(pt.s()).isEqualTo(1);
		assertThat(pt.cube()).isEqualTo(new eu.solven.matmul.catalog.Lineage.Atom("2x2x2-cube"));
		assertThat(pt.corner()).isEqualTo(new eu.solven.matmul.catalog.Lineage.Atom("2x1x1-corner"));
	}

	/** Trivial (naive) ⟨n,m,p⟩ matmul scheme, rank n·m·p. */
	private static eu.solven.matmul.NonCubicBilinearAlgorithm naive(int n, int m, int p) {
		int r = n * m * p;
		double[][] U = new double[n * m][r], V = new double[m * p][r], W = new double[n * p][r];
		int t = 0;
		for (int i = 0; i < n; i++) for (int j = 0; j < m; j++) for (int k = 0; k < p; k++) {
			U[i * m + j][t] = 1; V[j * p + k][t] = 1; W[i * p + k][t] = 1; t++;
		}
		return new eu.solven.matmul.NonCubicBilinearAlgorithm(n, m, p, U, V, W);
	}

	private static void addOuter(long[][][] t, int[] u, int[] v, int[] w) {
		for (int i = 0; i < u.length; i++)
			for (int j = 0; j < v.length; j++)
				for (int k = 0; k < w.length; k++)
					t[i][j][k] += (long) u[i] * v[j] * w[k];
	}

	@Test
	public void fused_rank_matches_FMM_2860() {
		// ⟨26,3,26⟩ ⊕ ⟨26,26,3⟩ via TA = 2028 + 676 + 78 + 78 = 2860 (FMM's TA term).
		assertThat(RectangularTrilinearAggregation.fusedRank(26, 3, 26)).isEqualTo(2860L);
		// general formula nrp + np + nr + rp
		assertThat(RectangularTrilinearAggregation.fusedRank(3, 2, 2)).isEqualTo(12 + 6 + 6 + 4L);
	}

	@Test
	public void build_computes_the_disjoint_sum_exactly_and_at_the_stated_rank() {
		int n = 3, r = 2, p = 2;  // distinct dims to catch index bugs
		var ta = RectangularTrilinearAggregation.build(n, r, p);
		assertThat(ta.terms()).hasSize((int) RectangularTrilinearAggregation.fusedRank(n, r, p));

		// expected disjoint-sum tensor: trivial LA of ⟨n,r,p⟩ (block 1) + ⟨p,n,r⟩ (block 2)
		long[][][] expected = new long[ta.dimU()][ta.dimV()][ta.dimW()];
		for (int i = 0; i < n; i++) for (int j = 0; j < r; j++) for (int k = 0; k < p; k++)
			expected[i * r + j][j * p + k][i * p + k] += 1;                    // A1⊗B1⊗C1
		for (int ip = 0; ip < p; ip++) for (int jp = 0; jp < n; jp++) for (int kp = 0; kp < r; kp++)
			expected[n * r + ip * n + jp][r * p + jp * r + kp][n * p + ip * r + kp] += 1; // A2⊗B2⊗C2

		long[][][] got = new long[ta.dimU()][ta.dimV()][ta.dimW()];
		for (var t : ta.terms()) addOuter(got, t.u(), t.v(), t.w());

		for (int u = 0; u < ta.dimU(); u++)
			for (int v = 0; v < ta.dimV(); v++)
				for (int w = 0; w < ta.dimW(); w++)
					assertThat(got[u][v][w]).as("(%d,%d,%d)", u, v, w).isEqualTo(expected[u][v][w]);
	}
}
