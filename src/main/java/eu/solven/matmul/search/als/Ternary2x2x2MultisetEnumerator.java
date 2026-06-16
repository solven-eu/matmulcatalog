package eu.solven.matmul.search.als;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.LinkedHashMap;
import java.util.List;
import java.util.Map;

import eu.solven.matmul.Verifier;

/**
 * Brute-force enumeration of rank-7 {@code ⟨2,2,2⟩} matmul schemes over the
 * ternary alphabet {@code {-1, 0, +1}}, grouped by the <b>recombination
 * multiset</b> each scheme induces.
 *
 * <h2>What "multiset" means here</h2>
 * Consider a 2-part block decomposition per axis {@code (n₁+n₂, m₁+m₂, p₁+p₂)}.
 * Recombining a {@code ⟨2,2,2⟩} base at that decomposition turns each of the 7
 * products into a smaller matmul {@code ⟨a,b,c⟩} whose sub-dimension on each
 * axis is the {@code max} of the block sizes that product touches (capped by the
 * {@code min} of the two relevant factor views — exactly
 * {@link eu.solven.matmul.search.AnalyticalMaskSearch#shapesAt}). The
 * <b>multiset of those 7 sub-shapes (with counts)</b> is the object of interest:
 * it is the complete <i>rank</i> invariant of plain additive recombination
 * (see {@code references/MULTISET_FRONTIER.md}), and it depends <b>only on the
 * zero/nonzero support pattern</b> of {@code U/V/W}, never on the coefficient
 * values. Strassen and Winograd induce <b>two different</b> such multisets — this
 * enumerator answers "how many distinct multisets does the whole ternary
 * {@code ⟨2,2,2⟩} space realise?".
 *
 * <h2>Why brute force is tractable</h2>
 * A naive sweep of 7 rank-1 terms over {@code ({-1,0,1}⁴)³} is {@code 3^84}.
 * We instead run a <b>covering DFS</b> on the residual tensor: at each step pick
 * the lex-first not-yet-cancelled tensor entry and branch only over rank-1 terms
 * that touch it, pruning by (a) the per-entry completion budget ({@code |R'|} ≤
 * remaining terms) and (b) matricization-rank of the residual once few terms
 * remain. Per-term sign gauge {@code (αβγ=1)} is fixed canonically so each term
 * value is visited once.
 *
 * <p>The canonical split used for reading the multiset is {@code (9,8)} on every
 * axis — i.e. the catalog's reference {@code ⟨17,17,17⟩} at {@code (9,8)³}, where
 * Strassen scores 2940 and Winograd 2930. Sub-dim {@code 8} = "touches only the
 * small block"; sub-dim {@code 9} = "touches the big block or both". The symbolic
 * rule generalises: {@code 8 ↦ N₂} (small block), {@code 9 ↦ max(N₁,N₂)}.
 */
public final class Ternary2x2x2MultisetEnumerator {

	private Ternary2x2x2MultisetEnumerator() {}

	private static final int D = 4; // (n·m) = (m·p) = (n·p) = 4 for ⟨2,2,2⟩
	private static final int R = 7; // target rank

	/** All 81 ternary 4-vectors. */
	private static final int[][] VEC = buildVecs();
	/** Vectors with first-nonzero entry positive (sign-gauge canonical for u, v). */
	private static final int[][] POS = buildPos();
	/** {@code NZ_POS[c]} = canonical-sign vectors with coordinate {@code c} nonzero. */
	private static final int[][][] NZ_POS = buildNonZeroAt(POS);
	/** {@code NZ_ALL[c]} = all 81 vectors with coordinate {@code c} nonzero (used for w). */
	private static final int[][][] NZ_ALL = buildNonZeroAt(VEC);

	private static int[][] buildVecs() {
		int[][] out = new int[81][D];
		for (int idx = 0; idx < 81; idx++) {
			int x = idx;
			for (int i = D - 1; i >= 0; i--) {
				out[idx][i] = (x % 3) - 1;
				x /= 3;
			}
		}
		return out;
	}

	private static int[][] buildPos() {
		List<int[]> pos = new ArrayList<>();
		for (int[] v : VEC) {
			int first = 0;
			for (int i = 0; i < D; i++) {
				if (v[i] != 0) { first = v[i]; break; }
			}
			if (first > 0) pos.add(v);
		}
		return pos.toArray(new int[0][]);
	}

	private static int[][][] buildNonZeroAt(int[][] pool) {
		int[][][] out = new int[D][][];
		for (int c = 0; c < D; c++) {
			List<int[]> list = new ArrayList<>();
			for (int[] v : pool) if (v[c] != 0) list.add(v);
			out[c] = list.toArray(new int[0][]);
		}
		return out;
	}

	/** Tunable run limits + collected results. */
	public static final class Config {
		public long maxNodes = Long.MAX_VALUE;
		public long maxMillis = Long.MAX_VALUE;
		public long progressEveryNodes = 50_000_000L;
		public boolean stopWhenStable = false; // stop early if no new multiset for a long while
	}

	public static final class Result {
		/** multiset key -> representative solution (7 terms, each {u,v,w}). */
		public final Map<String, int[][][]> representatives = new LinkedHashMap<>();
		/** multiset key -> count of distinct ternary solutions realising it. */
		public final Map<String, Long> solutionCounts = new LinkedHashMap<>();
		public long nodes;
		public long solutions;
		public long elapsedMillis;
		public boolean exhaustive; // true iff the DFS finished without hitting a cap
	}

	// ---- DFS state (single-threaded run) ------------------------------------

	private static int[][][] target() {
		return Verifier.intMatmulTensor(2); // 4×4×4, eight ones
	}

	public static Result enumerate(Config cfg) {
		Result res = new Result();
		int[][][] residual = target();
		int[][] terms = new int[R][]; // each = flattened [u(4) | v(4) | w(4)]
		long start = System.nanoTime();
		long[] counters = { 0L /*nodes*/, 0L /*solutions*/, 0L /*lastProgressNodes*/ };
		boolean exhaustive = dfs(residual, terms, 0, cfg, res, counters, start);
		res.nodes = counters[0];
		res.solutions = counters[1];
		res.elapsedMillis = (System.nanoTime() - start) / 1_000_000;
		res.exhaustive = exhaustive;
		return res;
	}

	/** @return true if the subtree was explored exhaustively (no cap hit). */
	private static boolean dfs(int[][][] residual, int[][] terms, int depth,
			Config cfg, Result res, long[] counters, long start) {
		counters[0]++; // node
		if (counters[0] - counters[2] >= cfg.progressEveryNodes) {
			counters[2] = counters[0];
			long ms = (System.nanoTime() - start) / 1_000_000;
			System.out.printf(
					"[progress] nodes=%,d solutions=%,d distinct-multisets=%d depth=%d %,d ms%n",
					counters[0], counters[1], res.representatives.size(), depth, ms);
			System.out.flush();
		}
		if (counters[0] % 4096 == 0) {
			long ms = (System.nanoTime() - start) / 1_000_000;
			if (counters[0] > cfg.maxNodes || ms > cfg.maxMillis) return false; // cap hit → not exhaustive
		}

		// Find lex-first non-zero residual entry.
		int ea = -1, eb = -1, ec = -1;
		outer:
		for (int a = 0; a < D; a++)
			for (int b = 0; b < D; b++)
				for (int c = 0; c < D; c++)
					if (residual[a][b][c] != 0) { ea = a; eb = b; ec = c; break outer; }

		if (ea == -1) {
			// Residual is zero. A genuine ⟨2,2,2⟩ decomposition needs exactly 7
			// essential terms; we only reach here at depth == R.
			if (depth == R) recordSolution(terms, res, counters);
			return true;
		}
		if (depth == R) return true; // no terms left but residual non-zero → dead, fully explored

		int remainingAfter = R - depth - 1;

		boolean exhaustive = true;
		// Covering: the next term must touch entry (ea,eb,ec).
		for (int[] u : NZ_POS[ea]) {
			int ua = u[ea];
			for (int[] v : NZ_POS[eb]) {
				int vb = v[eb];
				for (int[] w : NZ_ALL[ec]) {
					// Subtract the rank-1 term and check the per-entry budget.
					int[][][] next = subtractWithBudget(residual, u, v, w, remainingAfter);
					if (next == null) continue;
					if (remainingAfter <= 3 && exceedsRank(next, remainingAfter)) continue;

					terms[depth] = pack(u, v, w);
					boolean sub = dfs(next, terms, depth + 1, cfg, res, counters, start);
					exhaustive &= sub;
					if (!sub && (counters[0] > cfg.maxNodes
							|| (System.nanoTime() - start) / 1_000_000 > cfg.maxMillis)) {
						return false; // propagate the cap stop
					}
				}
			}
		}
		return exhaustive;
	}

	/**
	 * Return {@code residual - u⊗v⊗w}, or {@code null} if any entry's magnitude
	 * exceeds the number of terms still available to cancel it.
	 */
	private static int[][][] subtractWithBudget(int[][][] residual, int[] u, int[] v, int[] w, int budget) {
		int[][][] out = new int[D][D][D];
		for (int a = 0; a < D; a++) {
			int uaNon = u[a];
			for (int b = 0; b < D; b++) {
				int uv = uaNon * v[b];
				for (int c = 0; c < D; c++) {
					int val = residual[a][b][c] - uv * w[c];
					if (val > budget || val < -budget) return null;
					out[a][b][c] = val;
				}
			}
		}
		return out;
	}

	private static int[] pack(int[] u, int[] v, int[] w) {
		int[] t = new int[3 * D];
		System.arraycopy(u, 0, t, 0, D);
		System.arraycopy(v, 0, t, D, D);
		System.arraycopy(w, 0, t, 2 * D, D);
		return t;
	}

	private static void recordSolution(int[][] terms, Result res, long[] counters) {
		counters[1]++;
		String key = multisetKey(terms);
		res.solutionCounts.merge(key, 1L, Long::sum);
		res.representatives.computeIfAbsent(key, k -> deepCopy(terms));
	}

	private static int[][][] deepCopy(int[][] terms) {
		int[][][] out = new int[R][3][D];
		for (int k = 0; k < R; k++) {
			System.arraycopy(terms[k], 0, out[k][0], 0, D);
			System.arraycopy(terms[k], D, out[k][1], 0, D);
			System.arraycopy(terms[k], 2 * D, out[k][2], 0, D);
		}
		return out;
	}

	// ---- multiset computation (mirrors AnalyticalMaskSearch.shapesAt) --------

	// Canonical reference split: big block = 9, small block = 8 on every axis
	// → ⟨17,17,17⟩ at (9,8)³. block0 is the big one.
	private static final int BIG = 9, SMALL = 8;

	/** max over touched blocks: bit0 (big) dominates; else small. */
	private static int maxBlk(int bits) {
		return (bits & 1) != 0 ? BIG : SMALL;
	}

	/**
	 * The multiset of 7 sub-shapes at the (9,8)³ reference split, sorted and
	 * serialised. Index decode for ⟨2,2,2⟩: A-entry {@code a = i·2 + l}
	 * (i=row/n-block, l=col/m-block); B-entry {@code b = l·2 + j}; C-entry
	 * {@code c = i·2 + j}.
	 */
	static String multisetKey(int[][] terms) {
		int[][] shapes = new int[R][3];
		for (int k = 0; k < R; k++) {
			int[] t = terms[k];
			int uRowN = 0, uColM = 0, vRowM = 0, vColP = 0, wRowN = 0, wColP = 0;
			for (int a = 0; a < D; a++) if (t[a] != 0) { uRowN |= 1 << (a >> 1); uColM |= 1 << (a & 1); }
			for (int b = 0; b < D; b++) if (t[D + b] != 0) { vRowM |= 1 << (b >> 1); vColP |= 1 << (b & 1); }
			for (int c = 0; c < D; c++) if (t[2 * D + c] != 0) { wRowN |= 1 << (c >> 1); wColP |= 1 << (c & 1); }
			int subA = Math.min(maxBlk(uRowN), maxBlk(wRowN));
			int subB = Math.min(maxBlk(uColM), maxBlk(vRowM));
			int subC = Math.min(maxBlk(vColP), maxBlk(wColP));
			shapes[k] = new int[] { subA, subB, subC };
		}
		return renderMultiset(shapes);
	}

	/** Stable, human-readable multiset string: counts of each distinct ⟨a,b,c⟩. */
	static String renderMultiset(int[][] shapes) {
		String[] s = new String[shapes.length];
		for (int i = 0; i < shapes.length; i++) {
			s[i] = shapes[i][0] + "x" + shapes[i][1] + "x" + shapes[i][2];
		}
		Arrays.sort(s);
		// collapse to counted form: "3·8x9x9 + 1·9x9x9 ..."
		Map<String, Integer> counts = new LinkedHashMap<>();
		for (String x : s) counts.merge(x, 1, Integer::sum);
		StringBuilder sb = new StringBuilder();
		boolean first = true;
		for (var e : counts.entrySet()) {
			if (!first) sb.append(" + ");
			sb.append(e.getValue()).append("·⟨").append(e.getKey().replace('x', ',')).append('⟩');
			first = false;
		}
		return sb.toString();
	}

	// ---- residual matricization rank prune ----------------------------------

	/** True iff any of the three matricizations of {@code R} has rank > {@code maxRank}. */
	private static boolean exceedsRank(int[][][] r, int maxRank) {
		// flatten 1: rows a (4), cols (b·4+c) (16)
		double[][] m1 = new double[D][D * D];
		double[][] m2 = new double[D][D * D];
		double[][] m3 = new double[D][D * D];
		for (int a = 0; a < D; a++)
			for (int b = 0; b < D; b++)
				for (int c = 0; c < D; c++) {
					double v = r[a][b][c];
					m1[a][b * D + c] = v;
					m2[b][a * D + c] = v;
					m3[c][a * D + b] = v;
				}
		return matRank(m1) > maxRank || matRank(m2) > maxRank || matRank(m3) > maxRank;
	}

	/** Gaussian-elimination rank of a small real matrix. */
	private static int matRank(double[][] mat) {
		int rows = mat.length, cols = mat[0].length;
		int rank = 0;
		for (int col = 0; col < cols && rank < rows; col++) {
			int pivot = -1;
			double best = 1e-9;
			for (int row = rank; row < rows; row++) {
				double v = Math.abs(mat[row][col]);
				if (v > best) { best = v; pivot = row; }
			}
			if (pivot < 0) continue;
			double[] tmp = mat[rank]; mat[rank] = mat[pivot]; mat[pivot] = tmp;
			double pv = mat[rank][col];
			for (int row = 0; row < rows; row++) {
				if (row == rank) continue;
				double f = mat[row][col] / pv;
				if (f == 0.0) continue;
				for (int cc = col; cc < cols; cc++) mat[row][cc] -= f * mat[rank][cc];
			}
			rank++;
		}
		return rank;
	}
}
