package eu.solven.matmul;

import java.util.ArrayList;
import java.util.List;

/**
 * Analyses how a rank-r decomposition of a matmul tensor restricts to its
 * sub-block embeddings. For ⟨3,3,3⟩ there are 3³ = 27 distinct ⟨2,2,2⟩
 * sub-block embeddings, each characterized by a choice of 2 rows of A
 * (= rows of C), 2 cols of A (= rows of B), and 2 cols of B (= cols of C).
 *
 * <p>For each (rank-1 term, sub-block) pair, the term has a non-vanishing
 * restriction iff its u, v, w factors all have at least one non-zero entry
 * inside the sub-block's index set. By the Hopcroft–Kerr lower bound
 * `R(⟨2,2,2⟩) = 7`, every sub-block must have ≥ 7 non-vanishing terms in
 * any valid r-decomposition.</p>
 */
public final class SubblockAnalyzer {

	private SubblockAnalyzer() {}

	public static final class Subblock {
		public final int n; // outer matrix dimension (e.g., 3 for ⟨3,3,3⟩)
		public final int[] rowsA; // 2-subset of {0..n-1}
		public final int[] colsA; // 2-subset of {0..n-1} — same as rowsB
		public final int[] colsB; // 2-subset of {0..n-1}

		public Subblock(int n, int[] rowsA, int[] colsA, int[] colsB) {
			this.n = n;
			this.rowsA = rowsA;
			this.colsA = colsA;
			this.colsB = colsB;
		}

		public int[] aPositions() {
			int[] out = new int[rowsA.length * colsA.length];
			int idx = 0;
			for (int i : rowsA) for (int l : colsA) out[idx++] = i * n + l;
			return out;
		}

		public int[] bPositions() {
			int[] out = new int[colsA.length * colsB.length];
			int idx = 0;
			for (int l : colsA) for (int j : colsB) out[idx++] = l * n + j;
			return out;
		}

		public int[] cPositions() {
			int[] out = new int[rowsA.length * colsB.length];
			int idx = 0;
			for (int i : rowsA) for (int j : colsB) out[idx++] = i * n + j;
			return out;
		}

		public String label() {
			return "(rA=" + arr(rowsA) + " rB/cA=" + arr(colsA) + " cB=" + arr(colsB) + ")";
		}

		private static String arr(int[] a) {
			StringBuilder sb = new StringBuilder("{");
			for (int i = 0; i < a.length; i++) {
				if (i > 0) sb.append(",");
				sb.append(a[i]);
			}
			return sb.append("}").toString();
		}
	}

	/** All `C(n,2)³` distinct ⟨2,2,2⟩ sub-block embeddings of ⟨n,n,n⟩. */
	public static List<Subblock> enumerateAll222Subblocks(int n) {
		List<Subblock> out = new ArrayList<>();
		for (int[] rA : choose2(n)) {
			for (int[] cA : choose2(n)) {
				for (int[] cB : choose2(n)) {
					out.add(new Subblock(n, rA, cA, cB));
				}
			}
		}
		return out;
	}

	/** All ⟨2,2,3⟩-family embeddings of ⟨n,n,n⟩: union of ⟨2,2,3⟩, ⟨2,3,2⟩, ⟨3,2,2⟩. */
	public static List<Subblock> enumerateAll223FamilySubblocks(int n) {
		List<Subblock> out = new ArrayList<>();
		int[] allN = allOf(n);
		// ⟨2,2,3⟩: A 2×2, B 2×3, C 2×3 → rA size 2, cA size 2, cB size 3
		for (int[] rA : choose2(n))
			for (int[] cA : choose2(n))
				out.add(new Subblock(n, rA, cA, allN));
		// ⟨2,3,2⟩: A 2×3, B 3×2, C 2×2 → rA size 2, cA size 3, cB size 2
		for (int[] rA : choose2(n))
			for (int[] cB : choose2(n))
				out.add(new Subblock(n, rA, allN, cB));
		// ⟨3,2,2⟩: A 3×2, B 2×2, C 3×2 → rA size 3, cA size 2, cB size 2
		for (int[] cA : choose2(n))
			for (int[] cB : choose2(n))
				out.add(new Subblock(n, allN, cA, cB));
		return out;
	}

	/** All ⟨2,3,3⟩-family embeddings of ⟨n,n,n⟩: union of ⟨2,3,3⟩, ⟨3,2,3⟩, ⟨3,3,2⟩. */
	public static List<Subblock> enumerateAll233FamilySubblocks(int n) {
		List<Subblock> out = new ArrayList<>();
		int[] allN = allOf(n);
		// ⟨2,3,3⟩: only rA restricted to size 2
		for (int[] rA : choose2(n))
			out.add(new Subblock(n, rA, allN, allN));
		// ⟨3,2,3⟩: only cA restricted to size 2
		for (int[] cA : choose2(n))
			out.add(new Subblock(n, allN, cA, allN));
		// ⟨3,3,2⟩: only cB restricted to size 2
		for (int[] cB : choose2(n))
			out.add(new Subblock(n, allN, allN, cB));
		return out;
	}

	private static int[] allOf(int n) {
		int[] out = new int[n];
		for (int i = 0; i < n; i++) out[i] = i;
		return out;
	}

	/** Just the 4 "corner" sub-blocks of ⟨3,3,3⟩ (TL, TR, BL, BR). */
	public static List<Subblock> enumerateCornerSubblocks() {
		List<Subblock> out = new ArrayList<>();
		int[][] choices = { { 0, 1 }, { 1, 2 } };
		// TL = (rows {0,1}, cols/rows {0,1}, cols {0,1})
		// TR = (rows {0,1}, cols/rows {1,2}, cols {1,2})  — wait, this overlaps
		// Simplest: the 4 corners use consistent edge-touching:
		out.add(new Subblock(3, choices[0], choices[0], choices[0])); // TL
		out.add(new Subblock(3, choices[0], choices[0], choices[1])); // TR-ish
		out.add(new Subblock(3, choices[1], choices[1], choices[0])); // BL-ish
		out.add(new Subblock(3, choices[1], choices[1], choices[1])); // BR
		return out;
	}

	private static List<int[]> choose2(int n) {
		List<int[]> out = new ArrayList<>();
		for (int i = 0; i < n; i++) {
			for (int j = i + 1; j < n; j++) {
				out.add(new int[] { i, j });
			}
		}
		return out;
	}

	/**
	 * For a given rank-1 term k, returns true iff its restriction to the
	 * given sub-block is non-vanishing (i.e., at least one of u, v, w has a
	 * non-zero entry inside the sub-block's index set).
	 */
	public static boolean isTermActiveOnSubblock(BilinearAlgorithm alg, int k, Subblock sb) {
		return anyNonZero(alg.U, k, sb.aPositions())
				&& anyNonZero(alg.V, k, sb.bPositions())
				&& anyNonZero(alg.W, k, sb.cPositions());
	}

	private static boolean anyNonZero(double[][] factor, int k, int[] positions) {
		for (int p : positions) {
			if (factor[p][k] != 0.0) return true;
		}
		return false;
	}

	public static final class Analysis {
		public final List<Subblock> subblocks;
		public final boolean[][] contributionMatrix; // [term k][subblock s] → term k active on subblock s
		public final int[] termsPerSubblock; // length subblocks.size()
		public final int[] subblocksPerTerm; // length alg.r

		public Analysis(List<Subblock> sbs, boolean[][] mat, int[] tps, int[] spt) {
			this.subblocks = sbs;
			this.contributionMatrix = mat;
			this.termsPerSubblock = tps;
			this.subblocksPerTerm = spt;
		}
	}

	public static Analysis analyze(BilinearAlgorithm alg, List<Subblock> subblocks) {
		int r = alg.r;
		int s = subblocks.size();
		boolean[][] mat = new boolean[r][s];
		int[] tps = new int[s];
		int[] spt = new int[r];
		for (int k = 0; k < r; k++) {
			for (int j = 0; j < s; j++) {
				mat[k][j] = isTermActiveOnSubblock(alg, k, subblocks.get(j));
				if (mat[k][j]) {
					tps[j]++;
					spt[k]++;
				}
			}
		}
		return new Analysis(subblocks, mat, tps, spt);
	}
}
