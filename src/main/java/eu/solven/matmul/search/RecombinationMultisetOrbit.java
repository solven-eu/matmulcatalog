package eu.solven.matmul.search;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.LinkedHashMap;
import java.util.LinkedHashSet;
import java.util.List;
import java.util.Map;
import java.util.Set;
import java.util.TreeMap;

import eu.solven.matmul.NonCubicBilinearAlgorithm;

/**
 * Enumerate the distinct <b>recombination multisets</b> realisable by a base
 * matmul scheme {@code ⟨n,m,p⟩} over its entire {@code GL_n(ℚ)×GL_m(ℚ)×GL_p(ℚ)}
 * isotropy orbit — exactly, by symbolic enumeration rather than search.
 *
 * <h2>The object</h2>
 * Recombining a base {@code ⟨n,m,p⟩} at a block decomposition
 * {@code (n₁..n_n, m₁..m_m, p₁..p_p)} sends each of the {@code r} products to a
 * smaller matmul {@code ⟨a,b,c⟩}; the per-axis sub-dimension is the size of the
 * <i>largest</i> block that product touches (capped by the {@code min} of the two
 * relevant factor views — exactly {@link AnalyticalMaskSearch#shapesAt}). The
 * <b>multiset of those {@code r} sub-shapes</b> is the complete rank invariant of
 * plain additive recombination (see {@code references/MULTISET_FRONTIER.md}). It
 * is symbolic in the block sizes: with the blocks of an axis ordered descending
 * ({@code n₁≥n₂≥…}), each sub-dimension is one of the block sizes — so a product
 * is tagged, per axis, by the <b>block index</b> it lands on (0 = largest block).
 *
 * <h2>Why this is exact, not a search</h2>
 * By the isotropy action {@code U'_k = XᵀU_kYᵀ}, {@code V'_k = Y⁻ᵀV_kZᵀ},
 * {@code W'_k = X⁻¹W_kZ⁻¹}, the sub-dimension on each axis depends on <b>only
 * one</b> of {@code X,Y,Z} (an invertible factor on that axis kills no
 * rows/columns of the others), and only through that matrix's column/row
 * <b>directions</b>. So the three axes are independent and the realisable
 * multiset set is the product of three per-axis pattern sets, zipped by product
 * index. Each per-axis pattern set is piecewise-constant in the directions and
 * jumps only at the (finitely many) null-spaces of the fixed integer base
 * factors; a finite integer-direction sweep that is stable under refinement
 * therefore enumerates <i>every</i> realisable pattern. For a 2-part axis this is
 * provably complete (the null-spaces are the only critical directions, plus one
 * generic); for {@code ≥3}-part axes completeness is certified empirically by
 * direction-bound stability.
 *
 * <p>All arithmetic is exact integers (adjugates instead of inverses), so the
 * support (zero/non-zero) tests are exact. The base's factor coefficients must be
 * integers (the common case: Strassen, Winograd, Laderman, AlphaTensor-Z…).
 */
public final class RecombinationMultisetOrbit {

	private RecombinationMultisetOrbit() {}

	public static final class Result {
		/** Distinct canonical multiset keys (block-index encoding). */
		public final Set<String> canonicalMultisets = new LinkedHashSet<>();
		/** key → representative shape array {@code [r][3]} of block indices (0 = largest block). */
		public final Map<String, int[][]> representativeShapes = new LinkedHashMap<>();
		/** Per-axis count of distinct block-index patterns (diagnostics). */
		public int[] perAxisPatternCounts;
		public int dirBound;
		public long combinations;
		/** Base shape (number of parts per axis = axis dimension); fixes the canonicalising group. */
		public int n, m, p;

		/**
		 * The <b>dominance frontier</b> of {@link #canonicalMultisets}: the subset
		 * that is NOT sub-shape-dominated by another canonical multiset. Matmul rank
		 * is monotone (a bigger sub-multiply never costs less — pad the small one in)
		 * and the recombination cost is {@code Σ R(sub-shapeₖ)}; so if multiset A can
		 * be matched product-by-product to B with every A-block ≤ the paired B-block,
		 * then {@code cost(A) ≤ cost(B)} for <b>every</b> base recursion and B can
		 * never be rank-optimal. Pruning to this frontier is therefore <b>lossless</b>
		 * for the minimum and base/allocation-independent (block index 0 = largest =
		 * most expensive under any allocation).
		 *
		 * <p>Dominance is computed at the same {@code S}-quotient level as the keys
		 * themselves (the base's shape automorphisms): {@code A} dominates {@code B}
		 * iff <em>some</em> shape-symmetry image of {@code A} is product-by-product
		 * below {@code B}. This is the right level for symmetric / axis-permuted
		 * allocations and is matmul-rank exact there (rank is axis-permutation
		 * invariant). It is computed from the canonical KEYS (not the seed-dependent
		 * {@link #representativeShapes}), so it is seed-independent — e.g. ⟨2,2,2⟩
		 * yields 6 whether seeded from Strassen or Winograd. The finer
		 * <em>axis-tagged</em> frontier (no quotient, lossless for asymmetric
		 * allocations too) is strictly larger; see the paper's multiset section.
		 *
		 * @return frontier keys, a subset of {@link #canonicalMultisets}.
		 */
		public List<String> dominanceFrontier() {
			List<String> keys = new ArrayList<>(canonicalMultisets);
			int N = keys.size();
			int[][][] ms = new int[N][][];
			long[] cheapScore = new long[N];
			for (int i = 0; i < N; i++) {
				ms[i] = parseKey(keys.get(i)); // seed-independent: parse the canonical key
				long s = 0;
				for (int[] t : ms[i])
					for (int v : t)
						s += v; // larger block index = smaller block = cheaper
				cheapScore[i] = s;
			}
			int[][] stab = shapeStabilizer(n, m, p);
			// Process cheapest-first (highest score): a dominator always has score ≥
			// the dominated (axis-permutation preserves the score), so checking against
			// the kept frontier suffices — dominance is transitive (the stabiliser is a
			// group), so the chain's minimal element is always already kept.
			Integer[] order = new Integer[N];
			for (int i = 0; i < N; i++)
				order[i] = i;
			Arrays.sort(order, (a, b) -> Long.compare(cheapScore[b], cheapScore[a]));
			List<Integer> frontierIdx = new ArrayList<>();
			List<String> out = new ArrayList<>();
			for (int idx : order) {
				boolean dominated = false;
				outer: for (int f : frontierIdx) {
					for (int[] pi : stab) {
						if (dominatesBelow(permuteAxes(ms[f], pi), ms[idx])) {
							dominated = true;
							break outer;
						}
					}
				}
				if (!dominated) {
					frontierIdx.add(idx);
					out.add(keys.get(idx));
				}
			}
			return out;
		}

		/**
		 * The finer <b>axis-tagged</b> dominance frontier: dominance <em>without</em>
		 * the shape-symmetry quotient. It is lossless for asymmetric allocations too
		 * (where axis identity matters), and is strictly larger than
		 * {@link #dominanceFrontier()}. Returns axis-tagged multiset keys (the full
		 * shape-symmetry orbit of {@link #canonicalMultisets} is expanded first), so
		 * its members need not lie in {@link #canonicalMultisets}.
		 */
		public List<String> nonCanonicalFrontier() {
			int[][] stab = shapeStabilizer(n, m, p);
			java.util.LinkedHashSet<String> tagged = new java.util.LinkedHashSet<>();
			for (String k : canonicalMultisets)
				for (int[] pi : stab)
					tagged.add(normaliseKey(permuteAxes(parseKey(k), pi)));
			List<String> keys = new ArrayList<>(tagged);
			int N = keys.size();
			int[][][] ms = new int[N][][];
			long[] cheapScore = new long[N];
			for (int i = 0; i < N; i++) {
				ms[i] = parseKey(keys.get(i));
				long s = 0;
				for (int[] t : ms[i])
					for (int v : t)
						s += v;
				cheapScore[i] = s;
			}
			Integer[] order = new Integer[N];
			for (int i = 0; i < N; i++)
				order[i] = i;
			Arrays.sort(order, (a, b) -> Long.compare(cheapScore[b], cheapScore[a]));
			List<Integer> frontierIdx = new ArrayList<>();
			List<String> out = new ArrayList<>();
			for (int idx : order) {
				boolean dominated = false;
				for (int f : frontierIdx) {
					if (dominatesBelow(ms[f], ms[idx])) { // plain dominance, no axis relabel
						dominated = true;
						break;
					}
				}
				if (!dominated) {
					frontierIdx.add(idx);
					out.add(keys.get(idx));
				}
			}
			return out;
		}

		/** Count of distinct axis-tagged multisets (the shape-symmetry orbit of the canonical set). */
		public int nonCanonicalCount() {
			int[][] stab = shapeStabilizer(n, m, p);
			java.util.LinkedHashSet<String> tagged = new java.util.LinkedHashSet<>();
			for (String k : canonicalMultisets)
				for (int[] pi : stab)
					tagged.add(normaliseKey(permuteAxes(parseKey(k), pi)));
			return tagged.size();
		}
	}

	/** Parse a canonical key {@code "i,j,k|i,j,k|…"} into its {@code [r][3]} block-index shapes. */
	static int[][] parseKey(String key) {
		String[] parts = key.split("\\|");
		int[][] s = new int[parts.length][3];
		for (int k = 0; k < parts.length; k++) {
			String[] t = parts[k].split(",");
			for (int a = 0; a < 3; a++)
				s[k][a] = Integer.parseInt(t[a]);
		}
		return s;
	}

	/** Serialise an axis-tagged multiset to a sorted key {@code "i,j,k|…"} (order-independent). */
	static String normaliseKey(int[][] shapes) {
		String[] t = new String[shapes.length];
		for (int k = 0; k < shapes.length; k++)
			t[k] = shapes[k][0] + "," + shapes[k][1] + "," + shapes[k][2];
		Arrays.sort(t);
		return String.join("|", t);
	}

	/** Relabel the three axes of every product by permutation {@code pi}. */
	static int[][] permuteAxes(int[][] shapes, int[] pi) {
		int[][] out = new int[shapes.length][3];
		for (int k = 0; k < shapes.length; k++)
			for (int a = 0; a < 3; a++)
				out[k][a] = shapes[k][pi[a]];
		return out;
	}

	/**
	 * True iff multiset {@code a} is cheaper-or-equal to {@code b}: a perfect
	 * matching pairs each product of {@code a} to a distinct product of {@code b}
	 * with the {@code a}-block ≤ the paired {@code b}-block on every axis (i.e.
	 * {@code a}'s block index ≥ {@code b}'s, since a larger index is a smaller
	 * block). By matmul-rank monotonicity this implies {@code cost(a) ≤ cost(b)}.
	 */
	static boolean dominatesBelow(int[][] a, int[][] b) {
		int r = a.length;
		int[] matchOfB = new int[r];
		Arrays.fill(matchOfB, -1);
		for (int i = 0; i < r; i++) {
			boolean[] seen = new boolean[r];
			if (!augment(i, a, b, matchOfB, seen))
				return false;
		}
		return true;
	}

	private static boolean augment(int i, int[][] a, int[][] b, int[] matchOfB, boolean[] seen) {
		for (int j = 0; j < b.length; j++) {
			if (seen[j])
				continue;
			// a[i] lands on a block ≤ b[j] on every axis  ⇔  a-index ≥ b-index
			if (a[i][0] >= b[j][0] && a[i][1] >= b[j][1] && a[i][2] >= b[j][2]) {
				seen[j] = true;
				if (matchOfB[j] == -1 || augment(matchOfB[j], a, b, matchOfB, seen)) {
					matchOfB[j] = i;
					return true;
				}
			}
		}
		return false;
	}

	/**
	 * Enumerate all realisable canonical recombination multisets of {@code seed}.
	 *
	 * @param dirBound integer entry bound for the per-axis change-of-basis sweep
	 *                 (≥2 to include a generic direction). Use
	 *                 {@link #isStable} to certify completeness.
	 */
	public static Result enumerate(NonCubicBilinearAlgorithm seed, int dirBound) {
		int n = seed.n, m = seed.m, p = seed.p, r = seed.r;
		int[][][] U = reshape(seed.denseU(), r, n, m); // U[k] is n×m
		int[][][] V = reshape(seed.denseV(), r, m, p); // V[k] is m×p
		int[][][] W = reshape(seed.denseW(), r, n, p); // W[k] is n×p

		Set<int[]> nPat = axisPatterns(n, r, dirBound, (X, adjX, k) -> {
			int uMin = minRowNonZero_leftCols(X, U[k]); // rows of XᵀU_k
			int wMin = minRowNonZero_dualRows(adjX, W[k]); // rows of X⁻¹W_k
			return Math.max(uMin, wMin);
		});
		Set<int[]> mPat = axisPatterns(m, r, dirBound, (Y, adjY, k) -> {
			int uMin = minColNonZero_rightRows(U[k], Y); // cols of XᵀU_kYᵀ ← U_k·(rows of Y)
			int vMin = minRowNonZero_dualColsLeft(adjY, V[k]); // rows of Y⁻ᵀV_k
			return Math.max(uMin, vMin);
		});
		Set<int[]> pPat = axisPatterns(p, r, dirBound, (Z, adjZ, k) -> {
			int vMin = minColNonZero_rightRows(V[k], Z); // cols of …V_kZᵀ ← V_k·(rows of Z)
			int wMin = minColNonZero_dualCols(W[k], adjZ); // cols of …W_kZ⁻¹ ← W_k·(cols of adjZ)
			return Math.max(vMin, wMin);
		});

		int[][] stabilizer = shapeStabilizer(n, m, p);

		Result res = new Result();
		res.n = n;
		res.m = m;
		res.p = p;
		res.dirBound = dirBound;
		res.perAxisPatternCounts = new int[] { nPat.size(), mPat.size(), pPat.size() };
		for (int[] np : nPat)
			for (int[] mp : mPat)
				for (int[] pp : pPat) {
					res.combinations++;
					int[][] shapes = new int[r][3];
					for (int k = 0; k < r; k++) {
						shapes[k][0] = np[k];
						shapes[k][1] = mp[k];
						shapes[k][2] = pp[k];
					}
					String key = canonicalKey(shapes, stabilizer);
					if (res.canonicalMultisets.add(key)) res.representativeShapes.put(key, shapes);
				}
		return res;
	}

	/** True iff {@code enumerate} yields the same multiset set at {@code dirBound} and {@code dirBound+1}. */
	public static boolean isStable(NonCubicBilinearAlgorithm seed, int dirBound) {
		return enumerate(seed, dirBound).canonicalMultisets
				.equals(enumerate(seed, dirBound + 1).canonicalMultisets);
	}

	// ---- per-axis pattern enumeration ---------------------------------------

	private interface AxisIndex {
		/** block index (0 = largest) of product {@code k} given change-of-basis {@code M} and its adjugate. */
		int idx(int[][] M, int[][] adjM, int k);
	}

	/** Distinct block-index patterns (one {@code int[r]} per realisable pattern) over invertible M. */
	private static Set<int[]> axisPatterns(int dim, int r, int bound, AxisIndex fn) {
		Set<String> seenKeys = new java.util.HashSet<>();
		Set<int[]> patterns = new LinkedHashSet<>();
		int[][] M = new int[dim][dim];
		sweep(M, 0, 0, bound, () -> {
			if (determinant(M) == 0) return;
			int[][] adjM = adjugate(M);
			int[] pat = new int[r];
			for (int k = 0; k < r; k++) pat[k] = fn.idx(M, adjM, k);
			if (seenKeys.add(Arrays.toString(pat))) patterns.add(pat);
		});
		return patterns;
	}

	private interface Sink { void accept(); }

	/** Odometer over every {@code dim×dim} integer matrix with entries in {@code [-bound,bound]}. */
	private static void sweep(int[][] M, int i, int j, int bound, Sink sink) {
		int dim = M.length;
		if (i == dim) { sink.accept(); return; }
		int ni = (j + 1 == dim) ? i + 1 : i;
		int nj = (j + 1 == dim) ? 0 : j + 1;
		for (int v = -bound; v <= bound; v++) {
			M[i][j] = v;
			sweep(M, ni, nj, bound, sink);
		}
	}

	// minIndex helpers — block index = first index whose view is non-zero; the
	// recombination sub-dim is the max of the two relevant views' min-indices
	// (= the smaller block size, since a larger index is a smaller block).

	/** min row i with (col i of X)ᵀ·U_k ≠ 0  [rows of XᵀU_k]. */
	private static int minRowNonZero_leftCols(int[][] X, int[][] Uk) {
		int dim = X.length, cols = Uk[0].length;
		for (int i = 0; i < dim; i++) {
			for (int c = 0; c < cols; c++) {
				long s = 0;
				for (int a = 0; a < dim; a++) s += (long) X[a][i] * Uk[a][c];
				if (s != 0) return i;
			}
		}
		return dim - 1;
	}

	/** min row i with (row i of adjX)·W_k ≠ 0  [rows of X⁻¹W_k]. */
	private static int minRowNonZero_dualRows(int[][] adjX, int[][] Wk) {
		int dim = adjX.length, cols = Wk[0].length;
		for (int i = 0; i < dim; i++) {
			for (int c = 0; c < cols; c++) {
				long s = 0;
				for (int a = 0; a < dim; a++) s += (long) adjX[i][a] * Wk[a][c];
				if (s != 0) return i;
			}
		}
		return dim - 1;
	}

	/** min col j with U_k·(row j of Y) ≠ 0  [cols of XᵀU_kYᵀ]. */
	private static int minColNonZero_rightRows(int[][] Uk, int[][] Y) {
		int rows = Uk.length, dim = Y.length;
		for (int j = 0; j < dim; j++) {
			for (int rIdx = 0; rIdx < rows; rIdx++) {
				long s = 0;
				for (int a = 0; a < dim; a++) s += (long) Uk[rIdx][a] * Y[j][a];
				if (s != 0) return j;
			}
		}
		return dim - 1;
	}

	/** min row l with (col l of adjY)ᵀ·V_k ≠ 0  [rows of Y⁻ᵀV_k]. */
	private static int minRowNonZero_dualColsLeft(int[][] adjY, int[][] Vk) {
		int dim = adjY.length, cols = Vk[0].length;
		for (int l = 0; l < dim; l++) {
			for (int c = 0; c < cols; c++) {
				long s = 0;
				for (int a = 0; a < dim; a++) s += (long) adjY[a][l] * Vk[a][c];
				if (s != 0) return l;
			}
		}
		return dim - 1;
	}

	/** min col j with W_k·(col j of adjZ) ≠ 0  [cols of …W_kZ⁻¹]. */
	private static int minColNonZero_dualCols(int[][] Wk, int[][] adjZ) {
		int rows = Wk.length, dim = adjZ.length;
		for (int j = 0; j < dim; j++) {
			for (int rIdx = 0; rIdx < rows; rIdx++) {
				long s = 0;
				for (int a = 0; a < dim; a++) s += (long) Wk[rIdx][a] * adjZ[a][j];
				if (s != 0) return j;
			}
		}
		return dim - 1;
	}

	// ---- canonicalisation ----------------------------------------------------

	/** Axis permutations of {@code {0,1,2}} that fix the shape tuple (the base's automorphisms). */
	static int[][] shapeStabilizer(int n, int m, int p) {
		int[] s = { n, m, p };
		int[][] all = { { 0, 1, 2 }, { 0, 2, 1 }, { 1, 0, 2 }, { 1, 2, 0 }, { 2, 0, 1 }, { 2, 1, 0 } };
		List<int[]> keep = new ArrayList<>();
		for (int[] pi : all) {
			if (s[pi[0]] == s[0] && s[pi[1]] == s[1] && s[pi[2]] == s[2]) keep.add(pi);
		}
		return keep.toArray(new int[0][]);
	}

	/** Lex-min serialised multiset over the stabiliser group. */
	static String canonicalKey(int[][] shapes, int[][] stabilizer) {
		String best = null;
		for (int[] pi : stabilizer) {
			String[] s = new String[shapes.length];
			for (int k = 0; k < shapes.length; k++) {
				int[] sh = shapes[k];
				s[k] = sh[pi[0]] + "," + sh[pi[1]] + "," + sh[pi[2]];
			}
			Arrays.sort(s);
			String key = String.join("|", s);
			if (best == null || key.compareTo(best) < 0) best = key;
		}
		return best;
	}

	// ---- symbolic rendering --------------------------------------------------

	private static final String[] SUB = { "₁", "₂", "₃", "₄", "₅", "₆", "₇", "₈", "₉" };

	/**
	 * Render a canonical key symbolically in the block sizes, e.g.
	 * {@code 1·⟨n₁,m₁,p₁⟩ + 3·⟨n₂,m₁,p₁⟩ …} — block index {@code i} on an axis maps
	 * to that axis's {@code i+1}-th largest part. Pass axis letters (e.g.
	 * {@code "n","m","p"}); for cubic bases a single letter reads cleanest.
	 */
	public static String prettySymbolic(String key, String axisA, String axisB, String axisC) {
		String[] axes = { axisA, axisB, axisC };
		String[] parts = key.split("\\|");
		Map<String, Integer> counts = new TreeMap<>();
		for (String part : parts) {
			String[] idx = part.split(",");
			StringBuilder sb = new StringBuilder();
			for (int a = 0; a < 3; a++) {
				if (a > 0) sb.append(',');
				sb.append(axes[a]).append(SUB[Integer.parseInt(idx[a])]);
			}
			counts.merge(sb.toString(), 1, Integer::sum);
		}
		StringBuilder out = new StringBuilder();
		boolean first = true;
		for (var e : counts.entrySet()) {
			if (!first) out.append(" + ");
			out.append(e.getValue()).append("·⟨").append(e.getKey()).append('⟩');
			first = false;
		}
		return out.toString();
	}

	/** Render a canonical key with concrete block sizes, e.g. {@code sizes[axis] = {9,8}}. */
	public static String prettyConcrete(String key, int[][] sizes) {
		String[] parts = key.split("\\|");
		Map<String, Integer> counts = new TreeMap<>();
		for (String part : parts) {
			String[] idx = part.split(",");
			String sh = sizes[0][Integer.parseInt(idx[0])] + "," + sizes[1][Integer.parseInt(idx[1])]
					+ "," + sizes[2][Integer.parseInt(idx[2])];
			counts.merge(sh, 1, Integer::sum);
		}
		StringBuilder out = new StringBuilder();
		boolean first = true;
		for (var e : counts.entrySet()) {
			if (!first) out.append(" + ");
			out.append(e.getValue()).append("·⟨").append(e.getKey()).append('⟩');
			first = false;
		}
		return out.toString();
	}

	// ---- integer matrix utilities -------------------------------------------

	private static int[][][] reshape(double[][] factor, int r, int rows, int cols) {
		int[][][] out = new int[r][rows][cols];
		for (int k = 0; k < r; k++)
			for (int i = 0; i < rows; i++)
				for (int j = 0; j < cols; j++) {
					double v = factor[i * cols + j][k];
					long rv = Math.round(v);
					if (Math.abs(v - rv) > 1e-9)
						throw new IllegalArgumentException("base coefficient not integer: " + v);
					out[k][i][j] = (int) rv;
				}
		return out;
	}

	static long determinant(int[][] M) {
		int n = M.length;
		if (n == 1) return M[0][0];
		if (n == 2) return (long) M[0][0] * M[1][1] - (long) M[0][1] * M[1][0];
		long det = 0;
		for (int c = 0; c < n; c++) {
			long cof = (c % 2 == 0 ? 1 : -1) * M[0][c] * determinant(minor(M, 0, c));
			det += cof;
		}
		return det;
	}

	/** adjugate (transpose of the cofactor matrix); integer. */
	static int[][] adjugate(int[][] M) {
		int n = M.length;
		if (n == 1) return new int[][] { { 1 } };
		int[][] adj = new int[n][n];
		for (int i = 0; i < n; i++)
			for (int j = 0; j < n; j++) {
				long cof = ((i + j) % 2 == 0 ? 1 : -1) * determinant(minor(M, i, j));
				adj[j][i] = (int) cof; // transpose
			}
		return adj;
	}

	private static int[][] minor(int[][] M, int row, int col) {
		int n = M.length;
		int[][] out = new int[n - 1][n - 1];
		int ri = 0;
		for (int i = 0; i < n; i++) {
			if (i == row) continue;
			int ci = 0;
			for (int j = 0; j < n; j++) {
				if (j == col) continue;
				out[ri][ci++] = M[i][j];
			}
			ri++;
		}
		return out;
	}
}
