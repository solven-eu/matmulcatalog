package eu.solven.matmul.recombination;

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
		/**
		 * key → the integer GL transform {@code {X(n×n), Y(m×m), Z(p×p)}} that realises this
		 * multiset (only populated by the exact {@link #enumerate} path, not the structural one).
		 * {@link #materialise} applies it to the base to rebuild the concrete orbited scheme — the
		 * "1 scheme + transform config" so the frontier needn't persist a scheme per multiset.
		 */
		public final Map<String, int[][][]> representativeTransforms = new LinkedHashMap<>();
		/** Per-axis count of distinct block-index patterns (diagnostics). */
		public int[] perAxisPatternCounts;
		public int dirBound;
		public long combinations;
		/** Base shape (number of parts per axis = axis dimension); fixes the canonicalising group. */
		public int n, m, p;
		/**
		 * Whether {@link #dominanceFrontier()} is the EXACT, complete dominance antichain (true) or a
		 * possibly-INFLATED upper bound (false). True only for {@link #enumerate} (the complete
		 * canonical set). False for {@link #enumerateSampled}: dominance is computed over an INCOMPLETE
		 * sampled set, so a member whose true dominator was never sampled wrongly appears non-dominated
		 * — empirically 225 of 226 sampled ⟨2,3,3⟩ frontier members are actually dominated (exact
		 * frontier = 170). A sampled frontier is therefore a safe candidate SOURCE (every spurious
		 * member is still a valid GL support that can only tie, never beat, so build-verification stays
		 * sound) but must NOT be reported as "the frontier" or counted as exhaustive.
		 */
		public boolean frontierExact = false;

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

		// The three axes are independent (each sub-dimension depends on only one of
		// X,Y,Z — see class javadoc), so enumerate them concurrently. Each axis sweep
		// is itself internally parallel over the leading odometer entry.
		var fNPat = java.util.concurrent.CompletableFuture.supplyAsync(() -> axisPatterns(n, r, dirBound, (X, adjX, k) -> {
			int uMin = minRowNonZero_leftCols(X, U[k]); // rows of XᵀU_k
			int wMin = minRowNonZero_dualRows(adjX, W[k]); // rows of X⁻¹W_k
			return Math.max(uMin, wMin);
		}));
		var fMPat = java.util.concurrent.CompletableFuture.supplyAsync(() -> axisPatterns(m, r, dirBound, (Y, adjY, k) -> {
			int uMin = minColNonZero_rightRows(U[k], Y); // cols of XᵀU_kYᵀ ← U_k·(rows of Y)
			int vMin = minRowNonZero_dualColsLeft(adjY, V[k]); // rows of Y⁻ᵀV_k
			return Math.max(uMin, vMin);
		}));
		var fPPat = java.util.concurrent.CompletableFuture.supplyAsync(() -> axisPatterns(p, r, dirBound, (Z, adjZ, k) -> {
			int vMin = minColNonZero_rightRows(V[k], Z); // cols of …V_kZᵀ ← V_k·(rows of Z)
			int wMin = minColNonZero_dualCols(W[k], adjZ); // cols of …W_kZ⁻¹ ← W_k·(cols of adjZ)
			return Math.max(vMin, wMin);
		}));
		List<AxisRep> nPat = fNPat.join();
		List<AxisRep> mPat = fMPat.join();
		List<AxisRep> pPat = fPPat.join();

		int[][] stabilizer = shapeStabilizer(n, m, p);

		Result res = new Result();
		res.n = n;
		res.m = m;
		res.p = p;
		res.dirBound = dirBound;
		res.frontierExact = true; // complete canonical set ⇒ exact dominance antichain
		res.perAxisPatternCounts = new int[] { nPat.size(), mPat.size(), pPat.size() };
		for (AxisRep np : nPat)
			for (AxisRep mp : mPat)
				for (AxisRep pp : pPat) {
					res.combinations++;
					int[][] shapes = new int[r][3];
					for (int k = 0; k < r; k++) {
						shapes[k][0] = np.pattern[k];
						shapes[k][1] = mp.pattern[k];
						shapes[k][2] = pp.pattern[k];
					}
					String key = canonicalKey(shapes, stabilizer);
					if (res.canonicalMultisets.add(key)) {
						res.representativeShapes.put(key, shapes);
						// the (X,Y,Z) realising the RAW (np,mp,pp) → these shapes (key canonicalises them)
						res.representativeTransforms.put(key, new int[][][] { np.matrix, mp.matrix, pp.matrix });
					}
				}
		return res;
	}

	/**
	 * <b>Saturation-sampled enumeration</b> — the anytime path for bases where the exact GL odometer
	 * is intractable (any axis dim ≥ 4: {@code (2·bound+1)^(dim²)} blows up). Instead of sweeping
	 * every integer change-of-basis it SAMPLES random integer {@code (X,Y,Z)} per axis and keeps the
	 * distinct block-index patterns, stopping an axis once {@code maxTriesSinceNew} consecutive
	 * samples yield nothing new (saturation). Deterministic (fixed seed) for reproducibility. The
	 * frontier is therefore a <b>partial menu</b> — NOT certified complete (caller stamps
	 * {@code exhaustive=false}); but per [[Strassen-vs-Winograd]] it still exposes alternative
	 * supports the native does not.
	 *
	 * @param maxTriesSinceNew per-axis saturation budget (e.g. 100_000)
	 * @param entryBound       integer entry range {@code [-entryBound,entryBound]} for the random GL matrices
	 */
	/** Memory bound: stop sampling once this many distinct CANONICAL multisets are collected
	 *  (the dim-4 orbit is enormous; the dominanceFrontier is then computed over this capped set). */
	static int SAMPLED_MAX_CANONICAL = 6000;

	public static Result enumerateSampled(NonCubicBilinearAlgorithm seed, int maxTriesSinceNew, int entryBound) {
		int n = seed.n, m = seed.m, p = seed.p, r = seed.r;
		int[][][] U = reshape(seed.denseU(), r, n, m);
		int[][][] V = reshape(seed.denseV(), r, m, p);
		int[][][] W = reshape(seed.denseW(), r, n, p);
		int[][] stabilizer = shapeStabilizer(n, m, p);

		Result res = new Result();
		res.n = n; res.m = m; res.p = p; res.dirBound = entryBound;
		java.util.Random rng = new java.util.Random(917L + n * 31 + m * 7 + p);
		int since = 0;
		long total = 0;
		boolean first = true;
		while (since < maxTriesSinceNew && total < SAMPLED_MAX_TOTAL && res.canonicalMultisets.size() < SAMPLED_MAX_CANONICAL) {
			total++;
			// sample one GL point per axis (identity first → frontier ⊇ native)
			int[][] X = first ? eye(n) : randomInvertible(n, entryBound, rng);
			int[][] Y = first ? eye(m) : randomInvertible(m, entryBound, rng);
			int[][] Z = first ? eye(p) : randomInvertible(p, entryBound, rng);
			first = false;
			int[][] adjX = adjugate(X), adjY = adjugate(Y), adjZ = adjugate(Z);
			int[][] shapes = new int[r][3];
			for (int k = 0; k < r; k++) {
				shapes[k][0] = Math.max(minRowNonZero_leftCols(X, U[k]), minRowNonZero_dualRows(adjX, W[k]));
				shapes[k][1] = Math.max(minColNonZero_rightRows(U[k], Y), minRowNonZero_dualColsLeft(adjY, V[k]));
				shapes[k][2] = Math.max(minColNonZero_rightRows(V[k], Z), minColNonZero_dualCols(W[k], adjZ));
			}
			String key = canonicalKey(shapes, stabilizer);
			if (res.canonicalMultisets.add(key)) {
				res.representativeShapes.put(key, shapes);
				res.representativeTransforms.put(key, new int[][][] { deepCopy(X), deepCopy(Y), deepCopy(Z) });
				since = 0;
			} else {
				since++;
			}
		}
		res.combinations = total;
		res.perAxisPatternCounts = new int[] { res.canonicalMultisets.size(), 0, 0 };
		return res;
	}

	private static int[][] eye(int d) { int[][] I = new int[d][d]; for (int i = 0; i < d; i++) I[i][i] = 1; return I; }

	/** Random integer {@code d×d} matrix with entries in {@code [-bound,bound]} and nonzero det (rejection). */
	private static int[][] randomInvertible(int d, int bound, java.util.Random rng) {
		int[][] M = new int[d][d];
		for (int attempt = 0; attempt < 64; attempt++) {
			for (int i = 0; i < d; i++) for (int j = 0; j < d; j++) M[i][j] = rng.nextInt(2 * bound + 1) - bound;
			if (determinant(M) != 0) return M;
		}
		return eye(d); // degenerate fallback (vanishingly rare)
	}

	/** Hard per-axis total-sample backstop so high-dim sampling always terminates (saturation alone
	 *  may never trigger when the pattern set is huge, e.g. a dim-4 axis). */
	static long SAMPLED_MAX_TOTAL = 3_000_000L;
	/** Memory bound: stop an axis once this many distinct patterns are collected (a dim-4 axis can
	 *  have millions — storing them all OOMs). The frontier is then a partial menu (still ⊇ native). */
	static int SAMPLED_MAX_PATTERNS = 150_000;

	/** Random integer change-of-basis sampling for one axis, pruned ONLINE to the pointwise-maximal
	 *  antichain (a multiset using a dominated axis-pattern is itself dominated — sound for the
	 *  combined frontier and the only way to bound memory at dim-4). Stops at saturation
	 *  ({@code maxTriesSinceNew} no-new), the total-sample backstop, or the antichain cap. */
	private static List<AxisRep> axisPatternsSampled(int dim, int r, int maxTriesSinceNew, int entryBound, long seed, AxisIndex fn) {
		List<AxisRep> antichain = new ArrayList<>();
		java.util.Set<String> seenKeys = new java.util.HashSet<>();
		java.util.Random rng = new java.util.Random(seed * 1000003L + dim);
		int since = 0;
		long total = 0;
		int[][] M = new int[dim][dim];
		for (int i = 0; i < dim; i++) M[i][i] = 1; // seed with identity so the frontier always ⊇ native
		boolean first = true;
		while (since < maxTriesSinceNew && total < SAMPLED_MAX_TOTAL && antichain.size() < SAMPLED_MAX_PATTERNS) {
			total++;
			if (first) { first = false; } else
			for (int i = 0; i < dim; i++) for (int j = 0; j < dim; j++) M[i][j] = rng.nextInt(2 * entryBound + 1) - entryBound;
			if (determinant(M) == 0) { since++; continue; }
			int[][] adj = adjugate(M);
			int[] pat = new int[r];
			for (int k = 0; k < r; k++) pat[k] = fn.idx(M, adj, k);
			if (!seenKeys.add(Arrays.toString(pat))) { since++; continue; } // already considered this exact pattern
			boolean dominated = false;
			for (AxisRep q : antichain) if (dominatesPointwise(q.pattern, pat)) { dominated = true; break; }
			if (dominated) { since++; continue; }
			antichain.removeIf(q -> dominatesPointwise(pat, q.pattern));
			antichain.add(new AxisRep(pat, deepCopy(M)));
			since = 0;
		}
		return antichain;
	}

	/**
	 * Apply an integer GL transform {@code (X,Y,Z)} (from {@link Result#representativeTransforms})
	 * to a base, returning the orbited scheme — the isotropy action
	 * {@code U'_k = XᵀU_kYᵀ, V'_k = Y⁻ᵀV_kZᵀ, W'_k = X⁻¹W_kZ⁻¹}. The result computes the SAME
	 * ⟨n,m,p⟩ product at the same rank, but its support realises the target frontier multiset.
	 * Coefficients are rational (the inverses divide by det); {@code SchemeIO} stores them exactly.
	 */
	public static NonCubicBilinearAlgorithm materialise(NonCubicBilinearAlgorithm seed, int[][] X, int[][] Y, int[][] Z) {
		int n = seed.n, m = seed.m, p = seed.p, r = seed.r;
		int[][][] U = reshape(seed.denseU(), r, n, m);
		int[][][] V = reshape(seed.denseV(), r, m, p);
		int[][][] W = reshape(seed.denseW(), r, n, p);
		double[][] Xt = transposeI(X), Yt = transposeI(Y);
		double[][] Yinv_t = transpose(inverse(Y)), Xinv = inverse(X), Zt = transposeI(Z), Zinv = inverse(Z);
		double[][] Up = new double[n * m][r], Vp = new double[m * p][r], Wp = new double[n * p][r];
		for (int k = 0; k < r; k++) {
			double[][] uk = mul(mul(Xt, toD(U[k])), Yt);          // n×m
			double[][] vk = mul(mul(Yinv_t, toD(V[k])), Zt);      // m×p
			double[][] wk = mul(mul(Xinv, toD(W[k])), Zinv);      // n×p
			// snap floating residue to exact 0 so the support (nonzero pattern) is correct
			for (int i = 0; i < n; i++) for (int j = 0; j < m; j++) Up[i * m + j][k] = snap(uk[i][j]);
			for (int i = 0; i < m; i++) for (int j = 0; j < p; j++) Vp[i * p + j][k] = snap(vk[i][j]);
			for (int i = 0; i < n; i++) for (int j = 0; j < p; j++) Wp[i * p + j][k] = snap(wk[i][j]);
		}
		return new NonCubicBilinearAlgorithm(n, m, p, Up, Vp, Wp);
	}

	private static double snap(double x) { return Math.abs(x) < 1e-9 ? 0.0 : x; }
	private static double[][] toD(int[][] M) { double[][] d = new double[M.length][M[0].length]; for (int i = 0; i < M.length; i++) for (int j = 0; j < M[0].length; j++) d[i][j] = M[i][j]; return d; }
	private static double[][] transposeI(int[][] M) { int n = M.length, c = M[0].length; double[][] t = new double[c][n]; for (int i = 0; i < n; i++) for (int j = 0; j < c; j++) t[j][i] = M[i][j]; return t; }
	private static double[][] transpose(double[][] M) { int n = M.length, c = M[0].length; double[][] t = new double[c][n]; for (int i = 0; i < n; i++) for (int j = 0; j < c; j++) t[j][i] = M[i][j]; return t; }
	private static double[][] mul(double[][] A, double[][] B) { int n = A.length, kk = B.length, c = B[0].length; double[][] R = new double[n][c]; for (int i = 0; i < n; i++) for (int j = 0; j < c; j++) { double s = 0; for (int x = 0; x < kk; x++) s += A[i][x] * B[x][j]; R[i][j] = s; } return R; }
	/** Exact integer-matrix inverse via adjugate/det, as doubles (M small & integer ⇒ exact). */
	private static double[][] inverse(int[][] M) {
		long det = determinant(M);
		int[][] adj = adjugate(M);
		int d = M.length;
		double[][] inv = new double[d][d];
		for (int i = 0; i < d; i++) for (int j = 0; j < d; j++) inv[i][j] = adj[i][j] / (double) det;
		return inv;
	}

	/**
	 * <b>Structural (GL-free) enumeration.</b> Same {@link Result} as {@link #enumerate},
	 * but instead of an integer odometer over all of {@code GL(d)} it sweeps only flags
	 * whose columns are drawn from a finite <i>candidate-direction</i> set derived from the
	 * base's own factors: the footprint columns/rows that the index functions test against,
	 * plus the coordinate axes and small generic directions (for the "avoid the footprint"
	 * cells). The index pattern is piecewise-constant and jumps only when a flag step crosses
	 * a footprint subspace, so a representative flag per arrangement cell suffices.
	 *
	 * <p>Cost is {@code C·(C−1)···(C−d+1)} flags (C = #candidates ≈ O(r·d)) instead of
	 * {@code (2·bound+1)^(d²)} — e.g. ⟨2,4,4⟩ drops from {@code 9^16≈2e15} to ~1e8.
	 *
	 * <p><b>Optimality tier.</b> Provably complete for {@code d≤2} (footprint lines + one
	 * generic are the only critical directions). For {@code d≥3} it is
	 * <i>candidate-complete, not proven exhaustive</i> — treat as a <b>bound</b> until
	 * {@link #structuralMatchesGl} certifies it against the GL oracle for that base. The
	 * caller (and the test) cross-check the two on every shape where {@link #enumerate} is
	 * still tractable (⟨2,2,2⟩, ⟨2,2,3⟩, ⟨2,3,3⟩); a match there licenses trusting it where
	 * the odometer cannot run (⟨2,4,4⟩+).
	 */
	public static Result enumerateStructural(NonCubicBilinearAlgorithm seed, int genericBound) {
		int n = seed.n, m = seed.m, p = seed.p, r = seed.r;
		int[][][] U = reshape(seed.denseU(), r, n, m); // U[k] is n×m
		int[][][] V = reshape(seed.denseV(), r, m, p); // V[k] is m×p
		int[][][] W = reshape(seed.denseW(), r, n, p); // W[k] is n×p

		AxisCands cand = buildAxisCandidates(U, V, W, r, n, m, p, genericBound);

		// Flag orientation: the n-axis index reads COLUMNS of X (col_i = flag vector); the m/p
		// axes read ROWS of Y/Z (row_j = flag vector). So candidates are placed as columns for
		// n, as rows for m/p — otherwise the swept flag vectors are not the candidate directions.
		var fNPat = java.util.concurrent.CompletableFuture.supplyAsync(() -> axisPatternsStructural(n, r, cand.n, false, (X, adjX, k) ->
				Math.max(minRowNonZero_leftCols(X, U[k]), minRowNonZero_dualRows(adjX, W[k]))));
		var fMPat = java.util.concurrent.CompletableFuture.supplyAsync(() -> axisPatternsStructural(m, r, cand.m, true, (Y, adjY, k) ->
				Math.max(minColNonZero_rightRows(U[k], Y), minRowNonZero_dualColsLeft(adjY, V[k]))));
		var fPPat = java.util.concurrent.CompletableFuture.supplyAsync(() -> axisPatternsStructural(p, r, cand.p, true, (Z, adjZ, k) ->
				Math.max(minColNonZero_rightRows(V[k], Z), minColNonZero_dualCols(W[k], adjZ))));
		Set<int[]> nPat = fNPat.join();
		Set<int[]> mPat = fMPat.join();
		Set<int[]> pPat = fPPat.join();

		int[][] stabilizer = shapeStabilizer(n, m, p);
		Result res = new Result();
		res.n = n; res.m = m; res.p = p; res.dirBound = genericBound;
		res.perAxisPatternCounts = new int[] { nPat.size(), mPat.size(), pPat.size() };
		for (int[] np : nPat)
			for (int[] mp : mPat)
				for (int[] pp : pPat) {
					res.combinations++;
					int[][] shapes = new int[r][3];
					for (int k = 0; k < r; k++) { shapes[k][0] = np[k]; shapes[k][1] = mp[k]; shapes[k][2] = pp[k]; }
					String key = canonicalKey(shapes, stabilizer);
					if (res.canonicalMultisets.add(key)) res.representativeShapes.put(key, shapes);
				}
		return res;
	}

	/** True iff structural and GL enumeration yield the same canonical multiset set (the oracle check). */
	public static boolean structuralMatchesGl(NonCubicBilinearAlgorithm seed, int glBound, int genericBound) {
		return enumerate(seed, glBound).canonicalMultisets
				.equals(enumerateStructural(seed, genericBound).canonicalMultisets);
	}

	/**
	 * <b>Frontier-only structural enumeration — the {@code d=4}-capable path.</b> Prunes each
	 * axis to its pointwise-maximal antichain ({@link #axisFrontierStructural}) before combining,
	 * so the intractable full canonical set (the {@code |nPat|·|mPat|·|pPat|} triple product that
	 * OOMs at {@code d=4}) is never materialised. Returns the {@link Result#dominanceFrontier()}
	 * of the combined per-axis frontiers — identical to {@code enumerate(...).dominanceFrontier()}
	 * wherever the GL oracle is tractable (validated on ⟨2,2,2⟩, ⟨2,2,3⟩, ⟨2,3,3⟩).
	 *
	 * <p>The returned {@link Result#canonicalMultisets} holds the combined per-axis-frontier
	 * candidates (a frontier superset), NOT the full canonical set — read {@link Result#dominanceFrontier()}.
	 */
	public static Result enumerateStructuralFrontier(NonCubicBilinearAlgorithm seed, int genericBound) {
		int n = seed.n, m = seed.m, p = seed.p, r = seed.r;
		int[][][] U = reshape(seed.denseU(), r, n, m);
		int[][][] V = reshape(seed.denseV(), r, m, p);
		int[][][] W = reshape(seed.denseW(), r, n, p);

		AxisCands cand = buildAxisCandidates(U, V, W, r, n, m, p, genericBound);

		var fN = java.util.concurrent.CompletableFuture.supplyAsync(() -> axisFrontierStructural(n, r, cand.n, false, (X, adjX, k) ->
				Math.max(minRowNonZero_leftCols(X, U[k]), minRowNonZero_dualRows(adjX, W[k]))));
		var fM = java.util.concurrent.CompletableFuture.supplyAsync(() -> axisFrontierStructural(m, r, cand.m, true, (Y, adjY, k) ->
				Math.max(minColNonZero_rightRows(U[k], Y), minRowNonZero_dualColsLeft(adjY, V[k]))));
		var fP = java.util.concurrent.CompletableFuture.supplyAsync(() -> axisFrontierStructural(p, r, cand.p, true, (Z, adjZ, k) ->
				Math.max(minColNonZero_rightRows(V[k], Z), minColNonZero_dualCols(W[k], adjZ))));
		List<int[]> nPat = fN.join(), mPat = fM.join(), pPat = fP.join();

		int[][] stabilizer = shapeStabilizer(n, m, p);
		Result res = new Result();
		res.n = n; res.m = m; res.p = p; res.dirBound = genericBound;
		res.perAxisPatternCounts = new int[] { nPat.size(), mPat.size(), pPat.size() };
		for (int[] np : nPat)
			for (int[] mp : mPat)
				for (int[] pp : pPat) {
					res.combinations++;
					int[][] shapes = new int[r][3];
					for (int k = 0; k < r; k++) { shapes[k][0] = np[k]; shapes[k][1] = mp[k]; shapes[k][2] = pp[k]; }
					String key = canonicalKey(shapes, stabilizer);
					if (res.canonicalMultisets.add(key)) res.representativeShapes.put(key, shapes);
				}
		return res;
	}

	/** Column {@code c} of an {@code rows×cols} integer matrix, as a length-{@code rows} vector. */
	private static int[] col(int[][] M, int c, int rows) {
		int[] v = new int[rows];
		for (int i = 0; i < rows; i++) v[i] = M[i][c];
		return v;
	}

	/**
	 * Finite candidate flag-direction set for one axis: the footprint directions (deduped
	 * by ray), the coordinate axes, and every nonzero {@code {−g..g}^dim} generic (deduped
	 * by ray). Zero and parallel duplicates are dropped — only the ray matters to a flag.
	 */
	/** The three per-axis candidate-direction lists for a base. */
	private static final class AxisCands { List<int[]> n, m, p; }

	/** Build per-axis candidate flag directions (arrangement-aware) for both structural entry points. */
	private static AxisCands buildAxisCandidates(int[][][] U, int[][][] V, int[][][] W, int r, int n, int m, int p, int gen) {
		// Per axis, per product: the A-constraint set (its span's COMPLEMENT is the uMin subspace
		// A_k) and the S-constraint set (its SPAN is the wMin subspace S_k). See the axis table in
		// enumerate(): n-axis ← cols(U),cols(W); m-axis ← rows(U),cols(V); p-axis ← rows(V),rows(W).
		List<List<int[]>> nA = new ArrayList<>(), nS = new ArrayList<>();
		List<List<int[]>> mA = new ArrayList<>(), mS = new ArrayList<>();
		List<List<int[]>> pA = new ArrayList<>(), pS = new ArrayList<>();
		for (int k = 0; k < r; k++) {
			List<int[]> nAc = new ArrayList<>(), nSc = new ArrayList<>();
			for (int c = 0; c < m; c++) nAc.add(col(U[k], c, n)); // cols of U_k
			for (int c = 0; c < p; c++) nSc.add(col(W[k], c, n)); // cols of W_k
			nA.add(nAc); nS.add(nSc);

			List<int[]> mAc = new ArrayList<>(), mSc = new ArrayList<>();
			for (int i = 0; i < n; i++) mAc.add(U[k][i].clone()); // rows of U_k
			for (int c = 0; c < p; c++) mSc.add(col(V[k], c, m)); // cols of V_k
			mA.add(mAc); mS.add(mSc);

			List<int[]> pAc = new ArrayList<>(), pSc = new ArrayList<>();
			for (int l = 0; l < m; l++) pAc.add(V[k][l].clone()); // rows of V_k
			for (int i = 0; i < n; i++) pSc.add(W[k][i].clone()); // rows of W_k
			pA.add(pAc); pS.add(pSc);
		}
		AxisCands c = new AxisCands();
		c.n = candidateFlagDirections(n, gen, nA, nS);
		c.m = candidateFlagDirections(m, gen, mA, mS);
		c.p = candidateFlagDirections(p, gen, pA, pS);
		return c;
	}

	/**
	 * Candidate flag directions for one axis. Includes: the footprint columns themselves (land-ON
	 * cells), the coordinate axes, provably-generic Vandermonde directions (avoid every footprint),
	 * and — the part a dense cube only approximated — a generic basis of every member of the
	 * footprint <b>arrangement lattice</b>: the {@code A_k = constraintᗮ} closed under intersection
	 * (so a flag prefix can sit inside any ∩A_k, raising {@code uMin}) and the {@code S_k = span}
	 * closed under sum (so a flag tail can contain any ΣS_k, raising {@code wMin}).
	 */
	private static List<int[]> candidateFlagDirections(int dim, int genericBound,
			List<List<int[]>> aConstraints, List<List<int[]>> sConstraints) {
		Map<String, int[]> byRay = new LinkedHashMap<>();
		int maxCoef = 1;
		for (List<int[]> set : aConstraints) for (int[] v : set) { addRay(byRay, v); for (int x : v) maxCoef = Math.max(maxCoef, Math.abs(x)); }
		for (List<int[]> set : sConstraints) for (int[] v : set) { addRay(byRay, v); for (int x : v) maxCoef = Math.max(maxCoef, Math.abs(x)); }
		for (int i = 0; i < dim; i++) { int[] e = new int[dim]; e[i] = 1; addRay(byRay, e); }
		int t = maxCoef + 1 + Math.max(0, genericBound - 1);

		// A-side lattice: complements closed under intersection (raise uMin). Generic basis of each
		// member so a flag prefix can sit generically inside it. (A joint ∩+ closure of both sides
		// would capture the residual mixed cells too, but it ballooned the candidate set and
		// regressed even d=3 — reverted; per-side is fast and frontier-exact for AT-Z ⟨2,2,3⟩.)
		List<List<int[]>> aSubs = new ArrayList<>();
		for (List<int[]> set : aConstraints) { List<int[]> nul = SubspaceArrangement.nullspace(set, dim); if (!nul.isEmpty()) aSubs.add(nul); }
		addLatticeGenerics(byRay, SubspaceArrangement.closure(aSubs, dim, true), dim, t);
		// S-side lattice: spans closed under sum (raise wMin). Generic basis of each member for the tail.
		List<List<int[]>> sSubs = new ArrayList<>();
		for (List<int[]> set : sConstraints) { List<int[]> b = SubspaceArrangement.basisOf(set, dim); if (!b.isEmpty()) sSubs.add(b); }
		addLatticeGenerics(byRay, SubspaceArrangement.closure(sSubs, dim, false), dim, t);

		// Full-space Vandermonde generics (forward + reversed nodes) for the generic-position cells.
		for (int s = 0; s < 2 * dim + 2; s++) {
			int tt = t + s;
			int[] vandF = new int[dim], vandR = new int[dim];
			long pw = 1;
			for (int i = 0; i < dim; i++) { vandF[i] = (int) pw; vandR[dim - 1 - i] = (int) pw; pw *= tt; }
			addRay(byRay, vandF); addRay(byRay, vandR);
		}
		return new ArrayList<>(byRay.values());
	}

	/** For each lattice subspace add a generic basis (dim(member) independent generic-in vectors). */
	private static void addLatticeGenerics(Map<String, int[]> byRay, List<List<int[]>> members, int dim, int t) {
		for (List<int[]> sub : members) {
			if (sub.isEmpty()) continue;
			for (int[] b : sub) addRay(byRay, b); // the canonical basis vectors (land-on directions)
			for (int s = 0; s <= sub.size(); s++) addRay(byRay, SubspaceArrangement.genericIn(sub, dim, t + s)); // generic basis within
		}
	}

	/** Add a vector to the ray-map under its canonical (sign- and gcd-normalised) key; skip zero. */
	private static void addRay(Map<String, int[]> byRay, int[] v) {
		if (v == null) return;
		int gg = 0, firstSign = 0;
		for (int x : v) { gg = gcd(gg, Math.abs(x)); if (firstSign == 0 && x != 0) firstSign = x > 0 ? 1 : -1; }
		if (gg == 0) return; // zero vector
		int[] red = new int[v.length];
		for (int i = 0; i < v.length; i++) red[i] = firstSign * v[i] / gg;
		byRay.putIfAbsent(Arrays.toString(red), red);
	}

	private static int gcd(int a, int b) { while (b != 0) { int t = a % b; a = b; b = t; } return Math.abs(a); }

	/**
	 * Structural per-axis sweep: enumerate distinct index patterns over all ordered,
	 * independent {@code dim}-tuples of {@code candidates} (used as the columns of the
	 * change-of-basis M). No integer odometer over GL.
	 */
	private static Set<int[]> axisPatternsStructural(int dim, int r, List<int[]> candidates, boolean asRows, AxisIndex fn) {
		Map<String, int[]> patterns = new LinkedHashMap<>();
		int[][] cols = new int[dim][];
		structuralRec(0, dim, r, candidates, cols, asRows, fn, patterns, null);
		return new LinkedHashSet<>(patterns.values());
	}

	/**
	 * Like {@link #axisPatternsStructural} but stores ONLY the pointwise-maximal antichain of
	 * index patterns ({@code A} kept iff no other {@code B} has {@code B[k] ≥ A[k] ∀k} — a
	 * larger index is a smaller, cheaper block). Sound for the combined frontier: a multiset
	 * using a per-axis-dominated pattern is itself dominated (swap in the dominating axis pattern
	 * — identity matching makes it cheaper on this axis, equal on the others). This keeps memory
	 * bounded at {@code d=4}, where the full per-axis set is intractable to store.
	 */
	private static List<int[]> axisFrontierStructural(int dim, int r, List<int[]> candidates, boolean asRows, AxisIndex fn) {
		List<int[]> antichain = new ArrayList<>();
		int[][] cols = new int[dim][];
		structuralRec(0, dim, r, candidates, cols, asRows, fn, null, antichain);
		return antichain;
	}

	private static void structuralRec(int depth, int dim, int r, List<int[]> cand, int[][] cols, boolean asRows,
			AxisIndex fn, Map<String, int[]> patterns, List<int[]> antichain) {
		if (depth == dim) {
			int[][] M = new int[dim][dim];
			// cols[j] is the j-th chosen flag vector. Placed as columns (n-axis: index reads
			// col_i(M)) or as rows (m/p-axes: index reads row_j(M)).
			for (int j = 0; j < dim; j++) for (int i = 0; i < dim; i++) {
				if (asRows) M[j][i] = cols[j][i]; else M[i][j] = cols[j][i];
			}
			if (determinant(M) == 0) return;
			int[][] adjM = adjugate(M);
			int[] pat = new int[r];
			for (int k = 0; k < r; k++) pat[k] = fn.idx(M, adjM, k);
			if (patterns != null) {
				patterns.putIfAbsent(Arrays.toString(pat), pat);
			} else {
				insertMaximal(antichain, pat); // keep only the pointwise-maximal antichain
			}
			return;
		}
		for (int[] c : cand) {
			cols[depth] = c;
			if (independentPrefix(cols, depth + 1, dim)) structuralRec(depth + 1, dim, r, cand, cols, asRows, fn, patterns, antichain);
		}
	}

	/** Insert {@code p} into a pointwise-maximal antichain: skip if dominated, else drop those it dominates. */
	private static void insertMaximal(List<int[]> antichain, int[] p) {
		for (int[] q : antichain) if (dominatesPointwise(q, p)) return; // q ⪰ p → p not maximal
		antichain.removeIf(q -> dominatesPointwise(p, q)); // p ⪰ q → q no longer maximal
		antichain.add(p.clone());
	}

	/** True iff {@code a[k] ≥ b[k]} for all k (a is pointwise no-more-expensive than b). */
	private static boolean dominatesPointwise(int[] a, int[] b) {
		for (int k = 0; k < a.length; k++) if (a[k] < b[k]) return false;
		return true;
	}

	/** True iff the first {@code len} chosen columns are linearly independent (rational rank == len). */
	private static boolean independentPrefix(int[][] cols, int len, int dim) {
		double[][] A = new double[dim][len];
		for (int j = 0; j < len; j++) for (int i = 0; i < dim; i++) A[i][j] = cols[j][i];
		int rank = 0;
		for (int j = 0; j < len; j++) {
			int piv = -1;
			for (int i = rank; i < dim; i++) if (Math.abs(A[i][j]) > 1e-9) { piv = i; break; }
			if (piv < 0) return false; // column j dependent on earlier ones
			double[] tmp = A[piv]; A[piv] = A[rank]; A[rank] = tmp;
			for (int i = 0; i < dim; i++) if (i != rank) { double f = A[i][j] / A[rank][j]; for (int c = j; c < len; c++) A[i][c] -= f * A[rank][c]; }
			rank++;
		}
		return rank == len;
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

	/**
	 * Distinct block-index patterns (one {@code int[r]} per realisable pattern) over
	 * invertible M. The odometer over every {@code dim×dim} integer matrix is
	 * {@code (2·bound+1)^(dim²)} — embarrassingly parallel, so we split on the leading
	 * entry {@code M[0][0]} and sweep the rest in {@code (2·bound+1)} concurrent branches,
	 * each with a thread-local dedup map, merged at the end. (Constant — core-count —
	 * speedup only; it does <i>not</i> beat the {@code dim²} exponent: a dim-4 axis at
	 * bound 4 is {@code 9^16≈2e15} regardless. That needs critical-direction enumeration,
	 * not more threads.)
	 */
	/** A realisable axis block-index pattern together with one integer change-of-basis that realises it. */
	static final class AxisRep {
		final int[] pattern;
		final int[][] matrix; // the dim×dim integer M whose flag induces {@code pattern}
		AxisRep(int[] pattern, int[][] matrix) { this.pattern = pattern; this.matrix = matrix; }
	}

	private static List<AxisRep> axisPatterns(int dim, int r, int bound, AxisIndex fn) {
		Map<String, AxisRep> merged = java.util.stream.IntStream.rangeClosed(-bound, bound).parallel()
				.mapToObj(v0 -> {
					int[][] M = new int[dim][dim];
					M[0][0] = v0;
					Map<String, AxisRep> local = new java.util.HashMap<>();
					Sink eval = () -> {
						if (determinant(M) == 0) return;
						int[][] adjM = adjugate(M);
						int[] pat = new int[r];
						for (int k = 0; k < r; k++) pat[k] = fn.idx(M, adjM, k);
						String key = Arrays.toString(pat);
						// keep one representative matrix per pattern (deep-copy M, it is mutated by the sweep)
						local.computeIfAbsent(key, kk -> new AxisRep(pat, deepCopy(M)));
					};
					if (dim == 1) {
						eval.accept(); // 1×1: the single entry is already pinned
					} else {
						sweep(M, 0, 1, bound, eval); // continue the odometer past the pinned M[0][0]
					}
					return local;
				})
				.collect(java.util.HashMap::new, Map::putAll, Map::putAll);
		return new ArrayList<>(merged.values());
	}

	private static int[][] deepCopy(int[][] M) {
		int[][] c = new int[M.length][];
		for (int i = 0; i < M.length; i++) c[i] = M[i].clone();
		return c;
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
